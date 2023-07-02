#ifndef ABSTRACTTRANSPORTODESYSTEM_HPP_
#define ABSTRACTTRANSPORTODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractTransportReactionSystem.hpp"
#include "TransportOdeSystemInformation.hpp"
#include "StateVariableRegister.hpp"

/**
 * ODE system to model transport phenomena across a membrane seaprating two 
 * compartments. The compartments external (bulk) and internal (cell) are 
 * recorded separately.
 */
class AbstractTransportOdeSystem: public AbstractOdeSystem
{
private:

    AbstractTransportReactionSystem* mpTransportReactionSystem;

    StateVariableRegister* mPdeStateVariableRegister;

    std::vector<double> mReservedCellConcentration;

    std::vector<double> mChangeCellConc;

    unsigned mNumSpecies;

    unsigned mNumReactions;

    double mDeltaError;

    bool mIsCheckConcentration;

public:

    AbstractTransportOdeSystem(AbstractTransportReactionSystem* pReactionSystem = new AbstractTransportReactionSystem());

    virtual ~AbstractTransportOdeSystem()
    {
    };
    
    AbstractTransportOdeSystem(const AbstractTransportOdeSystem&);

    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    virtual void UpdateReservedCellConcentration(std::vector<double>);

    virtual void UpdateChangeCellConcentration(std::vector<double>);

    unsigned GetNumSpecies();

    unsigned GetNumReactions();

    AbstractTransportReactionSystem* GetReactionSystem();

    void SetNumSpecies(unsigned);

    void SetNumReactions(unsigned);

    void SetReactionSystem(AbstractTransportReactionSystem* );
    
    void SetReservedCellConcentration(std::vector<double>);
    
    void SetChangeCellConcentration(std::vector<double>);

    void UpdateReservedCellConcentrationByIndex(unsigned, double);

    void UpdateReservedCellConcentrationByName(std::string, double);

    void UpdateChangeCellConcentrationByIndex(unsigned, double);

    void UpdateChangeCellConcentrationByName(std::string, double);

    void SetDeltaError(double);

    double GetDeltaError();

    void SetPdeStateRegister(StateVariableRegister*);

    void CheckConcentration(const std::vector<double>&);

    void SetIsCheckConcentration(bool);

    bool GetIsCheckConcentration();

    std::vector<double> GetReservedCellConcentration();
   
    std::vector<double> GetChangeCellConcentration();

    StateVariableRegister* GetPdeStateVariableRegister();

    double RetrieveReservedCellConcentrationByIndex(unsigned);

    double RetrieveReservedCellConcentrationByName(std::string);

    double RetrieveChangeCellConcentrationByIndex(unsigned);

    double RetrieveChangeCellConcentrationByName(std::string);
};

AbstractTransportOdeSystem::AbstractTransportOdeSystem(AbstractTransportReactionSystem* pReactionSystem)
    : AbstractOdeSystem(pReactionSystem->GetNumBulkStates() + pReactionSystem->GetNumCellStates()),
      mpTransportReactionSystem(pReactionSystem)
{
    SetReactionSystem(mpTransportReactionSystem);
    SetNumSpecies(mpTransportReactionSystem->GetNumBulkStates()+mpTransportReactionSystem->GetNumCellStates());
    SetNumReactions(mpTransportReactionSystem->GetNumReactions());

    // Initialsie cell concentration vectors
    std::vector<double> cellConcentration(mpTransportReactionSystem->GetNumCellStates(), 1.0);
    SetReservedCellConcentration(cellConcentration);
    std::vector<double> cellConcentrationChange(mpTransportReactionSystem->GetNumCellStates(), 0.0);
    SetChangeCellConcentration(cellConcentrationChange);

    mpSystemInfo.reset(new TransportOdeSystemInformation<AbstractTransportOdeSystem>(mpTransportReactionSystem));

    SetStateVariables(GetInitialConditions()); // set up in the solver call; size of NumBulkStates + NumCellStates

    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractTransportOdeSystem::AbstractTransportOdeSystem(const AbstractTransportOdeSystem& existingSystem)
    : AbstractOdeSystem(existingSystem.mNumSpecies)
{
    mpTransportReactionSystem = existingSystem.mpTransportReactionSystem;
    mNumSpecies = existingSystem.mNumSpecies;
    mNumReactions = existingSystem.mNumReactions;
    mDeltaError = existingSystem.mDeltaError;
    mIsCheckConcentration = existingSystem.mIsCheckConcentration;
}

void AbstractTransportOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    unsigned num_bulk_states = mpTransportReactionSystem->GetNumBulkStates();
    unsigned num_cell_states = mpTransportReactionSystem->GetNumCellStates();

    std::vector<double> bulkConcentrations(num_bulk_states, 0.0);
    std::vector<double> changeBulkConcentrations(num_bulk_states, 0.0);

    std::vector<double> cellConcentrations(num_cell_states, 0.0);
    std::vector<double> changeCellConcentrations(num_cell_states, 0.0);
    
    // Parse the different concentration state vectors from the solver output to the different concentration systems
    std::string this_state = "";
    for (unsigned i=0; i<num_bulk_states; ++i)
    {
        bulkConcentrations[i] = rY.at(i);
    }

    for (unsigned i=0; i<num_cell_states; ++i)
    {
        cellConcentrations[i] = rY.at(i+num_bulk_states);
    }
    
    // Reset rDY
    for (unsigned i=0; i<rDY.size();++i)
    {
        rDY[i] = 0.0;
    }
    
    if (mIsCheckConcentration)
    {   
        CheckConcentration(bulkConcentrations);
        CheckConcentration(cellConcentrations);
    }

    mpTransportReactionSystem->ReactSystem(bulkConcentrations, cellConcentrations, changeBulkConcentrations, changeCellConcentrations);

    // Reform rDY for passing to the solver
    for (unsigned i=0; i<num_bulk_states; ++i)
    {
        rDY.at(i)= changeBulkConcentrations.at(i);
    }   

    for (unsigned i=0; i<num_cell_states; ++i)
    {
        rDY.at(i+num_bulk_states) = changeCellConcentrations.at(i);
    }
}  

void AbstractTransportOdeSystem::UpdateReservedCellConcentration(std::vector<double> reservedConcentration)
{
    SetReservedCellConcentration(reservedConcentration);
}

void AbstractTransportOdeSystem::UpdateChangeCellConcentration(std::vector<double> changeCellConc)
{
    mChangeCellConc = changeCellConc;
}

unsigned AbstractTransportOdeSystem::GetNumSpecies()
{
    return mNumSpecies;
}

unsigned AbstractTransportOdeSystem::GetNumReactions()
{
    return mNumReactions;
}

AbstractTransportReactionSystem* AbstractTransportOdeSystem::GetReactionSystem()
{
    return mpTransportReactionSystem;
}

void AbstractTransportOdeSystem::SetNumSpecies(unsigned numSpecies)
{
    mNumSpecies = numSpecies;
}

void AbstractTransportOdeSystem::SetNumReactions(unsigned numReactions)
{
    mNumReactions = numReactions;
}

void AbstractTransportOdeSystem::SetReactionSystem(AbstractTransportReactionSystem* p_reactionSystem)
{
    mpTransportReactionSystem = p_reactionSystem;
}

void AbstractTransportOdeSystem::SetReservedCellConcentration(std::vector<double> cellConcentration)
{
    mReservedCellConcentration = cellConcentration;
}
    
void AbstractTransportOdeSystem::SetChangeCellConcentration(std::vector<double> changeCellConcentration)
{
    mChangeCellConc = changeCellConcentration;
}

void AbstractTransportOdeSystem::UpdateReservedCellConcentrationByIndex(unsigned index, double concentration)
{
    mReservedCellConcentration[index] = concentration;
}

void AbstractTransportOdeSystem::UpdateReservedCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpTransportReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);

    mReservedCellConcentration[index] = concentration;
}

void AbstractTransportOdeSystem::UpdateChangeCellConcentrationByIndex(unsigned index, double concentration)
{
    mChangeCellConc[index] = concentration;
}

void AbstractTransportOdeSystem::UpdateChangeCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpTransportReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);

    mChangeCellConc[index] = concentration;
}

void AbstractTransportOdeSystem::SetDeltaError(double deltaError)
{
    mDeltaError = deltaError;
}

double AbstractTransportOdeSystem::GetDeltaError()
{
    return mDeltaError;
}

void AbstractTransportOdeSystem::SetPdeStateRegister(StateVariableRegister* pRegister)
{
    mPdeStateVariableRegister = pRegister;
}

void AbstractTransportOdeSystem::CheckConcentration(const std::vector<double>& rY)
{
    // if chemical concentration gets too low then NaNs can occur, concentration must be +ve
    for (unsigned i=0; i<rY.size(); ++i)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if (rY[i] < mDeltaError)
        {
            const_cast<double&>(rY[i]) = 0;
        }
    }
}

void AbstractTransportOdeSystem::SetIsCheckConcentration(bool isCheckConcentration)
{
    mIsCheckConcentration = isCheckConcentration;
}

bool AbstractTransportOdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

std::vector<double> AbstractTransportOdeSystem::GetReservedCellConcentration()
{
    return mReservedCellConcentration;
}

std::vector<double> AbstractTransportOdeSystem::GetChangeCellConcentration()
{
    return mChangeCellConc;
}

StateVariableRegister* AbstractTransportOdeSystem::GetPdeStateVariableRegister()
{
    return mPdeStateVariableRegister;
}

double AbstractTransportOdeSystem::RetrieveReservedCellConcentrationByIndex(unsigned index)
{
    return mReservedCellConcentration[index];
}

double AbstractTransportOdeSystem::RetrieveReservedCellConcentrationByName(std::string name)
{
    unsigned index = mpTransportReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);

    return mReservedCellConcentration[index];
}

double AbstractTransportOdeSystem::RetrieveChangeCellConcentrationByIndex(unsigned index)
{
    return mChangeCellConc[index];
}

double AbstractTransportOdeSystem::RetrieveChangeCellConcentrationByName(std::string name)
{
    unsigned index = mpTransportReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);

    return mChangeCellConc[index];
}

template<>
void TransportOdeSystemInformation<AbstractTransportOdeSystem>::Initialise()
{
    // Initialise the bulk state varibles
    for (unsigned i=0; i<mp_reaction_system->GetBulkChemistry()->GetNumberChemicals(); ++i)
    {
        this->mVariableNames.push_back(mp_reaction_system->GetBulkChemistry()->GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system->GetBulkChemistry()->GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    
    // Initialise the cell state varibles; appended onto the end of the bulk state vectors
    for (unsigned i=0; i<mp_reaction_system->GetCellChemistry()->GetNumberChemicals(); ++i)
    {
        this->mVariableNames.push_back(mp_reaction_system->GetCellChemistry()->GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system->GetCellChemistry()->GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    this->mInitialised = true;
}

#endif 