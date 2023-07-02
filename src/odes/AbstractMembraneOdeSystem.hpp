#ifndef ABSTRACTMEMBRANEODESYSTEM_HPP_
#define ABSTRACTMEMBRANEODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractMembraneReactionSystem.hpp"
#include "MembraneOdeSystemInformation.hpp"
#include "StateVariableRegister.hpp"

/**
 * Abstract ODE system for controllign the change in chemical concentrations 
 * (state variables) at the membrane of a cell. The class handles the different 
 * chemical concentrations at points eitherside of the infinitesmimal membrane, 
 * in this base class, parsing the differnet concentrations for the membrane 
 * reaction system. The class recombines the different concentrations to pass to 
 * handler methods and also performs concentration errror checking through set 
 * thresholds.
 */
class AbstractMembraneOdeSystem: public AbstractOdeSystem
{
private:

    AbstractMembraneReactionSystem* mpMembraneReactionSystem;

    StateVariableRegister* mPdeStateVariableRegister;

    std::vector<double> mReservedCellConcentration;

    std::vector<double> mChangeCellConc;

    unsigned mNumSpecies;

    unsigned mNumReactions;

    double mDeltaError;

    bool mIsCheckConcentration;

public:

    AbstractMembraneOdeSystem(AbstractMembraneReactionSystem* pReactionSystem = new AbstractMembraneReactionSystem());

    virtual ~AbstractMembraneOdeSystem()
    {
    };
    
    AbstractMembraneOdeSystem(const AbstractMembraneOdeSystem&);

    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    virtual void UpdateReservedCellConcentration(std::vector<double>);

    virtual void UpdateChangeCellConcentration(std::vector<double>);

    unsigned GetNumSpecies();

    unsigned GetNumReactions();

    AbstractMembraneReactionSystem* GetReactionSystem();

    void SetNumSpecies(unsigned);

    void SetNumReactions(unsigned);

    void SetReactionSystem(AbstractMembraneReactionSystem* );
    
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

AbstractMembraneOdeSystem::AbstractMembraneOdeSystem(AbstractMembraneReactionSystem* pReactionSystem)
    : AbstractOdeSystem(pReactionSystem->GetNumBulkStates() + pReactionSystem->GetNumCellStates()),
      mpMembraneReactionSystem(pReactionSystem)
{
    SetReactionSystem(mpMembraneReactionSystem);
    SetNumSpecies(mpMembraneReactionSystem->GetNumBulkStates()+mpMembraneReactionSystem->GetNumCellStates());
    SetNumReactions(mpMembraneReactionSystem->GetNumReactions());

    // Initialsie cell concentration vectors
    std::vector<double> cellConcentration(mpMembraneReactionSystem->GetNumCellStates(),1.0);
    SetReservedCellConcentration(cellConcentration);
    std::vector<double> cellConcentrationChange(mpMembraneReactionSystem->GetNumCellStates(),0.0);
    SetChangeCellConcentration(cellConcentrationChange);

    mpSystemInfo.reset(new MembraneOdeSystemInformation<AbstractMembraneOdeSystem>(mpMembraneReactionSystem));

    SetStateVariables(GetInitialConditions()); // set up in the solver call; size of NumBulkStates + NumCellStates

    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractMembraneOdeSystem::AbstractMembraneOdeSystem(const AbstractMembraneOdeSystem& existingSystem)
    : AbstractOdeSystem(existingSystem.mNumSpecies)
{
    mpMembraneReactionSystem = existingSystem.mpMembraneReactionSystem;
    mNumSpecies = existingSystem.mNumSpecies;
    mNumReactions = existingSystem.mNumReactions;
    mDeltaError = existingSystem.mDeltaError;
    mIsCheckConcentration = existingSystem.mIsCheckConcentration;
}

void AbstractMembraneOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    unsigned num_bulk_states = mpMembraneReactionSystem->GetNumBulkStates();
    unsigned num_cell_states = mpMembraneReactionSystem->GetNumCellStates();

    std::vector<double> bulkConcentrations(num_bulk_states, 0.0);
    std::vector<double> changeBulkConcentrations(num_bulk_states, 0.0);

    std::vector<double> cellConcentrations(num_cell_states, 0.0);
    std::vector<double> changeCellConcentrations(num_cell_states, 0.0);

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
 
    mpMembraneReactionSystem->ReactSystem(bulkConcentrations, cellConcentrations, changeBulkConcentrations, changeCellConcentrations);

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

void AbstractMembraneOdeSystem::UpdateReservedCellConcentration(std::vector<double> reservedConcentration)
{
    SetReservedCellConcentration(reservedConcentration);
}

void AbstractMembraneOdeSystem::UpdateChangeCellConcentration(std::vector<double> changeCellConc)
{
    mChangeCellConc = changeCellConc;
}

unsigned AbstractMembraneOdeSystem::GetNumSpecies()
{
    return mNumSpecies;
}

unsigned AbstractMembraneOdeSystem::GetNumReactions()
{
    return mNumReactions;
}

AbstractMembraneReactionSystem* AbstractMembraneOdeSystem::GetReactionSystem()
{
    return mpMembraneReactionSystem;
}

void AbstractMembraneOdeSystem::SetNumSpecies(unsigned numSpecies)
{
    mNumSpecies = numSpecies;
}

void AbstractMembraneOdeSystem::SetNumReactions(unsigned numReactions)
{
    mNumReactions = numReactions;
}

void AbstractMembraneOdeSystem::SetReactionSystem(AbstractMembraneReactionSystem* p_reactionSystem)
{
    mpMembraneReactionSystem = p_reactionSystem;
}

void AbstractMembraneOdeSystem::SetReservedCellConcentration(std::vector<double> cellConcentration)
{
    mReservedCellConcentration = cellConcentration;
}
    
void AbstractMembraneOdeSystem::SetChangeCellConcentration(std::vector<double> changeCellConcentration)
{
    mChangeCellConc = changeCellConcentration;
}

void AbstractMembraneOdeSystem::UpdateReservedCellConcentrationByIndex(unsigned index, double concentration)
{
    mReservedCellConcentration[index] = concentration;
}

void AbstractMembraneOdeSystem::UpdateReservedCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpMembraneReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);
    mReservedCellConcentration[index] = concentration;
}

void AbstractMembraneOdeSystem::UpdateChangeCellConcentrationByIndex(unsigned index, double concentration)
{
    mChangeCellConc[index] = concentration;
}

void AbstractMembraneOdeSystem::UpdateChangeCellConcentrationByName(std::string name, double concentration)
{
    unsigned index = mpMembraneReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);
    mChangeCellConc[index] = concentration;
}

void AbstractMembraneOdeSystem::SetDeltaError(double deltaError)
{
    mDeltaError = deltaError;
}

double AbstractMembraneOdeSystem::GetDeltaError()
{
    return mDeltaError;
}

void AbstractMembraneOdeSystem::SetPdeStateRegister(StateVariableRegister* pRegister)
{
    mPdeStateVariableRegister = pRegister;
}

void AbstractMembraneOdeSystem::CheckConcentration(const std::vector<double>& rY)
{
    // If chemical concentration gets too low then NaNs can occur, concentration must be +ve
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

void AbstractMembraneOdeSystem::SetIsCheckConcentration(bool isCheckConcentration)
{
    mIsCheckConcentration = isCheckConcentration;
}

bool AbstractMembraneOdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

std::vector<double> AbstractMembraneOdeSystem::GetReservedCellConcentration()
{
    return mReservedCellConcentration;
}

std::vector<double> AbstractMembraneOdeSystem::GetChangeCellConcentration()
{
    return mChangeCellConc;
}

StateVariableRegister* AbstractMembraneOdeSystem::GetPdeStateVariableRegister()
{
    return mPdeStateVariableRegister;
}

double AbstractMembraneOdeSystem::RetrieveReservedCellConcentrationByIndex(unsigned index)
{
    return mReservedCellConcentration[index];
}

double AbstractMembraneOdeSystem::RetrieveReservedCellConcentrationByName(std::string name)
{
    unsigned index = mpMembraneReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);

    return mReservedCellConcentration[index];
}

double AbstractMembraneOdeSystem::RetrieveChangeCellConcentrationByIndex(unsigned index)
{
    return mChangeCellConc[index];
}

double AbstractMembraneOdeSystem::RetrieveChangeCellConcentrationByName(std::string name)
{
    unsigned index = mpMembraneReactionSystem->GetCellChemistry()->GetChemicalIndexByName(name);
    return mChangeCellConc[index];
}

template<>
void MembraneOdeSystemInformation<AbstractMembraneOdeSystem>::Initialise()
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