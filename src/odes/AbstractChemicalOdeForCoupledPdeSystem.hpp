#ifndef AbstractChemicalOdeForCoupledPdeSystem_HPP_
#define AbstractChemicalOdeForCoupledPdeSystem_HPP_

#include "AbstractOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractReactionSystem.hpp"
#include "ChemicalOdeSystemInformation.hpp"

/**
 * Chemical ODE system to be coupled to a PDE system, contains the methods to 
 * evaluate a chemical reaction system. Also checks the output of the chemical 
 * reaction system to ensure non-negative concentrations (states) are removed.
 * Not used in the reaction-diffusion system as the chemical ODE is not required 
 * to be solved.
 */
class AbstractChemicalOdeForCoupledPdeSystem: public AbstractOdeSystemForCoupledPdeSystem
{
private:

    // reaction system object which contains the data structures and functions regarding the reactions
    AbstractReactionSystem* mpReactionSystem;

    // number of chemical species (state varibales) in the pde system 
    unsigned mNumSpecies;

    // numebr of the reacitons in the reaction system which need to be iterated through to fully react the system
    unsigned mNumReactions;
    
    // error tolerence for checking whether the concentrations are >0, i.e =0 if <tolerence
    // default 1e-6
    double mDeltaError; 
    
    // whether to check for negative concentrations
    // default true
    bool mIsCheckConcentration;

public:

    AbstractChemicalOdeForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem = new AbstractReactionSystem());

    virtual ~AbstractChemicalOdeForCoupledPdeSystem()
    {
    };
    
    AbstractChemicalOdeForCoupledPdeSystem(const AbstractChemicalOdeForCoupledPdeSystem&);

    // virtual methods

    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    // concrete methods

    void CheckConcentration(const std::vector<double>&);

    // set methods

    void SetNumSpecies(unsigned);

    void SetNumReactions(unsigned);

    void SetReactionSystem(AbstractReactionSystem* );

    void SetDeltaError(double);

    void SetIsCheckConcentration(bool);

    // get methods

    unsigned GetNumSpecies();

    unsigned GetNumReactions();

    AbstractReactionSystem* GetReactionSystem();

    double GetDeltaError();

    bool GetIsCheckConcentration();
};

AbstractChemicalOdeForCoupledPdeSystem::AbstractChemicalOdeForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem)
    : AbstractOdeSystemForCoupledPdeSystem(pReactionSystem->GetSystemChemistry()->GetNumberChemicals(), pReactionSystem->GetSystemChemistry()->GetNumberChemicals())            
{
    // Set the reaction insformation
    SetReactionSystem(pReactionSystem);
    SetNumSpecies(pReactionSystem->GetSystemChemistry()->GetNumberChemicals());
    SetNumReactions(pReactionSystem->GetNumReactions());

    // Fulfil the Chaste Ode system requirement of initialising the ode system information
    mpSystemInfo.reset(new ChemicalOdeSystemInformation<AbstractChemicalOdeForCoupledPdeSystem>(pReactionSystem));

    SetStateVariables(GetInitialConditions());

    // Setup the data needed for checking the concentration checking procedures
    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractChemicalOdeForCoupledPdeSystem::AbstractChemicalOdeForCoupledPdeSystem(const AbstractChemicalOdeForCoupledPdeSystem& existingSystem)
    : AbstractOdeSystemForCoupledPdeSystem(existingSystem.mNumSpecies,existingSystem.mNumSpecies)
{
    mpReactionSystem = existingSystem.mpReactionSystem;

    mNumSpecies = existingSystem.mNumSpecies;

    mNumReactions = existingSystem.mNumReactions;
}

void AbstractChemicalOdeForCoupledPdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double>& sol_pde  = this->rGetPdeSolution();

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        rDY[i] = 0.0;
    }
    
    // perform the reactions in the system updating rDY with the resultant change in concentrations
    mpReactionSystem->ReactSystem(rY, rDY);

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if (rDY[i]<mDeltaError && rDY[i]>mDeltaError)
        {
            rDY[i] = 0;
        }
    }

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        const_cast<double&>(rY[i]) = rDY[i]; // this can't be correct?
    }

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        rDY[i] = 0.0;
    }
}  

void AbstractChemicalOdeForCoupledPdeSystem::CheckConcentration(const std::vector<double>& rY)
{
    // if chemical concentration gets too low then nan can occur, concentration must be +ve
    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if (rY[i]<mDeltaError)
        {
            const_cast<double&>(rY[i]) = 0;
        }
    }
}

void AbstractChemicalOdeForCoupledPdeSystem::SetNumSpecies(unsigned numSpecies)
{
    mNumSpecies = numSpecies;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetNumReactions(unsigned numReactions)
{
    mNumReactions = numReactions;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetReactionSystem(AbstractReactionSystem* p_reactionSystem)
{
    mpReactionSystem = p_reactionSystem;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetDeltaError(double deltaError)
{
    mDeltaError = deltaError;
}

void AbstractChemicalOdeForCoupledPdeSystem::SetIsCheckConcentration(bool isCheckConcentration)
{
    mIsCheckConcentration = isCheckConcentration;
}

unsigned AbstractChemicalOdeForCoupledPdeSystem::GetNumSpecies()
{
    return mNumSpecies;
}

unsigned AbstractChemicalOdeForCoupledPdeSystem::GetNumReactions()
{
    return mNumReactions;
}

AbstractReactionSystem* AbstractChemicalOdeForCoupledPdeSystem::GetReactionSystem()
{
    return mpReactionSystem;
}

bool AbstractChemicalOdeForCoupledPdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

double AbstractChemicalOdeForCoupledPdeSystem::GetDeltaError()
{
    return mDeltaError;
}

template<>
void ChemicalOdeSystemInformation<AbstractChemicalOdeForCoupledPdeSystem>::Initialise()
{
    // initilaise the ode information required by Chaste for all ode systems
    // derive the states from the chemical vectors of the system chemistry
    for (unsigned i=0; i<mp_reaction_system->GetSystemChemistry()->GetNumberChemicals(); ++i)
    {
        this->mVariableNames.push_back(mp_reaction_system->GetSystemChemistry()->GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system->GetSystemChemistry()->GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    this->mInitialised = true;
}

#endif 