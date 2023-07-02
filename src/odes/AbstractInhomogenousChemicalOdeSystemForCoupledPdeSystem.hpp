#ifndef AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem_HPP_
#define AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem_HPP_

#include "AbstractInhomogenousOdeSystemForCoupledPdeSystem.hpp"
#include "AbstractReactionSystem.hpp"
#include "ChemicalOdeSystemInformation.hpp"
#include "StateVariableRegister.hpp"

/**
 * \todo Document class.
 */
class AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem: public AbstractInhomogenousOdeSystemForCoupledPdeSystem
{
private:

    // pointer to the reaction system which determines the chemical dynamics of the cell
    AbstractReactionSystem* mpReactionSystem;

    unsigned mNumSpecies;

    unsigned mNumReactions;

    double mDeltaError;

    bool mIsCheckConcentration;

public:

    AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem = new AbstractReactionSystem());

    virtual ~AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem()
    {
    };
    
    AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(const AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem&);

    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    unsigned GetNumSpecies();

    unsigned GetNumReactions();

    AbstractReactionSystem* GetReactionSystem();

    void SetNumSpecies(unsigned);

    void SetNumReactions(unsigned);

    void SetReactionSystem(AbstractReactionSystem* );

    void SetDeltaError(double);

    double GetDeltaError();

    void CheckConcentration(const std::vector<double>&);

    void SetIsCheckConcentration(bool);

    bool GetIsCheckConcentration();
};

AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(AbstractReactionSystem* pReactionSystem)
    : AbstractInhomogenousOdeSystemForCoupledPdeSystem(pReactionSystem->GetSystemChemistry()->GetNumberChemicals(),pReactionSystem->GetSystemChemistry()->GetNumberChemicals())
{
    SetReactionSystem(pReactionSystem);
    SetNumSpecies(pReactionSystem->GetSystemChemistry()->GetNumberChemicals());
    SetNumReactions(pReactionSystem->GetNumReactions());

    mpSystemInfo.reset(new ChemicalOdeSystemInformation<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem>(pReactionSystem));
    SetStateVariables(GetInitialConditions());

    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(const AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem& existingSystem)
    : AbstractInhomogenousOdeSystemForCoupledPdeSystem(existingSystem.mNumSpecies, existingSystem.mNumSpecies)
{
    mpReactionSystem = existingSystem.mpReactionSystem;
    mNumSpecies = existingSystem.mNumSpecies;
    mNumReactions = existingSystem.mNumReactions;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    if (mIsCheckConcentration)
    {
        CheckConcentration(rY);
    }

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        rDY[i] = 0.0;
    }

    mpReactionSystem->ReactSystem(rY, rDY);

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if (rDY[i] < mDeltaError && rDY[i] > mDeltaError)
        {
            rDY[i] = 0;
        }
    }

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        const_cast<double&>(rY[i]) = rDY[i];
    }

    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        rDY[i] = 0.0;
    }
}  

unsigned AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetNumSpecies()
{
    return mNumSpecies;
}

unsigned AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetNumReactions()
{
    return mNumReactions;
}

AbstractReactionSystem* AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetReactionSystem()
{
    return mpReactionSystem;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetNumSpecies(unsigned numSpecies)
{
    mNumSpecies = numSpecies;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetNumReactions(unsigned numReactions)
{
    mNumReactions = numReactions;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetReactionSystem(AbstractReactionSystem* p_reactionSystem)
{
    mpReactionSystem = p_reactionSystem;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetDeltaError(double deltaError)
{
    mDeltaError = deltaError;
}

double AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetDeltaError()
{
    return mDeltaError;
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::CheckConcentration(const std::vector<double>& rY)
{
    // if chemical concentration gets too low then nan can occur, concentration must be +ve
    for (unsigned i=0; i<mNumSpecies; ++i)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if (rY[i] < mDeltaError)
        {
            const_cast<double&>(rY[i]) = 0;
        }
    }
}

void AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::SetIsCheckConcentration(bool isCheckConcentration)
{
    mIsCheckConcentration = isCheckConcentration;
}

bool AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

template<>
void ChemicalOdeSystemInformation<AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem>::Initialise()
{
    for (unsigned i=0; i<mp_reaction_system->GetSystemChemistry()->GetNumberChemicals(); ++i)
    {
        this->mVariableNames.push_back(mp_reaction_system->GetSystemChemistry()->GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system->GetSystemChemistry()->GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    this->mInitialised = true;
}

#endif 