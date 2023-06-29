#ifndef ABSTRACTCHEMICALODESYSTEM_HPP_
#define ABSTRACTCHEMICALODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "AbstractReactionSystem.hpp"
#include "ChemicalOdeSystemInformation.hpp"

/**
 * ODE system for performing chemical reactions of templated chemical kinetic 
 * rate laws.
 */
class AbstractChemicalOdeSystem: public AbstractOdeSystem
{
private:

    /** The chemical reaction system to be modelled. */
    AbstractReactionSystem* mpReactionSystem;

    /** Number of chemical species within the chemical system. */
    unsigned mNumberOfSpecies;

    /** Number of chemical reactions within the reaction system to be solved. */
    unsigned mNumberOfReactions;

    /** Error threshold. Defaults to 1e-6. */
    double mDeltaError;

    /**
     * Whether to check the chemical concentrations (ODE states) for negative 
     * concentration values. Defaults to true.
     */
    bool mIsCheckConcentration;

public:

    /**
     * Default constructor.
     * 
     * @param pReactionSystem \todo document param
     */
    AbstractChemicalOdeSystem(AbstractReactionSystem* pReactionSystem = new AbstractReactionSystem());

    /**
     * Copy constructor.
     * 
     * @param rOtherReactionSystem \todo document param
     */
    AbstractChemicalOdeSystem(const AbstractChemicalOdeSystem& rOtherReactionSystem);

    /**
     * Destructor.
     */
    virtual ~AbstractChemicalOdeSystem();

    /**
     * Compute the RHS of the ODE system.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives
     */
    virtual void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    /**
     * \todo Document method.
     * 
     * @param rY value of the solution vector used to evaluate the RHS.
     */
    void CheckConcentration(const std::vector<double>& rY);

    /**
     * Set mpReactionSystem.
     * 
     * @param pReactionSystem The chemical reaction system to be modelled. 
     */
    void SetReactionSystem(AbstractReactionSystem* pReactionSystem);

    /**
     * Set mNumberOfSpecies.
     * 
     * @param numberOfSpecies Number of chemical species within the chemical system. 
     */
    void SetNumberOfSpecies(unsigned numberOfSpecies);

    /**
     * Set mNumberOfReactions.
     * 
     * @param numberOfSpecies Number of chemical reactions within the reaction system to be solved. 
     */
    void SetNumberOfReactions(unsigned numberOfReactions);

    /**
     * Set mDeltaError.
     * 
     * @param deltaError Error threshold. 
     */
    void SetDeltaError(double deltaError);

    /**
     * Set mIsCheckConcentration.
     * 
     * @param isCheckConcentration Whether to check the chemical concentrations 
     *                             (ODE states) for negative concentration values.
     */
    void SetIsCheckConcentration(bool isCheckConcentration);

    /**
     * @return mpReactionSystem.
     */
    AbstractReactionSystem* GetReactionSystem();

    /**
     * @return mNumberOfSpecies.
     */
    unsigned GetNumberOfSpecies();

    /**
     * @return mNumberOfReactions.
     */
    unsigned GetNumberOfReactions();

    /**
     * @return mDeltaError.
     */
    double GetDeltaError();

    /**
     * @return mIsCheckConcentration.
     */
    bool GetIsCheckConcentration();
};

AbstractChemicalOdeSystem::AbstractChemicalOdeSystem(AbstractReactionSystem* pReactionSystem)
    : AbstractOdeSystem(pReactionSystem->GetSystemChemistry()->GetNumberChemicals()) 
{
    // Set up the chemical reaction sytem
    SetReactionSystem(pReactionSystem);
    SetNumberOfSpecies(pReactionSystem->GetSystemChemistry()->GetNumberChemicals());
    SetNumberOfReactions(pReactionSystem->GetNumberOfReactions());

    // Fulfil the Chaste ODE information setup by initialising
    mpSystemInfo.reset(new ChemicalOdeSystemInformation<AbstractChemicalOdeSystem>(pReactionSystem));
    SetStateVariables(GetInitialConditions());

    // Set the chemical concentration error checking threhsold
    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

AbstractChemicalOdeSystem::AbstractChemicalOdeSystem(const AbstractChemicalOdeSystem& rOtherReactionSystem)
    : AbstractOdeSystem(rOtherReactionSystem.mNumberOfSpecies)
{
    mpReactionSystem = rOtherReactionSystem.mpReactionSystem;
    mNumberOfSpecies = rOtherReactionSystem.mNumberOfSpecies;
    mNumberOfReactions = rOtherReactionSystem.mNumberOfReactions;
    mDeltaError = rOtherReactionSystem.mDeltaError;
    mIsCheckConcentration = rOtherReactionSystem.mIsCheckConcentration;
}

AbstractChemicalOdeSystem::~AbstractChemicalOdeSystem()
{
}

void AbstractChemicalOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    if (mIsCheckConcentration)
    {
        CheckConcentration(rY);
    }

    // Reset rDY
    for (unsigned i = 0; i < mNumberOfSpecies; i++)
    {
        rDY[i] = 0.0;
    }

    // Perform the reaction, updating rDY with the resultant chemical change due to the reaction
    mpReactionSystem->ReactSystem(rY, rDY);
}  

void AbstractChemicalOdeSystem::CheckConcentration(const std::vector<double>& rY)
{
    // If chemical concentration gets too low then nan can occur, concentration must be +ve
    for (unsigned i = 0; i < mNumberOfSpecies; i++)
    {
        // Due to the discrete nature, rY can evaluate to < 0; ensure rY >= 0
        if (rY[i] < mDeltaError)
        {
            const_cast<double&>(rY[i]) = 0;
        }
    }
}

void AbstractChemicalOdeSystem::SetReactionSystem(AbstractReactionSystem* pReactionSystem)
{
    mpReactionSystem = pReactionSystem;
}

void AbstractChemicalOdeSystem::SetNumberOfSpecies(unsigned numberOfSpecies)
{
    mNumberOfSpecies = numberOfSpecies;
}

void AbstractChemicalOdeSystem::SetNumberOfReactions(unsigned numberOfReactions)
{
    mNumberOfReactions = numberOfReactions;
}

void AbstractChemicalOdeSystem::SetDeltaError(double deltaError)
{
    mDeltaError = deltaError;
}

void AbstractChemicalOdeSystem::SetIsCheckConcentration(bool isCheckConcentration)
{
    mIsCheckConcentration = isCheckConcentration;
}

AbstractReactionSystem* AbstractChemicalOdeSystem::GetReactionSystem()
{
    return mpReactionSystem;
}

unsigned AbstractChemicalOdeSystem::GetNumberOfSpecies()
{
    return mNumberOfSpecies;
}

unsigned AbstractChemicalOdeSystem::GetNumberOfReactions()
{
    return mNumberOfReactions;
}

bool AbstractChemicalOdeSystem::GetIsCheckConcentration()
{
    return mIsCheckConcentration;
}

double AbstractChemicalOdeSystem::GetDeltaError()
{
    return mDeltaError;
}

// System information template
template<>
void ChemicalOdeSystemInformation<AbstractChemicalOdeSystem>::Initialise()
{
    for (unsigned i = 0; i < mp_reaction_system->GetSystemChemistry()->GetNumberChemicals(); i++)
    {
        this->mVariableNames.push_back(mp_reaction_system->GetSystemChemistry()->GetChemicalNamesByIndex(i));
        this->mVariableUnits.push_back(mp_reaction_system->GetSystemChemistry()->GetChemicalDimensionsByIndex(i));
        this->mInitialConditions.push_back(1.0); // will be overridden
    }
    this->mInitialised = true;
}

#endif /* ABSTRACTCHEMICALODESYSTEM_HPP_ */