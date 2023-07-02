#include "AbstractReactionSystem.hpp"

AbstractReactionSystem::AbstractReactionSystem(AbstractChemistry* systemChemistry, 
                            std::vector<AbstractReaction*> reactionVector)
    : mpSystemChemistry(systemChemistry),
      mpReactionVector(reactionVector),
      mNumReactions(reactionVector.size())
{
}

AbstractReactionSystem::AbstractReactionSystem(const AbstractReactionSystem& existingReactionSystem)
{
    mpReactionVector = existingReactionSystem.mpReactionVector;
    mNumReactions = existingReactionSystem.mNumReactions;
    mpSystemChemistry = existingReactionSystem.mpSystemChemistry;
}

void AbstractReactionSystem::ReactSystem(const std::vector<double>& currentChemistryConc, std::vector<double>& changeChemistryConc)
{
    unsigned number_of_species = currentChemistryConc.size();

    // update the reaction system if any variables depend on the current system concentrations
    UpdateReactionSystem(currentChemistryConc);
    
    std::vector<double> deltaConcentration(number_of_species, 0.0);
    for (std::vector<AbstractReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        // iterate through the reactions in the system, modify a temporary change of concentration
        // update the system change in concentration
        
        AbstractReaction* p_system_reaction = dynamic_cast<AbstractReaction*>(*reaction_iter);
        
        p_system_reaction->React(mpSystemChemistry, currentChemistryConc, deltaConcentration);
       
        // update the change in concentration
        
        for (unsigned i=0; i<number_of_species; ++i)
        {
            changeChemistryConc[i] += deltaConcentration[i];

            deltaConcentration[i] = 0;
        }
    }
}

void AbstractReactionSystem::UpdateReactionSystem(const std::vector<double>& currentChemistryConc)
{
    // virtual
    return;
}

std::vector<AbstractReaction*> AbstractReactionSystem::GetReactionVector()
{
    return mpReactionVector;
}

AbstractReaction* AbstractReactionSystem::GetReactionByIndex(unsigned index)
{
    if (index < mNumReactions)
    {
        return mpReactionVector[index];
    }
    else
    {
        return new AbstractReaction();
    }
}

void AbstractReactionSystem::SetReactionVector(std::vector<AbstractReaction*> reactionVector)
{
    mpReactionVector = reactionVector;
}

void AbstractReactionSystem::SetReactionByIndex(AbstractReaction* reaction, unsigned index)
{
    if (index < mNumReactions)
    {
        mpReactionVector[index] = reaction;
    }
}

AbstractChemistry* AbstractReactionSystem::GetSystemChemistry()
{
    return mpSystemChemistry;
}

void AbstractReactionSystem::SetSystemChemistry(AbstractChemistry* systemChemistry)
{
    mpSystemChemistry = systemChemistry;
}

unsigned AbstractReactionSystem::GetNumReactions()
{
    return mNumReactions;
}

void AbstractReactionSystem::SetNumReactions(unsigned numReactions)
{
    mNumReactions = numReactions;
}