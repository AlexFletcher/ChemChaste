#ifndef ABSTRACTMEMBRANEREACTIONSYSTEM_HPP_
#define ABSTRACTMEMBRANEREACTIONSYSTEM_HPP_

#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include "AbstractMembraneReaction.hpp"
#include "AbstractChemistry.hpp"

/**
 * Class to hold a set of membrane reactions and handle the parsing and 
 * reforming of the different chemical concentration (state variable) vectors. 
 * These vectors being the concentrations external and internal to the 
 * infinitesimal thickness membrane.
 */
class AbstractMembraneReactionSystem
{
protected:

    AbstractChemistry* mpBulkChemistry;

    AbstractChemistry* mpCellChemistry;

    std::vector<AbstractMembraneReaction*> mpReactionVector;

    unsigned mNumReactions;

    unsigned mNumBulkStates=0;

    unsigned mNumCellStates=0;

public:
    AbstractMembraneReactionSystem( AbstractChemistry* bulkChemistry = new AbstractChemistry(), 
                                    AbstractChemistry* cellChemistry = new AbstractChemistry(), 
                                    std::vector<AbstractMembraneReaction*> reactionVector = std::vector<AbstractMembraneReaction*>());

    virtual ~AbstractMembraneReactionSystem()
    {
    };

    AbstractMembraneReactionSystem(const AbstractMembraneReactionSystem&);

    // virtual methods

    virtual void ReactSystem(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReactionSystem(const std::vector<double>&, const std::vector<double>&);

    // set methods

    void SetReactionVector(std::vector<AbstractMembraneReaction*>);

    void SetReactionByIndex(AbstractMembraneReaction*, unsigned);

    void SetBulkChemistry(AbstractChemistry*);

    void SetCellChemistry(AbstractChemistry*);

    void SetNumReactions(unsigned);


    // get methods 

    std::vector<AbstractMembraneReaction*> GetReactionVector();

    AbstractMembraneReaction* GetReactionByIndex(unsigned);

    AbstractChemistry* GetBulkChemistry();

    AbstractChemistry* GetCellChemistry();

    unsigned GetNumReactions();

    unsigned GetNumBulkStates();

    unsigned GetNumCellStates();

};

AbstractMembraneReactionSystem::AbstractMembraneReactionSystem( AbstractChemistry* bulkChemistry,
                                                                AbstractChemistry* cellChemistry, 
                                                                std::vector<AbstractMembraneReaction*> reactionVector)
    :   mpBulkChemistry(bulkChemistry),
        mpCellChemistry(cellChemistry),
        mpReactionVector(reactionVector),
        mNumReactions(reactionVector.size()),
        mNumBulkStates(bulkChemistry->GetNumberChemicals()),
        mNumCellStates(cellChemistry->GetNumberChemicals())
{
}

AbstractMembraneReactionSystem::AbstractMembraneReactionSystem(const AbstractMembraneReactionSystem& existingReactionSystem)
{
    mpBulkChemistry = existingReactionSystem.mpBulkChemistry;
    mpCellChemistry = existingReactionSystem.mpCellChemistry;
    mpReactionVector = existingReactionSystem.mpReactionVector;
    mNumReactions = existingReactionSystem.mNumReactions;
    mNumCellStates = existingReactionSystem.mNumCellStates;
    mNumBulkStates = existingReactionSystem.mNumBulkStates;
}

void AbstractMembraneReactionSystem::ReactSystem(const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc,std::vector<double>& changeCellConc)
{
    // update the reaction system if any variables depend on the current system concentrations
    UpdateReactionSystem(currentBulkConcentration,currentCellConcentration);

    std::vector<double> deltaBulkConcentration(mNumBulkStates, 0.0);
    std::vector<double> deltaCellConcentration(mNumCellStates, 0.0);


    for (std::vector<AbstractMembraneReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        // iterate through the reactions in the system, modify a temporary change of concentration
        // update the system change in concentration
        deltaBulkConcentration.resize(mNumBulkStates, 0.0);
        deltaCellConcentration.resize(mNumCellStates, 0.0);

        AbstractMembraneReaction *p_system_reaction = static_cast<AbstractMembraneReaction*>(*reaction_iter);

        p_system_reaction->React(mpBulkChemistry, mpCellChemistry, currentBulkConcentration, currentCellConcentration, deltaBulkConcentration, deltaCellConcentration);

        // update the change in concentrations
        
        for (unsigned i=0; i<mNumBulkStates; ++i)
        {
            changeBulkConc.at(i) += deltaBulkConcentration.at(i);
        }

        for (unsigned i=0; i<mNumCellStates; ++i)
        {
            changeCellConc.at(i) += deltaCellConcentration.at(i);
        }    
    }
}

void AbstractMembraneReactionSystem::UpdateReactionSystem(const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // virtual
    return;
}

void AbstractMembraneReactionSystem::SetReactionVector(std::vector<AbstractMembraneReaction*> reactionVector)
{
    mpReactionVector = reactionVector;
}

void AbstractMembraneReactionSystem::SetReactionByIndex(AbstractMembraneReaction* reaction, unsigned index)
{
    if (index < mNumReactions)
    {
        mpReactionVector[index] = reaction;
    }
}

void AbstractMembraneReactionSystem::SetBulkChemistry(AbstractChemistry* bulkChemistry)
{
    mpBulkChemistry = bulkChemistry;
    mNumBulkStates = bulkChemistry->GetNumberChemicals();
}

void AbstractMembraneReactionSystem::SetCellChemistry(AbstractChemistry* cellChemistry)
{
    mpCellChemistry = cellChemistry;
    mNumCellStates = cellChemistry->GetNumberChemicals();
}

void AbstractMembraneReactionSystem::SetNumReactions(unsigned numReactions)
{
    mNumReactions = numReactions;
}

std::vector<AbstractMembraneReaction*> AbstractMembraneReactionSystem::GetReactionVector()
{
    return mpReactionVector;
}

AbstractMembraneReaction* AbstractMembraneReactionSystem::GetReactionByIndex(unsigned index)
{
    if (index < mNumReactions)
    {
        return mpReactionVector[index];
    }else{
        return new AbstractMembraneReaction();
    }
}

AbstractChemistry* AbstractMembraneReactionSystem::GetBulkChemistry()
{
    return mpBulkChemistry;
}

AbstractChemistry* AbstractMembraneReactionSystem::GetCellChemistry()
{
    return mpCellChemistry;
}

unsigned AbstractMembraneReactionSystem::GetNumReactions()
{
    return mNumReactions;
}

unsigned AbstractMembraneReactionSystem::GetNumBulkStates()
{
    return mNumBulkStates;
}

unsigned AbstractMembraneReactionSystem::GetNumCellStates()
{
    return mNumCellStates;
}

#endif