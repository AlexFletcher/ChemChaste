#ifndef ABSTRACTTRANSPORTREACTIONSYSTEM_HPP_
#define ABSTRACTTRANSPORTREACTIONSYSTEM_HPP_

#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include "AbstractTransportReaction.hpp"
#include "AbstractChemistry.hpp"

/**
 * Class to contain and handle a system of transport reactions. Transport 
 * reactions are defined across a membrane where either side have their own 
 * concentration vectors.
 */
class AbstractTransportReactionSystem
{
protected:

    AbstractChemistry* mpBulkChemistry;

    AbstractChemistry* mpCellChemistry;

    std::vector<AbstractTransportReaction*> mpReactionVector;

    unsigned mNumReactions;

    unsigned mNumBulkStates=0;

    unsigned mNumCellStates=0;

public:

    AbstractTransportReactionSystem(AbstractChemistry* bulkChemistry = new AbstractChemistry(), 
                                    AbstractChemistry* cellChemistry = new AbstractChemistry(), 
                                    std::vector<AbstractTransportReaction*> reactionVector = std::vector<AbstractTransportReaction*>());

    virtual ~AbstractTransportReactionSystem()
    {
    };

    AbstractTransportReactionSystem(const AbstractTransportReactionSystem&);

    // virtual methods

    virtual void ReactSystem(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReactionSystem(const std::vector<double>&, const std::vector<double>&);

    // set methods

    void SetReactionVector(std::vector<AbstractTransportReaction*>);

    void SetReactionByIndex(AbstractTransportReaction*, unsigned);

    void SetBulkChemistry(AbstractChemistry*);

    void SetCellChemistry(AbstractChemistry*);

    void SetNumReactions(unsigned);


    // get methods 

    std::vector<AbstractTransportReaction*> GetReactionVector();

    AbstractTransportReaction* GetReactionByIndex(unsigned);

    AbstractChemistry* GetBulkChemistry();

    AbstractChemistry* GetCellChemistry();

    unsigned GetNumReactions();

    unsigned GetNumBulkStates();

    unsigned GetNumCellStates();
};

AbstractTransportReactionSystem::AbstractTransportReactionSystem(   AbstractChemistry* bulkChemistry,
                                                                    AbstractChemistry* cellChemistry, 
                                                                    std::vector<AbstractTransportReaction*> reactionVector)
    :   mpBulkChemistry(bulkChemistry),
        mpCellChemistry(cellChemistry),
        mpReactionVector(reactionVector),
        mNumReactions(reactionVector.size()),
        mNumBulkStates(bulkChemistry->GetNumberChemicals()),
        mNumCellStates(cellChemistry
        GetNumberChemicals())
{
}

AbstractTransportReactionSystem::AbstractTransportReactionSystem(const AbstractTransportReactionSystem& existingReactionSystem)
{
    mpBulkChemistry = existingReactionSystem.mpBulkChemistry;
    mpCellChemistry = existingReactionSystem.mpCellChemistry;
    mpReactionVector = existingReactionSystem.mpReactionVector;
    mNumReactions = existingReactionSystem.mNumReactions;
    mNumCellStates = existingReactionSystem.mNumCellStates;
    mNumBulkStates = existingReactionSystem.mNumBulkStates;
}

void AbstractTransportReactionSystem::ReactSystem(const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc,std::vector<double>& changeCellConc)
{
    // update the reaction system if any variables depend on the current system concentrations
    UpdateReactionSystem(currentBulkConcentration,currentCellConcentration);

    std::vector<double> deltaBulkConcentration(mNumBulkStates, 0.0);
    std::vector<double> deltaCellConcentration(mNumCellStates, 0.0);

    for (std::vector<AbstractTransportReaction*>::iterator reaction_iter = mpReactionVector.begin();
            reaction_iter != mpReactionVector.end();
            ++reaction_iter)
    {
        // iterate through the reactions in the system, modify a temporary change of concentration
        // update the system change in concentration
        deltaBulkConcentration.resize(mNumBulkStates, 0.0);
        deltaCellConcentration.resize(mNumCellStates, 0.0);

        AbstractTransportReaction *p_system_reaction = static_cast<AbstractTransportReaction*>(*reaction_iter);

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

void AbstractTransportReactionSystem::UpdateReactionSystem(const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    return;
}

void AbstractTransportReactionSystem::SetReactionVector(std::vector<AbstractTransportReaction*> reactionVector)
{
    mpReactionVector = reactionVector;
}

void AbstractTransportReactionSystem::SetReactionByIndex(AbstractTransportReaction* reaction, unsigned index)
{
    if (index < mNumReactions)
    {
        mpReactionVector[index] = reaction;
    }
}

void AbstractTransportReactionSystem::SetBulkChemistry(AbstractChemistry* bulkChemistry)
{
    mpBulkChemistry = bulkChemistry;
    mNumBulkStates = bulkChemistry->GetNumberChemicals();
}

void AbstractTransportReactionSystem::SetCellChemistry(AbstractChemistry* cellChemistry)
{
    mpCellChemistry = cellChemistry;
    mNumCellStates = cellChemistry->GetNumberChemicals();
}

void AbstractTransportReactionSystem::SetNumReactions(unsigned numReactions)
{
    mNumReactions = numReactions;
}

std::vector<AbstractTransportReaction*> AbstractTransportReactionSystem::GetReactionVector()
{
    return mpReactionVector;
}

AbstractTransportReaction* AbstractTransportReactionSystem::GetReactionByIndex(unsigned index)
{
    if (index < mNumReactions)
    {
        return mpReactionVector[index];
    }else{
        return new AbstractTransportReaction();
    }
}

AbstractChemistry* AbstractTransportReactionSystem::GetBulkChemistry()
{
    return mpBulkChemistry;
}

AbstractChemistry* AbstractTransportReactionSystem::GetCellChemistry()
{
    return mpCellChemistry;
}

unsigned AbstractTransportReactionSystem::GetNumReactions()
{
    return mNumReactions;
}

unsigned AbstractTransportReactionSystem::GetNumBulkStates()
{
    return mNumBulkStates;
}

unsigned AbstractTransportReactionSystem::GetNumCellStates()
{
    return mNumCellStates;
}

#endif