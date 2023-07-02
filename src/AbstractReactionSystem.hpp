#ifndef ABSTRACTREACTIONSYSTEM_HPP_
#define ABSTRACTREACTIONSYSTEM_HPP_

#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include "AbstractReaction.hpp"
#include "AbstractChemistry.hpp"

/**
 * Class to hold the information and data structures defining a system of 
 * chemical reactions.
 */
class AbstractReactionSystem
{
protected:

    AbstractChemistry* mpSystemChemistry;

    std::vector<AbstractReaction*> mpReactionVector;

    unsigned mNumReactions;

public:
    AbstractReactionSystem( AbstractChemistry* systemChemistry = new AbstractChemistry(), 
                            std::vector<AbstractReaction*> reactionVector = std::vector<AbstractReaction*>());

    virtual ~AbstractReactionSystem()
    {
    };

    AbstractReactionSystem(const AbstractReactionSystem&);

    virtual void ReactSystem(const std::vector<double>&, std::vector<double>&);

    virtual void UpdateReactionSystem(const std::vector<double>&);

    std::vector<AbstractReaction*> GetReactionVector();

    AbstractReaction* GetReactionByIndex(unsigned);

    void SetReactionVector(std::vector<AbstractReaction*>);

    void SetReactionByIndex(AbstractReaction*, unsigned);

    AbstractChemistry* GetSystemChemistry();

    void SetSystemChemistry(AbstractChemistry*);

    unsigned GetNumReactions();

    void SetNumReactions(unsigned);
};

#endif /* ABSTRACTREACTIONSYSTEM_HPP_ */