#ifndef CHEMICALSRNFROMFILE_HPP_
#define CHEMICALSRNFROMFILE_HPP_

#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include "AbstractReactionSystem.hpp"
#include "AbstractReactionSystemFromFile.hpp"
#include "AbstractReaction.hpp"
#include "AbstractChemistry.hpp"
#include "ReactionTypeDatabase.hpp"
#include "ChemicalSrnModel.hpp"
 
/**
 * \todo Document class.
 */
class ChemicalSrnFromFile
{
protected:

    std::string mSrnReactionFilename;

    ChemicalSrnModel* mpChemicalSrnModel;

public:

    ChemicalSrnFromFile(std::string);

    virtual ~ChemicalSrnFromFile()
    {
    };

    void SetChemicalSrnModel(ChemicalSrnModel*);

    ChemicalSrnModel* GetChemicalSrnModel();

    std::string GetReactionFileName();

};

#endif /* CHEMICALSRNFROMFILE_HPP_ */