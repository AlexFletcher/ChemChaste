#ifndef MEMBRANECELLPROPERTYFROMFILE_HPP
#define MEMBRANECELLPROPERTYFROMFILE_HPP

//general includes
#include <vector>
#include <boost/shared_ptr.hpp>

#include "MembraneCellProperty.hpp"
#include "ChemicalCell.hpp"
#include "AbstractMembraneReactionSystemFromFile.hpp"

/**
 * \todo Document this class.
 */
class MembraneCellPropertyFromFile
{
protected:

    std::string mMembraneFilename;

    boost::shared_ptr<MembraneCellProperty> mpMembraneProperty;

public:

    MembraneCellPropertyFromFile(std::string filename="");

    ~MembraneCellPropertyFromFile()
    {
    };

    void SetUpMembraneProperty(CellPtr pCell);

    void SetMembraneFilename(std::string);

    void SetMembraneProperty(boost::shared_ptr<MembraneCellProperty>);

    std::string GetMembraneFilename();

    boost::shared_ptr<MembraneCellProperty> GetMembraneProperty();
};

#endif