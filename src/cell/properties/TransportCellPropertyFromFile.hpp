#ifndef TRANSPORTCELLPROPERTYFROMFILE_HPP
#define TRANSPORTCELLPROPERTYFROMFILE_HPP

#include <vector>
#include <boost/shared_ptr.hpp>

#include "TransportCellProperty.hpp"
#include "ChemicalCell.hpp"
#include "AbstractTransportReactionSystemFromFile.hpp"

/**
 * \todo Document class.
 */
class TransportCellPropertyFromFile
{
protected:

    std::string mTransportFilename;

    boost::shared_ptr<TransportCellProperty> mpTransportProperty;

public:

    TransportCellPropertyFromFile(std::string filename ="");

    ~TransportCellPropertyFromFile()
    {
    };

    void SetUpTransportProperty(CellPtr pCell);

    void SetTransportFilename(std::string);

    void SetTransportProperty(boost::shared_ptr<TransportCellProperty>);

    std::string GetTransportFilename();

    boost::shared_ptr<TransportCellProperty> GetTransportProperty();
};

#endif