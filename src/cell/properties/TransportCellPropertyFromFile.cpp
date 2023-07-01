#include "TransportCellPropertyFromFile.hpp"

TransportCellPropertyFromFile::TransportCellPropertyFromFile(std::string filename)
    : mTransportFilename(filename)
{
    if (filename != "")
    {
        SetTransportFilename(filename);
    }
}

void TransportCellPropertyFromFile::SetUpTransportProperty(CellPtr pCell)
{
    boost::shared_ptr<TransportCellProperty> p_transport = boost::static_pointer_cast<TransportCellProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());
    AbstractTransportReactionSystemFromFile* p_transport_system_from_file = new AbstractTransportReactionSystemFromFile(mTransportFilename);

    p_transport->SetUp(p_transport_system_from_file, pCell);

    SetTransportProperty(p_transport);
}

void TransportCellPropertyFromFile::SetTransportFilename(std::string filename)
{
    mTransportFilename = filename;
}

void TransportCellPropertyFromFile::SetTransportProperty(boost::shared_ptr<TransportCellProperty> pTransportProperty)
{
    mpTransportProperty = pTransportProperty;
}

std::string TransportCellPropertyFromFile::GetTransportFilename()
{
    return mTransportFilename;
}

boost::shared_ptr<TransportCellProperty> TransportCellPropertyFromFile::GetTransportProperty()
{
    return mpTransportProperty;
}