#include "MembraneCellPropertyFromFile.hpp"

MembraneCellPropertyFromFile::MembraneCellPropertyFromFile(std::string filename)
    : mMembraneFilename(filename)
{
    if (filename != "")
    {
        SetMembraneFilename(filename);
    }
}

void MembraneCellPropertyFromFile::SetUpMembraneProperty(CellPtr pCell)
{
    boost::shared_ptr<MembraneCellProperty> p_membrane = boost::static_pointer_cast<MembraneCellProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());
    AbstractMembraneReactionSystemFromFile* p_membrane_system_from_file = new AbstractMembraneReactionSystemFromFile(mMembraneFilename);

    p_membrane->SetUp(p_membrane_system_from_file, pCell);
    SetMembraneProperty(p_membrane);
}

void MembraneCellPropertyFromFile::SetMembraneFilename(std::string filename)
{
    mMembraneFilename = filename;
}

void MembraneCellPropertyFromFile::SetMembraneProperty(boost::shared_ptr<MembraneCellProperty> pMembraneProperty)
{
    mpMembraneProperty = pMembraneProperty;
}

std::string MembraneCellPropertyFromFile::GetMembraneFilename()
{
    return mMembraneFilename;
}

boost::shared_ptr<MembraneCellProperty> MembraneCellPropertyFromFile::GetMembraneProperty()
{
    return mpMembraneProperty;
}
