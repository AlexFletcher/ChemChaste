#include "ChemicalSrnFromFile.hpp"

ChemicalSrnFromFile::ChemicalSrnFromFile(std::string reactionFilename) 
    : mSrnReactionFilename(reactionFilename)
{
    AbstractReactionSystemFromFile* p_file_reaction_system = new  AbstractReactionSystemFromFile(reactionFilename);

    AbstractChemistry* this_cell_chemistry = p_file_reaction_system->GetSystemChemistry();
    unsigned numChemicals = this_cell_chemistry->GetNumberChemicals();

    ChemicalSrnModel* pChemicalSrnModel = new ChemicalSrnModel(p_file_reaction_system);
    SetChemicalSrnModel(pChemicalSrnModel);

    mpChemicalSrnModel->SetReactionSystem(p_file_reaction_system);
    mpChemicalSrnModel->Initialise();
}

void ChemicalSrnFromFile::SetChemicalSrnModel(ChemicalSrnModel* pChemicalSrnModel)
{
    mpChemicalSrnModel = pChemicalSrnModel;
}

ChemicalSrnModel* ChemicalSrnFromFile::GetChemicalSrnModel()
{
    return mpChemicalSrnModel;
}

std::string ChemicalSrnFromFile::GetReactionFileName()
{
    return mSrnReactionFilename;
}