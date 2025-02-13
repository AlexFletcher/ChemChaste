#include "ChemicalSrnModel.hpp"

ChemicalSrnModel::ChemicalSrnModel(AbstractReactionSystem* pReactionSystem,boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(pReactionSystem->GetSystemChemistry()->GetNumberChemicals(), pOdeSolver),
    mpReactionSystem(pReactionSystem),
    mpCellChemistry(pReactionSystem->GetSystemChemistry())
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<ChemicalSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<ChemicalSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

ChemicalSrnModel::ChemicalSrnModel(const ChemicalSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    assert(rModel.GetOdeSystem());

    mpReactionSystem = rModel.mpReactionSystem;
    mpCellChemistry = rModel.mpCellChemistry;

    SetOdeSystem(new AbstractChemicalOdeSystem(mpReactionSystem));

    std::vector<double> stateVector = rModel.GetOdeSystem()->rGetStateVariables();
    mpOdeSystem->SetStateVariables(stateVector);
}

AbstractSrnModel* ChemicalSrnModel::CreateSrnModel()
{
    return new ChemicalSrnModel(*this);
}

void ChemicalSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise( new AbstractChemicalOdeSystem(mpReactionSystem));
}

void ChemicalSrnModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateOdeParameters();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void ChemicalSrnModel::SetReactionSystem(AbstractReactionSystem* reactionSystem)
{
    mpReactionSystem = reactionSystem;
}

AbstractReactionSystem* ChemicalSrnModel::GetReactionSystem()
{
    return mpReactionSystem;
}

void ChemicalSrnModel::SetCellChemistry(AbstractChemistry* chemistry)
{
    mpCellChemistry = chemistry;
}

AbstractChemistry* ChemicalSrnModel::GetCellChemistry()
{
    return mpCellChemistry;
}

void ChemicalSrnModel::UpdateOdeStatesFromCellData()
{
    unsigned numSpecies = mpCellChemistry->GetNumberChemicals();

    double species_value = 0.0;
    std::string species_name = "";
    std::vector<double> current_state_values(numSpecies, 0.0);
    for (unsigned i = 0; i < numSpecies; ++i)
    {
        species_name = mpCellChemistry->GetChemicalNamesByIndex(i);
        species_value = mpCell->GetCellData()->GetItem(species_name);
        current_state_values[i] = species_value;
    }

    mpOdeSystem->SetStateVariables(current_state_values);
}

void ChemicalSrnModel::UpdateOdeParameters()
{
}

double ChemicalSrnModel::GetStateValueByIndex(unsigned index)
{
    return mpOdeSystem->rGetStateVariables()[index];
}

double ChemicalSrnModel::GetStateValueByName(std::string name)
{
    unsigned index = mpCellChemistry->GetChemicalIndexByName(name);
    return GetStateValueByIndex(index);
}

void ChemicalSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}