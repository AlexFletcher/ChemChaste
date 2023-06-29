#include "SimpleChemicalThresholdCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "ApoptoticCellProperty.hpp"

SimpleChemicalThresholdCellCycleModel::SimpleChemicalThresholdCellCycleModel()
{ 
}

SimpleChemicalThresholdCellCycleModel::SimpleChemicalThresholdCellCycleModel(const SimpleChemicalThresholdCellCycleModel& rModel)
    : AbstractSimplePhaseBasedCellCycleModel(rModel),
      mCurrentStarvationDuration(rModel.mCurrentStarvationDuration),
      mCurrentStarvationOnsetTime(rModel.mCurrentStarvationOnsetTime),
      mCriticalStarvationDuration(rModel.mCriticalStarvationDuration),
      mpThresholdChemistry(rModel.mpThresholdChemistry),
      mNumberThresholdSpecies(rModel.mNumberThresholdSpecies),
      mMaxThresholdConcentrationValues(rModel.mMaxThresholdConcentrationValues),
      mIsMaximumThresholdSet(rModel.mIsMaximumThresholdSet),
      mMinThresholdConcentrationValues(rModel.mMinThresholdConcentrationValues),
      mIsMinimumThresholdSet(rModel.mIsMinimumThresholdSet)
{
}

SimpleChemicalThresholdCellCycleModel::~SimpleChemicalThresholdCellCycleModel()
{
}

void SimpleChemicalThresholdCellCycleModel::SetUp(AbstractChemistry* pChemistry)
{
    SetThresholdChemistry(pChemistry);
    SetNumberThresholdSpecies(pChemistry->GetNumberChemicals());

    std::vector<double> maxThresholdConcentrationValues(mNumberThresholdSpecies, 0.0);
    SetMaximumSpeciesThreshold(maxThresholdConcentrationValues);

    std::vector<bool> isMaximumThresholdSet(mNumberThresholdSpecies, false);
    SetMaximumThresholdCheck(isMaximumThresholdSet);

    std::vector<double> minThresholdConcentrationValues(mNumberThresholdSpecies, 0.0);
    SetMinimumSpeciesThreshold(minThresholdConcentrationValues);

    std::vector<bool> isMinimumThresholdSet(mNumberThresholdSpecies, false);
    SetMinimumThresholdCheck(isMinimumThresholdSet);

    std::vector<double> mSpeciesConcentrations(mNumberThresholdSpecies, 0.0);
    SetSpeciesConcentrations(mSpeciesConcentrations);

    mCurrentStarvationOnsetTime = SimulationTime::Instance()->GetTime();
}

void SimpleChemicalThresholdCellCycleModel::Initialise()
{
    AbstractSimplePhaseBasedCellCycleModel::Initialise();
}

void SimpleChemicalThresholdCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}

double SimpleChemicalThresholdCellCycleModel::GetCurrentStarvationDuration() const
{
    return mCurrentStarvationDuration;
}

double SimpleChemicalThresholdCellCycleModel::GetCurrentStarvationOnsetTime() const
{
    return mCurrentStarvationOnsetTime;
}

void SimpleChemicalThresholdCellCycleModel::UpdateCellCyclePhase()
{
    bool cell_is_apoptotic = mpCell->HasCellProperty<ApoptoticCellProperty>();

    if (!cell_is_apoptotic)
    {
        UpdateStarvationDuration();

        AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();

        if (mCurrentCellCyclePhase == G_ONE_PHASE)
        {
            double dt = SimulationTime::Instance()->GetTimeStep();
            mG1Duration += dt;
        }
    }
}

bool SimpleChemicalThresholdCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);
    
    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();

        if ((mCurrentCellCyclePhase != G_ZERO_PHASE) && ConcentrationAboveMaxThreshold())
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
}

AbstractCellCycleModel* SimpleChemicalThresholdCellCycleModel::CreateCellCycleModel()
{
    return new SimpleChemicalThresholdCellCycleModel(*this);
}

void SimpleChemicalThresholdCellCycleModel::UpdateStarvationDuration()
{
    assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());

    std::vector<std::string> chemicalNames = mpThresholdChemistry->GetChemicalNames();

    // Get concentration vector
    for (unsigned species=0; species<mNumberThresholdSpecies; species++)
    {
        SetSpeciesConcentrationsByIndex(mpCell->GetCellData()->GetItem(mpThresholdChemistry->GetChemicalNamesByIndex(species)), species);
    }

    // Check each species for starvation
    bool is_starving = false;
    for (unsigned species = 0; species < mNumberThresholdSpecies; species++)
    {
        if (GetMinimumThresholdCheckByIndex(species) && 
            GetSpeciesConcentrationsByIndex(species) < GetMinimumThresholdByIndex(species))
        {
            mCurrentStarvationDuration = SimulationTime::Instance()->GetTime() - mCurrentStarvationOnsetTime;
 
            if (mCurrentStarvationDuration > mCriticalStarvationDuration && 
                RandomNumberGenerator::Instance()->ranf() < CellDeathProbability({GetSpeciesConcentrationsByIndex(species), GetMinimumThresholdByIndex(species)}))
            {
                boost::shared_ptr<AbstractCellProperty> p_apoptotic_property = mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
                mpCell->AddCellProperty(p_apoptotic_property);
                is_starving = true;
                break;
            }
        }
    }
    if (!is_starving)
    {
        // Reset the cell's starvation duration and update the time at which the onset of Starvation occurs
        mCurrentStarvationDuration = 0.0;
        mCurrentStarvationOnsetTime = SimulationTime::Instance()->GetTime();
    }
}

double SimpleChemicalThresholdCellCycleModel::CellDeathProbability(std::vector<double> parameters)
{
    double prob_of_death = 0.9;
    return prob_of_death;
}

bool SimpleChemicalThresholdCellCycleModel::ConcentrationAboveMaxThreshold()
{
    // Search through all the threhsold species for whether nay are above the maximum threshold
    for (unsigned species = 0; species < mNumberThresholdSpecies; species++)
    {
        if (GetMaximumThresholdCheckByIndex(species) && 
            GetSpeciesConcentrationsByIndex(species) > GetMaximumThresholdByIndex(species))
        {
            return true;
        }
    }
    return false;
}

bool SimpleChemicalThresholdCellCycleModel::ConcentrationBelowMinThreshold()
{
    // search through all the threhsold species for whether nay are above the maximum threshold
    for (unsigned species=0; species<mNumberThresholdSpecies; species++)
    {
        if (GetMinimumThresholdCheckByIndex(species) && GetSpeciesConcentrationsByIndex(species)<GetMinimumThresholdByIndex(species))
        {
            return true;
        }
    }
    return false;
}

void SimpleChemicalThresholdCellCycleModel::PrepareForDivision()
{   
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumSpeciesThreshold(std::vector<double> maxThreshold)
{
    mMaxThresholdConcentrationValues = maxThreshold;
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumSpeciesThreshold(std::vector<double> minThreshold)
{
    mMinThresholdConcentrationValues = minThreshold;
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdByIndex(double threshold, unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        mMaxThresholdConcentrationValues[index] = threshold;
    }
    else
    {
        ///\todo Replace with EXCEPTION
        std::cout << "Error: SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdByIndex, index out of bounds" << std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdByIndex(double threshold, unsigned index)
{
    if (index<mNumberThresholdSpecies)
    {
        mMinThresholdConcentrationValues[index] = threshold;
    }
    else
    {
        ///\todo Replace with EXCEPTION
        std::cout << "Error: SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdByIndex, index out of bounds" << std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetSpeciesConcentrations(std::vector<double> concentrations)
{
    mSpeciesConcentrations = concentrations;
}

void SimpleChemicalThresholdCellCycleModel::SetSpeciesConcentrationsByIndex(double threshold, unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        mSpeciesConcentrations[index] = threshold;
    }
    else
    {
        ///\todo Replace with EXCEPTION
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetSpeciesConcentrationsByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetNumberThresholdSpecies(unsigned speciesNumber)
{
    mNumberThresholdSpecies = speciesNumber;
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdCheck(std::vector<bool> thresholdCheck)
{
    mIsMaximumThresholdSet = thresholdCheck;
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdCheck(std::vector<bool> thresholdCheck)
{
    mIsMinimumThresholdSet = thresholdCheck;
}

void SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdCheckByIndex(bool thresholdCheck, unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        mIsMaximumThresholdSet[index] = thresholdCheck;
    }
    else
    {
        ///\todo Replace with EXCEPTION
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetMaximumThresholdCheckByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdCheckByIndex(bool thresholdCheck, unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        mIsMinimumThresholdSet[index] = thresholdCheck;
    }
    else
    {
        ///\todo Replace with EXCEPTION
        std::cout<<"Error: SimpleChemicalThresholdCellCycleModel::SetMinimumThresholdCheckByIndex, index out of bounds"<<std::endl;
    }
}

void SimpleChemicalThresholdCellCycleModel::SetThresholdChemistry(AbstractChemistry* p_chemistry)
{
    mpThresholdChemistry = p_chemistry;
}

AbstractChemistry* SimpleChemicalThresholdCellCycleModel::GetThresholdChemistry()
{
    return mpThresholdChemistry;
}

std::vector<double> SimpleChemicalThresholdCellCycleModel::GetMaximumSpeciesThreshold()
{
    return mMaxThresholdConcentrationValues;
}

std::vector<double> SimpleChemicalThresholdCellCycleModel::GetMinimumSpeciesThreshold()
{
    return mMinThresholdConcentrationValues;
}

double SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdByIndex(unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        return mMaxThresholdConcentrationValues[index];
    }

    ///\todo Replace with EXCEPTION
    std::cout << "Error: SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdByIndex, index out of bounds" << std::endl;
    return 0.0;
}

double SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdByIndex(unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        return mMinThresholdConcentrationValues[index];
    }

    ///\todo Replace with EXCEPTION
    std::cout << "Error: SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdByIndex, index out of bounds" << std::endl;
    return 0.0;
}

std::vector<double> SimpleChemicalThresholdCellCycleModel::GetSpeciesConcentrations()
{
    return mSpeciesConcentrations;
}

double SimpleChemicalThresholdCellCycleModel::GetSpeciesConcentrationsByIndex(unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        return mSpeciesConcentrations[index];
    }
    std::cout << "Error: SimpleChemicalThresholdCellCycleModel::GetSpeciesConcentrationsByIndex, index out of bounds" << std::endl;
    return 0.0;
}

unsigned SimpleChemicalThresholdCellCycleModel::GetNumberThresholdSpecies()
{
    return mNumberThresholdSpecies;
}

std::vector<bool> SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdCheck()
{
    return mIsMaximumThresholdSet;
}

std::vector<bool> SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdCheck()
{
    return mIsMinimumThresholdSet;
}

bool SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdCheckByIndex(unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        return mIsMaximumThresholdSet[index];
    }
    std::cout << "Error: SimpleChemicalThresholdCellCycleModel::GetMaximumThresholdCheckByIndex, index out of bounds" << std::endl;
    return false;
}

bool SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdCheckByIndex(unsigned index)
{
    if (index < mNumberThresholdSpecies)
    {
        return mIsMinimumThresholdSet[index];
    }

    ///\todo Replace with EXCEPTION
    std::cout << "Error: SimpleChemicalThresholdCellCycleModel::GetMinimumThresholdCheckByIndex, index out of bounds" << std::endl;
    return false;
}

double SimpleChemicalThresholdCellCycleModel::GetCriticalStarvationDuration() const
{
    return mCriticalStarvationDuration;
}

void SimpleChemicalThresholdCellCycleModel::SetCriticalStarvationDuration(double criticalStarvationDuration)
{
    assert(criticalStarvationDuration >= 0.0);
    mCriticalStarvationDuration = criticalStarvationDuration;
}

void SimpleChemicalThresholdCellCycleModel::SetCurrentStarvationOnsetTime(double currentStarvationOnsetTime)
{
    assert(currentStarvationOnsetTime >= 0.0);
    mCurrentStarvationOnsetTime = currentStarvationOnsetTime;
}