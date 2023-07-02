#include "AbstractChemistry.hpp"

AbstractChemistry::AbstractChemistry()
{
    mChemicalVector = std::vector<AbstractChemical*>();
    mChemicalNames = std::vector<std::string>();
    mChemicalDimensions = std::vector<std::string>();
    mNumberChemicals = 0;
}

AbstractChemistry::~AbstractChemistry()
{
}

void AbstractChemistry::AddChemistry(AbstractChemistry* pChemistry)
{
    /*
     * Function to combine the chemistries in the effort to form an overall 
     * union of chemistries, for use as domain. Domain chemistry determines 
     * the solver state variable ordering. Chemistry classes need not be of 
     * the same type, hence derived classes will need to determine the 
     * "highest class" define the highest class to be the parent chemistry. 
     * Want to add the chemical vectors while preventing duplicates.
     */
    std::vector<AbstractChemical*> chemical_vector = pChemistry->rGetChemicalVector();
    for (auto chemical_iter = chemical_vector.begin(); 
         chemical_iter != chemical_vector.end();
         ++chemical_iter)
    {
        /* 
         * Add each of the chemicals from the additional chemistry in turn 
         * (this checks for duplicates implicitly).
         */
        AddChemical(dynamic_cast<AbstractChemical*>(*chemical_iter));
    }
}

void AbstractChemistry::AddChemical(AbstractChemical* pChemical)
{
    bool new_chem = false;
    if (mChemicalVector.empty())
    {
        new_chem = true;
    }
    else
    {
        new_chem = CheckChemical(chemical);
    }

    if (new_chem)
    {    
        UpdateChemicalVectors(chemical);
    }
}

bool AbstractChemistry::CheckChemical(AbstractChemical* pChemical)
{
    bool new_chem  = true;
    for (auto chemical_iter = mChemicalVector.begin(); 
         chemical_iter != mChemicalVector.end();
         ++chemical_iter)
    {    
        if (pChemical->GetChemicalName() == dynamic_cast<AbstractChemical*>(*chemical_iter)->GetChemicalName())
        {
            new_chem = false;
            break;
        }
    }
    return new_chem;
}

void AbstractChemistry::UpdateChemicalVectors(AbstractChemical* pChemical)
{
    mChemicalVector.push_back(pChemical);

    // Properties needed for general function
    mChemicalNames.push_back(pChemical->GetChemicalName()); // may remove this, will lead to duplicates
    mChemicalDimensions.push_back(pChemical->GetChemicalDimensions()); // may remove this, will lead to duplicates

    // Further properties, to be overridded in derived classes
    mNumberChemicals++;
}

void AbstractChemistry::SetChemicalVector(std::vector<AbstractChemical*> chemicalVector)
{
    mChemicalVector = chemicalVector;
}

void AbstractChemistry::SetChemicalNames(std::vector<std::string> chemicalNames)
{
    mChemicalNames = chemicalNames;
}

void AbstractChemistry::SetChemicalDimensions(std::vector<std::string> chemicalDimensions)
{
    mChemicalDimensions = chemicalDimensions;
}

std::vector<AbstractChemical*> AbstractChemistry::rGetChemicalVector()
{
    return mChemicalVector;
}

std::vector<std::string> AbstractChemistry::GetChemicalNames()
{
    return mChemicalNames;
}

std::string AbstractChemistry::GetChemicalNamesByIndex(unsigned index)
{
    if (index < mNumberChemicals)
    {
        return mChemicalNames[index];
    }
    else
    {
        ///\todo replace with EXCEPTION
        std::cout << "Error: AbstractChemistry::GetChemicalNamesByIndex(unsigned index), index out of bounds" << std::endl;
        return "Null";
    } 
}

unsigned AbstractChemistry::GetChemicalIndexByName(std::string name)
{
    unsigned index = 0;
    for (auto chemical_iter = mChemicalVector.begin(); 
         chemical_iter != mChemicalVector.end();
         ++chemical_iter)
    {
        if (name == dynamic_cast<AbstractChemical*>(*chemical_iter)->GetChemicalName())
        {
            break;
        }
        index++;
    }

    // If name is not one of the chemicals then index == mNumberChemicals; would throw error later on
    return index;
}

std::vector<std::string> AbstractChemistry::GetChemicalDimensions()
{
    return mChemicalDimensions;
}

std::string AbstractChemistry::GetChemicalDimensionsByIndex(unsigned index)
{
    if (index < mNumberChemicals)
    {
        return mChemicalDimensions[index];
    }
    else
    {
        ///\todo replace with EXCEPTION
        std::cout << "Error: AbstractChemistry::GetChemicalDimensionsByIndex(unsigned index), index out of bounds" << std::endl;
        return "Null";
    } 
}

std::string AbstractChemistry::GetChemistryType()
{
    return "AbstractChemistry";
}

unsigned AbstractChemistry::GetNumberChemicals()
{
    return mNumberChemicals;
}
