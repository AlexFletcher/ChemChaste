#ifndef ABSTRACTTRANSPORTREACTION_HPP_
#define ABSTRACTTRANSPORTREACTION_HPP_

#include <string>
#include <stdlib.h> 
#include <vector>
#include <iostream>

#include "AbstractChemistry.hpp"
#include "AbstractChemical.hpp"

/**
 * Abstract property to contain information about the interactions of chemical 
 * species at a cell boundary:
 * bulk->cell
 * 
 * Reaction substrates denote the species in the bulk while products denote 
 * species in the cell.
 */
class AbstractTransportReaction 
{
protected:

    // vector with AbstractChemical of bulk species, in the forward reaction sense.  To be used to track system concentrations
    std::vector<AbstractChemical*> mpBulkReactionSpecies;

    // vector with the AbstractChemical cell species, in the forward reaction sense.  To be used to track system concentrations
    std::vector<AbstractChemical*> mpCellReactionSpecies;

    unsigned mNumBulkSpecies;

    unsigned mNumCellSpecies;


    // vector containing the stoichmetry of the substrates, bulk species.  Will be "-ve" in reaction
    std::vector<unsigned> mStoichBulk;

    // vector containing the stoichmetry of the products, cell species.  Will be "+ve" in reaction
    std::vector<unsigned> mStoichCell;

    // container for the reaction rate constant.  '+ve' value. May be overriden in derived classes to account for functional rates
    double mReactionRate;

    bool mIsRateCheck;

    double mDeltaRateMin;

    double mDeltaRateMax;

    std::string mIrreversibleDelimiter = "->";
  
    std::string mIrreversibleRateName = "kf =";

public:

    // constructor
    AbstractTransportReaction(  std::vector<AbstractChemical*> bulkReactionSpecies = std::vector<AbstractChemical*>(),
                                std::vector<AbstractChemical*> cellReactionSpecies = std::vector<AbstractChemical*>(),
                                std::vector<unsigned> stoichBulk = std::vector<unsigned>(),
                                std::vector<unsigned> stoichCell = std::vector<unsigned>(),
                                double reactionRate = 1.0
    );
    
    // copy constructor
    AbstractTransportReaction(const AbstractTransportReaction&);

    // destructor
    virtual ~AbstractTransportReaction()
    {
    };


    // function to take in pointer to current concentration state vector of the state vector for change in cocnentration.  Basic update via multiplying by constant reaction rate
    virtual void React(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual void SetReactionRate(double);

    virtual double GetReactionRate();

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    std::vector<AbstractChemical*> GetBulkSpecies();

    AbstractChemical* GetBulkSpeciesByIndex(unsigned);

    std::vector<AbstractChemical*> GetCellSpecies();

    AbstractChemical* GetCellSpeciesByIndex(unsigned);

    void SetBulkSpecies(std::vector<AbstractChemical*>);

    void SetCellSpecies(std::vector<AbstractChemical*>);

    std::vector<unsigned> GetStoichBulk();

    unsigned GetStoichBulkByIndex(unsigned);

    void SetStoichBulk(std::vector<unsigned>);

    std::vector<unsigned> GetStoichCell();

    unsigned GetStoichCellByIndex(unsigned);

    void SetStoichCell(std::vector<unsigned>);

    void SetNumBulkSpecies(unsigned);

    unsigned GetNumBulkSpecies();

    void SetNumCellSpecies(unsigned);

    unsigned GetNumCellSpecies();

    void SetIsRateCheck(bool);

    bool GetIsRateCheck();

    void SetDeltaErrorRateMin(double);

    double GetDeltaErrorRateMin();

    void SetDeltaErrorRateMax(double);

    double GetDeltaErrorRateMax();

    double CheckRate(double);

    void SetIrreversibleDelimiter(std::string);

    std::string GetIrreversibleDelimiter();

    void SetIrreversibleRateName(std::string);

    std::string GetIrreversibleRateName();
};

AbstractTransportReaction::AbstractTransportReaction(   std::vector<AbstractChemical*> bulkReactionSpecies,
                                                        std::vector<AbstractChemical*> cellReactionSpecies,
                                                        std::vector<unsigned> stoichBulk,
                                                        std::vector<unsigned> stoichCell,
                                                        double reactionRate)
        :   mpBulkReactionSpecies(bulkReactionSpecies),
            mpCellReactionSpecies(cellReactionSpecies),
            mStoichBulk(stoichBulk),
            mStoichCell(stoichCell),
            mReactionRate(reactionRate)
{
    mNumCellSpecies = cellReactionSpecies.size();
    mNumBulkSpecies = bulkReactionSpecies.size();

    SetIsRateCheck(true);
    SetDeltaErrorRateMin(1e-8);
    SetDeltaErrorRateMax(1e8);
}

AbstractTransportReaction::AbstractTransportReaction(const AbstractTransportReaction& existingReaction)
{
    
    mpBulkReactionSpecies = existingReaction.mpBulkReactionSpecies;
    mpCellReactionSpecies = existingReaction.mpCellReactionSpecies;
    mStoichBulk = existingReaction.mStoichBulk;
    mStoichCell = existingReaction.mStoichCell;
    mReactionRate = existingReaction.mReactionRate;
    mNumCellSpecies = existingReaction.mNumCellSpecies;
    mNumBulkSpecies = existingReaction.mNumBulkSpecies;
    mIsRateCheck = existingReaction.mIsRateCheck;
    mDeltaRateMin = existingReaction.mDeltaRateMin;
    mDeltaRateMax = existingReaction.mDeltaRateMax;
}

void AbstractTransportReaction::React(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
{   
    std::vector<AbstractChemical*> p_bulk_chemical_vector = bulkChemistry->rGetChemicalVector();
    
    std::vector<AbstractChemical*> p_cell_chemical_vector = cellChemistry->rGetChemicalVector();
    
    UpdateReactionRate(bulkChemistry, cellChemistry, currentBulkConcentration, currentCellConcentration);
    
    // perform the reaction

    // run through the bulk species
    unsigned index=0;
    for (std::vector<AbstractChemical*>::iterator chem_iter = p_bulk_chemical_vector.begin();
         chem_iter != p_bulk_chemical_vector.end();
         ++chem_iter, ++index)
    {
        AbstractChemical *p_bulk_chemical = static_cast<AbstractChemical*>(*chem_iter);

        // for each bulk chemical, parse whether it is involved in this reaction.
        for (unsigned j=0; j<mNumBulkSpecies; j++)
        {
            if (mpBulkReactionSpecies.at(j)->GetChemicalName()==p_bulk_chemical->GetChemicalName())
            {
                changeBulkConc.at(index)  -= mStoichBulk.at(j)*GetReactionRate();
                break;
            }
        }
    }

    // run through the cell species
    index = 0;
    for (std::vector<AbstractChemical*>::iterator chem_iter = p_cell_chemical_vector.begin();
         chem_iter != p_cell_chemical_vector.end();
         ++chem_iter, ++index)
    {
        AbstractChemical *p_cell_chemical = static_cast<AbstractChemical*>(*chem_iter);
   
        for (unsigned j=0; j<mNumCellSpecies; j++)
        {
            if (mpCellReactionSpecies.at(j)->GetChemicalName()==p_cell_chemical->GetChemicalName())
            {
                changeCellConc.at(index) += mStoichCell.at(j) *GetReactionRate();
                break;
            }
        }
    }
}

void AbstractTransportReaction::UpdateReaction()
{
    return;
}

void AbstractTransportReaction::UpdateReactionRate(AbstractChemistry* bulkChemistry, AbstractChemistry* cellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // Default to zeroth order reaction
    SetReactionRate(mReactionRate);
}

void AbstractTransportReaction::SetReactionRate(double reactionRate)
{
    if (mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    
    mReactionRate = reactionRate;
}

double AbstractTransportReaction::GetReactionRate()
{
    return mReactionRate;
}

std::string AbstractTransportReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderTransportIntoCell";
}

void AbstractTransportReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    if (reaction_information.find(mIrreversibleRateName) != std::string::npos)
    {
        size_t pos = reaction_information.find(mIrreversibleRateName);
        SetReactionRate(atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1, std::string::npos).c_str()));
    }
}

std::vector<AbstractChemical*> AbstractTransportReaction::GetBulkSpecies()
{
    return mpBulkReactionSpecies;
}

AbstractChemical* AbstractTransportReaction::GetBulkSpeciesByIndex(unsigned index)
{
    if (index < mNumBulkSpecies)
    {
        return mpBulkReactionSpecies[index];
    }
    else
    {
        std::cout << "Error: AbstractTransportReaction::GetBulkSpeciesByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractTransportReaction::GetCellSpecies()
{
    return mpCellReactionSpecies;
}

AbstractChemical* AbstractTransportReaction::GetCellSpeciesByIndex(unsigned index)
{
    if (index < mNumCellSpecies)
    {
        return mpCellReactionSpecies[index];
    }
    else
    {
        std::cout << "Error: AbstractTransportReaction::GetCellSpeciesByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

void AbstractTransportReaction::SetBulkSpecies(std::vector<AbstractChemical*> bulkReactionSpecies)
{
    mpBulkReactionSpecies = bulkReactionSpecies;
}

void AbstractTransportReaction::SetCellSpecies(std::vector<AbstractChemical*> cellReactionSpecies)
{
    mpCellReactionSpecies = cellReactionSpecies;
}

std::vector<unsigned> AbstractTransportReaction::GetStoichBulk()
{
    return mStoichBulk;
}

unsigned AbstractTransportReaction::GetStoichBulkByIndex(unsigned index)
{
    if (index < mNumBulkSpecies)
    {
        return mStoichBulk[index];
    }
    else
    {
        std::cout << "Error: AbstractTransportReaction::GetStoichBulkByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractTransportReaction::SetStoichBulk(std::vector<unsigned> stoichBulk)
{
    mStoichBulk = stoichBulk;
}

std::vector<unsigned> AbstractTransportReaction::GetStoichCell()
{
    return mStoichCell;
}

unsigned AbstractTransportReaction::GetStoichCellByIndex(unsigned index)
{
    if (index < mNumCellSpecies)
    {
        return mStoichCell[index];
    }
    else
    {
        std::cout << "Error: AbstractTransportReaction::GetStoichCellByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractTransportReaction::SetStoichCell(std::vector<unsigned> StoichCell)
{
    mStoichCell = StoichCell;
}

void AbstractTransportReaction::SetNumBulkSpecies(unsigned numBulkSpecies)
{
    mNumBulkSpecies = numBulkSpecies;
}

unsigned AbstractTransportReaction::GetNumBulkSpecies()
{
    return mNumBulkSpecies;
}

void AbstractTransportReaction::SetNumCellSpecies(unsigned numCellSpecies)
{
    mNumCellSpecies = numCellSpecies;
}

unsigned AbstractTransportReaction::GetNumCellSpecies()
{
    return mNumCellSpecies;
}

void AbstractTransportReaction::SetIsRateCheck(bool IsRateCheck)
{
    mIsRateCheck = IsRateCheck;
}

bool AbstractTransportReaction::GetIsRateCheck()
{
    return mIsRateCheck;
}

void AbstractTransportReaction::SetDeltaErrorRateMin(double delta_rate_min)
{
    mDeltaRateMin = delta_rate_min;
}

double AbstractTransportReaction::GetDeltaErrorRateMin()
{
    return mDeltaRateMin;
}

void AbstractTransportReaction::SetDeltaErrorRateMax(double delta_rate_max)
{
    mDeltaRateMax = delta_rate_max;
}

double AbstractTransportReaction::GetDeltaErrorRateMax()
{
    return mDeltaRateMax;
}

double AbstractTransportReaction::CheckRate(double rate)
{
    // if reaction rate gets too low or high then undefined behaviour can occur

    if (abs(rate) < mDeltaRateMin)
    {
        rate = 0.0;
    }else if (abs(rate) > mDeltaRateMax)
    {
        rate = mDeltaRateMax;
    }
    return rate;
}

void AbstractTransportReaction::SetIrreversibleDelimiter(std::string delim)
{
    mIrreversibleDelimiter = delim;
}

std::string AbstractTransportReaction::GetIrreversibleDelimiter()
{
    return mIrreversibleDelimiter;
}

void AbstractTransportReaction::SetIrreversibleRateName(std::string rateName)
{
    mIrreversibleRateName = rateName;
}

std::string AbstractTransportReaction::GetIrreversibleRateName()
{
    return mIrreversibleRateName;
}

#endif