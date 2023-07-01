#ifndef ABSTRACTREACTION_HPP_
#define ABSTRACTREACTION_HPP_

#include <string>
#include <stdlib.h> 
#include <vector>
#include <iostream>

#include "AbstractChemistry.hpp"
#include "AbstractChemical.hpp"

/**
 * Abstract property to contain information about the interactions of chemical 
 * species. Base class for further reaction types.
 */
class AbstractReaction 
{
protected:
    // vector with AbstractChemical substrates, in the forward reaction sense.  To be used to track system concentrations
    std::vector<AbstractChemical*> mpSubstrates;

    // vector with the AbstractChemical products, in the forward reaction sense.  To be used to track system concentrations
    std::vector<AbstractChemical*> mpProducts;

    unsigned mNumProducts;

    unsigned mNumSubstrates;

    // vector containing the stoichmetry of the substrates.  Will be "-ve" in reaction
    std::vector<unsigned> mStoichSubstrates;

    // vector containing the stoichmetry of the products.  Will be "+ve" in reaction
    std::vector<unsigned> mStoichProducts;

    // container for the reaction rate constant.  '+ve' value. May be overriden in derived classes to account for functional rates
    double mReactionRate;

    bool mIsRateCheck;

    double mDeltaRateMin;

    double mDeltaRateMax;

    std::string mIrreversibleDelimiter = "->";
  
    std::string mIrreversibleRateName = "kf =";

public:

    // constructor
    AbstractReaction(   std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                        std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                        std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                        double reactionRate = 1.0
    );

    // copy constructor
    AbstractReaction(const AbstractReaction&);

    // destructor
    virtual ~AbstractReaction()
    {
    };

    // function to take in pointer to current concentration state vector of the state vector for change in concentration.  Basic update via multiplying by constant reaction rate
    virtual void React(AbstractChemistry*, const std::vector<double>&, std::vector<double>&);

    // pure vitual function to update member variables and conditions of derived classes
    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual void SetReactionRate(double);

    virtual double GetReactionRate();

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    std::vector<AbstractChemical*> GetSubstrates();

    AbstractChemical* GetSubstratesByIndex(unsigned);

    std::vector<AbstractChemical*> GetProducts();

    AbstractChemical* GetProductsByIndex(unsigned);

    void SetSubstrates(std::vector<AbstractChemical*>);

    void SetProducts(std::vector<AbstractChemical*>);

    std::vector<unsigned> GetStoichSubstrates();

    unsigned GetStoichSubstratesByIndex(unsigned);

    void SetStoichSubstrates(std::vector<unsigned>);

    std::vector<unsigned> GetStoichProducts();

    unsigned GetStoichProductsByIndex(unsigned);

    void SetStoichProducts(std::vector<unsigned>);

    void SetNumSubstrates(unsigned);

    unsigned GetNumSubstrates();

    void SetNumProducts(unsigned);

    unsigned GetNumProducts();

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

    unsigned FindIndexOfLastDelimiterPosition(std::vector<std::string>, std::string);

    unsigned FindIndexOfFirstDelimiterPosition(std::vector<std::string>, std::string);
};

AbstractReaction::AbstractReaction( std::vector<AbstractChemical*> substrates,
                                    std::vector<AbstractChemical*> products,
                                    std::vector<unsigned> stoichSubstrates,
                                    std::vector<unsigned> stoichProducts,
                                    double reactionRate)
        :   mpSubstrates(substrates),
            mpProducts(products),
            mStoichSubstrates(stoichSubstrates),
            mStoichProducts(stoichProducts),
            mReactionRate(reactionRate)
{
    mNumProducts = products.size();
    mNumSubstrates  =substrates.size();

    SetIsRateCheck(true);
    SetDeltaErrorRateMin(1e-8);
    SetDeltaErrorRateMax(1e8);
}

AbstractReaction::AbstractReaction(const AbstractReaction& existingReaction)
{
    
    mpSubstrates = existingReaction.mpSubstrates;
    mpProducts = existingReaction.mpProducts;
    mStoichSubstrates = existingReaction.mStoichSubstrates;
    mStoichProducts = existingReaction.mStoichProducts;
    mReactionRate = existingReaction.mReactionRate;
    mNumProducts = existingReaction.mNumProducts;
    mNumSubstrates = existingReaction.mNumSubstrates;
}

void AbstractReaction::React(AbstractChemistry* systemChemistry, const std::vector<double>& currentChemistryConc, std::vector<double>& changeChemistryConc)
{

    std::vector<AbstractChemical*> p_chemical_vector = systemChemistry->rGetChemicalVector();

    UpdateReactionRate(systemChemistry, currentChemistryConc);
    
    // perform the reaction
    unsigned index=0;
    for (std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // for each system chemical, parse whether it is involved in this reaction.
        for (unsigned j=0; j<mNumSubstrates; j++)
        {
            if (mpSubstrates[j]->GetChemicalName()==p_system_chemical->GetChemicalName())
            {
                changeChemistryConc[index] -= mStoichSubstrates[j]*GetReactionRate();
                break;
            }
        }
        // a reactant may be present on both sides of the reaction and may convert at different functional rates
        for (unsigned j=0; j<mNumProducts; j++)
        {
            if (mpProducts[j]->GetChemicalName()==p_system_chemical->GetChemicalName())
            {
                changeChemistryConc[index] += mStoichProducts[j]*GetReactionRate();
                break;
            }
        }
    }
}

void AbstractReaction::UpdateReaction()
{
    return;
}

void AbstractReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentChemistryConc)
{
    // default to zeroth order reaction
    SetReactionRate(mReactionRate);
}

void AbstractReaction::SetReactionRate(double reactionRate)
{
    if (mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    
    mReactionRate = reactionRate;
}

double AbstractReaction::GetReactionRate()
{
    return mReactionRate;
}

std::string AbstractReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderReaction";
}

void AbstractReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    if (reaction_information.find(mIrreversibleRateName) != std::string::npos)
    {

        size_t pos= reaction_information.find(mIrreversibleRateName);

        SetReactionRate(atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1,std::string::npos).c_str()));
    }
}

std::vector<AbstractChemical*> AbstractReaction::GetSubstrates()
{
    return mpSubstrates;
}

AbstractChemical* AbstractReaction::GetSubstratesByIndex(unsigned index)
{
    if (index < mNumSubstrates)
    {
        return mpSubstrates[index];
    }
    else
    {
        std::cout << "Error: AbstractReaction::GetSubstratesByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractReaction::GetProducts()
{
    return mpProducts;
}

AbstractChemical* AbstractReaction::GetProductsByIndex(unsigned index)
{
    if (index < mNumProducts)
    {
        return mpProducts[index];
    }
    else
    {
        std::cout << "Error: AbstractReaction::GetProductsByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

void AbstractReaction::SetSubstrates(std::vector<AbstractChemical*> substrates)
{
    mpSubstrates = substrates;
}

void AbstractReaction::SetProducts(std::vector<AbstractChemical*> products)
{
    mpProducts = products;
}

std::vector<unsigned> AbstractReaction::GetStoichSubstrates()
{
    return mStoichSubstrates;
}

unsigned AbstractReaction::GetStoichSubstratesByIndex(unsigned index)
{
    if (index < mNumSubstrates)
    {
        return mStoichSubstrates[index];
    }
    else
    {
        std::cout << "Error: AbstractReaction::GetStoichSubstratesByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractReaction::SetStoichSubstrates(std::vector<unsigned> stoichStustrates)
{
    mStoichSubstrates = stoichStustrates;
}

std::vector<unsigned> AbstractReaction::GetStoichProducts()
{
    return mStoichProducts;
}

unsigned AbstractReaction::GetStoichProductsByIndex(unsigned index)
{
    if (index < mNumProducts)
    {
        return mStoichProducts[index];
    }
    else
    {
        std::cout << "Error: AbstractReaction::GetStoichProductsByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractReaction::SetStoichProducts(std::vector<unsigned> StoichProducts)
{
    mStoichProducts = StoichProducts;
}

void AbstractReaction::SetNumSubstrates(unsigned numSubstrates)
{
    mNumSubstrates = numSubstrates;
}

unsigned AbstractReaction::GetNumSubstrates()
{
    return mNumSubstrates;
}

void AbstractReaction::SetNumProducts(unsigned numProducts)
{
    mNumProducts = numProducts;
}

unsigned AbstractReaction::GetNumProducts()
{
    return mNumProducts;
}

void AbstractReaction::SetIsRateCheck(bool IsRateCheck)
{
    mIsRateCheck = IsRateCheck;
}

bool AbstractReaction::GetIsRateCheck()
{
    return mIsRateCheck;
}

void AbstractReaction::SetDeltaErrorRateMin(double delta_rate_min)
{
    mDeltaRateMin = delta_rate_min;
}

double AbstractReaction::GetDeltaErrorRateMin()
{
    return mDeltaRateMin;
}

void AbstractReaction::SetDeltaErrorRateMax(double delta_rate_max)
{
    mDeltaRateMax = delta_rate_max;
}

double AbstractReaction::GetDeltaErrorRateMax()
{
    return mDeltaRateMax;
}

double AbstractReaction::CheckRate(double rate)
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

void AbstractReaction::SetIrreversibleDelimiter(std::string delim)
{
    mIrreversibleDelimiter = delim;
}

std::string AbstractReaction::GetIrreversibleDelimiter()
{
    return mIrreversibleDelimiter;
}

void AbstractReaction::SetIrreversibleRateName(std::string rateName)
{
    mIrreversibleRateName = rateName;
}

std::string AbstractReaction::GetIrreversibleRateName()
{
    return mIrreversibleRateName;
}

unsigned AbstractReaction::FindIndexOfLastDelimiterPosition(std::vector<std::string> delimiterVector, std::string textString)
{
    // function to find the location of the delimiter strings in "delimiterVector" within the string textString
    // function return the index which occurs last. There may be multiple occurances of each delimiter 
    // within the textString. If all delimiters occur at npos then no delimiter occurs within the string
    // and index delimiterVector.size() is returned. 

    unsigned index=0; 
    unsigned numDelimiters = delimiterVector.size();
    std::vector<size_t> delimiterPositions(numDelimiters,0);
    bool areAllNPOS = true; 
    size_t delimiterPosition=0;

    for (unsigned i=0; i<numDelimiters; ++i)
    {
        // find the positions of the last occurance of each delimiter
        delimiterPosition = textString.rfind(delimiterVector[i]);
    
        if (delimiterPosition != std::string::npos)
        {
            areAllNPOS = false;
        }
        delimiterPositions[i] = delimiterPosition;
    }

    // if at least one of the delimiters 
    bool IsNotFound=true;
    if (!areAllNPOS)
    {
        for (unsigned i=0; i<numDelimiters; ++i)
        {
            if (delimiterPositions[i] != std::string::npos)
            {
                if (IsNotFound)
                {
                    index = i;
                    IsNotFound = false;
                }
                else if (delimiterPositions[i]>delimiterPositions[index])
                {
                    index = i;
                }
            }
        }
    }
    else
    {
        // no delimiter occurs within the textString so return the number of delimiters
        index = numDelimiters;
    }

    return index;
}

unsigned AbstractReaction::FindIndexOfFirstDelimiterPosition(std::vector<std::string> delimiterVector, std::string textString)
{
    // function to find the location of the delimiter strings in "delimiterVector" within the string textString
    // function return the index which occurs first. There may be multiple occurances of each delimiter 
    // within the textString. If all delimiters occur at npos then no delimiter occurs within the string
    // and index delimiterVector.size() is returned. 

    unsigned index=0; 
    unsigned numDelimiters = delimiterVector.size();
    std::vector<size_t> delimiterPositions(numDelimiters,0);
    bool areAllNPOS = true; 
    size_t delimiterPosition=0;
    
    for (unsigned i=0; i<numDelimiters; ++i)
    {
        // find the positions of the first occurance of each delimiter
        delimiterPosition = textString.find(delimiterVector[i]);

        if (delimiterPosition != std::string::npos)
        {
            areAllNPOS = false;
        }
        delimiterPositions[i] = delimiterPosition;
    }

    // if at least one of the delimiters 
    if (!areAllNPOS)
    {
        for (unsigned i=1; i<numDelimiters; ++i)
        {
            if (delimiterPositions[index] == std::string::npos || delimiterPositions[i]<delimiterPositions[index])
            {
                index = i;
            }
        }
    }
    else
    {
        // no delimiter occurs within the textString so return the number of delimiters
        index = numDelimiters;
    }

    return index;
}

#endif