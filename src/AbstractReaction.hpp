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
    AbstractReaction(std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                     std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                     std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                     std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                     double reactionRate = 1.0);

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

#endif /* ABSTRACTREACTION_HPP_ */