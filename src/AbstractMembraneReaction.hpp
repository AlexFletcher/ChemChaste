#ifndef ABSTRACTMEMBRANEREACTION_HPP_
#define ABSTRACTMEMBRANEREACTION_HPP_

#include <string>
#include <stdlib.h> 
#include <vector>
#include <iostream>

#include "AbstractChemistry.hpp"
#include "AbstractChemical.hpp"

/**
 * Abstract property to contain information about the interactions of chemical 
 * species at a cell boundary. Define two coupled reactions; first the reaction 
 * in the bulk, second the reaction on the cell side: bulk->bulk | cell->cell
 */
class AbstractMembraneReaction 
{
protected:
    
    std::vector<AbstractChemical*> mpBulkSubstrates;

    std::vector<AbstractChemical*> mpBulkProducts;

    std::vector<AbstractChemical*> mpCellSubstrates;

    std::vector<AbstractChemical*> mpCellProducts;

    unsigned mNumBulkProducts;

    unsigned mNumBulkSubstrates;

    unsigned mNumCellProducts;

    unsigned mNumCellSubstrates;

    // vector containing the stoichmetry of the bulk substrates.  Will be "-ve" in reaction
    std::vector<unsigned> mStoichBulkSubstrates;

    // vector containing the stoichmetry of the bulk products.  Will be "+ve" in reaction
    std::vector<unsigned> mStoichBulkProducts;

    // vector containing the stoichmetry of the cell substrates.  Will be "-ve" in reaction
    std::vector<unsigned> mStoichCellSubstrates;

    // vector containing the stoichmetry of the cell products.  Will be "+ve" in reaction
    std::vector<unsigned> mStoichCellProducts;


    // container for the reaction rate constant.  '+ve' value. May be overriden in derived classes to account for functional rates
    double mReactionRate;

    bool mIsRateCheck;

    double mDeltaRateMin;

    double mDeltaRateMax;

    std::string mIrreversibleDelimiter = "->";
  
    std::string mIrreversibleRateName = "kf =";

public:

    // constructor
    AbstractMembraneReaction(std::vector<AbstractChemical*> bulkSubstrates = std::vector<AbstractChemical*>(),
                             std::vector<AbstractChemical*> bulkProducts = std::vector<AbstractChemical*>(),
                             std::vector<AbstractChemical*> cellSubstrates = std::vector<AbstractChemical*>(),
                             std::vector<AbstractChemical*> cellProducts = std::vector<AbstractChemical*>(),
                             std::vector<unsigned> stoichBulkSubstrates = std::vector<unsigned>(),
                             std::vector<unsigned> stoichBulkProducts = std::vector<unsigned>(),
                             std::vector<unsigned> stoichCellSubstrates = std::vector<unsigned>(),
                             std::vector<unsigned> stoichCellProducts = std::vector<unsigned>(),
                             double reactionRate = 1.0);

    // copy constructor
    AbstractMembraneReaction(const AbstractMembraneReaction&);

    // destructor
    virtual ~AbstractMembraneReaction()
    {
    };

    // function to take in pointer to current concentration state vector of the state vector for change in concentration.  Basic update via multiplying by constant reaction rate
    virtual void React(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    // pure vitual function to update member variables and conditions of derived classes
    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, AbstractChemistry*, const std::vector<double>&, const std::vector<double>&);

    virtual void SetReactionRate(double);

    virtual double GetReactionRate();

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    // Chemical handeling functions
    std::vector<AbstractChemical*> GetBulkSubstrates();

    AbstractChemical* GetBulkSubstratesByIndex(unsigned);

    std::vector<AbstractChemical*> GetCellSubstrates();

    AbstractChemical* GetCellSubstratesByIndex(unsigned);

    std::vector<AbstractChemical*> GetBulkProducts();

    AbstractChemical* GetBulkProductsByIndex(unsigned);

    std::vector<AbstractChemical*> GetCellProducts();

    AbstractChemical* GetCellProductsByIndex(unsigned);

    void SetBulkSubstrates(std::vector<AbstractChemical*>);

    void SetCellSubstrates(std::vector<AbstractChemical*>);

    void SetBulkProducts(std::vector<AbstractChemical*>);

    void SetCellProducts(std::vector<AbstractChemical*>);

    std::vector<unsigned> GetBulkStoichSubstrates();

    unsigned GetBulkStoichSubstratesByIndex(unsigned);

    void SetBulkStoichSubstrates(std::vector<unsigned>);

    std::vector<unsigned> GetCellStoichSubstrates();

    unsigned GetCellStoichSubstratesByIndex(unsigned);

    void SetCellStoichSubstrates(std::vector<unsigned>);

    std::vector<unsigned> GetBulkStoichProducts();

    unsigned GetBulkStoichProductsByIndex(unsigned);

    void SetBulkStoichProducts(std::vector<unsigned>);

    std::vector<unsigned> GetCellStoichProducts();

    unsigned GetCellStoichProductsByIndex(unsigned);

    void SetCellStoichProducts(std::vector<unsigned>);

    void SetNumBulkSubstrates(unsigned);

    unsigned GetNumBulkSubstrates();

    void SetNumCellSubstrates(unsigned);

    unsigned GetNumCellSubstrates();

    void SetNumBulkProducts(unsigned);

    unsigned GetNumBulkProducts();

    void SetNumCellProducts(unsigned);

    unsigned GetNumCellProducts();

    // reaction concentration checking functions
    void SetIsRateCheck(bool);

    bool GetIsRateCheck();

    void SetDeltaErrorRateMin(double);

    double GetDeltaErrorRateMin();

    void SetDeltaErrorRateMax(double);

    double GetDeltaErrorRateMax();

    double CheckRate(double);

    // file read functions
    void SetIrreversibleDelimiter(std::string);

    std::string GetIrreversibleDelimiter();

    void SetIrreversibleRateName(std::string);

    std::string GetIrreversibleRateName();
};

#endif /* ABSTRACTMEMBRANEREACTION_HPP_ */