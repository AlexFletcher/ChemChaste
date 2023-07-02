#include "AbstractMembraneReaction.hpp"

AbstractMembraneReaction::AbstractMembraneReaction(   
    std::vector<AbstractChemical*> bulkSubstrates,
    std::vector<AbstractChemical*> bulkProducts,
    std::vector<AbstractChemical*> cellSubstrates,
    std::vector<AbstractChemical*> cellProducts,
    std::vector<unsigned> stoichBulkSubstrates,
    std::vector<unsigned> stoichBulkProducts,
    std::vector<unsigned> stoichCellSubstrates,
    std::vector<unsigned> stoichCellProducts,
    double reactionRate)
    : mpBulkSubstrates(bulkSubstrates),
      mpBulkProducts(bulkProducts),
      mpCellSubstrates(cellSubstrates),
      mpCellProducts(cellProducts),
      mStoichBulkSubstrates(stoichBulkSubstrates),
      mStoichBulkProducts(stoichBulkProducts),
      mStoichCellSubstrates(stoichCellSubstrates),
      mStoichCellProducts(stoichCellProducts),
      mReactionRate(reactionRate)
{
    mNumBulkProducts = bulkProducts.size();
    mNumBulkSubstrates = bulkSubstrates.size();
    mNumCellProducts = cellProducts.size();
    mNumCellSubstrates = cellSubstrates.size();

    SetIsRateCheck(true);
    SetDeltaErrorRateMin(1e-8);
    SetDeltaErrorRateMax(1e8);
}

AbstractMembraneReaction::AbstractMembraneReaction(const AbstractMembraneReaction& rOtherReaction)
{
    mpBulkSubstrates = rOtherReaction.mpBulkSubstrates;
    mpBulkProducts = rOtherReaction.mpBulkProducts;
    mpCellSubstrates = rOtherReaction.mpCellSubstrates;
    mpCellProducts = rOtherReaction.mpCellProducts;
    mNumBulkProducts = rOtherReaction.mNumBulkProducts;
    mNumBulkSubstrates = rOtherReaction.mNumBulkSubstrates;
    mNumCellProducts = rOtherReaction.mNumCellProducts;
    mNumCellSubstrates = rOtherReaction.mNumCellSubstrates;
    mStoichBulkSubstrates = rOtherReaction.mStoichBulkSubstrates;
    mStoichBulkProducts = rOtherReaction.mStoichBulkProducts;
    mStoichCellSubstrates = rOtherReaction.mStoichCellSubstrates;
    mStoichCellProducts = rOtherReaction.mStoichCellProducts;
    mReactionRate = rOtherReaction.mReactionRate;
    mIsRateCheck = rOtherReaction.mIsRateCheck;
    mDeltaRateMin = rOtherReaction.mDeltaRateMin;
    mDeltaRateMax = rOtherReaction.mDeltaRateMax;
}

void AbstractMembraneReaction::React(AbstractChemistry* pBulkChemistry, 
                                     AbstractChemistry* pCellChemistry, 
                                     const std::vector<double>& rCurrentBulkConcentration, 
                                     const std::vector<double>& rCurrentCellConcentration, 
                                     std::vector<double>& rChangeBulkConc, 
                                     std::vector<double>& rChangeCellConc)
{
    std::vector<AbstractChemical*> p_bulk_chemical_vector = pBulkChemistry->rGetChemicalVector();  
    std::vector<AbstractChemical*> p_cell_chemical_vector = pCellChemistry->rGetChemicalVector();

    UpdateReactionRate(pBulkChemistry, pCellChemistry, rCurrentBulkConcentration, rCurrentCellConcentration);
    
    // Perform the reaction

    // Run through the bulk species
    unsigned index = 0;
    for (std::vector<AbstractChemical*>::iterator chem_iter = p_bulk_chemical_vector.begin();
         chem_iter != p_bulk_chemical_vector.end();
         ++chem_iter, ++index)
    {
        AbstractChemical* p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // For each bulk chemical, parse whether it is involved in this reaction.
        for (unsigned j = 0; j < mNumBulkSubstrates; ++j)
        {
            if (mpBulkSubstrates[j]->GetChemicalName() == p_system_chemical->GetChemicalName())
            {
                rChangeBulkConc[index] -= mStoichBulkSubstrates[j]*GetReactionRate();
                break;
            }
        }

        // A reactant may be present on both sides of the reaction and may convert at different functional rates
        for (unsigned j = 0; j < mNumBulkProducts; ++j)
        {
            if (mpBulkProducts[j]->GetChemicalName() == p_system_chemical->GetChemicalName())
            {
                rChangeBulkConc[index] += mStoichBulkProducts[j]*GetReactionRate();
                break;
            }
        }
    }

    // Run through the cell species
    index = 0;
    for (std::vector<AbstractChemical*>::iterator chem_iter = p_cell_chemical_vector.begin();
         chem_iter != p_cell_chemical_vector.end();
         ++chem_iter, ++index)
    {
        AbstractChemical* p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        // For each bulk chemical, parse whether it is involved in this reaction.
        for (unsigned j = 0; j < mNumCellSubstrates; ++j)
        {
            if (mpCellSubstrates[j]->GetChemicalName() == p_system_chemical->GetChemicalName())
            {
                rChangeCellConc[index] -= mStoichCellSubstrates[j]*GetReactionRate();
                break;
            }
        }

        // A reactant may be present on both sides of the reaction and may convert at different functional rates
        for (unsigned j = 0; j < mNumCellProducts; ++j)
        {
            if (mpCellProducts[j]->GetChemicalName() == p_system_chemical->GetChemicalName())
            {
                rChangeCellConc[index] += mStoichCellProducts[j]*GetReactionRate();
                break;
            }
        }
    }
}

// pure vitual function to update member variables and conditions of derived classes
void AbstractMembraneReaction::UpdateReaction()
{
    return;
}

void AbstractMembraneReaction::UpdateReactionRate(AbstractChemistry* pBulkChemistry, AbstractChemistry* pCellChemistry, const std::vector<double>& currentBulkConc, const std::vector<double>& currentCellConc)
{
    // default to zeroth order reaction
    SetReactionRate(mReactionRate);
}

void AbstractMembraneReaction::SetReactionRate(double reactionRate)
{
    if (mIsRateCheck)
    {
        reactionRate = CheckRate(reactionRate);
    }
    
    mReactionRate = reactionRate;
}

double AbstractMembraneReaction::GetReactionRate()
{
    return mReactionRate;
}

std::string AbstractMembraneReaction::GetReactionType()
{
    // virtual function to be overriden in derived classes, used to idnetify reaction inheritance 
    return "ZerothOrderCoupledMembrane";
}

void AbstractMembraneReaction::ParseReactionInformation(std::string reaction_information, bool isReversible=false)
{
    if (reaction_information.find(mIrreversibleRateName) != std::string::npos)
    {
        if (reaction_information.find(";") != std::string::npos)
        {
            size_t pos= reaction_information.find(";");
            reaction_information = reaction_information.substr(0,pos);
        }
        size_t pos= reaction_information.find(mIrreversibleRateName);

        SetReactionRate(atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1,std::string::npos).c_str()));
    }
}

// Chemical handeling functions
std::vector<AbstractChemical*> AbstractMembraneReaction::GetBulkSubstrates()
{
    return mpBulkSubstrates;
}

AbstractChemical* AbstractMembraneReaction::GetBulkSubstratesByIndex(unsigned index)
{
    if (index < mNumBulkSubstrates)
    {
        return mpBulkSubstrates[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetBulkSubstratesByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractMembraneReaction::GetCellSubstrates()
{
    return mpCellSubstrates;
}

AbstractChemical* AbstractMembraneReaction::GetCellSubstratesByIndex(unsigned index)
{
    if (index < mNumCellSubstrates)
    {
        return mpCellSubstrates[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetCellSubstratesByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractMembraneReaction::GetBulkProducts()
{
    return mpBulkProducts;
}

AbstractChemical* AbstractMembraneReaction::GetBulkProductsByIndex(unsigned index)
{
    if (index < mNumBulkProducts)
    {
        return mpBulkProducts[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetBulkProductsByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

std::vector<AbstractChemical*> AbstractMembraneReaction::GetCellProducts()
{
    return mpCellProducts;
}

AbstractChemical* AbstractMembraneReaction::GetCellProductsByIndex(unsigned index)
{
    if (index < mNumCellProducts)
    {
        return mpCellProducts[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetCellProductsByIndex(unsigned index), index out of bounds" << std::endl;
        return new AbstractChemical();
    } 
}

void AbstractMembraneReaction::SetBulkSubstrates(std::vector<AbstractChemical*> substrates)
{
    mpBulkSubstrates = substrates;
}

void AbstractMembraneReaction::SetCellSubstrates(std::vector<AbstractChemical*> substrates)
{
    mpCellSubstrates = substrates;
}

void AbstractMembraneReaction::SetBulkProducts(std::vector<AbstractChemical*> products)
{
    mpBulkProducts = products;
}

void AbstractMembraneReaction::SetCellProducts(std::vector<AbstractChemical*> products)
{
    mpCellProducts = products;
}

std::vector<unsigned> AbstractMembraneReaction::GetBulkStoichSubstrates()
{
    return mStoichBulkSubstrates;
}

unsigned AbstractMembraneReaction::GetBulkStoichSubstratesByIndex(unsigned index)
{
    if (index < mNumBulkSubstrates)
    {
        return mStoichBulkSubstrates[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetBulkStoichSubstratesByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetBulkStoichSubstrates(std::vector<unsigned> stoich)
{
    mStoichBulkSubstrates = stoich;
}

std::vector<unsigned> AbstractMembraneReaction::GetCellStoichSubstrates()
{
    return mStoichCellSubstrates;
}

unsigned AbstractMembraneReaction::GetCellStoichSubstratesByIndex(unsigned index)
{
    if (index < mNumCellSubstrates)
    {
        return mStoichCellSubstrates[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetCellStoichSubstratesByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetCellStoichSubstrates(std::vector<unsigned> stoich)
{
    mStoichCellSubstrates = stoich;
}

std::vector<unsigned> AbstractMembraneReaction::GetBulkStoichProducts()
{
    return mStoichBulkProducts;
}

unsigned AbstractMembraneReaction::GetBulkStoichProductsByIndex(unsigned index)
{
    if (index < mNumBulkProducts)
    {
        return mStoichBulkProducts[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetBulkStoichProductsByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetBulkStoichProducts(std::vector<unsigned> stoich)
{
    mStoichBulkProducts = stoich;
}

std::vector<unsigned> AbstractMembraneReaction::GetCellStoichProducts()
{
    return mStoichCellProducts;
}

unsigned AbstractMembraneReaction::GetCellStoichProductsByIndex(unsigned index)
{
    if (index < mNumCellProducts)
    {
        return mStoichCellProducts[index];
    }
    else
    {
        std::cout << "Error: AbstractMembraneReaction::GetCellStoichProductsByIndex(unsigned index), index out of bounds" << std::endl;
        return 0;
    } 
}

void AbstractMembraneReaction::SetCellStoichProducts(std::vector<unsigned> stoich)
{
    mStoichCellProducts = stoich;
}

void AbstractMembraneReaction::SetNumBulkSubstrates(unsigned numberSubstrates)
{
    mNumBulkSubstrates = numberSubstrates;
}

unsigned AbstractMembraneReaction::GetNumBulkSubstrates()
{
    return mNumBulkSubstrates;
}

void AbstractMembraneReaction::SetNumCellSubstrates(unsigned numberSubstrates)
{
    mNumCellSubstrates = numberSubstrates;
}

unsigned AbstractMembraneReaction::GetNumCellSubstrates()
{
    return mNumCellSubstrates;
}

void AbstractMembraneReaction::SetNumBulkProducts(unsigned numberProducts)
{
    mNumBulkProducts = numberProducts;
}

unsigned AbstractMembraneReaction::GetNumBulkProducts()
{
    return mNumBulkProducts;
}

void AbstractMembraneReaction::SetNumCellProducts(unsigned numberProducts)
{
    mNumCellProducts = numberProducts;
}

unsigned AbstractMembraneReaction::GetNumCellProducts()
{
    return mNumCellProducts;
}

void AbstractMembraneReaction::SetIsRateCheck(bool isRateCheck)
{
    mIsRateCheck = isRateCheck;
}

bool AbstractMembraneReaction::GetIsRateCheck()
{
    return mIsRateCheck;
}

void AbstractMembraneReaction::SetDeltaErrorRateMin(double delta_rate_min)
{
    mDeltaRateMin = delta_rate_min;
}

double AbstractMembraneReaction::GetDeltaErrorRateMin()
{
    return mDeltaRateMin;
}

void AbstractMembraneReaction::SetDeltaErrorRateMax(double delta_rate_max)
{
    mDeltaRateMax = delta_rate_max;
}

double AbstractMembraneReaction::GetDeltaErrorRateMax()
{
    return mDeltaRateMax;
}

double AbstractMembraneReaction::CheckRate(double rate)
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

void AbstractMembraneReaction::SetIrreversibleDelimiter(std::string delim)
{
    mIrreversibleDelimiter = delim;
}

std::string AbstractMembraneReaction::GetIrreversibleDelimiter()
{
    return mIrreversibleDelimiter;
}

void AbstractMembraneReaction::SetIrreversibleRateName(std::string rateName)
{
    mIrreversibleRateName = rateName;
}

std::string AbstractMembraneReaction::GetIrreversibleRateName()
{
    return mIrreversibleRateName;
}

#endif