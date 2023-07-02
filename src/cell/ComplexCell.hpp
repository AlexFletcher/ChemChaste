#ifndef COMPLEXCELL_HPP_
#define COMPLEXCELL_HPP_

#include <vector>
#include <string>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "Cell.hpp"

#include "AbstractChemistry.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "CellAnalyticsProperty.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).
class AbstractSrnModel; // Circular definition (cells need to know about subcellular reaction network models and vice-versa).
class Cell;

/**
 * Version of Cell.hpp in which the cell division properties may be overridden
 */
class ComplexCell : public Cell
{
protected:
    using Cell::Divide;
    //using Cell::ReadyToDivide;

    double mSplitRatio =0.5; // proportion of parent cell volume retained, the rest goes to daughter 

    std::vector<std::string> mShareKey = {"share","Share","shared","Shared","split"};

    std::vector<std::string> mChemicalNames;

    std::vector<std::string> mChemicalDivsionRules;

public:

    ComplexCell(boost::shared_ptr<AbstractCellProperty> pMutationState,
         AbstractCellCycleModel* pCellCycleModel,
         AbstractSrnModel* pSrnModel=nullptr,
         bool archiving=false,
         CellPropertyCollection cellPropertyCollection=CellPropertyCollection());

    virtual ~ComplexCell()
    {
    };

    virtual CellPtr Divide();

    virtual void DetermineSplitRatio();

    virtual double SplitParentCellData(double);

    bool IsChemicalShared(std::string);

    double GetSplitRatio();

    void SetSplitRatio(double);

    std::vector<std::string> GetShareKey();

    void SetShareKey(std::vector<std::string>);

    std::vector<std::string> GetChemicalNames();

    void SetChemicalNames(std::vector<std::string>);

    std::vector<std::string> GetChemicalDivsionRules();

    void SetChemicalDivsionRules(std::vector<std::string>);
};

#endif /* COMPLEXCELL_HPP_ */