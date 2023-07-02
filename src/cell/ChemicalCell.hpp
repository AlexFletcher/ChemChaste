#ifndef CHEMICALCELL_HPP_
#define CHEMICALCELL_HPP_

#include "Cell.hpp"
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellMutationState.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "CellData.hpp"
#include "CellVecData.hpp"
#include "Cell.hpp"

#include "AbstractCellCycleModel.hpp"
#include "AbstractSrnModel.hpp"
#include "CellPropertyCollection.hpp"

#include "AbstractChemistry.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "CellAnalyticsProperty.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).
class AbstractSrnModel; // Circular definition (cells need to know about subcellular reaction network models and vice-versa).
class Cell;

/**
 * Subclass of Cell in which the cell division properties may be overridden.
 */
class ChemicalCell : public Cell
{
protected:

    using Cell::Divide;

    double mSplitRatio = 0.5; // proportion of parent cell volume retained, rest goes to daughter 

public:

    ChemicalCell(boost::shared_ptr<AbstractCellProperty> pMutationState,
         AbstractCellCycleModel* pCellCycleModel,
         AbstractSrnModel* pSrnModel=nullptr,
         bool archiving=false,
         CellPropertyCollection cellPropertyCollection=CellPropertyCollection());

    virtual ~ChemicalCell()
    {
    };

    virtual CellPtr Divide();

    virtual void DetermineSplitRatio();

    virtual double SplitParentCellData(double);

    double GetSplitRatio();

    void SetSplitRatio(double);
};

#endif /* CHEMICALCELL_HPP_ */