#ifndef CELLANALYTICSPROPERTY_HPP_
#define CELLANALYTICSPROPERTY_HPP_

#include "Cell.hpp"

/**
 * \todo Document class.
 */
class CellAnalyticsProperty : public AbstractCellProperty
{
protected:

    // CellPtr to a given cell for access to cell data and for properties
    CellPtr mpCell;

    unsigned mCellId;

public:

    /**
     * Constructor.
     */
    CellAnalyticsProperty();

    /**
     * Destructor.
     */
    virtual ~CellAnalyticsProperty();

    CellAnalyticsProperty(const CellAnalyticsProperty& rOtherProperty);

    virtual void SetUp(CellPtr pCell, unsigned cellId);

    virtual void PreparePostDivisionParent(double splitRatio);
    
    virtual void PreparePostDivisionDaughter(const CellAnalyticsProperty& rParentProperty, double splitRatio);

    void SetCellPtr(CellPtr pCell);

    void SetCellId(unsigned cellId);

    CellPtr GetCellPtr();

    unsigned GetCellId();
};

#endif /* CELLANALYTICSPROPERTY_HPP_ */