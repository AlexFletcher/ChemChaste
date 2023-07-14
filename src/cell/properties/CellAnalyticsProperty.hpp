#ifndef CELLANALYTICSPROPERTY_HPP_
#define CELLANALYTICSPROPERTY_HPP_

#include "Cell.hpp"

/**
 * \todo Document class.
 */
class CellAnalyticsProperty : public AbstractCellProperty
{
protected:

    /** CellPtr to a given cell for access to cell data and for properties. */
    CellPtr mpCell;

    /** \todo Document member variable */
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

    /**
     * \todo Document method.
     * 
     * @param rOtherProperty
     */
    CellAnalyticsProperty(const CellAnalyticsProperty& rOtherProperty);

    /**
     * \todo Document method.
     * 
     * @param pCell
     * @param cellId
     */
    virtual void SetUp(CellPtr pCell, unsigned cellId);

    /**
     * \todo Document method.
     * 
     * @param splitRatio
     */
    virtual void PreparePostDivisionParent(double splitRatio);
    
    /**
     * \todo Document method.
     * 
     * @param rParentProperty
     * @param splitRatio
     */
    virtual void PreparePostDivisionDaughter(
        const CellAnalyticsProperty& rParentProperty,
        double splitRatio);
    
    /**
     * Set mpCell.
     * 
     * @param pCell the new pointer to the cell
     */
    void SetCellPtr(CellPtr pCell);
    
    /**
     * Set mCellId.
     * 
     * @param cellId the new cell ID
     */
    void SetCellId(unsigned cellId);

    /**
     * @return mpCell
     */
    CellPtr GetCellPtr();

    /**
     * @return mCellId
     */
    unsigned GetCellId();
};

#endif /* CELLANALYTICSPROPERTY_HPP_ */