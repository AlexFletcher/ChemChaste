#ifndef CELLANALYTICSPROPERTYFROMCELLID_HPP_
#define CELLANALYTICSPROPERTYFROMCELLID_HPP_

#include <boost/shared_ptr.hpp>

#include "CellAnalyticsProperty.hpp"

/**
 * \todo Document class.
 */
class CellAnalyticsPropertyFromCellId
{
protected:

    /** The cell ID. */
    unsigned mCellId;

    /** \todo document member variable */
    boost::shared_ptr<CellAnalyticsProperty> mpCellAnalyticsProperty;

public:

    /**
     * \todo document method
     * 
     * @param cellId
     */
    CellAnalyticsPropertyFromCellId(unsigned cellId);

     /**
      * Destructor.
      */
    ~CellAnalyticsPropertyFromCellId();

    /**
     * \todo document method
     * 
     * @param pCell
     */
    void SetUpCellAnalyticsProperty(CellPtr pCell);

    /**
     * Set mCellId.
     * 
     * @param cellID the new cell ID
     */
    void SetCellId(unsigned cellId);

    /**
     * Set mpCellAnalyticsProperty.
     * 
     * @param pProperty \todo
     */
    void SetCellAnalyticsProperty(
        boost::shared_ptr<CellAnalyticsProperty> pProperty);

    /**
     * @return mCellId
     */
    unsigned GetCellId();

    /**
     * @return mpCellAnalyticsProperty
     */
    boost::shared_ptr<CellAnalyticsProperty> GetCellAnalyticsProperty();
};

#endif /* CELLANALYTICPROPERTYFROMCELLID_HPP_ */