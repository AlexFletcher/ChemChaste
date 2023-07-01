#ifndef CELLANALYTICSPROPERTYFROMCELLID_HPP_
#define CELLANALYTICSPROPERTYFROMCELLID_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>

#include "CellAnalyticsProperty.hpp"
#include "ChemicalCell.hpp"
#include "AbstractTransportReactionSystemFromFile.hpp"

/**
 * \todo Document class.
 */
class CellAnalyticsPropertyFromCellId
{
protected:

    unsigned mCellId;

    boost::shared_ptr<CellAnalyticsProperty> mpCellAnalyticsProperty;

public:

    CellAnalyticsPropertyFromCellId(unsigned);

    ~CellAnalyticsPropertyFromCellId()
    {
    };

    void SetUpCellAnalyticsProperty(CellPtr);

    void SetCellId(unsigned);

    void SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty>);

    unsigned GetCellId();

    boost::shared_ptr<CellAnalyticsProperty> GetCellAnalyticsProperty();
};

#endif /* CELLANALYTICPROPERTYFROMCELLID_HPP_ */