#include "CellAnalyticsProperty.hpp"

CellAnalyticsProperty::CellAnalyticsProperty()
    : AbstractCellProperty()
{
}

CellAnalyticsProperty::~CellAnalyticsProperty()
{
}

CellAnalyticsProperty::CellAnalyticsProperty(const CellAnalyticsProperty& rOtherProperty)
{
    mCellId = rOtherProperty.mCellId;
}

void CellAnalyticsProperty::SetUp(CellPtr pCell, unsigned cellId)
{
    SetCellId(cellId);
    SetCellPtr(pCell);
}

void CellAnalyticsProperty::PreparePostDivisionParent(double splitRatio)
{
}
    
void CellAnalyticsProperty::PreparePostDivisionDaughter(const CellAnalyticsProperty& rParentProperty, double splitRatio)
{
}

void CellAnalyticsProperty::SetCellPtr(CellPtr pCell)
{
    mpCell = pCell;
}

void CellAnalyticsProperty::SetCellId(unsigned cellId)
{
    mCellId = cellId;
}

CellPtr CellAnalyticsProperty::GetCellPtr()
{
    return mpCell;
}

unsigned CellAnalyticsProperty::GetCellId()
{
    return mCellId;
}