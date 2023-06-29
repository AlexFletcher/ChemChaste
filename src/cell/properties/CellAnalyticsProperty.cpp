#include "CellAnalyticsProperty.hpp"

CellAnalyticsProperty::CellAnalyticsProperty()
    : AbstractCellProperty()
{
}

CellAnalyticsProperty::~CellAnalyticsProperty()
{
}

CellAnalyticsProperty::CellAnalyticsProperty(const CellAnalyticsProperty& existingProperty)
{
    mCellID = existingProperty.mCellID;
}

void CellAnalyticsProperty::SetUp(CellPtr pCell,unsigned cellID)
{
    SetCellID(cellID);
    SetCellPtr(pCell);
}

void CellAnalyticsProperty::PreparePostDivisionParent(double splitRatio)
{
}
    
void CellAnalyticsProperty::PreparePostDivisionDaughter(const CellAnalyticsProperty& parentProperty, double splitRatio)
{
}

void CellAnalyticsProperty::SetCellPtr(CellPtr pCell)
{
    mpCell = pCell;
}

void CellAnalyticsProperty::SetCellID(unsigned cellID)
{
    mCellID = cellID;
}

CellPtr CellAnalyticsProperty::GetCellPtr()
{
    return mpCell;
}

unsigned CellAnalyticsProperty::GetCellID()
{
    return mCellID;
}