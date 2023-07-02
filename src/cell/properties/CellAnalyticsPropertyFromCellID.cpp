#include "CellAnalyticsPropertyFromCellId.hpp"

CellAnalyticsPropertyFromCellId::CellAnalyticsPropertyFromCellId(unsigned cellId)
    : mCellId(cellId)
{
    SetCellId(cellId);
}

void CellAnalyticsPropertyFromCellId::SetUpCellAnalyticsProperty(CellPtr p_cell)
{
    auto p_analytics = boost::static_pointer_cast<CellAnalyticsProperty>(p_cell->rGetCellPropertyCollection().GetPropertiesType<CellAnalyticsProperty>().GetProperty());
    p_analytics->SetUp(p_cell, mCellId);
    SetCellAnalyticsProperty(p_analytics);
}

void CellAnalyticsPropertyFromCellId::SetCellId(unsigned cellId)
{
    mCellId = cellId;
}

void CellAnalyticsPropertyFromCellId::SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty> pCellAnalyticsProperty)
{
    mpCellAnalyticsProperty = pCellAnalyticsProperty;
}

unsigned CellAnalyticsPropertyFromCellId::GetCellId()
{
    return mCellId;
}

boost::shared_ptr<CellAnalyticsProperty> CellAnalyticsPropertyFromCellId::GetCellAnalyticsProperty()
{
    return mpCellAnalyticsProperty;
}