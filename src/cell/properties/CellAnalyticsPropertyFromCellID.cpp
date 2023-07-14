#include "CellAnalyticsPropertyFromCellId.hpp"

CellAnalyticsPropertyFromCellId::CellAnalyticsPropertyFromCellId(unsigned cellId)
    : mCellId(cellId)
{
    SetCellId(cellId);
}

CellAnalyticsPropertyFromCellId::~CellAnalyticsPropertyFromCellId()
{
}

void CellAnalyticsPropertyFromCellId::SetUpCellAnalyticsProperty(CellPtr pCell)
{
    auto p_analytics = boost::static_pointer_cast<CellAnalyticsProperty>(pCell->rGetCellPropertyCollection().GetPropertiesType<CellAnalyticsProperty>().GetProperty());
    p_analytics->SetUp(pCell, mCellId);
    SetCellAnalyticsProperty(p_analytics);
}

void CellAnalyticsPropertyFromCellId::SetCellId(unsigned cellId)
{
    mCellId = cellId;
}

void CellAnalyticsPropertyFromCellId::SetCellAnalyticsProperty(
    boost::shared_ptr<CellAnalyticsProperty> pProperty)
{
    mpCellAnalyticsProperty = pProperty;
}

unsigned CellAnalyticsPropertyFromCellId::GetCellId()
{
    return mCellId;
}

boost::shared_ptr<CellAnalyticsProperty> CellAnalyticsPropertyFromCellId::GetCellAnalyticsProperty()
{
    return mpCellAnalyticsProperty;
}