#include "AveragedSourceParabolicPde.hpp"
#include "ApoptoticCellProperty.hpp"

template<unsigned DIM>
AveragedSourceParabolicPde<DIM>::AveragedSourceParabolicPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                            double duDtCoefficient,
                                                            double diffusionCoefficient,
                                                            double sourceCoefficient)
    : mrCellPopulation(rCellPopulation),
      mDuDtCoefficient(duDtCoefficient),
      mDiffusionCoefficient(diffusionCoefficient),
      mSourceCoefficient(sourceCoefficient)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& AveragedSourceParabolicPde<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
void AveragedSourceParabolicPde<DIM>::SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap) // must be called before solve
{
    // Allocate memory
    mCellDensityOnCoarseElements.resize(rCoarseMesh.GetNumElements());
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        mCellDensityOnCoarseElements[elem_index] = 0.0;
    }

    // Loop over cells, find which coarse element it is in, and add 1 to mSourceTermOnCoarseElements[elem_index]
    for (auto cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        unsigned elem_index = 0;
        const ChastePoint<DIM>& r_position_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

        if (pCellPdeElementMap != nullptr)
        {
            elem_index = (*pCellPdeElementMap)[*cell_iter];
        }
        else
        {
            elem_index = rCoarseMesh.GetContainingElementIndex(r_position_of_cell);
        }

        // Update element map if cell has moved
        bool cell_is_apoptotic = cell_iter->template HasCellProperty<ApoptoticCellProperty>();

        if (!cell_is_apoptotic)
        {
            mCellDensityOnCoarseElements[elem_index] += 1.0;
        }
    }

    // Then divide each entry of mSourceTermOnCoarseElements by the element's area
    c_matrix<double, DIM, DIM> jacobian;
    double det;
    for (unsigned elem_index=0; elem_index<mCellDensityOnCoarseElements.size(); elem_index++)
    {
        rCoarseMesh.GetElement(elem_index)->CalculateJacobian(jacobian, det);
        mCellDensityOnCoarseElements[elem_index] /= rCoarseMesh.GetElement(elem_index)->GetVolume(det);
    }
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return mDuDtCoefficient;
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    assert(!mCellDensityOnCoarseElements.empty());
    double coefficient = mSourceCoefficient * mCellDensityOnCoarseElements[pElement->GetIndex()];

    // The source term is C*u
    return coefficient*u;
}

// LCOV_EXCL_START
template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    NEVER_REACHED;
    return 0.0;
}
// LCOV_EXCL_STOP

template<unsigned DIM>
c_matrix<double,DIM,DIM> AveragedSourceParabolicPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{ 
    return mDiffusionCoefficient*identity_matrix<double>(DIM);
}

template<unsigned DIM>
double AveragedSourceParabolicPde<DIM>::GetUptakeRateForElement(unsigned elementIndex)
{
    return this->mCellDensityOnCoarseElements[elementIndex];
}

// Explicit instantiation
template class AveragedSourceParabolicPde<1>;
template class AveragedSourceParabolicPde<2>;
template class AveragedSourceParabolicPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AveragedSourceParabolicPde)
