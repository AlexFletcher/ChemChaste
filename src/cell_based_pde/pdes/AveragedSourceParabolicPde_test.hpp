#ifndef AVERAGESOURCEPARABOLICPDE_TEST_HPP_
#define AVERAGESOURCEPARABOLICPDE_TEST_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractLinearParabolicPde.hpp"

/**
 * A parabolic PDE to be solved numerically using the finite element method, for
 * coupling to a cell-based simulation.
 *
 * The PDE takes the form
 *
 * c*du/dt = Grad.(D*Grad(u)) + k*u*rho(x),
 *
 * where the scalars c, D and k are specified by the members mDuDtCoefficient,
 * mDiffusionCoefficient and mSourceCoefficient, respectively. Their values must
 * be set in the constructor.
 *
 * The function rho(x) denotes the local density of non-apoptotic cells. This
 * quantity is computed for each element of a 'coarse' finite element mesh that is
 * passed to the method SetupSourceTerms() and stored in the member mCellDensityOnCoarseElements.
 * For a point x, rho(x) is defined to be the number of non-apoptotic cells whose
 * centres lie in each finite element containing that point, scaled by the area of
 * that element.
 */
template<unsigned DIM>
class AveragedSourceParabolicPde : public AbstractLinearParabolicPde<DIM,DIM>
{
    friend class TestCellBasedParabolicPdes;

private:

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the PDE and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractLinearParabolicPde<DIM, DIM> >(*this);
       archive & mDuDtCoefficient;
       archive & mDiffusionCoefficient;
       archive & mSourceCoefficient;
       archive & mCellDensityOnCoarseElements;
    }

protected:

    /** The cell population member. */
    AbstractCellPopulation<DIM, DIM>& mrCellPopulation;

    /** Coefficient of rate of change term.  */
    double mDuDtCoefficient;

    /** Diffusion coefficient. */
    double mDiffusionCoefficient;

    /** Coefficient of the rate of uptake of the dependent variable by non-apoptotic cells. */
    double mSourceCoefficient;

    /** Vector of averaged cell densities on elements of the coarse mesh. */
    std::vector<double> mCellDensityOnCoarseElements;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param duDtCoefficient rate of reaction (defaults to 1.0)
     * @param diffusionCoefficient rate of diffusion (defaults to 1.0)
     * @param sourceCoefficient the source term coefficient (defaults to 0.0)
     */
    AveragedSourceParabolicPde(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                               double duDtCoefficient=1.0,
                               double diffusionCoefficient=1.0,
                               double sourceCoefficient=0.0);

    /**
     * @return const reference to the cell population (used in archiving).
     */
    const AbstractCellPopulation<DIM>& rGetCellPopulation() const;

    /**
     * Set up the source terms.
     *
     * \todo this is identical to the one in AveragedSourceEllipticPde so refactor.
     *
     * @param rCoarseMesh reference to the coarse mesh
     * @param pCellPdeElementMap optional pointer to the map from cells to coarse elements
     */
    void virtual SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=nullptr);

    /**
     * Overridden ComputeDuDtCoefficientFunction() method.
     *
     * @return the function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
     *
     * @param rX the point in space at which the function c is computed
     */
    virtual double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX);

    /**
     * Overridden ComputeSourceTerm() method.
     *
     * @return computed source term.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     */
    virtual double ComputeSourceTerm(const ChastePoint<DIM>& rX,
                                     double u,
                                     Element<DIM,DIM>* pElement=NULL);

    /**
     * Overridden ComputeSourceTermAtNode() method. That is never called.
     *
     * @return computed source term at a node.
     *
     * @param rNode the node at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the node
     */
    virtual double ComputeSourceTermAtNode(const Node<DIM>& rNode, double u);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @param rX the point in space at which the diffusion term is computed
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     *
     * @return a matrix.
     */
    virtual c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement=NULL);

    /**
     * @return the uptake rate.
     *
     * @param elementIndex the element we wish to return the uptake rate for
     */
    double GetUptakeRateForElement(unsigned elementIndex);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AveragedSourceParabolicPde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a AveragedSourceParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const AveragedSourceParabolicPde<DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM, DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a AveragedSourceParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, AveragedSourceParabolicPde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM, DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)AveragedSourceParabolicPde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*AVERAGESOURCEPARABOLICPDE_HPP_*/
