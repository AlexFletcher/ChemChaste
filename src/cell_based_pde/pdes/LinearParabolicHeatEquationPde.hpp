#ifndef LINEARPARABOLICHEATPDE_HPP_
#define LINEARPARABOLICHEATPDE_HPP_

#include "AbstractLinearParabolicPdeSystemForCoupledOdeSystem.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"

/**
 * The heat diffusion equation with constant source term.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class LinearParabolicHeatEquationPde : public AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    // diffusion rate
    std::vector<double> mDiffusionRates;      

public:

    LinearParabolicHeatEquationPde(std::vector<double> diffusionRates)
        : AbstractLinearParabolicPdeSystemForCoupledOdeSystem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(),
          mDiffusionRates(diffusionRates)
    {
    }

    double ComputeDuDtCoefficientFunction(const ChastePoint<SPACE_DIM>& rX, unsigned index)
    {
        return 1.0;
    }

    double ComputeSourceTerm(const ChastePoint<SPACE_DIM>& rX, c_vector<double,PROBLEM_DIM>& rU, std::vector<double>& rOdeSolution, unsigned pdeIndex)
    {
        return 0.0;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>& rX, unsigned pdeIndex, Element<ELEMENT_DIM,SPACE_DIM>* pElement=NULL)
    {
        assert(pdeIndex<PROBLEM_DIM);

        c_matrix<double, SPACE_DIM, SPACE_DIM> diffusion_term;

        diffusion_term = mDiffusionRates[pdeIndex]*identity_matrix<double>(SPACE_DIM);
        
        return diffusion_term;
    }
};

#endif /* LINEARPARABOLICHEATPDE_HPP_ */