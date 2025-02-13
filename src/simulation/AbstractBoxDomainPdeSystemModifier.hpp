#ifndef ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_
#define ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_

#include "AbstractPdeSystemModifier.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"
#include "SmartPointers.hpp"
#include "AbstractDomainField.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "ExtendedCellProperty.hpp"

/**
 * \todo Document class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractBoxDomainPdeSystemModifier : public AbstractPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
protected:

    /** Map between cells and the elements of the FE mesh containing them. */
    std::map<CellPtr, unsigned> mCellPdeElementMap;

    /**
     * Pointer to a ChasteCuboid storing the outer boundary for the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    boost::shared_ptr<ChasteCuboid<SPACE_DIM> > mpMeshCuboid;

    /**
     * The step size to be used in the FE mesh.
     * Stored as a member to facilitate archiving.
     */
    double mStepSize;

    /**
     * Whether to set the boundary condition on the edge of the box domain 
     * rather than the cell population. Defaults to true.
     */
    bool mSetBcsOnBoxBoundary;
    
public:

    /**
     * Default constructor.
     * 
     * @param pDomainField \todo document param
     * @param pMeshCuboid A shared pointer to a ChasteCuboid specifying the 
     *                    outer boundary for the FE mesh (defaults to NULL)
     * @param stepSize step size to be used in the FE mesh (defaults to 1.0, 
     *                 i.e. the default cell size)
     * @param solution solution vector (defaults to NULL)
    */
    AbstractBoxDomainPdeSystemModifier(ChemicalDomainFieldForCellCoupling<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pDomainField,
                                       boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<SPACE_DIM> >(),
                                       double stepSize=1.0,
                                       Vec solution=nullptr);

    /**
     * Destructor.
    */
    virtual ~AbstractBoxDomainPdeSystemModifier();

    /**
     * @return mStepSize.
     */
    double GetStepSize();

    /**
     * Set mSetBcsOnCoarseBoundary.
     *
     * @param setBcsOnBoxBoundary whether to set the boundary condition on the 
     * edge of the box domain rather than the cell population
     */
    void SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary);

    /**
     * @return mSetBcsOnCoarseBoundary.
     */
    bool AreBcsSetOnBoxBoundary();

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * Here we just initialize the cell PDE element map.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste 
     *                        output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, 
                            std::string outputDirectory);

    /**
     * Helper method to generate the mesh.
     *
     * @param pMeshCuboid the outer boundary for the FE mesh.
     * @param stepSize the step size to be used in the FE mesh.
     */
    void GenerateFeMesh(boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid, 
                        double stepSize);

    /**
     * Helper method to copy the PDE solution to CellData.
     *
     * Here we need to interpolate from the FE mesh onto the cells.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Initialise mCellPdeElementMap.
     *
     * @param rCellPopulation reference to the cell population
     */
    void InitialiseCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Update the mCellPdeElementMap.
     *
     * This method should be called before sending the element map to a PDE 
     * class to ensure the map is up to date.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractBoxDomainPdeSystemModifier(
    ChemicalDomainFieldForCellCoupling<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* pDomainField,
    boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid,
    double stepSize,
    Vec solution) 
    : AbstractPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pDomainField, solution), 
      mpMeshCuboid(pMeshCuboid),
      mStepSize(stepSize),
      mSetBcsOnBoxBoundary(true)
{
    if (pMeshCuboid)
    {
        // We only need to generate mpFeMesh once, as it does not vary with time
        this->GenerateFeMesh(mpMeshCuboid, mStepSize);
        this->mDeleteFeMesh = true;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~AbstractBoxDomainPdeSystemModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetStepSize()
{
     return mStepSize;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary)
{
    mSetBcsOnBoxBoundary = setBcsOnBoxBoundary;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AreBcsSetOnBoxBoundary()
{
    return mSetBcsOnBoxBoundary;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetupSolve(rCellPopulation, outputDirectory);
    InitialiseCellPdeElementMap(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GenerateFeMesh(boost::shared_ptr<ChasteCuboid<SPACE_DIM> > pMeshCuboid, double stepSize)
{
    // Generate mesh and PDE system from the Abstract domain and the mesh cuboid

    // Create a regular coarse tetrahedral mesh

    // Get centroid of meshCuboid
    ChastePoint<SPACE_DIM> upper = pMeshCuboid->rGetUpperCorner();
    ChastePoint<SPACE_DIM> lower = pMeshCuboid->rGetLowerCorner();
    c_vector<double, SPACE_DIM> centre_of_cuboid = 0.5*(upper.rGetLocation() + lower.rGetLocation());

    this->mpCoupledDomainField->GenerateFeMesh();
    this->mpCoupledDomainField->SetUpDomainFromFiles();
    
    this->mpFeMesh = this->mpCoupledDomainField->rGetDomainFeMesh();

    // Set the lower point of the cuboid as the origin of the mesh
    std::vector<double> origin(SPACE_DIM, 0.0);
    c_vector<double, SPACE_DIM> offset = zero_vector<double>(SPACE_DIM);
    for (unsigned dim = 0; dim < SPACE_DIM; dim++)
    {
        offset[dim] = lower[dim];
        origin[dim] = lower[dim];
    }

    // Now move the mesh to the correct location
    this->mpFeMesh->Translate(offset);
    this->mpCoupledDomainField->SetLabelOrigin(origin);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution); // nodal solution from pdeSolver

    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        /*
         * In both cases, with or without the extended cell property, need to 
         * refresh the cell concentration vector avaliable for transport from 
         * cell data to the appropriate concentration containers.
         * /
        CellPropertyCollection& prop_collection = cell_iter->rGetCellPropertyCollection();
        
        unsigned elem_index = mCellPdeElementMap[*cell_iter];
        
        if (prop_collection.HasProperty<ExtendedCellProperty<SPACE_DIM>>())
        {
            ///\todo Remove commented code?
            /*

            // the cells have a spatial extent where the state varibales are sorted with location keys in the extended cell property
            // therefore no need to interpolate form the nodes but need to recover and set the CellData
            auto extended_cell_property = boost::static_pointer_cast<ExtendedCellProperty<SPACE_DIM>>(prop_collection.GetPropertiesType<ExtendedCellProperty<SPACE_DIM>>().GetProperty());
                    
            // for each state variable in the cell data, retrieve the stored cell data value and add the total internal cell concetration
            
            std::vector<double> this_cell_next_total_boundary_concentration = extended_cell_property->GetNextTimestepConcentrationVector();
            unsigned number_cell_states =this_cell_next_total_boundary_concentration.size();
            std::vector<double> this_cell_total_concentration(number_cell_states,0.0);
            for (unsigned i=0; i<number_cell_states;++i)
            {
                // add the change in cell concentration due to the effects of the internal end of the cell boundary 
                double next_state_value = cell_iter->GetCellData()->GetItem(extended_cell_property->GetStateVariableRegister()->RetrieveStateVariableName(i)) + this_cell_next_total_concentration_change[i];
                cell_iter->GetCellData()->SetItem(extended_cell_property->GetStateVariableRegister()->RetrieveStateVariableName(i), next_state_value);
                this_cell_total_concentration[i] = next_state_value;
            }

            // next determine how much of the cell data may be utilised by the next solver call
            // this is a property of the cell so store in extended cell property
            extended_cell_property->UpdateCellConcentrationVector(this_cell_total_concentration);

            // remove this reserved concentration vector from the cell data
            for (unsigned i=0; i<number_cell_states;++i)
            {
                double cell_value = cell_iter->GetCellData()->GetItem(extended_cell_property->GetStateVariableRegister()->RetrieveStateVariableName(i))
                cell_iter->GetCellData()->SetItem(extended_cell_property->GetStateVariableRegister()->RetrieveStateVariableName(i), cell_value - this_cell_total_concentration[i]);
            }
            */
        }
        else
        {
            /*
             * The nodes are not extended so have no saved data regarding the 
             * simulation concentrations. Interpolate the nodal values to the 
             * location of the cells within their associated FeMesh elements.
             * 
             * For each cell in the population, find the element it belongs to 
             * then for each node in the elemnt (SPACE_DIM+1) sum the weighted 
             * nodal values for the pde solutions, then save sum to the cell.
             */

            // The cells are not nodes of the mesh, so we must interpolate
            std::vector<double> solution_vector_at_cell(PROBLEM_DIM,0.0);

            Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->mpFeMesh->GetElement(elem_index);

            const ChastePoint<SPACE_DIM>& node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            // Find weights of the nodal values based on the location of the cell relative to nodes
            c_vector<double, SPACE_DIM+1> weights = p_element->CalculateInterpolationWeights(node_location);

            for (unsigned pd = 0; pd < PROBLEM_DIM; pd++)
            {         
                // Start with whatever was left over from the transport ode call
                for (unsigned i = 0; i < SPACE_DIM+1; ++i)
                {                    
                    double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(i) + pd]; // serialised [node number * pde_index]
                    solution_vector_at_cell[pd] += nodal_value * weights(i);
                }
            }
       
            StateVariableRegister* p_bulk_register_pde = this->mpCoupledDomainField->GetDomainStateVariableRegister();

            if (prop_collection.HasProperty<TransportCellProperty>())
            { 
                auto transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(prop_collection.GetPropertiesType<TransportCellProperty>().GetProperty());
        
                // Only a subset of the solution vector at the cell are used by the transport property
                StateVariableRegister* p_bulk_register_cell = transport_cell_property->GetBulkStateVariableRegister();
                std::vector<double> subset_vector_at_cell(p_bulk_register_cell->GetNumStateVariables(), 0.0);
                unsigned domain_index = 0;
                for (unsigned i = 0; i < p_bulk_register_cell->GetNumStateVariables(); ++i)
                {
                    if (p_bulk_register_pde->IsStateVariablePresent(p_bulk_register_cell->RetrieveStateVariableName(i)))
                    {
                        domain_index =  p_bulk_register_pde->RetrieveStateVariableIndex(p_bulk_register_cell->RetrieveStateVariableName(i));
                        subset_vector_at_cell[i] = solution_vector_at_cell[domain_index];
                    }
                    // If state not found in domain register then automatically has 0.0 value
                }

                transport_cell_property->UpdateBulkConcentrationVector(subset_vector_at_cell);
            }

            if (prop_collection.HasProperty<MembraneCellProperty>())
            {
                auto membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(prop_collection.GetPropertiesType<MembraneCellProperty>().GetProperty());

                // Only a subset of the solution vector at the cell are used by the transport property
                StateVariableRegister* p_bulk_register_cell = membrane_cell_property->GetBulkStateVariableRegister();
                std::vector<double> subset_vector_at_cell(p_bulk_register_cell->GetNumStateVariables(), 0.0);
                unsigned domain_index = 0;
                for (unsigned i = 0; i < p_bulk_register_cell->GetNumStateVariables(); ++i)
                {
                    if (p_bulk_register_pde->IsStateVariablePresent(p_bulk_register_cell->RetrieveStateVariableName(i)))
                    {
                        domain_index =  p_bulk_register_pde->RetrieveStateVariableIndex(p_bulk_register_cell->RetrieveStateVariableName(i));
                        subset_vector_at_cell[i] = solution_vector_at_cell[domain_index];
                    }
                    // If state not found in domain register then automatically has 0.0 value
                }

                membrane_cell_property->UpdateBulkConcentrationVector(subset_vector_at_cell);
            }
        }

        if (this->mOutputGradient)
        {
            // Now calculate the gradient of the solution and store this in CellVecData
            for (unsigned pd = 0; pd < PROBLEM_DIM; pd++)
            {
                c_vector<double, SPACE_DIM> solution_gradient = zero_vector<double>(SPACE_DIM); // change?
                Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->mpFeMesh->GetElement(elem_index);

                // Calculate the basis functions at any point (e.g. zero) in the element
                c_matrix<double, SPACE_DIM, SPACE_DIM> jacobian, inverse_jacobian; // change? lloks like it's for only PROBLEM_DIM=1 but not nodal
                double jacobian_det;
                this->mpFeMesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
                const ChastePoint<SPACE_DIM> zero_point;
                c_matrix<double, SPACE_DIM, SPACE_DIM+1> grad_phi;
                LinearBasisFunction<SPACE_DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

                for (unsigned node_index=0; node_index<SPACE_DIM+1; node_index++)
                {
                    double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(node_index) + pd];
                    for (unsigned j = 0; j < SPACE_DIM; j++)
                    {
                        solution_gradient(j) += nodal_value* grad_phi(j, node_index); 
                    }
                }

                cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_x", solution_gradient(0));
                if (SPACE_DIM > 1)
                {
                    cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_y", solution_gradient(1));
                }
                if (SPACE_DIM > 2)
                {
                    cell_iter->GetCellData()->SetItem("var_node_"+this->mpPdeSystem->GetStateVariableRegister()->RetrieveStateVariableName(pd)+"_grad_z", solution_gradient(2));
                }
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    mCellPdeElementMap.clear();

    // Find the element of mpFeMesh that contains each cell and populate mCellPdeElementMap
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<SPACE_DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndex(r_position_of_cell);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::UpdateCellPdeElementMap(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Find the element of mpCoarsePdeMesh that contains each cell and populate mCellPdeElementMap
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<SPACE_DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndexWithInitialGuess(r_position_of_cell, mCellPdeElementMap[*cell_iter]);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractBoxDomainPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPdeSystemModifier<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*ABSTRACTBOXDOMAINPDESYSTEMMODIFIER_HPP_*/