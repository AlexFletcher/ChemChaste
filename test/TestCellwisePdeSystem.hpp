#ifndef TESTCELLWISEPDESYSTEM_HPP_
#define TESTCELLWISEPDESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
 
#include "AbstractDomainField_templated.hpp"
#include "AbstractTransportOutReaction.hpp"
#include "AbstractTransportReaction.hpp"
#include "AbstractTransportReactionSystem.hpp"
#include "AbstractTransportReactionSystemFromFile.hpp"
#include "AbstractReversibleTransportReaction.hpp"
#include "ApoptoticCellProperty.hpp"
#include "ArchiveOpener.hpp"
#include "AveragedSourceParabolicPde_test.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ChemicalCellProperty.hpp"
#include "ChemicalDomainField_templated.hpp"
#include "ChemicalDomainFieldForCellCoupling.hpp"
#include "ConstBoundaryCondition.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "MembraneCellProperty.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ReplicatableVector.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransportCellProperty.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "WildTypeCellMutationState.hpp"

class TestCellwisePdeSystem : public AbstractCellBasedTestSuite
{
public:
    void TestModifiedTumourSpheroid()
    {
        EXIT_IF_PARALLEL;

        // Create mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<WildTypeCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<StemCellProliferativeType> p_stem_type(new StemCellProliferativeType);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*(p_model->GetStemCellG1Duration() + p_model->GetSG2MDuration());
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Define rank 1 PDE system with knowledge of the cell population and a parameter
        boost::shared_ptr<CellwiseSourceEllipticPde<2> > p_pde(new CellwiseSourceEllipticPde<2>(cell_population, -0.03));

        // Dirichlet boundary conditions held at 1.0 for the diffusiong species 
        boost::shared_ptr<ConstBoundaryCondition<2> > p_bc(new ConstBoundaryCondition<2>(1.0));
        bool is_neumann_bc = false;

        // Create domain box
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        boost::shared_ptr<ChasteCuboid<2> > p_cuboid (new ChasteCuboid<2>(lower, upper));

        // mCreate a PDE modifier for the elliptic pde with a box domain
        boost::shared_ptr<EllipticGrowingDomainPdeModifier<2> > p_pde_modifier(new EllipticGrowingDomainPdeModifier<2>(p_pde, p_bc, is_neumann_bc));
        p_pde_modifier->SetDependentVariableName("oxygen");

        // Create simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.AddSimulationModifier(p_pde_modifier);
        simulator.SetOutputDirectory("TestBoxTumourSpheroid_test");
        simulator.SetEndTime(1.0);

        // Add linear force to the simulation mesh
        boost::shared_ptr<GeneralisedLinearSpringForce<2> > p_linear_force (new GeneralisedLinearSpringForce<2>());
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();
    }

    void TestMeshBasedSquareMonolayerWithNeumanBcs()
    {
        // Create mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term), assume cell at each node in cell layer mesh
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property = cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); ++i)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 3.0 || cell_location(0) < 3.6)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, 0.1, 1.0, -10.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetOutputSolutionAtPdeNodes(true);
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithNeumannWithMeshOnSquare_test");

        // Run for 10 time steps
        for (unsigned i=0; i<10; ++i)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }
    }

    void TestChemicalDomainFieldForCellCoupling()
    {
        // Variables for the user modify
        std::string data_file_root = "/home/chaste/projects/ChemChaste/src/Data/CellwiseSchnakenberg/";
        std::string domain_file_name = "CellwiseDomain.csv";
        std::string domain_key_file_name = "CellwiseDomainLabelKey.csv";
        std::string ode_label_file_name = "CellwiseDomain.csv";
        std::string ode_key_file_name = "CellwiseOdeKey.csv";
        std::string diffusion_file_name = "Cellwise_DiffusionDatabase.csv";
        std::string initial_conditions_file_name = "CellwiseInitialConditions.csv";
        std::string boundary_conditions_file_name = "Cellwise_Boundary_Conditions.csv";

        // System properties
        const unsigned prob_dim = 2; // need to set manually to the number of diffusive variables for the pde solver to solve
        const unsigned space_dim = 2;
        const unsigned element_dim = 2;

        // solver properties
        double t_end = 10;
        double simulation_time_step = 1e-2;
        double sampling_time_step = 1e-2;
        std::string output_filename = "TestChemicalDomainFieldForCellCoupling";
    
        // generate domain
        // run the domain field set up and parse files
        ChemicalDomainFieldForCellCoupling<element_dim,space_dim,prob_dim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<element_dim,space_dim,prob_dim>(data_file_root,data_file_root+domain_file_name, data_file_root+domain_key_file_name, data_file_root+domain_file_name, data_file_root+domain_key_file_name, data_file_root+ode_label_file_name, data_file_root+ode_key_file_name, data_file_root+diffusion_file_name);

        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<element_dim,space_dim,prob_dim> solver(p_Pde_field->GetMeshGenerator()->GetMesh(), p_Pde_field->ReturnSharedPtrPdeSystem().get(), p_Pde_field->ReturnSharedPtrBoundaryConditionsContainer().get(), p_Pde_field->GetNodalOdeSystems(), p_Pde_field->GetNodalOdeSolvers());

        // solver properties
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(simulation_time_step);
        solver.SetSamplingTimeStep(sampling_time_step);
        solver.SetOutputDirectory(output_filename);
        Vec initial_condition = PetscTools::CreateVec(p_Pde_field->GetInitialNodeConditions());
        solver.SetInitialCondition(initial_condition);
    
        // solve
        solver.SolveAndWriteResultsToFile();

        // clean
        PetscTools::Destroy(initial_condition);
    }

    void TestCellTransportProperty()
    {
        AbstractChemistry* p_bulk_system_chemistry = new AbstractChemistry();
        AbstractChemistry* p_cell_system_chemistry = new AbstractChemistry();

        std::vector<AbstractChemical*> p_bulk_1 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_bulk_2 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_bulk_3 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_cell_1 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_cell_2 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_cell_3 = std::vector<AbstractChemical*>();

        std::vector<unsigned> stoich_bulk_1 = std::vector<unsigned>();
        std::vector<unsigned> stoich_cell_1 = std::vector<unsigned>();
        std::vector<unsigned> stoich_bulk_2 = std::vector<unsigned>();
        std::vector<unsigned> stoich_cell_2 = std::vector<unsigned>();
        std::vector<unsigned> stoich_bulk_3 = std::vector<unsigned>();
        std::vector<unsigned> stoich_cell_3 = std::vector<unsigned>();

        // r1: 2U + V->3U     forwardRate = 0.1
        // r2: U <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: V <- V           forwardRate = 0.3

        AbstractChemical* p_chemical_U = new AbstractChemical("U");
        p_bulk_system_chemistry->AddChemical(p_chemical_U);
        p_cell_system_chemistry->AddChemical(p_chemical_U);

        // add U to reactions
        p_bulk_1.push_back(p_chemical_U);
        stoich_bulk_1.push_back(2);
        p_cell_1.push_back(p_chemical_U);
        stoich_cell_1.push_back(3);
        p_bulk_2.push_back(p_chemical_U);
        stoich_bulk_2.push_back(1);
        p_cell_2.push_back(p_chemical_U);
        stoich_cell_2.push_back(1);

        AbstractChemical* p_chemical_V = new AbstractChemical("V");
        p_bulk_system_chemistry->AddChemical(p_chemical_V);
        p_cell_system_chemistry->AddChemical(p_chemical_V);
        // add U to reactions
        p_bulk_1.push_back(p_chemical_V);
        stoich_bulk_1.push_back(1);
        p_bulk_3.push_back(p_chemical_V);
        stoich_bulk_3.push_back(1);
        p_cell_3.push_back(p_chemical_V);
        stoich_cell_3.push_back(1);

        // for the sake of testing multiple species in the bulk which do not take part in the transport
        AbstractChemical* p_chemical_C = new AbstractChemical("C");
        p_bulk_system_chemistry->AddChemical(p_chemical_C);

        double reaction_1_rate = 0.1;
        double reaction_2_forward_rate = 0.1;
        double reaction_2_reverse_rate = 0.2;
        double reaction_3_rate = 0.3;

        AbstractTransportReaction* p_reaction_1 = new AbstractTransportReaction(p_bulk_1, p_cell_1, stoich_bulk_1, stoich_cell_1,reaction_1_rate);
        AbstractReversibleTransportReaction* p_reaction_2 = new AbstractReversibleTransportReaction(p_bulk_2, p_cell_2, stoich_bulk_2, stoich_cell_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractTransportOutReaction* p_reaction_3 = new AbstractTransportOutReaction(p_bulk_3, p_cell_3, stoich_bulk_3, stoich_cell_3,reaction_3_rate);

        std::vector<AbstractTransportReaction*> p_reaction_vector_2;
        p_reaction_vector_2.push_back(p_reaction_1);
        p_reaction_vector_2.push_back(p_reaction_2);
        p_reaction_vector_2.push_back(p_reaction_3);

        AbstractTransportReactionSystem* p_transport_reaction_system = new AbstractTransportReactionSystem(p_system_chemistry, p_reaction_vector_2);

        // states in the bulk environemnt
        std::vector<std::string> environment_states = {"U", "V", "C"};
        std::vector<double> environment_concentration_at_point = {1.0, 1.0, 1.0};
        std::vector<double> change_environment_concentration_vector = {0.0, 0.0, 0.0};

        // states bound within the cell
        std::vector<std::string> cell_states = {"U", "V"};
        std::vector<double> cell_concentration = {1.0, 1.0};
        std::vector<double> change_cell_concentration_vector = {0.0, 0.0};

        p_transport_reaction_system->ReactSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);

        // Make a single cell object        
        boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());
        boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
        p_cell_transport->SetUp(p_transport_reaction_system);
        change_environment_concentration_vector = {0.0, 0.0, 0.0};
        change_cell_concentration_vector = {0.0, 0.0, 0.0};
        p_cell_transport->PerformTransportSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);

        boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
        p_cell_membrane->SetMembraneThickness(1.0);

        AbstractTransportOdeSystem transportOde(p_transport_reaction_system);
        EulerIvpOdeSolver euler_solver;
        std::vector<double> initial_condition = {1.0, 1.0};

        OdeSolution solutions = euler_solver.Solve(&chemicalOde, initial_condition, 0, 1, 0.01, 0.1);

        /*
        std::cout << "Transport cell property" << std::endl;

        boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
        p_cell_transport->SetUp(p_transport_reaction_system);
        change_environment_concentration_vector= {0.0,0.0,0.0};
        change_cell_concentration_vector= {0.0,0.0,0.0};
        p_cell_transport->PerformTransportSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);

        for (unsigned i=0; i<environment_concentration_at_point.size();++i)
        {
            std::cout << "Bulk State "<<environment_states[i]<<": initial: "<<environment_concentration_at_point[i]<<" change: "<<change_environment_concentration_vector[i] << std::endl;
        }

        for (unsigned i=0; i<cell_concentration.size();++i)
        {
            std::cout << "Cell State "<<cell_states[i]<<": initial: "<<cell_concentration[i]<<" change: "<<change_cell_concentration_vector[i] << std::endl;
        }

        std::cout << "Membrane cell property" << std::endl;

        boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
        
        p_cell_membrane->SetMembraneThickness(1.0);

        std::cout << "Membrane thickness: "<<p_cell_membrane->GetMembraneThickness() << std::endl;

        // form the cells on a mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;

        std::vector<std::vector<double>> cell_concentration_vector;
        std::vector<std::vector<double>> change_cell_concentration_cell_vector;

        change_environment_concentration_vector= {0.0,0.0,0.0};

        for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
        {
            std::vector<std::string> membrane_states = {"U", "V"};
            std::vector<double> cell_concentration = {1.0,1.0};
            std::vector<double> change_cell_concentration= {0.0,0.0};
            cell_concentration_vector.push_back(cell_concentration);
            change_cell_concentration_cell_vector.push_back(change_cell_concentration);

            CellPropertyCollection collection;
            boost::shared_ptr<MembraneCellProperty> p_cell_membrane_iter(new MembraneCellProperty());

            p_cell_membrane_iter->InitialseCell(membrane_states, cell_concentration);

            AbstractTransportReactionSystem* p_reaction_system_2_iter = new AbstractTransportReactionSystem(p_system_chemistry, p_reaction_vector_2);

            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
            p_cell_transport->SetUp(p_reaction_system_2_iter);

            NoCellCycleModel* p_cell_cycle = new NoCellCycleModel();
            p_cell_membrane_iter->SetMembraneThickness((double)i);
            collection.AddProperty(p_cell_transport);
            collection.AddProperty(p_cell_membrane_iter);

            CellPtr p_cell(new Cell(p_state, p_cell_cycle, NULL, false, collection));
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // iterate through the cells
        unsigned cell_count=0;
        for (auto cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (cell_iter->HasCellProperty<MembraneCellProperty>())
            {
                auto property = boost::static_pointer_cast<MembraneCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<MembraneCellProperty>().GetProperty());

                std::cout << "Cell "<<cell_count<<": Membrane thickness: "<<property->GetMembraneThickness() << std::endl;
                for (unsigned i=0; i<property->GetStateVariableRegister()->GetNumStateVariables();++i)
                {
                    std::cout << "State: "<<property->GetStateVariableRegister()->RetrieveStateVariableName(i)<<" Concentration: "<< property->GetCellConcentrationByIndex(i) << std::endl;
                }
            }

            if (cell_iter->HasCellProperty<TransportCellProperty>())
            {
                auto p_cell_transport = boost::static_pointer_cast<TransportCellProperty>(cell_iter->rGetCellPropertyCollection().GetPropertiesType<TransportCellProperty>().GetProperty());

                p_cell_transport->PerformTransportSystem(environment_concentration_at_point,cell_concentration_vector[cell_count],change_environment_concentration_vector,change_cell_concentration_cell_vector[cell_count]);
            }
            cell_count++;
        }
        */
    }

    void TestAbstractCellCoupling()
    {
        /*
        AbstractChemistry* p_system_chemistry = new AbstractChemistry();

        std::vector<AbstractChemical*> p_bulk_1 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_bulk_2 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_bulk_3 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_cell_1 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_cell_2 = std::vector<AbstractChemical*>();
        std::vector<AbstractChemical*> p_cell_3 = std::vector<AbstractChemical*>();

        std::vector<unsigned> stoich_bulk_1 = std::vector<unsigned>();
        std::vector<unsigned> stoich_cell_1 = std::vector<unsigned>();
        std::vector<unsigned> stoich_bulk_2 = std::vector<unsigned>();
        std::vector<unsigned> stoich_cell_2 = std::vector<unsigned>();
        std::vector<unsigned> stoich_bulk_3 = std::vector<unsigned>();
        std::vector<unsigned> stoich_cell_3 = std::vector<unsigned>();

        // r1: 2U + V->3U     forwardRate = 0.1
        // r2: U <-> U          forwardRate = 0.1 reverseRate = 0.2
        // r3: V <- V           forwardRate = 0.3

        AbstractChemical* p_chemical_U = new AbstractChemical("U");
        p_system_chemistry->AddChemical(p_chemical_U);
        // add U to reactions
        p_bulk_1.push_back(p_chemical_U);
        stoich_bulk_1.push_back(2);
        p_cell_1.push_back(p_chemical_U);
        stoich_cell_1.push_back(3);
        p_bulk_2.push_back(p_chemical_U);
        stoich_bulk_2.push_back(1);
        p_cell_2.push_back(p_chemical_U);
        stoich_cell_2.push_back(1);

        AbstractChemical* p_chemical_V = new AbstractChemical("V");
        p_system_chemistry->AddChemical(p_chemical_V);
        // add U to reactions
        p_bulk_1.push_back(p_chemical_V);
        stoich_bulk_1.push_back(1);
        p_bulk_3.push_back(p_chemical_V);
        stoich_bulk_3.push_back(1);
        p_cell_3.push_back(p_chemical_V);
        stoich_cell_3.push_back(1);

        double reaction_1_rate = 0.1;
        double reaction_2_forward_rate = 0.1;
        double reaction_2_reverse_rate = 0.2;
        double reaction_3_rate = 0.3;

        AbstractTransportReaction* p_reaction_1 = new AbstractTransportReaction(p_bulk_1, p_cell_1, stoich_bulk_1, stoich_cell_1,reaction_1_rate);
        AbstractReversibleTransportReaction* p_reaction_2 = new AbstractReversibleTransportReaction(p_bulk_2, p_cell_2, stoich_bulk_2, stoich_cell_2,reaction_2_forward_rate, reaction_2_reverse_rate);
        AbstractTransportOutReaction* p_reaction_3 = new AbstractTransportOutReaction(p_bulk_3, p_cell_3, stoich_bulk_3, stoich_cell_3,reaction_3_rate);

        std::vector<AbstractTransportReaction*> p_reaction_vector_2;
        p_reaction_vector_2.push_back(p_reaction_1);
        p_reaction_vector_2.push_back(p_reaction_2);
        p_reaction_vector_2.push_back(p_reaction_3);

        AbstractTransportReactionSystem* p_transport_reaction_system = new AbstractTransportReactionSystem(p_system_chemistry, p_reaction_vector_2);

        // staes in the bulk environemnt
        std::vector<std::string> environment_states = {"U", "V"};
        std::vector<double> environment_concentration_at_point = {1.0,1.0};
        std::vector<double> change_environment_concentration_vector= {0.0,0.0};
        //StateVariableRegister* environment_register  = new StateVariableRegister(environment_states); 

        // states bound within the cell
        std::vector<std::string> cell_states = {"U", "V"};
        std::vector<double> cell_concentration = {1.0,1.0};
        std::vector<double> change_cell_concentration_vector= {0.0,0.0};

        p_transport_reaction_system->ReactSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);

        boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
        p_cell_transport->SetUp(p_transport_reaction_system);
        change_environment_concentration_vector= {0.0,0.0};
        change_cell_concentration_vector= {0.0,0.0};
        p_cell_transport->PerformTransportSystem(environment_concentration_at_point,cell_concentration,change_environment_concentration_vector,change_cell_concentration_vector);

        boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());
        p_cell_membrane->SetMembraneThickness(1.0);
        
        // form the cells on a mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;

        std::vector<std::vector<double>> cell_concentration_vector;
        std::vector<std::vector<double>> change_cell_concentration_cell_vector;

        change_environment_concentration_vector= {0.0,0.0};

        for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
        {
            CellPropertyCollection collection;
            boost::shared_ptr<MembraneCellProperty> p_cell_membrane_iter(new MembraneCellProperty());

            AbstractTransportReactionSystem* p_reaction_system_2_iter = new AbstractTransportReactionSystem(p_system_chemistry, p_reaction_vector_2);

            boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());
            p_cell_transport->SetUp(p_reaction_system_2_iter);

            NoCellCycleModel* p_cell_cycle = new NoCellCycleModel();
            p_cell_membrane_iter->SetMembraneThickness((double)i);
            collection.AddProperty(p_cell_transport);
            collection.AddProperty(p_cell_membrane_iter);

            CellPtr p_cell(new Cell(p_state, p_cell_cycle, NULL, false, collection));
            
            std::vector<double> cell_concentration = {1.0,1.0};
            std::vector<double> change_cell_concentration= {0.0,0.0};
            cell_concentration_vector.push_back(cell_concentration);
            change_cell_concentration_cell_vector.push_back(change_cell_concentration);
            
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        std::string data_file_root = "/home/chaste/projects/ChemChaste/src/Data/CellwiseSchnakenberg/";
        std::string domain_file_name = "CellwiseDomain.csv";
        std::string domain_key_file_name = "CellwiseDomainLabelKey.csv";
        std::string ode_label_file_name = "CellwiseDomain.csv";
        std::string ode_key_file_name = "CellwiseOdeKey.csv";
        std::string diffusion_file_name = "Cellwise_DiffusionDatabase.csv";
        std::string initial_conditions_file_name = "CellwiseInitialConditions.csv";
        std::string boundary_conditions_file_name = "Cellwise_Boundary_Conditions.csv";

        // System properties
        const unsigned prob_dim =2; // need to set manually to the number of diffusive variables for the pde solver to solve
        const unsigned space_dim = 2;
        const unsigned element_dim = 2;

        // solver properties
        double t_end = 10;
        double simulation_time_step = 1e-2;
        double sampling_time_step = 1e-2;
        std::string output_filename = "TestAbstractCellCoupling";        
        
        // generate domain
        // run the domain field set up and parse files
        ChemicalDomainFieldForCellCoupling<element_dim,space_dim,prob_dim>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<element_dim,space_dim,prob_dim>(data_file_root,data_file_root+domain_file_name, data_file_root+domain_key_file_name, data_file_root+domain_file_name, data_file_root+domain_key_file_name, data_file_root+ode_label_file_name, data_file_root+ode_key_file_name, data_file_root+diffusion_file_name);

        // create domain box
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        boost::shared_ptr<ChasteCuboid<2> > p_cuboid (new ChasteCuboid<2>(lower, upper));

        boost::shared_ptr<AbstractPdeSystemModifier<element_dim,space_dim,prob_dim>> p_pde_modifier(new AbstractPdeSystemModifier<element_dim,space_dim,prob_dim>(p_Pde_field));
    */

        /*
        // Make a PDE modifier for the elliptic pde with a box domain
        boost::shared_ptr<EllipticGrowingDomainPdeModifier<2> > p_pde_modifier(new EllipticGrowingDomainPdeModifier<2>(p_pde, p_bc, is_neumann_bc));
        p_pde_modifier->SetDependentVariableName("oxygen");

        // run simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.SetOutputDirectory("output_filename");
        simulator.SetEndTime(t_end);

        // add linear force to the simulation mesh
        boost::shared_ptr<GeneralisedLinearSpringForce<2> > p_linear_force (new GeneralisedLinearSpringForce<2>());
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
        
        // solver
        InhomogenousCoupledPdeOdeSolverTemplated<element_dim,space_dim,prob_dim> solver(p_Pde_field->GetMeshGenerator()->GetMesh(), p_Pde_field->ReturnSharedPtrPdeSystem().get(), p_Pde_field->ReturnSharedPtrBoundaryConditionsContainer().get(), p_Pde_field->GetNodalOdeSystems(), p_Pde_field->GetNodalOdeSolvers());

        // solver properties
        solver.SetTimes(0, t_end);
        solver.SetTimeStep(simulation_time_step);
        solver.SetSamplingTimeStep(sampling_time_step);
        solver.SetOutputDirectory(output_filename);
        Vec initial_condition = PetscTools::CreateVec(p_Pde_field->GetInitialNodeConditions());
        solver.SetInitialCondition(initial_condition);
        
        // solve
        solver.SolveAndWriteResultsToFile();
        
        // clean
        PetscTools::Destroy(initial_condition);
        */
    }

    void TestCellBasedPdesWithKnownSolution()
    {
        /*
        // PDE variables

        // Variables for the user modify
        std::string data_file_root = "/home/chaste/projects/ChemChaste/src/Data/CellwiseSchnakenberg/";
        std::string domain_file_name = "CellwiseDomain.csv";
        std::string domain_key_file_name = "CellwiseDomainLabelKey.csv";
        std::string ode_label_file_name = "CellwiseOdeSelector.csv";
        std::string ode_key_file_name = "CellwiseOdeKey.csv";
        std::string diffusion_file_name = "Cellwise_DiffusionDatabase.csv";
        std::string initial_conditions_file_name = "CellwiseInitialConditions.csv";
        std::string boundary_conditions_file_name = "Cellwise_Boundary_Conditions.csv";
        
        // System properties
        const unsigned prob_dim = 2; // need to set manually to the number of diffusive variables for the pde solver to solve
        const unsigned space_dim = 2;
        const unsigned element_dim = 2;
        
        // solver properties
        double t_end = 10;
        double simulation_time_step = 1e-2;
        double sampling_time_step = 1e-2;
        std::string output_filename = "TestCellBasedPdesWithKnownSolution";

        // cell properties
        std::vector<double> initialCellValuesOfStateVariables = {2.0, 0.75};
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;

        // experiment domain
        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<space_dim>, p_cuboid, (lower, upper));

        // generate domain
        // run the domain field set up and parse files
        ChemicalDomainField* p_Pde_field = new ChemicalDomainField(data_file_root,data_file_root+domain_file_name, data_file_root+domain_key_file_name, data_file_root+ode_label_file_name, data_file_root+ode_key_file_name, data_file_root+diffusion_file_name);

        // check that the file input problem dimension is the same as the user defined problem dimension
        std::cout << "File prob_dim: "<<p_Pde_field->GetProblemDimensions() << std::endl;
        std::cout << "User prob_dim: "<<prob_dim << std::endl;

        std::vector<std::string> stateVariableNames = p_Pde_field->GetStateVariableVector()->GetStateVariableRegisterVector();

        // cell mesh generator
        HoneycombMeshGenerator generator(5,5,0);
        MutableMesh<element_dim,space_dim>* p_cell_mesh = generator.GetMesh();
        cells_generator.GenerateBasicRandom(cells, p_cell_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property = cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); ++i)
        {
            c_vector<double,2> cell_location;
            cell_location = p_cell_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for PDE
            for (unsigned i=0; i<prob_dim; ++i)
            {
                cells[i]->GetCellData()->SetItem(stateVariableNames[0],initialCellValuesOfStateVariables[0]);
                cells[i]->GetCellData()->SetItem(stateVariableNames[1],initialCellValuesOfStateVariables[1]);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        MeshBasedCellPopulation<space_dim> cell_population(*p_cell_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<element_dim,space_dim,prob_dim>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<element_dim,space_dim,prob_dim>(p_Pde_field, p_cuboid));
        p_pde_modifier->SetDependentVariableNameVector(stateVariableNames); //needs to happen before SetupSolve

        // dont change below
        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupSolve(cell_population,output_filename);

        // Run for 10 time steps
        for (unsigned i=0; i<10; ++i)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }
        */
    }
};

#endif
