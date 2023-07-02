#ifndef TESTCHEMCHASTE_HPP_
#define TESTCHEMCHASTE_HPP_

#include <cxxtest/TestSuite.h>

#include "BoundaryConditionsContainer_extended.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OffLatticeSimulation.hpp"
#include "ApoptoticCellProperty.hpp"
#include "ArchiveOpener.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ReplicatableVector.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "BoundaryConditionsContainer_extended.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "UblasIncludes.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

#include "AbstractBoxDomainPdeSystemModifier.hpp"
#include "AbstractChemical.hpp"
#include "AbstractChemicalOdeForCoupledPdeSystem.hpp"
#include "AbstractChemicalOdeSystem.hpp"
#include "AbstractChemistry.hpp"
#include "AbstractDiffusiveChemical.hpp"
#include "AbstractDiffusiveChemistry.hpp" 
#include "AbstractMembraneOdeSystem.hpp"
#include "AbstractMembraneReaction.hpp"
#include "AbstractMembraneReactionSystem.hpp"
#include "AbstractPdeSystemModifier.hpp"
#include "AbstractReaction.hpp"
#include "AbstractReactionSystem.hpp"
#include "AbstractReactionSystemFromFile.hpp"
#include "AbstractReversibleMembraneReaction.hpp"
#include "AbstractReversibleReaction.hpp"
#include "AbstractReversibleTransportReaction.hpp"
#include "AbstractTransportOdeSystem.hpp"
#include "AbstractTransportOutReaction.hpp"
#include "AbstractTransportReaction.hpp"
#include "AbstractTransportReactionSystem.hpp"
#include "ChemicalCell.hpp"
#include "ChemicalCellProperty.hpp"
#include "ChemicalDomainField_templated.hpp"
#include "ChemicalDomainFieldForCellCoupling.hpp"
#include "ChemicalSrnModel.hpp"
#include "ChemicalTrackingModifier.hpp"
#include "ExtendedCellProperty.hpp"
#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystem_templated.hpp"
#include "InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusion_templated.hpp"
#include "MassActionCoupledMembraneReaction.hpp"
#include "MassActionReaction.hpp"
#include "MassActionTransportReaction.hpp"
#include "MembraneCellProperty.hpp"
#include "NullSrnModel.hpp"
#include "OdeConsumerProducer.hpp"
#include "OdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "ParabolicBoxDomainPdeSystemModifier.hpp"
#include "PdeConsumerProducer.hpp"
#include "PdeSchnackenbergCoupledPdeOdeSystem.hpp"
#include "ReactionTypeDatabase.hpp"
#include "SchnackenbergCoupledPdeSystem.hpp"
#include "SchnackenbergSrnModel.hpp"
#include "SimpleChemicalThresholdCellCycleModel.hpp"
#include "TransportCellProperty.hpp"

#include "ChemicalCellFromFile.hpp"
#include "ChemicalStructuresForTests.hpp"
#include "InhomogenousFisherPde.hpp"
#include "InhomogenousFisherDiffusiveInhibitionPde.hpp"

struct ControlStruct
{
    bool ReactionSystemWithoutCells = false;
    bool Fisher = true;
    bool FisherDiffusiveInhibition = false;
    bool ReactionSystemWithoutCellsUsingInhomogenousSolver = false;
    bool ReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver = false;
    bool ReactionSystemWithCells = false;
    bool ReactionSystemWithNodeBasedCells = false;
} control;

class TestChemChaste : public AbstractCellBasedTestSuite
{
public:

    void TestReactionSystemWithoutCells()
    {
        if (control.ReactionSystemWithoutCells)
        {
            std::cout << "Mass action kinetics of Schnackenberg using the AbstractChemicalOdeSystem Class" << std::endl;
            std::cout << "-----------------------------" << std::endl;

            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
            chemical_structure->SetUpPdeChemicalReactionSystem();
            AbstractReactionSystem* chemicalReactionSystem = chemical_structure->rGetPtrPdeChemicalReactionSystem();

            // form ode system
            std::cout << "SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem" << std::endl;

            // system properties
            std::vector<double> initValues = {2.0, 0.75};
            std::vector<double> bcValues = {0.0, 0.0};

            // mesh
            HoneycombMeshGenerator generator(100, 100, 0);
            MutableMesh<2, 2>* p_mesh = generator.GetMesh();


            std::cout << "number mesh nodes: "<<p_mesh->GetNumNodes() << std::endl;

            // Process Boundary Conditions
            BoundaryConditionsContainer<2, 2, 2> bcc;
            std::vector<bool> areNeumannBoundaryConditions(2, true);
            std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
            for (unsigned pdeDim=0; pdeDim < 2; pdeDim++)
            {
                vectorConstBCs.push_back(new ConstBoundaryCondition<2>(bcValues[pdeDim]));
            }
            for (unsigned pdeDim=0; pdeDim < 2; pdeDim++)
            {
                if (areNeumannBoundaryConditions[pdeDim] == false)
                {
                    for (TetrahedralMesh<2, 2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                         node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                         ++node_iter)
                    {
                        bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
                else
                {
                    for (TetrahedralMesh<2, 2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                         boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                         boundary_iter++)
                    {
                        bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
            }
            std::cout << "Initial conditions" << std::endl;
            // initial conditions
            std::vector<double> init_conds(2*p_mesh->GetNumNodes());
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {
                // set as being a random perturbation about the boundary values
                for (unsigned pdeDim=0; pdeDim < 2; pdeDim++)
                {
                       // serialised for nodes
                    init_conds[2*i + pdeDim] = fabs(initValues[pdeDim]+ RandomNumberGenerator::Instance()->ranf());
                }
            }
            // PETSc Vec
            std::cout << "PETSc Vec" << std::endl;
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            // coupled ode system
            std::cout << "Ode loop" << std::endl;
            std::vector<AbstractOdeSystemForCoupledPdeSystem*> odeSystem;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {
                // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
                odeSystem.push_back(new AbstractChemicalOdeForCoupledPdeSystem(chemicalReactionSystem));
            }

            // PDE system
            PdeSchnackenbergCoupledPdeOdeSystem<2, 2, 2> pde(odeSystem[0], 1e-4, 1e-2);
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

            // solver
            LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2, 2, 2> solver(p_mesh, &pde, &bcc,odeSystem,p_solver);

            // solver properties        
            solver.SetTimes(0, 30);
            solver.SetTimeStep(1e-2);
            solver.SetSamplingTimeStep(1e-2);
            solver.SetOutputDirectory("TestReactionSystemWithoutCells_ver2");
            solver.SetInitialCondition(initial_condition);

            // solve
            solver.SolveAndWriteResultsToFile();
        }
    }

    void TestFisher()
    {
        if (control.Fisher)
        {
            // system properties
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValuesHigh = {0.5};
            std::vector<double> initValuesLow = {0.0};
            std::vector<double> bcValues = {0.0};
            std::vector<double> areBCsNeumann = {true};
            std::vector<double> diffusionRates = {1.0};
            std::vector<double> growthRates = {1.0};
            std::vector<double> carryingCapacities = {1.0};

            // mesh
            TetrahedralMesh<2, 2>* p_mesh = new TetrahedralMesh<2, 2>();
            p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);

            // Process Boundary Conditions
            std::cout << "Process Boundary Conditions" << std::endl;
            BoundaryConditionsContainer<2, 2, 1> bcc;
            std::vector<bool> areNeumannBoundaryConditions(1, true);
            std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
            
            for (unsigned pdeDim=0; pdeDim<1; pdeDim++)
{
                vectorConstBCs.push_back(new ConstBoundaryCondition<2>(bcValues[pdeDim]));
            }
            
            for (unsigned pdeDim=0; pdeDim<1; pdeDim++)
            {
                if (areNeumannBoundaryConditions[pdeDim] == false)
                {
                    for (TetrahedralMesh<2, 2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                         node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                         ++node_iter)
                    {

                        bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
                else
                {
                    for (TetrahedralMesh<2, 2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                         boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                         boundary_iter++)
                    {
                        bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
            }

            // initial conditions
            std::vector<double> init_conds(1*p_mesh->GetNumNodes(), 0.0);
            unsigned column_num = 0;
            unsigned row_num = 0;
            for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
            {   
                // Set as being a random perturbation about the boundary values
                column_num = 0;
                row_num = 0;
                while (i >= row_num*(MeshDimensions[0]+1))
                {
                    row_num = row_num + 1;                    
                }
                
                column_num = i - (row_num-1)*(MeshDimensions[0]+1);
                
                if (column_num < 10)
                {
                    // Serialised for nodes
                    init_conds[i] = fabs(initValuesHigh[0]);
                }
                else
                {
                    // Serialised for nodes
                    init_conds[i] = fabs(initValuesLow[0]);
                }
            }
            // PETSc Vec
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;
    
            // PDE system
            InhomogenousFisherPde<2, 2, 1> pde(diffusionRates,growthRates,carryingCapacities);

            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<2, 2,1> solver(p_mesh, &pde, &bcc);

            // solver properties        
            solver.SetTimes(0, 200);
            solver.SetTimeStep(1e-1);
            solver.SetSamplingTimeStep(1e-1);
            solver.SetOutputDirectory("TestFisherPaper");
            solver.SetInitialCondition(initial_condition);

            // solve
            solver.SolveAndWriteResultsToFile();
        }
    }

    void TestFisherDiffusiveInhibition()
    {
        if (control.FisherDiffusiveInhibition)
        {
            // system properties
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValuesHigh = {0.5};
            std::vector<double> initValuesLow = {0.0};
            std::vector<double> bcValues = {0.0};
            std::vector<double> areBCsNeumann = {true};
            std::vector<double> diffusionRates = {1.0};
            std::vector<double> growthRates = {1.0};
            std::vector<double> carryingCapacities = {1.0};

            // mesh
            TetrahedralMesh<2, 2>* p_mesh = new TetrahedralMesh<2, 2>();
             p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);

            // Process Boundary Conditions
            std::cout << "Process Boundary Conditions" << std::endl;
            BoundaryConditionsContainer<2, 2, 1> bcc;
            std::vector<bool> areNeumannBoundaryConditions(1, true);
            std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
            
            for (unsigned pdeDim=0; pdeDim<1; pdeDim++)
            {
                vectorConstBCs.push_back(new ConstBoundaryCondition<2>(bcValues[pdeDim]));
            }
            
            for (unsigned pdeDim=0; pdeDim<1; pdeDim++)
            {
                if (areNeumannBoundaryConditions[pdeDim] == false)
                {
                    for (TetrahedralMesh<2, 2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                         node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                         ++node_iter)
                    {

                        bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
                else
                {
                    for (TetrahedralMesh<2, 2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                         boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                         boundary_iter++)
                    {
                        bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
            }

            // initial conditions
            std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
            unsigned column_num = 0;
            unsigned row_num = 0;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {
                // set as being a random perturbation about the boundary values
                column_num = 0;
                row_num = 0;
                while (i >= row_num*(MeshDimensions[0]+1))
                {
                    row_num = row_num + 1;
                }
                
                column_num = i - (row_num-1)*(MeshDimensions[0]+1);
                if (column_num<10)
                {
                    for (unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {
                        // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesHigh[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }
                else
                {
                    for (unsigned pdeDim=0; pdeDim<1; pdeDim++)
                    {   // serialised for nodes
                        init_conds[1*i + pdeDim] = fabs(initValuesLow[pdeDim]);// + RandomNumberGenerator::Instance()->ranf());
                    }
                }

            }
            // PETSc Vec
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;

            // PDE system
            InhomogenousFisherDiffusiveInhibitionPde<2, 2, 1> pde(diffusionRates,growthRates,carryingCapacities);

            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<2, 2, 1> solver(p_mesh, &pde, &bcc);

            // solver properties
            solver.SetTimes(0,      00);
            solver.SetTimeStep(1);
            solver.SetSamplingTimeStep(1e-1);
            solver.SetOutputDirectory("TestFisherPaperDiffusionInhibition_2_close_regions");
            solver.SetInitialCondition(initial_condition);

            // solve
            solver.SolveAndWriteResultsToFile();
        }
    }

    void TestReactionSystemWithoutCellsInhomogenousSolver()
    {
        if (control.ReactionSystemWithoutCellsUsingInhomogenousSolver)
        {            
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
            chemical_structure->SetUpPdeChemicalReactionSystem();
            chemical_structure->SetUpChemicalDomainField();
            AbstractReactionSystem* chemicalReactionSystem = chemical_structure->rGetPtrPdeChemicalReactionSystem();

            // form ode system
            std::cout << "SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem" << std::endl;
        
            // system properties
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValues = {2.0, 0.75};
            std::vector<double> bcValues = {0.0, 0.0};
            std::vector<double> diffusionRates = {0.05, 0.05};

            // mesh
            TetrahedralMesh<2, 2>* p_mesh = new TetrahedralMesh<2, 2>();
            p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);


            // Process Boundary Conditions
            std::cout << "Process Boundary Conditions" << std::endl;
            BoundaryConditionsContainer<2, 2, 2> bcc;
            std::vector<bool> areNeumannBoundaryConditions(2, true);
            std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
            
            for (unsigned pdeDim=0; pdeDim<2; pdeDim++)
            {
                vectorConstBCs.push_back(new ConstBoundaryCondition<2>(bcValues[pdeDim]));
            }
            
            for (unsigned pdeDim=0; pdeDim<2; pdeDim++)
            {
                if (areNeumannBoundaryConditions[pdeDim]==false)
                {
                    for (TetrahedralMesh<2, 2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                         node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                         ++node_iter)
                    {

                        bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }else{
                    for (TetrahedralMesh<2, >::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                         boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                         boundary_iter++)
                    {
                        bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
            }

            // Initial conditions
            std::vector<double> init_conds(2*p_mesh->GetNumNodes());
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {
                // set as being a random perturbation about the boundary values
                for (unsigned pdeDim=0; pdeDim<2; pdeDim++)
                {
                    // serialised for nodes
                    init_conds[2*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }
            // PETSc Vec
            Vec initial_condition = PetscTools::CreateVec(init_conds);
            std::cout << "Here" << std::endl;
            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {
                // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
                odeSystem.push_back(new AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(chemicalReactionSystem));
                boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
                solverSystem.push_back(p_solver);//std::dynamic_pointer_cast<AbstractIvpOdeSolver>(p_solver));
            }

            // PDE system
            InhomogenousParabolicPdeForCoupledOdeSystemTemplated<2, 2, 2> pde(chemical_structure->rGetPtrChemicalDomain());

            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<2, 2, 2> solver(p_mesh, &pde, &bcc,odeSystem,solverSystem);

            // solver properties
            solver.SetTimes(0, 20);
            solver.SetTimeStep(1e-2);
            solver.SetSamplingTimeStep(1e-2);
            solver.SetOutputDirectory("TestReactionSystemWithoutCellsUsingInhomogenousSolver_single_larger_only_reaction_within");
            solver.SetInitialCondition(initial_condition);

            // solve
            solver.SolveAndWriteResultsToFile();
        }
    }

    void TestReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver()
    {
        if (control.ReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver)
        {            
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
            chemical_structure->SetUpPdeChemicalReactionSystem();
            chemical_structure->SetUpChemicalDomainField();
            AbstractReactionSystem* chemicalReactionSystem = chemical_structure->rGetPtrPdeChemicalReactionSystem();

            // form ode system
            std::cout << "SchnackenbergCoupledPdeOdeSystem  -As AbstractChemicalODeSystem" << std::endl;

            // system properties
            double MeshStepSize = 1.0;
            std::vector<unsigned> MeshDimensions = {100,100,0};
            std::vector<double> initValues = {2.0, 0.75};
            std::vector<double> bcValues = {0.0, 0.0};
            std::vector<double> diffusionRates = {0.05, 0.05};

            // mesh
            TetrahedralMesh<2, 2>* p_mesh = new TetrahedralMesh<2, 2>();
            p_mesh->ConstructRegularSlabMesh(MeshStepSize, MeshDimensions[0], MeshDimensions[1]);
            
            // Process Boundary Conditions
            BoundaryConditionsContainer<2, 2, 2> bcc;
            std::vector<bool> areNeumannBoundaryConditions(2, true);
            std::vector<ConstBoundaryCondition<2>*> vectorConstBCs;
            
            for (unsigned pdeDim=0; pdeDim< 2; pdeDim++)
            {
                vectorConstBCs.push_back(new ConstBoundaryCondition<2>(bcValues[pdeDim]));
            }
            
            for (unsigned pdeDim=0; pdeDim< 2; pdeDim++)
            {
                if (areNeumannBoundaryConditions[pdeDim] == false)
                {
                    for (TetrahedralMesh<2, 2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                         node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                         ++node_iter)
                    {

                        bcc.AddDirichletBoundaryCondition(*node_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
                else
                {
                    for (TetrahedralMesh<2, 2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
                         boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
                         boundary_iter++)
                    {
                        bcc.AddNeumannBoundaryCondition(*boundary_iter, vectorConstBCs[pdeDim], pdeDim);
                    }
                }
            }

            // initial conditions
            std::vector<double> init_conds(2*p_mesh->GetNumNodes());
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {   
                // set as being a random perturbation about the boundary values
                for (unsigned pdeDim=0; pdeDim<2; pdeDim++)
                {   
                    // serialised for nodes
                    init_conds[2*i + pdeDim] = fabs(initValues[pdeDim] + RandomNumberGenerator::Instance()->ranf());
                }
            }
            // PETSc Vec
            Vec initial_condition = PetscTools::CreateVec(init_conds);

            std::vector<AbstractInhomogenousOdeSystemForCoupledPdeSystem*> odeSystem;
            std::vector<boost::shared_ptr<AbstractIvpOdeSolver> > solverSystem;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
{
                // number of ode system objects must match the number of nodes, i.e the individual odes may be multi-dimensional
                odeSystem.push_back(new AbstractInhomogenousChemicalOdeSystemForCoupledPdeSystem(chemicalReactionSystem));
                boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
                solverSystem.push_back(p_solver);//std::dynamic_pointer_cast<AbstractIvpOdeSolver>(p_solver));
            }

            // PDE system
            InhomogenousParabolicPdeForCoupledOdeSystemInhibitedDiffusionTemplated<2, 2, 2> pde(chemical_structure->rGetPtrChemicalDomain());

            // solver
            InhomogenousCoupledPdeOdeSolverTemplated<2, 2, 2> solver(p_mesh, &pde, &bcc,odeSystem,solverSystem);

            // solver properties        
            solver.SetTimes(0, 20);
            solver.SetTimeStep(1e-2);
            solver.SetSamplingTimeStep(1e-2);
            solver.SetOutputDirectory("TestReactionSystemInhibitedDiffusionWithoutCellsUsingInhomogenousSolver_single_larger_only_reaction_within");
            solver.SetInitialCondition(initial_condition);

            // solve
            solver.SolveAndWriteResultsToFile();
        }
    }

    void TestReactionSystemWithCells()
    {    
        if (control.ReactionSystemWithCells)
        {
            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
            MutableMesh<2,2>* p_mesh = generator.GetMesh();

            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states
                ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
                
                if (i==0)
                {
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();
                    cells.push_back(p_cell);
                }
                else
                {
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();

                    cells.push_back(p_cell);
                }                
            }
        
            MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<2> lower(-4.0, -4.0);
            ChastePoint<2> upper(7.0, 7.0);
            MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

            chemical_structure->SetUpChemicalDomainFieldForCellCoupling();
            
            ChemicalDomainFieldForCellCoupling<2,2,2>*& p_Pde_field = chemical_structure->rGetPtrChemicalDomainFieldForCellCoupling();

            std::cout << "initial conditions" << std::endl;
            std::vector<double> init_conditions = p_Pde_field->GetInitialNodeConditions();
            for (unsigned i=0; i<init_conditions.size(); ++i)
            {
                std::cout << init_conditions[i] << std::endl;
            }
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<2,2,2>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<2,2,2>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<2,2>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<2,2>()); //= chemical_structure->rGetPtrChemicalTrackingModifier();
            
            // writers
            cell_population.SetWriteVtkAsPoints(false);
            cell_population.AddPopulationWriter<VoronoiDataWriter>();

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            cell_population.AddCellWriter<CellLocationIndexWriter>();

            std::vector<std::string> chemicalCellASpeciesNames = chemical_structure->GetChemicalCellASpeciesNames();

            for (unsigned i=0; i<chemicalCellASpeciesNames.size(); ++i)
            {
                boost::shared_ptr<CellDataItemWriter<2,2>> dataWriter(new CellDataItemWriter<2,2>(chemicalCellASpeciesNames[i]));
                cell_population.AddCellWriter(dataWriter);
            }
        
            OffLatticeSimulation<2> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
        
            simulator.SetOutputDirectory("TestReactionSystemWithCells_OscillatingCase_meshWriters");
            simulator.SetEndTime(4.0);

            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
            p_linear_force->SetCutOffLength(1.5);
            simulator.AddForce(p_linear_force);

            simulator.Solve();
        }
    }

    void TestReactionSystemWithNodeBasedCells()
    {
        if (control.ReactionSystemWithNodeBasedCells)
        {
            // make mesh of cell with associated mesh based cell population
            HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
            MutableMesh<2,2>* p_mesh = generator.GetMesh();

            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(*p_mesh, 1.5);

            std::vector<CellPtr> cells;
            // assume cell at each node in cell layer mesh
            for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
            {
                // provide each cell with a transport cell property and membrane property, cell cycle, wild type states
                ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();
                
                if (i==0)
                {
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();
                    cells.push_back(p_cell);
                }
                else
                {
                    CellPtr p_cell = chemical_structure->SetUpChemicalCellObjectA();

                    cells.push_back(p_cell);
                }                
            }

            NodeBasedCellPopulation<2> cell_population(mesh, cells);
        
            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<2> lower(-4.0, -4.0);
            ChastePoint<2> upper(7.0, 7.0);
            MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

            // set up chemical domain field
            ChemicalStructuresForTests* chemical_structure = new ChemicalStructuresForTests();

            chemical_structure->SetUpChemicalDomainFieldForCellCoupling();
            
            ChemicalDomainFieldForCellCoupling<2,2,2>*& p_Pde_field = chemical_structure->rGetPtrChemicalDomainFieldForCellCoupling();

            std::cout << "initial conditions" << std::endl;
            std::vector<double> init_conditions = p_Pde_field->GetInitialNodeConditions();
            for (unsigned i=0; i<init_conditions.size(); ++i)
            {
                std::cout << init_conditions[i] << std::endl;
            }
            
            boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<2,2,2>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<2,2,2>(p_Pde_field, p_cuboid));
            
            boost::shared_ptr<ChemicalTrackingModifier<2,2>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<2,2>()); //= chemical_structure->rGetPtrChemicalTrackingModifier();

            OffLatticeSimulation<2> simulator(cell_population);

            simulator.AddSimulationModifier(p_pde_modifier);

            simulator.AddSimulationModifier(p_chemical_tracking_modifier);
        
            simulator.SetOutputDirectory("TestReactionSystemWithCells_OscillatingCase_meshWriters_nodeBased");
            simulator.SetEndTime(4.0);

            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
            p_linear_force->SetCutOffLength(1.5);
            simulator.AddForce(p_linear_force);

            simulator.Solve();
        }
    }

    void TestCellLayerRead()
    {
        std::string dataFileRoot = "/home/chaste/projects/ChemChaste/DataInput/Data/MulticellCase/DomainField/";
        std::string cellFileRoot = "/home/chaste/projects/ChemChaste/DataInput/Data/MulticellCase/Cell/";
        std::string cellLabelFilename = "CellLayerTopology.csv";
        std::string cellKeyFilename = "CellLayerKey.csv";
        std::string domainFilename = "Domain.csv";
        std::string domainKeyFilename = "DomainKey.csv";
        std::string odeLabelFilename = "NodeSelector.csv";
        std::string odeKeyFilename = "OdeReactionFileKey.csv";
        std::string diffusionFilename = "DiffusionDatabaseFile.csv";
        std::string initialConditionsFilename = "InitialConditionFile.csv";
        std::string boundaryConditionsFilename = "BoundaryConditionFile.csv";
        
        // generate domain
        // run the domain field set up and parse files
        ChemicalDomainFieldForCellCoupling<2, 2, 4>* p_Pde_field = new ChemicalDomainFieldForCellCoupling<2, 2, 4>(dataFileRoot,cellFileRoot+cellLabelFilename,cellFileRoot+cellKeyFilename,dataFileRoot+domainFilename, dataFileRoot+domainKeyFilename, dataFileRoot+odeLabelFilename, dataFileRoot+odeKeyFilename, dataFileRoot+diffusionFilename, dataFileRoot+initialConditionsFilename, dataFileRoot+boundaryConditionsFilename);

        TetrahedralMesh<2, 2>* p_cell_mesh = p_Pde_field->rGetCellMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_cell_mesh, 1.5);

        std::vector<CellPtr> cells;

        // assume cell at each node in cell layer mesh
        std::string cell_label;
        std::string cell_key;
        std::string given_cell_root;
        for (unsigned i=0; i<p_cell_mesh->GetNumNodes(); ++i)
        {
            cell_label = p_Pde_field->GetCellLabelByIndex(i);
            cell_key = p_Pde_field->ReturnCellKeyFromCellLabel(cell_label);
            std::cout << "Cell i: "<<i<<" label: "<<cell_label << std::endl;
            std::cout << "Cell i: "<<i<<" label: "<<cell_key << std::endl;
            given_cell_root = cellFileRoot+cell_key+"/";

            ChemicalCellFromFile* p_cell_reader = new ChemicalCellFromFile(
                                given_cell_root+"SpeciesThreshold.csv", 
                                given_cell_root+"Srn.txt",
                                given_cell_root+"InitialCellConcentrations.csv",
                                given_cell_root+"TransportReactions.txt",
                                given_cell_root+"MembraneReactions.txt");

            cells.push_back(p_cell_reader->GetCellPtr());
        }   

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        
        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-4.0, -4.0);
        ChastePoint<2> upper(7.0, 7.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        std::cout << "initial conditions" << std::endl;
        std::vector<double> init_conditions = p_Pde_field->GetInitialNodeConditions();
        for (unsigned i=0; i<init_conditions.size(); ++i)
        {
            std::cout << init_conditions[i] << std::endl;
        }
        
        boost::shared_ptr<ParabolicBoxDomainPdeSystemModifier<2, 2, 4>> p_pde_modifier(new ParabolicBoxDomainPdeSystemModifier<2, 2, 4>(p_Pde_field, p_cuboid));
        
        boost::shared_ptr<ChemicalTrackingModifier<2,2>> p_chemical_tracking_modifier(new ChemicalTrackingModifier<2,2>()); //= chemical_structure->rGetPtrChemicalTrackingModifier();
        
        OffLatticeSimulation<2> simulator(cell_population);

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.AddSimulationModifier(p_chemical_tracking_modifier);
    
        simulator.SetOutputDirectory("TestNodeBasedCellsFromFile");
        simulator.SetEndTime(10.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
    }
};

#endif