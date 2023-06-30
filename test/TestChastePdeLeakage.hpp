#ifndef TESTCHASTEPDELEAKAGE_HPP_
#define TESTCHASTEPDELEAKAGE_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UblasIncludes.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "LinearParabolicHeatEquationPde.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "SimpleHeatEquation.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "TetrahedralMesh.hpp"

struct ParameterStruct
{
    // Standard values
    double MeshStepSize = 1.0;
    std::vector<unsigned> MeshDimensions = {100, 10};
    std::vector<double> initValuesHigh = {1.0};
    std::vector<double> initValuesLow = {0.0};
    std::vector<double> bcValuesDir = {0.0};
    std::vector<double> bcValuesNeu = {0.0};
    std::vector<double> diffusionRates = {100.0};

    // Solver properties
    double startTime = 0.0;
    double endTime = 10.0;
    double timestep = 0.01;
    double samplingTimestep = 0.1;
    std::string outputDirName = "ChasteLeakage/";
} params;

class TestChastePdeLeakage : public AbstractCellBasedTestSuite
{
public:

    void TestSimpleHeatDiffusionWithoutSource()
    {
        // Create mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();
        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);

        // Process boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        std::vector<ConstBoundaryCondition<2>*> vector_const_bcs;
        for (unsigned pde_dim = 0; pde_dim < 1; pde_dim++)
        {
            vector_const_bcs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[pde_dim]));
        }
        for (TetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
             boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
             boundary_iter++)
        {
            bcc.AddNeumannBoundaryCondition(*boundary_iter, vector_const_bcs[0], 0);
        }

        // Create initial conditions
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(), 0.0);
        unsigned column_num = 0;
        unsigned row_num = 0;
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            column_num = 0;
            row_num = 0;
            while (i >= row_num*(params.MeshDimensions[0]+1))
            {
                row_num++;
            }
            
            column_num = i - (row_num - 1)*(params.MeshDimensions[0] + 1);
            if ((column_num==3 || column_num==4 || column_num==5 || column_num==6) && (row_num ==3 || row_num ==4 || row_num ==5 || row_num ==6 ))
            {
                // Serialised for nodes
                init_conds[i] = fabs(params.initValuesHigh[0]);
            }
            else
            {
                // Serialised for nodes
                init_conds[i] = fabs(params.initValuesLow[0]);
            }
        }

        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        SimpleHeatEquation<2> pde;

        // Solver
        SimpleLinearParabolicSolver<2,2> solver(p_mesh, &pde, &bcc);
        solver.SetTimes(params.startTime, params.endTime);
        solver.SetTimeStep(params.timestep);
        solver.SetInitialCondition(initial_condition);
        solver.SetOutputDirectoryAndPrefix(params.outputDirName+"ChasteLeakageSimpleLinearParabolicSolver","results");
        solver.SetOutputToVtk(true);
        solver.SetPrintingTimestepMultiple(10);

        // Solve
        Vec solution = solver.Solve();
        ReplicatableVector solution_repl(solution);
    }

    void TestHeatDiffusionWithoutSource()
    {
        // Create mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();
        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
    
        // Process boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        
        std::vector<ConstBoundaryCondition<2>*> vector_const_bcs;
        vector_const_bcs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vector_const_bcs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0])); 

        for (TetrahedralMesh<2,2>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
             node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            bcc.AddDirichletBoundaryCondition(*node_iter, vector_const_bcs[0]);
        }

        for (TetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
        boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
        boundary_iter++)
        {
            bcc.AddNeumannBoundaryCondition(*boundary_iter, vector_const_bcs[1]);
        }
    
        // Create initial conditions
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned column_num = 0;
        unsigned row_num = 0;
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            // set as being a random perturbation about the boundary values
            column_num = 0;
            row_num = 0;

            while (i >= row_num*(params.MeshDimensions[0]+1))
            {
                row_num++;            
            }
            
            column_num = i - (row_num-1)*(params.MeshDimensions[0]+1);
            if ((column_num==3 || column_num==4 || column_num==5 || column_num==6) && (row_num ==3 || row_num ==4 || row_num ==5 || row_num ==6 ))
            {
                for (unsigned pde_dim=0; pde_dim<1; pde_dim++)
                {   // serialised for nodes
                    init_conds[1*i + pde_dim] = fabs(params.initValuesHigh[pde_dim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }
            else
            {
                for (unsigned pde_dim=0; pde_dim<1; pde_dim++)
                {
                    // serialised for nodes
                    init_conds[1*i + pde_dim] = fabs(params.initValuesLow[pde_dim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }
        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        LinearParabolicHeatEquationPde<2, 2, 1> pde(params.diffusionRates);

        // solver
        LinearParabolicPdeSystemWithCoupledOdeSystemSolver<2,2,1> solver(p_mesh, &pde, &bcc);
        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep);
        if (params.timestep > params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
        solver.SetOutputDirectory(params.outputDirName+"LinearParabolicPdeSystemWithCoupledOdeSystemSolver/");
        solver.SetInitialCondition(initial_condition);

        // solve
        solver.SolveAndWriteResultsToFile();
    }
};

#endif
