#ifndef TESTTOTALMASS_HPP_
#define TESTTOTALMASS_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ChemChasteVolumeAssembler.hpp"
#include "ChemChasteSurfaceAssembler.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "InhomogenousCoupledPdeOdeSolver_templated.hpp"
#include "InhomogenousHeatEquationPde.hpp"
#include "LinearParabolicPdeSystemWithCoupledOdeSystemSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleHeatEquation.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "UblasIncludes.hpp"

struct ParameterStruct
{
    // Standard values
    double MeshStepSize = 1.0;
    std::vector<unsigned> MeshDimensions = {100,10};
    std::vector<double> initValuesHigh = {1.0};
    std::vector<double> initValuesLow = {0.0};
    std::vector<double> bcValuesDir = {0.0};
    std::vector<double> bcValuesNeu = {0.0};
    std::vector<double> diffusionRates = {100.0};

    // Solver properties
    double startTime= 0.0;
    double endTime = 10.0;
    double timestep = 0.01;
    double samplingTimestep = 0.1;
    std::string outputDirName = "TotalMass";

} params;

class TestTotalMass : public AbstractCellBasedTestSuite
{
public:

    void TestTotalMassHeatDiffusionWithoutSource()
    {
        ofstream myfile("/home/chaste/testoutput/totalMassExample.txt");
        if (myfile.is_open())
        {
            myfile << "Writing this to a file.\n";
            myfile.close();
        }
        else
        {
            std::cout << "not open" << std::endl;
        }    

        // Mesh
        TetrahedralMesh<2,2>* p_mesh = new TetrahedralMesh<2,2>();
        p_mesh->ConstructRegularSlabMesh(params.MeshStepSize, params.MeshDimensions[0], params.MeshDimensions[1]);
    
        // Process boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        std::vector<ConstBoundaryCondition<2>*> vector_const_bcs;        
        vector_const_bcs.push_back(new ConstBoundaryCondition<2>(params.bcValuesDir[0]));
        vector_const_bcs.push_back(new ConstBoundaryCondition<2>(params.bcValuesNeu[0]));

        for (TetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter = p_mesh->GetBoundaryElementIteratorBegin();
             boundary_iter != p_mesh->GetBoundaryElementIteratorEnd();
             boundary_iter++)
        {
            bcc.AddNeumannBoundaryCondition(*boundary_iter, vector_const_bcs[1]);
        }
    
        // Initial conditions
        std::vector<double> init_conds(1*p_mesh->GetNumNodes(),0.0);
        unsigned column_num = 0;
        unsigned row_num = 0;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); ++i)
        {
            // Set as being a random perturbation about the boundary values            
            column_num = 0;
            row_num = 0;
            while (i >= row_num*(params.MeshDimensions[0]+1))
            {
                row_num = row_num + 1;            
            }
            
            column_num = i - (row_num-1)*(params.MeshDimensions[0]+1);
            if ((column_num==3 || column_num==4 || column_num==5 || column_num==6) && (row_num ==3 || row_num ==4 || row_num ==5 || row_num ==6 ))
            {
                for (unsigned pde_dim=0; pde_dim<1; pde_dim++)
                {
                    // Serialised for nodes
                    init_conds[1*i + pde_dim] = fabs(params.initValuesHigh[pde_dim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }
            else
            {
                for (unsigned pde_dim=0; pde_dim<1; pde_dim++)
                {
                    // Serialised for nodes
                    init_conds[1*i + pde_dim] = fabs(params.initValuesLow[pde_dim]);// + RandomNumberGenerator::Instance()->ranf());
                }
            }

        }
        // PETSc Vec
        Vec initial_condition = PetscTools::CreateVec(init_conds);

        InhomogenousHeatEquationPde<2, 2, 1> pde(params.diffusionRates);

        // Solver
        InhomogenousCoupledPdeOdeSolverTemplated<2,2,1> solver(p_mesh, &pde, &bcc);

        solver.SetTimes(params.startTime, params.endTime);

        solver.SetTimeStep(params.timestep);
        if (params.timestep>params.samplingTimestep)
        {
            solver.SetSamplingTimeStep(params.timestep);   
        }
        else
        {
            solver.SetSamplingTimeStep(params.samplingTimestep);
        }
        solver.SetOutputDirectory(params.outputDirName);
        solver.SetInitialCondition(initial_condition);

        // Solve
        solver.SolveAndWriteResultsToFile();
    }
};

#endif