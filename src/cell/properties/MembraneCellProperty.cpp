#include "MembraneCellProperty.hpp"

MembraneCellProperty::MembraneCellProperty()
    : AbstractCellProperty()
{
}

MembraneCellProperty::~MembraneCellProperty()
{
}

MembraneCellProperty::MembraneCellProperty(const MembraneCellProperty& existingProperty)
{
    mIsExtent = existingProperty.mIsExtent;
    mMembraneThickness = existingProperty.mMembraneThickness;
    mpMembraneStateVariableRegister = existingProperty.mpMembraneStateVariableRegister;
    mMembraneConcentrationVector = existingProperty.mMembraneConcentrationVector;
    mpBulkStateVariableRegister = existingProperty.mpBulkStateVariableRegister;
    mpCellStateVariableRegister = existingProperty.mpCellStateVariableRegister;
    mpMembraneReactionSystem = existingProperty.mpMembraneReactionSystem;
    mpMembraneOdeSystem = existingProperty.mpMembraneOdeSystem;
    mpMembraneOdeSolver = existingProperty.mpMembraneOdeSolver;
    mIsConstantCellConcentration = existingProperty.mIsConstantCellConcentration;
    mBulkBoundaryConcentrationVector = existingProperty.mBulkBoundaryConcentrationVector;
    mChangeBulkBoundaryConcentrationVector = existingProperty.mChangeBulkBoundaryConcentrationVector;
    mCellBoundaryConcentrationVector = existingProperty.mCellBoundaryConcentrationVector;
    mChangeCellBoundaryConcentrationVector = existingProperty.mChangeCellBoundaryConcentrationVector;
    mIncludeMembraneOdeInterpolationOnBoundary = existingProperty.mIncludeMembraneOdeInterpolationOnBoundary;
}

void MembraneCellProperty::SetUp(AbstractMembraneReactionSystem* membraneReactionSystem, CellPtr pCell)
{
    UpdateMembraneReactionSystem(membraneReactionSystem);
    SetCellPtr(pCell);

    AbstractMembraneOdeSystem* p_ode_system = new AbstractMembraneOdeSystem(membraneReactionSystem);
    UpdateMembraneOdeSystem(p_ode_system);

    StateVariableRegister* p_bulk_state_register = new StateVariableRegister(membraneReactionSystem->GetBulkChemistry()->GetChemicalNames());
    SetBulkStateVariableRegister(p_bulk_state_register);

    StateVariableRegister* p_cell_state_register = new StateVariableRegister(membraneReactionSystem->GetCellChemistry()->GetChemicalNames());
    SetCellStateVariableRegister(p_cell_state_register);

    SetUpCellConcentrationVector(p_cell_state_register->GetNumStateVariables());
    SetUpBulkConcentrationVector(p_bulk_state_register->GetNumStateVariables());
    SetUpChangeCellConcentrationVector(p_cell_state_register->GetNumStateVariables());
    SetUpChangeBulkConcentrationVector(p_bulk_state_register->GetNumStateVariables());
}

void MembraneCellProperty::UpdateMembraneConcentrationVector(std::vector<double> cellBoundaryConcentrationVector)
{
    mMembraneConcentrationVector = cellBoundaryConcentrationVector;
}

void MembraneCellProperty::UpdateCellConcentrationVector(std::vector<double> cellBoundaryConcentrationVector)
{
    mCellBoundaryConcentrationVector = cellBoundaryConcentrationVector;
}

void MembraneCellProperty::UpdateBulkConcentrationVector(std::vector<double> bulkBoundaryConcentrationVector)
{
    mBulkBoundaryConcentrationVector = bulkBoundaryConcentrationVector;
}

void MembraneCellProperty::SetUpCellConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables,0.0);
    mCellBoundaryConcentrationVector = reset;
}

void MembraneCellProperty::SetUpBulkConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables,0.0);
    mBulkBoundaryConcentrationVector = reset;
}

void MembraneCellProperty::SetUpChangeCellConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables,0.0);
    mChangeCellBoundaryConcentrationVector = reset;
}

void MembraneCellProperty::SetUpChangeBulkConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables,0.0);
    mChangeBulkBoundaryConcentrationVector = reset;
}

void MembraneCellProperty::PerformMembraneSystem(const std::vector<double>& currentBulkConcentration, const std::vector<double>& currentCellConcentration, std::vector<double>& changeBulkConc, std::vector<double>& changeCellConc)
{
    mpMembraneReactionSystem->ReactSystem(currentBulkConcentration, currentCellConcentration, changeBulkConc, changeCellConc);    
}

void MembraneCellProperty::UpdateMembraneReactionSystem(AbstractMembraneReactionSystem* p_system)
{
    mpMembraneReactionSystem = p_system;
    boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
    SetMembraneOdeSolver(p_solver);
}

void MembraneCellProperty::UpdateMembraneOdeSystem(AbstractMembraneOdeSystem* p_ode_system)
{
    mpMembraneOdeSystem = p_ode_system;
}

void MembraneCellProperty::PreparePostDivisionParent(double splitRatio)
{
    // Split any properties that are shared
    for (unsigned i = 0; i < mMembraneConcentrationVector.size(); ++i)
    {
        mMembraneConcentrationVector[i] = splitRatio*mMembraneConcentrationVector[i];
    }
}

void MembraneCellProperty::PreparePostDivisionDaughter(const MembraneCellProperty& parentProperty ,double splitRatio)
{
    // Split any properties that are shared
    std::vector<double> parentConcentrationVector = parentProperty.mMembraneConcentrationVector;

    for (unsigned i = 0; i < parentConcentrationVector.size(); ++i)
    {
        mMembraneConcentrationVector[i] = (1 - splitRatio)*parentConcentrationVector[i];
    }
}

void MembraneCellProperty::InitialiseMembrane(std::vector<std::string> stateNameVector, std::vector<double> concentrationVector)
{
    SetMembraneStateVariableRegister(new StateVariableRegister(stateNameVector));

    SetMembraneConcentrationVector(concentrationVector);
}

void MembraneCellProperty::InitialiseMembrane(StateVariableRegister* p_StateVariableRegister, std::vector<double> concentrationVector)
{
    SetMembraneStateVariableRegister(p_StateVariableRegister);

    SetMembraneConcentrationVector(concentrationVector);
}

double MembraneCellProperty::RetrieveBoundarySourceByStateName(std::string stateName)
{
    // If cell is on boundary add the result of the transport Ode system
    unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);
    return GetExternalCellBoundaryConcentrationByIndex(index);// - mInitBulkBoundaryConcentrationVector[index];
}

double MembraneCellProperty::RetrieveChangeBoundarySourceByStateName(std::string stateName)
{
    // If cell is on boundary add the result of the transport ODE system

    // Check is named state is present on the bulk boundary side
    if (mpBulkStateVariableRegister->IsStateVariablePresent(stateName))
    {
        unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);
        mNumCallsThisReactionStep++;
        return GetChangeExternalCellBoundaryConcentrationByIndex(index);
    }

    // Else return 0 change
    return 0.0;
}

void MembraneCellProperty::AppendInternalCellBoundaryConcentrations(std::vector<double>& rY)
{
    rY.insert(rY.end(), mCellBoundaryConcentrationVector.begin(), mCellBoundaryConcentrationVector.end());
}

void MembraneCellProperty::ReplaceBoundaryStateVariables(std::vector<double>& rY)
{
    // Partition the ODE state variable vector into the constitient internal and bulk state vectors
    unsigned num_bulk_states = mBulkBoundaryConcentrationVector.size();
    for (unsigned i = 0; i < num_bulk_states; ++i)
    {
        mBulkBoundaryConcentrationVector[i] = rY[i];
    }

    unsigned num_cell_states = mCellBoundaryConcentrationVector.size();
    for (unsigned i = 0; i < num_cell_states; ++i)
    {
        mCellBoundaryConcentrationVector[i] = rY[i + num_bulk_states];
    }
}

void MembraneCellProperty::ReplaceChangeBoundaryStateVariables(std::vector<double>& rDY)
{
    // Partition the ODE state variable vector into the constitient internal and bulk state vectors
    unsigned num_cell_states = mCellBoundaryConcentrationVector.size();
 
    AbstractChemistry* p_bulk_chemistry = mpMembraneReactionSystem->GetBulkChemistry();
    std::string this_state = "";
    unsigned num_bulk_reaction_states = p_bulk_chemistry->GetNumberChemicals();

    for (unsigned i = 0; i < num_bulk_reaction_states; ++i)
    {
        mChangeBulkBoundaryConcentrationVector[i] = rDY[i];
    }

    for (unsigned i = 0; i < num_cell_states; ++i)
    {
        mChangeCellBoundaryConcentrationVector[i] = rDY[i + num_bulk_reaction_states];
    }
}

void MembraneCellProperty::ResetReactionCalls()
{
    mNumCallsThisReactionStep = 0;
}

unsigned MembraneCellProperty::GetReactionCalls()
{
    return mNumCallsThisReactionStep;
}

void MembraneCellProperty::SetMembraneExtentBool(bool isExtent)
{
    mIsExtent = isExtent;
}

void MembraneCellProperty::SetMembraneThickness(double thickness)
{
    mMembraneThickness = thickness;
}

void MembraneCellProperty::SetDoubleMembraneBool(bool isDouble)
{
    mIsDoubleMembrane = isDouble;
}

void MembraneCellProperty::SetBulkStateVariableRegister(StateVariableRegister* pStateRegister)
{
    mpBulkStateVariableRegister = pStateRegister;
}

void MembraneCellProperty::SetCellStateVariableRegister(StateVariableRegister* pStateRegister)
{
    mpCellStateVariableRegister = pStateRegister;
}

void MembraneCellProperty::SetConstantCellConcentrationBool(bool isConstantCellConcentration)
{
    mIsConstantCellConcentration = isConstantCellConcentration;
}

void MembraneCellProperty::SetCellBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mCellBoundaryConcentrationVector = concentrationVector;
}

void MembraneCellProperty::SetBulkBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mBulkBoundaryConcentrationVector = concentrationVector;
}

void MembraneCellProperty::SetInitBulkBoundaryConcentrationVector(std::vector<double> concentrationVector)
{
    mInitBulkBoundaryConcentrationVector = concentrationVector;
}

void MembraneCellProperty::SetMembraneOdeSolver(boost::shared_ptr<AbstractIvpOdeSolver> pSolver)
{
    mpMembraneOdeSolver = pSolver;
}
 
void MembraneCellProperty::SetIncludeMembraneOdeInterpolationOnBoundary(bool includeOdeInterpolationOnBoundary)
{
    mIncludeMembraneOdeInterpolationOnBoundary = includeOdeInterpolationOnBoundary;
}

void MembraneCellProperty::SetCellPtr(CellPtr pCell)
{
    mpCell = pCell;
}

bool MembraneCellProperty::GetMembraneExtentBool()
{
    return mIsExtent;
}

double MembraneCellProperty::GetMembraneThickness()
{
    return mMembraneThickness;
}

bool MembraneCellProperty::GetDoubleMembraneBool()
{
    return mIsDoubleMembrane;
}

StateVariableRegister* MembraneCellProperty::GetBulkStateVariableRegister()
{
    return mpBulkStateVariableRegister;
}

StateVariableRegister* MembraneCellProperty::GetCellStateVariableRegister()
{
    return mpCellStateVariableRegister;
}

bool MembraneCellProperty::GetConstantCellConcentrationBool()
{
    return mIsConstantCellConcentration;
}

std::vector<double> MembraneCellProperty::GetInternalCellBoundaryConcentrationVector()
{
    return mCellBoundaryConcentrationVector;
}

double MembraneCellProperty::GetInternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if (index<mpCellStateVariableRegister->GetNumStateVariables())
    {
        return mCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if (index < mpCellStateVariableRegister->GetNumStateVariables())
    {
        return mCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetChangeInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);
    if (index < mpCellStateVariableRegister->GetNumStateVariables())
    {
        return mChangeCellBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

std::vector<double> MembraneCellProperty::GetExternalCellBoundaryConcentrationVector()
{
    return mBulkBoundaryConcentrationVector;
}

double MembraneCellProperty::GetExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if (index<mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetChangeExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if (index<mChangeBulkBoundaryConcentrationVector.size())
    {
        return mChangeBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

double MembraneCellProperty::GetExternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);

    if (index<mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }

    return 0.0; 
}

AbstractMembraneReactionSystem* MembraneCellProperty::GetMembraneReactionSystem()
{
    return mpMembraneReactionSystem;
}

AbstractMembraneOdeSystem* MembraneCellProperty::GetMembraneOdeSystem()
{
    return mpMembraneOdeSystem;
}

boost::shared_ptr<AbstractIvpOdeSolver> MembraneCellProperty::GetMembraneOdeSolver()
{
    return mpMembraneOdeSolver;
}

bool MembraneCellProperty::GetIncludeMembraneOdeInterpolationOnBoundary()
{
    return mIncludeMembraneOdeInterpolationOnBoundary;
}

CellPtr MembraneCellProperty::GetCellPtr()
{
    return mpCell;
}

void MembraneCellProperty::SetMembraneStateVariableRegister(StateVariableRegister* pRegister)
{
    mpMembraneStateVariableRegister = pRegister;
}

void MembraneCellProperty::SetMembraneConcentrationVector(std::vector<double> concentrationVector)
{
    mMembraneConcentrationVector = concentrationVector;
}

StateVariableRegister* MembraneCellProperty::GetMembraneStateVariableRegister()
{
    return mpMembraneStateVariableRegister;
}

std::vector<double> MembraneCellProperty::GetMembraneConcentrationVector()
{
    return mMembraneConcentrationVector;
}

double MembraneCellProperty::GetMembraneConcentrationByIndex(unsigned index)
{
    if (index < mMembraneConcentrationVector.size())
    {
        return mMembraneConcentrationVector[index];
    }
    
    return 0.0;
}

double MembraneCellProperty::GetMembraneConcentrationByName(std::string name)
{
    if (mpMembraneStateVariableRegister->IsStateVariablePresent(name))
    {
        return mMembraneConcentrationVector[mpMembraneStateVariableRegister->RetrieveStateVariableIndex(name)];
    }
    
    return 0.0;
}