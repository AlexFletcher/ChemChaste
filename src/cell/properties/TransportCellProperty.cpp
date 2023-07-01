#include "TransportCellProperty.hpp"

TransportCellProperty::TransportCellProperty()
    : AbstractCellProperty()
{
}

TransportCellProperty::~TransportCellProperty()
{
}

TransportCellProperty::TransportCellProperty(const TransportCellProperty& rOtherTransportCellProperty)
{
    mpBulkStateVariableRegister = rOtherTransportCellProperty.mpBulkStateVariableRegister;
    mpCellStateVariableRegister = rOtherTransportCellProperty.mpCellStateVariableRegister;
    mpTransportReactionSystem = rOtherTransportCellProperty.mpTransportReactionSystem;
    mpTransportOdeSystem = rOtherTransportCellProperty.mpTransportOdeSystem;
    mpSolver = rOtherTransportCellProperty.mpSolver;
    mIsConstantCellConcentration = rOtherTransportCellProperty.mIsConstantCellConcentration;
    mBulkBoundaryConcentrationVector = rOtherTransportCellProperty.mBulkBoundaryConcentrationVector;
    mChangeBulkBoundaryConcentrationVector = rOtherTransportCellProperty.mChangeBulkBoundaryConcentrationVector;
    mCellBoundaryConcentrationVector = rOtherTransportCellProperty.mCellBoundaryConcentrationVector;
    mChangeCellBoundaryConcentrationVector = rOtherTransportCellProperty.mChangeCellBoundaryConcentrationVector;
    mIncludeOdeInterpolationOnBoundary = rOtherTransportCellProperty.mIncludeOdeInterpolationOnBoundary;
}

void TransportCellProperty::SetUp(AbstractTransportReactionSystem* pTransportReactionSystem, CellPtr pCell)
{ 
    UpdateTransportReactionSystem(pTransportReactionSystem);
    SetCellPtr(pCell);

    AbstractTransportOdeSystem* p_ode_system = new AbstractTransportOdeSystem(pTransportReactionSystem);
    UpdateTransportOdeSystem(p_ode_system);

    StateVariableRegister* p_bulk_state_register = new StateVariableRegister(pTransportReactionSystem->GetBulkChemistry()->GetChemicalNames());
    SetBulkStateVariableRegister(p_bulk_state_register);

    StateVariableRegister* p_cell_state_register = new StateVariableRegister(pTransportReactionSystem->GetCellChemistry()->GetChemicalNames());  
    SetCellStateVariableRegister(p_cell_state_register);

    SetUpCellConcentrationVector(p_cell_state_register->GetNumStateVariables());
    SetUpBulkConcentrationVector(p_bulk_state_register->GetNumStateVariables());
    SetUpChangeCellConcentrationVector(p_cell_state_register->GetNumStateVariables());
    SetUpChangeBulkConcentrationVector(p_bulk_state_register->GetNumStateVariables());
}

void TransportCellProperty::UpdateCellConcentrationVector(std::vector<double> cellBoundaryConcentrationVector)
{
    mCellBoundaryConcentrationVector = cellBoundaryConcentrationVector;
}

void TransportCellProperty::UpdateBulkConcentrationVector(std::vector<double> bulkBoundaryConcentrationVector)
{
    mBulkBoundaryConcentrationVector = bulkBoundaryConcentrationVector;
}

void TransportCellProperty::SetUpCellConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables, 0.0);
    mCellBoundaryConcentrationVector = reset;
}

void TransportCellProperty::SetUpBulkConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables, 0.0);
    mBulkBoundaryConcentrationVector = reset;
}

void TransportCellProperty::SetUpChangeCellConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables, 0.0);
    mChangeCellBoundaryConcentrationVector = reset;
}

void TransportCellProperty::SetUpChangeBulkConcentrationVector(unsigned numStateVariables)
{
    std::vector<double> reset(numStateVariables, 0.0);
    mChangeBulkBoundaryConcentrationVector = reset;
}

void TransportCellProperty::PerformTransportSystem(
    const std::vector<double>& rCurrentBulkConcentration, 
    const std::vector<double>& rCurrentCellConcentration, 
    std::vector<double>& rChangeBulkConcentration, 
    std::vector<double>& rChangeCellConcentration)
{
    mpTransportReactionSystem->ReactSystem(rCurrentBulkConcentration, 
                                           rCurrentCellConcentration, 
                                           rChangeBulkConcentration, 
                                           rChangeCellConcentration);
}

void TransportCellProperty::UpdateTransportReactionSystem(AbstractTransportReactionSystem* pReactionSystem)
{
    mpTransportReactionSystem = pReactionSystem;
    boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
    SetTransportOdeSolver(p_solver);
}   

void TransportCellProperty::UpdateTransportOdeSystem(AbstractTransportOdeSystem* pOdeSystem)
{
    mpTransportOdeSystem = pOdeSystem;
}

void TransportCellProperty::PreparePostDivisionParent(double splitRatio)
{
}
    
void TransportCellProperty::PreparePostDivisionDaughter(const TransportCellProperty& rParentProperty, double splitRatio)
{
}

double TransportCellProperty::RetrieveBoundarySourceByStateName(std::string stateName)
{
    unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);
    return GetExternalCellBoundaryConcentrationByIndex(index);
}

double TransportCellProperty::RetrieveChangeBoundarySourceByStateName(std::string stateName)
{
    if (mpBulkStateVariableRegister->IsStateVariablePresent(stateName))
    {
        unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);
        mNumCallsThisReactionStep++;
        return GetChangeExternalCellBoundaryConcentrationByIndex(index);
    }
    return 0.0;
}

void TransportCellProperty::AppendInternalCellBoundaryConcentrations(std::vector<double>& rY)
{
    rY.insert(rY.end(), mCellBoundaryConcentrationVector.begin(), mCellBoundaryConcentrationVector.end());
}

void TransportCellProperty::ReplaceBoundaryStateVariables(std::vector<double>& rY)
{
    /*
     * Partition the ODE state variable vector into the constituent internal and 
     * bulk state vectors
     */
    unsigned num_bulk_states = mBulkBoundaryConcentrationVector.size();
    unsigned num_cell_states = mCellBoundaryConcentrationVector.size();

    for (unsigned i = 0; i < num_bulk_states; ++i)
    {
        mBulkBoundaryConcentrationVector[i] = rY[i];
    }

    for (unsigned i = 0; i < num_cell_states; ++i)
    {
        mCellBoundaryConcentrationVector[i] = rY[i + num_bulk_states];
    }

}

void TransportCellProperty::ReplaceChangeBoundaryStateVariables(std::vector<double>& rDY)
{
    /*
     * Partition the ODE state variable vector into the constituent internal and 
     * bulk state vectors
     */
    unsigned num_bulk_states = mBulkBoundaryConcentrationVector.size();
    unsigned num_cell_states = mCellBoundaryConcentrationVector.size();

    for (unsigned i = 0; i < num_bulk_states; ++i)
    {
        mChangeBulkBoundaryConcentrationVector[i] = rDY[i];
    }

    for (unsigned i=0; i<num_cell_states; ++i)
    {

        mChangeCellBoundaryConcentrationVector[i] = rDY[i + num_bulk_states];
    }
}

void TransportCellProperty::ResetReactionCalls()
{
    mNumCallsThisReactionStep = 0;
}

unsigned TransportCellProperty::GetReactionCalls()
{
    return mNumCallsThisReactionStep;
}

void TransportCellProperty::SetBulkStateVariableRegister(StateVariableRegister* pStateRegister)
{
    mpBulkStateVariableRegister = pStateRegister;
}

void TransportCellProperty::SetCellStateVariableRegister(StateVariableRegister* pStateRegister)
{
    mpCellStateVariableRegister = pStateRegister;
}

void TransportCellProperty::SetConstantCellConcentrationBool(bool isConstantCellConcentration)
{
    mIsConstantCellConcentration = isConstantCellConcentration;
}

void TransportCellProperty::SetCellBoundaryConcentrationVector(std::vector<double> vector)
{
    mCellBoundaryConcentrationVector = vector;
}

void TransportCellProperty::SetBulkBoundaryConcentrationVector(std::vector<double> vector)
{
    mBulkBoundaryConcentrationVector = vector;
}

void TransportCellProperty::SetInitBulkBoundaryConcentrationVector(std::vector<double> vector)
{
    mInitBulkBoundaryConcentrationVector = vector;
}

void TransportCellProperty::SetTransportOdeSolver(boost::shared_ptr<AbstractIvpOdeSolver> pSolver)
{
    mpSolver = pSolver;
}

void TransportCellProperty::SetIncludeOdeInterpolationOnBoundary(bool includeOdeInterpolationOnBoundary)
{
    mIncludeOdeInterpolationOnBoundary = includeOdeInterpolationOnBoundary;
}

void TransportCellProperty::SetCellPtr(CellPtr pCell)
{
    mpCell = pCell;
}

StateVariableRegister* TransportCellProperty::GetBulkStateVariableRegister()
{
    return mpBulkStateVariableRegister;
}

StateVariableRegister* TransportCellProperty::GetCellStateVariableRegister()
{
    return mpCellStateVariableRegister;
}

bool TransportCellProperty::GetConstantCellConcentrationBool()
{
    return mIsConstantCellConcentration;
}

std::vector<double> TransportCellProperty::GetInternalCellBoundaryConcentrationVector()
{
    return mCellBoundaryConcentrationVector;
}

double TransportCellProperty::GetInternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if (index <mpCellStateVariableRegister->GetNumStateVariables())
    { 
        return mCellBoundaryConcentrationVector[index];
    }
    return 0.0; 
}

double TransportCellProperty::GetInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);
    if (index < mpCellStateVariableRegister->GetNumStateVariables())
    {
        return mCellBoundaryConcentrationVector[index];
    }
    return 0.0; 
}

double TransportCellProperty::GetChangeInternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpCellStateVariableRegister->RetrieveStateVariableIndex(stateName);
    if (index < mpCellStateVariableRegister->GetNumStateVariables())
    {
        return mChangeCellBoundaryConcentrationVector[index];
    }
    return 0.0; 
}

std::vector<double> TransportCellProperty::GetExternalCellBoundaryConcentrationVector()
{
    return mBulkBoundaryConcentrationVector;
}

double TransportCellProperty::GetExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if (index < mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }
    return 0.0; 
}

double TransportCellProperty::GetChangeExternalCellBoundaryConcentrationByIndex(unsigned index)
{
    if (index < mChangeBulkBoundaryConcentrationVector.size())
    {
        return mChangeBulkBoundaryConcentrationVector[index];
    }
    return 0.0; 
}

double TransportCellProperty::GetExternalCellBoundaryConcentrationByName(std::string stateName)
{
    unsigned index = mpBulkStateVariableRegister->RetrieveStateVariableIndex(stateName);
    if (index < mBulkBoundaryConcentrationVector.size())
    {
        return mBulkBoundaryConcentrationVector[index];
    }
    return 0.0; 
}

AbstractTransportReactionSystem* TransportCellProperty::GetTransportReactionSystem()
{
    return mpTransportReactionSystem;
}

AbstractTransportOdeSystem* TransportCellProperty::GetTransportOdeSystem()
{
    return mpTransportOdeSystem;
}

boost::shared_ptr<AbstractIvpOdeSolver> TransportCellProperty::GetTransportOdeSolver()
{
    return mpSolver;
}

bool TransportCellProperty::GetIncludeOdeInterpolationOnBoundary()
{
    return mIncludeOdeInterpolationOnBoundary;
}

CellPtr TransportCellProperty::GetCellPtr()
{
    return mpCell;
}