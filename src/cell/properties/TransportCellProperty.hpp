#ifndef TRANSPORTCELLPROPERTY_HPP
#define TRANSPORTCELLPROPERTY_HPP

#include <vector>
#include <boost/shared_ptr.hpp>
#include "Cell.hpp"
#include "AbstractCellProperty.hpp"
#include "StateVariableRegister.hpp"
#include "AbstractTransportReactionSystem.hpp"
#include "AbstractTransportOdeSystem.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ChastePoint.hpp"

/**
 * \todo Document class.
 */
class TransportCellProperty : public AbstractCellProperty
{
protected:

    // CellPtr to a given cell for access to cell data and for properties that relate ot the transport ode
    CellPtr mpCell;

    // register for variables that cross the membrane; either side
    StateVariableRegister* mpBulkStateVariableRegister;
    StateVariableRegister* mpCellStateVariableRegister;

    // transport reaction system for the passage of state variables across the boundary
    AbstractTransportReactionSystem* mpTransportReactionSystem;

    // transport reaction system Ode corresponding to the mpTransportReactionSystem
    AbstractTransportOdeSystem* mpTransportOdeSystem;

    // transport Ode solver
    boost::shared_ptr<AbstractIvpOdeSolver> mpSolver;

    // concentration of state variables at the internal boundary of the cell if this is to remain constant
    bool mIsConstantCellConcentration = false; // for use in averaged source for example?

    std::vector<double> mBulkBoundaryConcentrationVector;

    std::vector<double> mChangeBulkBoundaryConcentrationVector;

    std::vector<double> mCellBoundaryConcentrationVector;

    std::vector<double> mChangeCellBoundaryConcentrationVector;

    bool mIncludeOdeInterpolationOnBoundary = false;
    
    std::vector<double> mInitBulkBoundaryConcentrationVector;

    unsigned mNumberOCalls_this_reaction_step = 0;

public:

    /**
     * Default constructor.
     */
    TransportCellProperty();

    /**
     * Destructor.
     */
    virtual ~TransportCellProperty();

    /**
     * Copy constructor.
     * \todo Document method.
     * 
     * @param rOtherTransportCellProperty
     */
    TransportCellProperty(
        const TransportCellProperty& rOtherTransportCellProperty);

    /**
     * \todo Document method.
     * 
     * @param pTransportReactionSystem
     * @param pCell
     */
    virtual void SetUp(
        AbstractTransportReactionSystem* pTransportReactionSystem, 
        CellPtr pCell);

    /**
     * \todo Document method.
     * 
     * @param cellBoundaryConcentrationVector
     */
    virtual void UpdateCellConcentrationVector(
        std::vector<double> cellBoundaryConcentrationVector);

    /**
     * \todo Document method.
     * 
     * @param bulkBoundaryConcentrationVector
     */
    virtual void UpdateBulkConcentrationVector(
        std::vector<double> bulkBoundaryConcentrationVector);

    /**
     * \todo Document method.
     * 
     * @param numStateVariables
     */
    void SetUpCellConcentrationVector(unsigned numStateVariables);

    /**
     * \todo Document method.
     * 
     * @param numStateVariables
     */
    void SetUpBulkConcentrationVector(unsigned numStateVariables);

    /**
     * \todo Document method.
     * 
     * @param numStateVariables
     */
    void SetUpChangeCellConcentrationVector(unsigned numStateVariables);

    /**
     * \todo Document method.
     * 
     * @param numStateVariables
     */
    void SetUpChangeBulkConcentrationVector(unsigned numStateVariables);

    /**
     * \todo Document method.
     * 
     * @param rCurrentBulkConcentration
     * @param rCurrentCellConcentration
     * @param rChangeBulkConcentration
     * @param rChangeCellConcentration
     */
    virtual void PerformTransportSystem(
        const std::vector<double>& rCurrentBulkConcentration, 
        const std::vector<double>& rCurrentCellConcentration, 
        std::vector<double>& rChangeBulkConcentration, 
        std::vector<double>& rChangeCellConcentration);

    /**
     * \todo Document method.
     * 
     * @param pReactionSystem
     */
    virtual void UpdateTransportReactionSystem(
        AbstractTransportReactionSystem* pReactionSystem);

    /**
     * \todo Document method.
     * 
     * @param pOdeSystem
     */
    virtual void UpdateTransportOdeSystem(
        AbstractTransportOdeSystem* pOdeSystem);

    /**
     * \todo Document method.
     * 
     * @param splitRatio
     */
    virtual void PreparePostDivisionParent(double splitRatio);

    /**
     * \todo Document method.
     * 
     * @param rParentProperty
     * @param splitRatio
     */
    virtual void PreparePostDivisionDaughter(
        const TransportCellProperty& rParentProperty, 
        double splitRatio);

    /**
     * \todo Document method.
     * 
     * @param stateName
     * @return 
     */
    double RetrieveBoundarySourceByStateName(std::string stateName);

    /**
     * \todo Document method.
     * 
     * @param stateName
     * @return 
     */
    double RetrieveChangeBoundarySourceByStateName(std::string stateName);

    /**
     * \todo Document method.
     * 
     * @param rY
     */    
    void AppendInternalCellBoundaryConcentrations(std::vector<double>& rY);
    
    /**
     * \todo Document method.
     * 
     * @param rY
     */    
    void ReplaceBoundaryStateVariables(std::vector<double>& rY);
    
    /**
     * \todo Document method.
     * 
     * @param rDY
     */    
    void ReplaceChangeBoundaryStateVariables(std::vector<double>& rDY);
    
    /**
     * \todo Document method.
     * 
     * @param pStateRegister
     */    
    void SetBulkStateVariableRegister(StateVariableRegister* pStateRegister);
    
    /**
     * \todo Document method.
     * 
     * @param pStateRegister
     */    
    void SetCellStateVariableRegister(StateVariableRegister* pStateRegister);
    
    /**
     * \todo Document method.
     * 
     * @param isConstantCellConcentration
     */    
    void SetConstantCellConcentrationBool(bool isConstantCellConcentration);

    /**
     * \todo Document method.
     * 
     * @param vector
     */    
    void SetCellBoundaryConcentrationVector(std::vector<double> vector);

    /**
     * \todo Document method.
     * 
     * @param vector
     */    
    void SetBulkBoundaryConcentrationVector(std::vector<double> vector);
    
    /**
     * \todo Document method.
     * 
     * @param vector
     */    
    void SetInitBulkBoundaryConcentrationVector(std::vector<double> vector);
    
    /**
     * \todo Document method.
     * 
     * @param pSolver
     */    
    void SetTransportOdeSolver(boost::shared_ptr<AbstractIvpOdeSolver> pSolver);
    
    /**
     * \todo Document method.
     * 
     * @param includeOdeInterpolationOnBoundary
     */    
    void SetIncludeOdeInterpolationOnBoundary(
        bool includeOdeInterpolationOnBoundary);

    /**
     * \todo Document method.
     * 
     * @param pCell
     */    
    void SetCellPtr(CellPtr pCell);
    
    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    double GetInternalCellBoundaryConcentrationByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @param stateName
     * 
     * @return 
     */
    double GetInternalCellBoundaryConcentrationByName(std::string stateName);

    /**
     * \todo Document method.
     * 
     * @param stateName
     * 
     * @return 
     */
    double GetChangeInternalCellBoundaryConcentrationByName(
        std::string stateName);

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    double GetExternalCellBoundaryConcentrationByIndex(unsigned index); 

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    double GetChangeExternalCellBoundaryConcentrationByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @param stateName
     * 
     * @return 
     */
    double GetExternalCellBoundaryConcentrationByName(std::string stateName);
};

#endif /* TRANSPORTCELLPROPERTY_HPP_ */