#ifndef MEMBRANECELLPROPERTY_HPP
#define MEMBRANECELLPROPERTY_HPP

//general includes
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

// chaste includes
#include "Cell.hpp"
#include "AbstractCellProperty.hpp"
#include "StateVariableRegister.hpp"
#include "AbstractMembraneReactionSystem.hpp"
#include "AbstractMembraneOdeSystem.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"

/**
 * \todo Document this class.
 */
class MembraneCellProperty : public AbstractCellProperty
{
protected:

    // CellPtr to a given cell for access to cell data and for properties that relate to the membrane ode
    CellPtr mpCell;

    // for the calculation of gradients
    bool mIsExtent = false;

    double mMembraneThickness = 0.0;

    bool mIsDoubleMembrane = false;

    // membrane constituent ; lipids, proteins, channels etc
    // register for the state variables present in the membrane region
    StateVariableRegister* mpMembraneStateVariableRegister;

    // concentration vector of the membrane, if double memebrane
    std::vector<double> mMembraneConcentrationVector; 

    // register for variables that cross the membrane; either side
    StateVariableRegister* mpBulkStateVariableRegister;

    StateVariableRegister* mpCellStateVariableRegister;

    AbstractMembraneReactionSystem* mpMembraneReactionSystem;

    // membrane reaction system Ode corresponding to the mpMembraneReactionSystem
    AbstractMembraneOdeSystem* mpMembraneOdeSystem;

    boost::shared_ptr<AbstractIvpOdeSolver> mpMembraneOdeSolver;

    // concentration of state variables at the internal boundary of the cell if this is to remain constant
    bool mIsConstantCellConcentration = false; // for use in averaged source for example?

    std::vector<double> mBulkBoundaryConcentrationVector;

    std::vector<double> mChangeBulkBoundaryConcentrationVector;

    std::vector<double> mCellBoundaryConcentrationVector;

    std::vector<double> mChangeCellBoundaryConcentrationVector;

    bool mIncludeMembraneOdeInterpolationOnBoundary=false;

    std::vector<double> mInitBulkBoundaryConcentrationVector;

    unsigned mNumCallsThisReactionStep = 0;

public:

    MembraneCellProperty();

    virtual ~MembraneCellProperty();

    MembraneCellProperty(const MembraneCellProperty&);

    virtual void SetUp(AbstractMembraneReactionSystem*, CellPtr);

    virtual void UpdateMembraneConcentrationVector(std::vector<double>);

    virtual void UpdateCellConcentrationVector(std::vector<double>);

    virtual void UpdateBulkConcentrationVector(std::vector<double>);

    virtual void PerformMembraneSystem(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, std::vector<double>&);

    virtual void UpdateMembraneReactionSystem(AbstractMembraneReactionSystem*);

    virtual void UpdateMembraneOdeSystem(AbstractMembraneOdeSystem*);

    virtual void PreparePostDivisionParent(double);
    
    virtual void PreparePostDivisionDaughter(const MembraneCellProperty& ,double);

    virtual void InitialiseMembrane(std::vector<std::string>, std::vector<double>);

    virtual void InitialiseMembrane(StateVariableRegister*, std::vector<double>);

    void SetUpCellConcentrationVector(unsigned);

    void SetUpBulkConcentrationVector(unsigned);

    void SetUpChangeCellConcentrationVector(unsigned);

    void SetUpChangeBulkConcentrationVector(unsigned);

    double RetrieveBoundarySourceByStateName(std::string);

    double RetrieveChangeBoundarySourceByStateName(std::string);
    
    void AppendInternalCellBoundaryConcentrations(std::vector<double>&);

    void ReplaceBoundaryStateVariables(std::vector<double>&);

    void ReplaceChangeBoundaryStateVariables(std::vector<double>&);

    void ResetReactionCalls();

    unsigned GetReactionCalls();

    void SetMembraneExtentBool(bool);

    void SetMembraneThickness(double);

    void SetDoubleMembraneBool(bool);

    void SetBulkStateVariableRegister(StateVariableRegister*);

    void SetCellStateVariableRegister(StateVariableRegister*);

    void SetConstantCellConcentrationBool(bool);

    void SetCellBoundaryConcentrationVector(std::vector<double>);

    void SetBulkBoundaryConcentrationVector(std::vector<double>);

    void SetInitBulkBoundaryConcentrationVector(std::vector<double>);

    void SetMembraneOdeSolver(boost::shared_ptr<AbstractIvpOdeSolver>);

    void SetIncludeMembraneOdeInterpolationOnBoundary(bool);

    void SetCellPtr(CellPtr);

    bool GetMembraneExtentBool();

    double GetMembraneThickness();

    bool GetDoubleMembraneBool();

    StateVariableRegister* GetBulkStateVariableRegister();

    StateVariableRegister* GetCellStateVariableRegister();

    bool GetConstantCellConcentrationBool();

    std::vector<double> GetInternalCellBoundaryConcentrationVector();

    double GetInternalCellBoundaryConcentrationByIndex(unsigned);

    double GetInternalCellBoundaryConcentrationByName(std::string);

    double GetChangeInternalCellBoundaryConcentrationByName(std::string);

    std::vector<double> GetExternalCellBoundaryConcentrationVector();

    double GetExternalCellBoundaryConcentrationByIndex(unsigned);

    double GetChangeExternalCellBoundaryConcentrationByIndex(unsigned);

    double GetExternalCellBoundaryConcentrationByName(std::string);

    AbstractMembraneReactionSystem* GetMembraneReactionSystem();

    AbstractMembraneOdeSystem* GetMembraneOdeSystem();

    boost::shared_ptr<AbstractIvpOdeSolver> GetMembraneOdeSolver();

    bool GetIncludeMembraneOdeInterpolationOnBoundary();

    CellPtr GetCellPtr();

    void SetMembraneStateVariableRegister(StateVariableRegister*);

    void SetMembraneConcentrationVector(std::vector<double>);

    StateVariableRegister* GetMembraneStateVariableRegister();

    std::vector<double> GetMembraneConcentrationVector();

    double GetMembraneConcentrationByIndex(unsigned);

    double GetMembraneConcentrationByName(std::string);
};

#endif