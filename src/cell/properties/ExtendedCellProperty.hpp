#ifndef EXTENDEDCELLPROPERTY_HPP_
#define EXTENDEDCELLPROPERTY_HPP_

#include <vector>
#include <string>
#include <math.h>
#include <boost/shared_ptr.hpp>

#include "ChemicalCellProperty.hpp"

// base case assume spherical geometry

/**
 * \todo Document class.
 */
template<unsigned SPACE_DIM>
class ExtendedCellProperty : public ChemicalCellProperty
{
private:

    // inherited:
    //StateVariableRegister* mpStateVariableRegister;
    //std::vector<double> mConcentrationVector;

    using ChemicalCellProperty::UpdateCellConcentrationVector;
    using ChemicalCellProperty::InitialiseCell;
protected:

    // base case; one value equal to the radius
    std::vector<double> mCellDimensions = {1.0}; 
    
    // vector of minimal dimensions before dimension is assumed to be null.
    std::vector<double> mMinimalDimensions = {1e-5};

    double mCellVolume = 0.0;

    double mVolumePerFineMeshElement = 0.0;

    bool mIsCellVolumeDefined = false;

    bool mIsMeshElementVolumeDefined = false;

    unsigned mTotalNumberMeshVoxels = 0;

    unsigned mNumberMeshVoxelsInCell = 0;

    unsigned mNumberMeshVoxelsOnBoundary = 0;

    // is the cell permitted to grow and/or change in dimensions
    bool mIsCellDynamic = false;

    std::vector<double> mChangeInCellDimensions = {0.0};

    std::vector<double> mMeshDomainScale = std::vector<double>();

    // option switchs for interpolating into the bulk FeMesh
    bool mIncludeOdeInterpolationInCell = false;

    bool mAverageInternalCellStates = true;

    // member variables for the cell boundary voxels

    std::vector<std::vector<double>> mVectorOfExternalBoundaryStateVariables;

    std::vector<std::vector<double>> mVectorOfInternalBoundaryStateVariables;

    std::vector<ChastePoint<SPACE_DIM>> mVectorOfBoundaryLocations;

    // double in [0.0,1.0] relating the proportion of cell data avaliable for transport processes
    double mTransportConcentrationAvaliability = 1.0; 

    // resultant cell boundary concentration vector; total contribution to the cell concentration for all the boundary voxels
    std::vector<double> mNextTimestepConcentrationVector;

public:

    ExtendedCellProperty();

    virtual ~ExtendedCellProperty();

    virtual void SetUpExtendedCell(StateVariableRegister*, std::vector<double>, std::vector<double>);

    virtual void SetUpExtendedCell(StateVariableRegister*, std::vector<double>, double radius=0.0);

    virtual void UpdateExtendedCell();

    virtual void InitialiseCell(StateVariableRegister*, std::vector<double>);

    virtual void CellGrowth();

    virtual bool IsPointInCell(const ChastePoint<SPACE_DIM>&, ChastePoint<SPACE_DIM>&);

    virtual bool IsPointOnCellBoundary(const ChastePoint<SPACE_DIM>&, ChastePoint<SPACE_DIM>&);

    virtual bool CheckCellBoundary(const ChastePoint<SPACE_DIM>&, ChastePoint<SPACE_DIM>&);

    virtual void CalculateVolumeOfCell();

    virtual void CalculateNumCellMeshElements(ChastePoint<SPACE_DIM>);

    virtual void UpdateCellConcentrationVector(std::vector<double>&);

    virtual double TransportConcentrationAvaliabilityThisState(std::string state_name="", unsigned state_index=0);

    void StepThroughDirectionalXBoundaryVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughDirectionalXInternalVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughXVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&);

    void StepThroughDirectionalXYBoundaryVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughDirectionalXYInternalVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughXYVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&);

    void StepThroughDirectionalXYZBoundaryVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughDirectionalXYZInternalVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&, unsigned);

    void StepThroughXYZVoxels(ChastePoint<SPACE_DIM>&,ChastePoint<SPACE_DIM>&);

    double RetrieveBoundarySourceByStateName(std::string, ChastePoint<SPACE_DIM>);

    double RetrieveInternalCellSourceByStateName(std::string);

    double RetrieveInternalCellSourceByStateName(std::string, ChastePoint<SPACE_DIM>);

    double RetrieveBoundaryCellSourceByStateNameAndLocation(std::string, ChastePoint<SPACE_DIM>);

    double RetrieveBoundaryInternalCellSourceByStateNameAndLocation(std::string, ChastePoint<SPACE_DIM>);

    std::vector<double>& RetrieveBoundaryCellSourceByLocation(ChastePoint<SPACE_DIM>);

    std::vector<double>& RetrieveInternalBoundaryCellSourceByLocation(ChastePoint<SPACE_DIM>);

    bool CheckChastePointsForEquality(ChastePoint<SPACE_DIM>, ChastePoint<SPACE_DIM>);

    void ResetNextTimestepConcentrationVector();

    void ResetNextTimestepConcentrationVector(unsigned);

    void ResetVectorOfBoundaryStateVariables();

    void ResetVectorOfBoundaryLocations();

    void ResetVectorOfInternalBoundaryStateVariables();

    void AddToVectorOfBoundaryStateVariables(std::vector<double>);

    void AddToVectorOfInternalBoundaryStateVariables(std::vector<double>);

    void AddToVectorOfBoundaryLocations(ChastePoint<SPACE_DIM>);

    void AppendInternalCellBoundaryConcentrations(std::vector<double>&, unsigned);

    void ReplaceBoundaryStateVariables(unsigned, std::vector<double>&);

    void RecordLocationAndStateVariable(ChastePoint<SPACE_DIM>, std::vector<double>);

    void SetCellRadius(double);

    void SetCellDimensions(std::vector<double>);

    void SetMinimalCellDimensions(std::vector<double>);

    void SetCellVolume(double);

    void SetIsCellDynamic(bool);

    void SetMeshDomainScale(std::vector<double>);

    void SetTotalNumberMeshVoxels(unsigned);

    void SetNumberMeshVoxelsInCell(unsigned);

    void SetNumberMeshVoxelsOnBoundary(unsigned);

    void SetChangeInCellDimensions(std::vector<double>);

    void SetIncludeOdeInterpolationInCell(bool);

    void SetAverageInternalCellStates(bool);

    void SetVectorOfBoundaryStateVariables(std::vector<std::vector<double> >);

    void SetVectorOfBoundaryStateVariablesByIndex(unsigned,std::vector<double>);

    void SetVectorOfInternalBoundaryStateVariables(std::vector<std::vector<double> >);

    void SetVectorOfInternalBoundaryStateVariablesByIndex(unsigned,std::vector<double>);

    void SetVectorOfBoundaryLocations(std::vector<ChastePoint<SPACE_DIM> >);

    void SetVectorOfBoundaryLocationsByIndex(unsigned,ChastePoint<SPACE_DIM>);

    void SetNextTimestepConcentrationVector(std::vector<double>);

    double GetCellRadius();

    std::vector<double> GetCellDimensions();

    std::vector<double> GetMinimalCellDimensions();

    double GetCellVolume();

    bool GetIsCellVolumeDefined();

    bool GetIsMeshElementVolumeDefined();

    bool GetIsCellDynamic();

    std::vector<double> GetMeshDomainScale();

    unsigned GetTotalNumberMeshVoxels();

    unsigned GetNumberMeshVoxelsInCell();

    unsigned GetNumberMeshVoxelsOnBoundary();
    
    std::vector<double> GetChangeInCellDimensions();

    bool GetIncludeOdeInterpolationInCell();

    bool GetAverageInternalCellStates();

    std::vector<std::vector<double>> GetVectorOfBoundaryStateVariables();

    std::vector<double> GetVectorOfBoundaryStateVariablesByLocationIndex(unsigned);

    std::vector<std::vector<double>> GetVectorOfInternalBoundaryStateVariables();

    std::vector<double> GetVectorOfInternalBoundaryStateVariablesByLocationIndex(unsigned);

    std::vector<ChastePoint<SPACE_DIM>> GetVectorOfBoundaryLocations();

    ChastePoint<SPACE_DIM> GetVectorOfBoundaryLocationsByIndex(unsigned);

    std::vector<double> GetNextTimestepConcentrationVector();
};

#endif /* EXTENDEDCELLPROPERTY_HPP_ */