#ifndef CHEMICALCELLFROMFILE_HPP_
#define CHEMICALCELLFROMFILE_HPP_

#include "ChemicalCell.hpp"
#include "SimpleChemicalThresholdCellCycleFromFile.hpp"
#include "AbstractSrnModel.hpp"
#include "ChemicalSrnFromFile.hpp"
#include "AbstractCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "InitialCellConditionsFromFile.hpp"
#include "TransportCellPropertyFromFile.hpp"
#include "MembraneCellPropertyFromFile.hpp"
#include "CellAnalyticsPropertyFromCellId.hpp"

/**
 * \todo Document class.
 */
class ChemicalCellFromFile
{
protected:

    CellPtr mpCell;

    CellPropertyCollection mPropertyCollection;

    boost::shared_ptr<ChemicalCellProperty> mpCellChemicalProperty;

    boost::shared_ptr<MembraneCellProperty> mpCellMembraneProperty;

    boost::shared_ptr<TransportCellProperty> mpCellTransportProperty;

    boost::shared_ptr<CellAnalyticsProperty> mpCellAnalyticsProperty;

    std::string mCellCycleFilename;

    bool mIsCellCycleSet = false;

    std::string mSrnFilename;

    bool mIsSRNSet = false;

    std::string mInitialConditionsFilename;

    bool mIsInitConditionsSet = false;

    std::string mTransportPropertyFilename;

    bool mIsTransportPropertySet = false;

    std::string mMembranePropertyFilename;

    bool mIsMembranePropertySet = false;

    unsigned mCellId;

    bool mIsCellIdSet = false;

    StateVariableRegister* mpFullChemicalStateRegister; 

    SimpleChemicalThresholdCellCycleModel*  mpSimpleChemicalThresholdCellCycleModel;

    ChemicalSrnModel*   mpChemicalSrnModel;

public:

    ChemicalCellFromFile(std::string cellCycleFilename="", 
                         std::string srnFilename="",
                         std::string initialConditionFilename="",
                         std::string transportPropertyFilename="",
                         std::string membranePropertyFilename="",
                         unsigned cellId=0,
                         bool isCellIdSet=false);

    virtual ~ChemicalCellFromFile()
    {
    };

    void SetUpSRNandCellCycle();

    void SetUpCellObject();

    void SetUpCellProperties();

    void SetUpCellInitialConditions(CellPtr, std::vector<std::string>, std::vector<double>);

    CellPtr GetCellPtr();

    CellPropertyCollection GetCellPropertyCollection();
    
    boost::shared_ptr<ChemicalCellProperty> GetChemicalCellProperty();

    boost::shared_ptr<MembraneCellProperty> GetMembraneCellProperty();

    boost::shared_ptr<TransportCellProperty> GetTransportCellProperty();
    
    boost::shared_ptr<CellAnalyticsProperty> GetCellAnalyticsProperty();

    ChemicalSrnModel* GetChemicalSrnModel();

    SimpleChemicalThresholdCellCycleModel* GetChemicalCellCycleModel();

    std::string GetCellCycleFilename();

    bool GetIsCellCycleSet();

    std::string GetSrnFilename();

    bool GetIsSrnSet();

    std::string GetInitialConditionsFilename();

    bool GetInitConditionsSet();

    std::string GetTransportPropertyFilename();

    bool GetIsTransportPropertySet();

    std::string GetMembranePropertyFilename();

    bool GetIsMembranePropertySet();

    unsigned GetCellId();

    bool GetIsCellIdSet();

    StateVariableRegister* GetFullChemicalStateRegister();

    std::vector<std::string> GetFullChemicalNamesVector();

    void SetCellPtr(CellPtr);

    void SetCellPropertyCollection(CellPropertyCollection);

    void SetChemicalCellProperty(boost::shared_ptr<ChemicalCellProperty>);

    void SetMembraneCellProperty(boost::shared_ptr<MembraneCellProperty>);

    void SetTransportCellProperty(boost::shared_ptr<TransportCellProperty>);

    void SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty>);

    void SetChemicalSrnModel(ChemicalSrnModel*);

    void SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel*);

    void SetCellCycleFilename(std::string);

    void SetSrnFilename(std::string);

    void SetInitialConditionsFilename(std::string);

    void SetTransportPropertyFilename(std::string);

    void SetMembranePropertyFilename(std::string);

    void SetFullChemicalStateRegister(StateVariableRegister*);
};

#endif /* CHEMICALCELLFROMFILE_HPP_ */