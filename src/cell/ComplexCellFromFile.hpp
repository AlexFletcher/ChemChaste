#ifndef COMPLEXCELLFROMFILE_HPP_
#define COMPLEXCELLFROMFILE_HPP_

#include "ComplexCell.hpp"
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
class ComplexCellFromFile
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

    std::string mDivisionRulesFilename;

    bool mIsDivisionRulesSet = false;

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

    ComplexCellFromFile(    std::string cellCycleFilename="",
                            std::string divisionRulesFilename="", 
                            std::string srnFilename="",
                            std::string initialConditionFilename="",
                            std::string transportPropertyFilename="",
                            std::string membranePropertyFilename="",
                            unsigned cellId =0,
                            bool isCellIdSet = false);

    virtual ~ComplexCellFromFile()
    {
    };

    void SetUpSRNandCellCycle();

    void SetUpCellObject();

    void SetUpCellProperties();

    void SetUpCellInitialConditions(CellPtr, std::vector<std::string>, std::vector<double>);

    void SetUpCellDivisionRules(CellPtr);

    std::vector<std::string> parseMatrixLineString(std::string);

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

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

    std::string GetDivisionRulesFilename();

    bool GetIsDivisionRulesSet();

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

    void SetDivisionRulesFilename(std::string);

    void SetSrnFilename(std::string);

    void SetInitialConditionsFilename(std::string);

    void SetTransportPropertyFilename(std::string);

    void SetMembranePropertyFilename(std::string);

    void SetFullChemicalStateRegister(StateVariableRegister*);
};

ComplexCellFromFile::ComplexCellFromFile(
                        std::string cellCycleFilename,
                        std::string divisionRulesFilename,
                        std::string srnFilename,
                        std::string initialConditionFilename,
                        std::string transportPropertyFilename,
                        std::string membranePropertyFilename,
                        unsigned cellId,
                        bool isCellIdSet)
    : mCellCycleFilename(cellCycleFilename),
      mDivisionRulesFilename(divisionRulesFilename),
      mSrnFilename(srnFilename),
      mInitialConditionsFilename(initialConditionFilename),
      mTransportPropertyFilename(transportPropertyFilename),
      mMembranePropertyFilename(membranePropertyFilename),
      mCellId(cellId),
      mIsCellIdSet(isCellIdSet)
{
    if (mCellCycleFilename != "")
    {
        mIsCellCycleSet = true;
    }
    
    if (mDivisionRulesFilename != "")
    {
        mIsDivisionRulesSet = true;
    }

    if (mSrnFilename != "")
    {
        mIsSRNSet = true;
    }

    if (mInitialConditionsFilename != "")
    {
        mIsInitConditionsSet = true;
        boost::shared_ptr<ChemicalCellProperty> p_cell_chemical(new ChemicalCellProperty());

        SetChemicalCellProperty(p_cell_chemical);
    }

    if (mTransportPropertyFilename != "")
    {
        mIsTransportPropertySet = true;
        boost::shared_ptr<TransportCellProperty> p_cell_transport(new TransportCellProperty());

        SetTransportCellProperty(p_cell_transport);
    }

    if (mMembranePropertyFilename != "")
    {
        mIsMembranePropertySet = true;
        boost::shared_ptr<MembraneCellProperty> p_cell_membrane(new MembraneCellProperty());

        SetMembraneCellProperty(p_cell_membrane);
    }
    
    if (mIsCellIdSet)
    {
        boost::shared_ptr<CellAnalyticsProperty> p_cell_analytics(new CellAnalyticsProperty());

        SetCellAnalyticsProperty(p_cell_analytics);
    }

    SetUpSRNandCellCycle();

    SetUpCellObject();

    SetUpCellProperties();

    // update the initial conditons of the cell

    // superset of all internal cell chemical species, if not defined in intial conditions then set concentration to 0

    // chemical species from initial conditions
    
    std::vector<std::string> cellChemicalNames = mpCellChemicalProperty->GetStateVariableRegister()->GetStateVariableRegisterVector();
    std::vector<double> cellConcentrationVector = mpCellChemicalProperty->GetCellConcentrationVector();

    StateVariableRegister* pFullStateVariableRegister = new StateVariableRegister(cellChemicalNames);

    // chemical species from srn
    if (mIsSRNSet)
    {
        std::vector<std::string> srn_chemical_names = mpChemicalSrnModel->GetCellChemistry()->GetChemicalNames();
        pFullStateVariableRegister->AddStateVariableVector(srn_chemical_names);
    }  

    // chemical species from cell cycle
    if (mIsCellCycleSet)
    {
        std::vector<std::string> cell_cycle_chemical_names = mpSimpleChemicalThresholdCellCycleModel->GetThresholdChemistry()->GetChemicalNames();
        pFullStateVariableRegister->AddStateVariableVector(cell_cycle_chemical_names);
    }    

    // chemical species from membrane property
    if (mIsMembranePropertySet)
    {
        std::vector<std::string> cell_membrane_chemical_names = mpCellMembraneProperty->GetCellStateVariableRegister()->GetStateVariableRegisterVector();
        pFullStateVariableRegister->AddStateVariableVector(cell_membrane_chemical_names);
    }    

    // chemical species from transport property
    if (mIsTransportPropertySet)
    {
        std::vector<std::string> cell_transport_chemical_names = mpCellTransportProperty->GetCellStateVariableRegister()->GetStateVariableRegisterVector();
        pFullStateVariableRegister->AddStateVariableVector(cell_transport_chemical_names);
    } 

    // ensure all species are accounted for in initial conditions

    std::vector<std::string> fullChemicalNames = pFullStateVariableRegister->GetStateVariableRegisterVector();
    if (cellChemicalNames.size() != fullChemicalNames.size())
    {
        for (unsigned i=cellChemicalNames.size(); i<fullChemicalNames.size(); ++i)
        {
            cellConcentrationVector.push_back(0.0);
        }
    }
    
    SetFullChemicalStateRegister(pFullStateVariableRegister);

    SetUpCellInitialConditions(mpCell, cellChemicalNames, cellConcentrationVector);
}

void ComplexCellFromFile::SetUpCellProperties()
{
    if (mIsInitConditionsSet)
    {
        InitialCellConditionsFromFile* p_init_conditions_from_file = new InitialCellConditionsFromFile(mInitialConditionsFilename);
        StateVariableRegister* p_state_register = new StateVariableRegister(p_init_conditions_from_file->GetChemicalNamesVector());
        mpCellChemicalProperty->InitialiseCell(p_state_register,p_init_conditions_from_file->GetConcentrationVector());
        SetUpCellInitialConditions(mpCell, p_state_register->GetStateVariableRegisterVector(), p_init_conditions_from_file->GetConcentrationVector());
    }

    if (mIsTransportPropertySet)
    {
        TransportCellPropertyFromFile* p_transport_property_from_file = new TransportCellPropertyFromFile(mTransportPropertyFilename);
        p_transport_property_from_file->SetUpTransportProperty(mpCell);
        mpCellTransportProperty = p_transport_property_from_file->GetTransportProperty();
    }

    if (mIsMembranePropertySet)
    {
        MembraneCellPropertyFromFile* p_membrane_property_from_file = new MembraneCellPropertyFromFile(mMembranePropertyFilename);

        p_membrane_property_from_file->SetUpMembraneProperty(mpCell);

        mpCellMembraneProperty = p_membrane_property_from_file->GetMembraneProperty();

        mpCellMembraneProperty->SetMembraneThickness(5.0);
    }

    if (mIsCellIdSet)
    {
        CellAnalyticsPropertyFromCellId* p_cell_analytics_property_from_cellId = new CellAnalyticsPropertyFromCellId(mCellId);

        p_cell_analytics_property_from_cellId->SetUpCellAnalyticsProperty(mpCell);

        mpCellAnalyticsProperty = p_cell_analytics_property_from_cellId->GetCellAnalyticsProperty();
    }
}

void ComplexCellFromFile::SetUpSRNandCellCycle()
{
    if (mIsSRNSet)
    {
         ChemicalSrnFromFile* p_srn_reaction_system_from_file = new ChemicalSrnFromFile(mSrnFilename);
        SetChemicalSrnModel(p_srn_reaction_system_from_file->GetChemicalSrnModel());
        AbstractChemistry* this_cell_chemistry = mpChemicalSrnModel->GetCellChemistry();
        unsigned numChemicals = this_cell_chemistry->GetNumberChemicals();
    }

    if (mIsCellCycleSet)
    {
        SimpleChemicalThresholdCellCycleFromFile* pCellCycle = new SimpleChemicalThresholdCellCycleFromFile(mCellCycleFilename);

        pCellCycle->SetUp();
        SetChemicalCellCycleModel(pCellCycle);
    }
}

void ComplexCellFromFile::SetUpCellObject()
{
    // Form cell
    CellPropertyCollection collection;
    MAKE_PTR(WildTypeCellMutationState, p_state);
    if (mIsInitConditionsSet)
    {
        collection.AddProperty(mpCellChemicalProperty);
    }

    if (mIsTransportPropertySet)
    {
        collection.AddProperty(mpCellTransportProperty);
    }
    
    if (mIsMembranePropertySet)
    {
        collection.AddProperty(mpCellMembraneProperty);
    }

    if (mIsCellIdSet)
    {
        collection.AddProperty(mpCellAnalyticsProperty);
    }

    // Check if SRN and cell cycle are set, if not provide default models
    AbstractSrnModel* pSrnModel=nullptr;

    if (mIsSRNSet)
    {
        pSrnModel = GetChemicalSrnModel();
    }
    
    AbstractCellCycleModel* pCellCycleModel = new NoCellCycleModel();

    if (mIsCellCycleSet)
    {
        pCellCycleModel = GetChemicalCellCycleModel();
    }

    // Call cell constructor
    CellPtr p_cell(new ComplexCell(p_state, pCellCycleModel, pSrnModel, false, collection));

    // At present cell state and proliferation type is default
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    p_cell->SetCellProliferativeType(p_stem_type);
    
    if (mIsDivisionRulesSet)
    {
        SetUpCellDivisionRules(p_cell);
    }
    
    SetCellPtr(p_cell);
}

void ComplexCellFromFile::SetUpCellInitialConditions(CellPtr p_cell, std::vector<std::string> speciesNames, std::vector<double> initValue)
{
    for (unsigned i=0; i<speciesNames.size(); ++i)
    {
        p_cell->GetCellData()->SetItem(speciesNames[i],initValue[i]);
    }
}

void ComplexCellFromFile::SetUpCellDivisionRules(CellPtr p_cell)
{
    // function to read the division rules from a file and supply them to the cell model
    std::vector<std::string> chemicalNames;
    std::vector<std::string> chemicalDivisionRules;

    std::vector<std::vector<std::string>> divisionDatabase;

    divisionDatabase = ReadMatrix(mDivisionRulesFilename);

    unsigned numChemicals = divisionDatabase.size();

    for (unsigned row=0; row<numChemicals;row++)
    {
        chemicalNames.push_back(divisionDatabase[row][0]);
        chemicalDivisionRules.push_back(divisionDatabase[row][1]);
    }

    auto p_cell_new = boost::static_pointer_cast<ComplexCell>(p_cell);
    p_cell_new->SetChemicalNames(chemicalNames);

    p_cell_new->SetChemicalDivsionRules(chemicalDivisionRules);
}

std::vector<std::vector<std::string>> ComplexCellFromFile::ReadMatrix(std::string filename)
{
    // parse a matrix file (.csv) line by line, ignore escape line,s containing file information
    // that is lines starting with '#' 
    std::string line;
    std::ifstream inputFile(filename);
    
    // read all data types as std::string therefore return the matrix of strings for personalised
    // methods down stream
    std::vector<std::vector<std::string>> outputMatrix = std::vector<std::vector<std::string>>();

    // check file exists and is openable
    if (inputFile.is_open())
    {
        // open the matrix file
        while (getline(inputFile,line))
        {
            // while the file still has lines not read.
            // read line left to right, top to bottom.
            if (!line.empty())
            {
                if (line.at(0) != '#')
                {
                    outputMatrix.push_back(parseMatrixLineString(line));
                }   
            }
        }
        inputFile.close();

        return outputMatrix;
    }
    else
    {
        ///\todo replace with EXCEPTION
        std::cout << "Error: Unable to open file: "<<filename << std::endl;
        return outputMatrix;
    }
}

std::vector<std::string> ComplexCellFromFile::parseMatrixLineString(std::string line)
{
    // For a line string in the matrix read, parse into vector data entries based on delimiters ','
    std::vector<std::string> rowVector = std::vector<std::string>();

    // Delimiter, may be modified by further methods
    std::string delim = ",";
    std::string matrixCell;

    // Determine the position of the delimiter
    size_t delimiter_pos = line.find(delim);

    bool is_end_of_line = false;
    
    while (!is_end_of_line)
    {
        // while not at the end of the file, sample sub strings from the posiiton of the delimiter
        if (delimiter_pos == std::string::npos)
        {
            is_end_of_line = true;
        }
        
        // sample substring from begining of the string to the delimiter positioon, store as data entry
        matrixCell = line.substr(0,delimiter_pos);

        // remove the sampled entry from the string
        line = line.substr(delimiter_pos+1, std::string::npos);

        rowVector.push_back(matrixCell);

        // update delimiter position
        delimiter_pos = line.find(delim);
    }
    return rowVector;
}

CellPtr ComplexCellFromFile::GetCellPtr()
{
    return mpCell;
}

CellPropertyCollection ComplexCellFromFile::GetCellPropertyCollection()
{
    return mPropertyCollection;
}

boost::shared_ptr<ChemicalCellProperty> ComplexCellFromFile::GetChemicalCellProperty()
{
    return mpCellChemicalProperty;
}

boost::shared_ptr<MembraneCellProperty> ComplexCellFromFile::GetMembraneCellProperty()
{
    return mpCellMembraneProperty;
}

boost::shared_ptr<TransportCellProperty> ComplexCellFromFile::GetTransportCellProperty()
{
    return mpCellTransportProperty;
}

boost::shared_ptr<CellAnalyticsProperty> ComplexCellFromFile::GetCellAnalyticsProperty()
{
    return mpCellAnalyticsProperty;
}

ChemicalSrnModel* ComplexCellFromFile::GetChemicalSrnModel()
{
    return mpChemicalSrnModel;
}

SimpleChemicalThresholdCellCycleModel* ComplexCellFromFile::GetChemicalCellCycleModel()
{
    return mpSimpleChemicalThresholdCellCycleModel;
}

std::string ComplexCellFromFile::GetCellCycleFilename()
{
    return mCellCycleFilename;
}

bool ComplexCellFromFile::GetIsCellCycleSet()
{
    return mIsCellCycleSet;
}

std::string ComplexCellFromFile::GetSrnFilename()
{
    return mSrnFilename;
}

bool ComplexCellFromFile::GetIsSrnSet()
{
    return mIsSRNSet;
}

std::string ComplexCellFromFile::GetInitialConditionsFilename()
{
    return mInitialConditionsFilename;
}

bool ComplexCellFromFile::GetInitConditionsSet()
{
    return mIsInitConditionsSet;
}

std::string ComplexCellFromFile::GetTransportPropertyFilename()
{
    return mTransportPropertyFilename;
}

bool ComplexCellFromFile::GetIsTransportPropertySet()
{
    return mIsTransportPropertySet;
}

std::string ComplexCellFromFile::GetMembranePropertyFilename()
{
    return mMembranePropertyFilename;
}

bool ComplexCellFromFile::GetIsMembranePropertySet()
{
    return mIsMembranePropertySet;
}

unsigned ComplexCellFromFile::GetCellId()
{
    return mCellId;
}

bool ComplexCellFromFile::GetIsCellIdSet()
{
    return mIsCellIdSet;
}

StateVariableRegister* ComplexCellFromFile::GetFullChemicalStateRegister()
{
    return mpFullChemicalStateRegister;
}

std::vector<std::string> ComplexCellFromFile::GetFullChemicalNamesVector()
{
    return mpFullChemicalStateRegister->GetStateVariableRegisterVector();
}

void ComplexCellFromFile::SetCellPtr(CellPtr p_cell)
{
    mpCell = p_cell;
}

void ComplexCellFromFile::SetCellPropertyCollection(CellPropertyCollection propertyCollection)
{
    mPropertyCollection = propertyCollection;
}

void ComplexCellFromFile::SetChemicalCellProperty(boost::shared_ptr<ChemicalCellProperty> p_chemical)
{
    mpCellChemicalProperty = p_chemical;
}

void ComplexCellFromFile::SetMembraneCellProperty(boost::shared_ptr<MembraneCellProperty> p_membrane)
{
    mpCellMembraneProperty = p_membrane;
}

void ComplexCellFromFile::SetTransportCellProperty(boost::shared_ptr<TransportCellProperty> p_transport)
{
    mpCellTransportProperty = p_transport;
}

void ComplexCellFromFile::SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty> p_cellAnalytics)
{
    mpCellAnalyticsProperty = p_cellAnalytics;
}

void ComplexCellFromFile::SetChemicalSrnModel(ChemicalSrnModel* p_chemicalSrn)
{
    mpChemicalSrnModel = p_chemicalSrn;
}

void ComplexCellFromFile::SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel* p_chemicalCellCycleModel)
{
    mpSimpleChemicalThresholdCellCycleModel = p_chemicalCellCycleModel;
}

void ComplexCellFromFile::SetCellCycleFilename(std::string filename)
{
    mCellCycleFilename = filename;
    mIsCellCycleSet = true;
}

void ComplexCellFromFile::SetSrnFilename(std::string filename)
{
    mSrnFilename = filename;
    mIsSRNSet = true;
}

void ComplexCellFromFile::SetInitialConditionsFilename(std::string filename)
{
    mInitialConditionsFilename = filename;
    mIsInitConditionsSet = true;
}

void ComplexCellFromFile::SetTransportPropertyFilename(std::string filename)
{
    mTransportPropertyFilename = filename;
    mIsTransportPropertySet = true;
}

void ComplexCellFromFile::SetMembranePropertyFilename(std::string filename)
{
    mMembranePropertyFilename = filename;
    mIsMembranePropertySet = true;
}

void ComplexCellFromFile::SetFullChemicalStateRegister(StateVariableRegister* pRegister)
{
    mpFullChemicalStateRegister = pRegister;
}

#endif /* COMPLEXCELLFROMFILE_HPP_ */