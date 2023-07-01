#include "ChemicalCellFromFile.hpp"

ChemicalCellFromFile::ChemicalCellFromFile(
                        std::string cellCycleFilename,
                        std::string srnFilename,
                        std::string initialConditionFilename,
                        std::string transportPropertyFilename,
                        std::string membranePropertyFilename,
                        unsigned cellId,
                        bool isCellIdSet)
     : mCellCycleFilename(cellCycleFilename),
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
    std::vector<std::string> full_chemical_names = pFullStateVariableRegister->GetStateVariableRegisterVector();
    if (cellChemicalNames.size() != full_chemical_names.size())
    {
        for (unsigned i=cellChemicalNames.size(); i<full_chemical_names.size(); ++i)
        {
            cellConcentrationVector.push_back(0.0);
        }
    }
    
    SetFullChemicalStateRegister(pFullStateVariableRegister);

    SetUpCellInitialConditions(mpCell, cellChemicalNames, cellConcentrationVector);
}

void ChemicalCellFromFile::SetUpCellProperties()
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

void ChemicalCellFromFile::SetUpSRNandCellCycle()
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

void ChemicalCellFromFile::SetUpCellObject()
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
    CellPtr p_cell(new ChemicalCell(p_state, pCellCycleModel, pSrnModel, false, collection));

    // At present cell state and proliferation type is default
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    p_cell->SetCellProliferativeType(p_stem_type);
    
    SetCellPtr(p_cell);
}

void ChemicalCellFromFile::SetUpCellInitialConditions(CellPtr p_cell, std::vector<std::string> speciesNames, std::vector<double> initValue)
{
    for (unsigned i = 0; i < speciesNames.size(); ++i)
    {
        p_cell->GetCellData()->SetItem(speciesNames[i],initValue[i]);
    }
}

CellPtr ChemicalCellFromFile::GetCellPtr()
{
    return mpCell;
}

CellPropertyCollection ChemicalCellFromFile::GetCellPropertyCollection()
{
    return mPropertyCollection;
}

boost::shared_ptr<ChemicalCellProperty> ChemicalCellFromFile::GetChemicalCellProperty()
{
    return mpCellChemicalProperty;
}

boost::shared_ptr<MembraneCellProperty> ChemicalCellFromFile::GetMembraneCellProperty()
{
    return mpCellMembraneProperty;
}

boost::shared_ptr<TransportCellProperty> ChemicalCellFromFile::GetTransportCellProperty()
{
    return mpCellTransportProperty;
}

boost::shared_ptr<CellAnalyticsProperty> ChemicalCellFromFile::GetCellAnalyticsProperty()
{
    return mpCellAnalyticsProperty;
}

ChemicalSrnModel* ChemicalCellFromFile::GetChemicalSrnModel()
{
    return mpChemicalSrnModel;
}

SimpleChemicalThresholdCellCycleModel* ChemicalCellFromFile::GetChemicalCellCycleModel()
{
    return mpSimpleChemicalThresholdCellCycleModel;
}

std::string ChemicalCellFromFile::GetCellCycleFilename()
{
    return mCellCycleFilename;
}

bool ChemicalCellFromFile::GetIsCellCycleSet()
{
    return mIsCellCycleSet;
}

std::string ChemicalCellFromFile::GetSrnFilename()
{
    return mSrnFilename;
}

bool ChemicalCellFromFile::GetIsSrnSet()
{
    return mIsSRNSet;
}

std::string ChemicalCellFromFile::GetInitialConditionsFilename()
{
    return mInitialConditionsFilename;
}

bool ChemicalCellFromFile::GetInitConditionsSet()
{
    return mIsInitConditionsSet;
}

std::string ChemicalCellFromFile::GetTransportPropertyFilename()
{
    return mTransportPropertyFilename;
}

bool ChemicalCellFromFile::GetIsTransportPropertySet()
{
    return mIsTransportPropertySet;
}

std::string ChemicalCellFromFile::GetMembranePropertyFilename()
{
    return mMembranePropertyFilename;
}

bool ChemicalCellFromFile::GetIsMembranePropertySet()
{
    return mIsMembranePropertySet;
}

unsigned ChemicalCellFromFile::GetCellId()
{
    return mCellId;
}

bool ChemicalCellFromFile::GetIsCellIdSet()
{
    return mIsCellIdSet;
}

StateVariableRegister* ChemicalCellFromFile::GetFullChemicalStateRegister()
{
    return mpFullChemicalStateRegister;
}

std::vector<std::string> ChemicalCellFromFile::GetFullChemicalNamesVector()
{
    return mpFullChemicalStateRegister->GetStateVariableRegisterVector();
}

void ChemicalCellFromFile::SetCellPtr(CellPtr pCell)
{
    mpCell = pCell;
}

void ChemicalCellFromFile::SetCellPropertyCollection(CellPropertyCollection propertyCollection)
{
    mPropertyCollection = propertyCollection;
}

void ChemicalCellFromFile::SetChemicalCellProperty(boost::shared_ptr<ChemicalCellProperty> pChemical)
{
    mpCellChemicalProperty = pChemical;
}

void ChemicalCellFromFile::SetMembraneCellProperty(boost::shared_ptr<MembraneCellProperty> pMembrane)
{
    mpCellMembraneProperty = pMembrane;
}

void ChemicalCellFromFile::SetTransportCellProperty(boost::shared_ptr<TransportCellProperty> pTransport)
{
    mpCellTransportProperty = pTransport;
}

void ChemicalCellFromFile::SetCellAnalyticsProperty(boost::shared_ptr<CellAnalyticsProperty> pCellAnalytics)
{
    mpCellAnalyticsProperty = pCellAnalytics;
}

void ChemicalCellFromFile::SetChemicalSrnModel(ChemicalSrnModel* pChemicalSrn)
{
    mpChemicalSrnModel = pChemicalSrn;
}

void ChemicalCellFromFile::SetChemicalCellCycleModel(SimpleChemicalThresholdCellCycleModel* pChemicalCellCycleModel)
{
    mpSimpleChemicalThresholdCellCycleModel = pChemicalCellCycleModel;
}

void ChemicalCellFromFile::SetCellCycleFilename(std::string filename)
{
    mCellCycleFilename = filename;
    mIsCellCycleSet = true;
}

void ChemicalCellFromFile::SetSrnFilename(std::string filename)
{
    mSrnFilename = filename;
    mIsSRNSet = true;
}

void ChemicalCellFromFile::SetInitialConditionsFilename(std::string filename)
{
    mInitialConditionsFilename = filename;
    mIsInitConditionsSet = true;
}

void ChemicalCellFromFile::SetTransportPropertyFilename(std::string filename)
{
    mTransportPropertyFilename = filename;
    mIsTransportPropertySet = true;
}

void ChemicalCellFromFile::SetMembranePropertyFilename(std::string filename)
{
    mMembranePropertyFilename = filename;
    mIsMembranePropertySet = true;
}

void ChemicalCellFromFile::SetFullChemicalStateRegister(StateVariableRegister* pRegister)
{
    mpFullChemicalStateRegister = pRegister;
}