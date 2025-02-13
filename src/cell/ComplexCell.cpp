#include "ComplexCell.hpp"

ComplexCell::ComplexCell(boost::shared_ptr<AbstractCellProperty> pMutationState,
                         AbstractCellCycleModel* pCellCycleModel,
                         AbstractSrnModel* pSrnModel,
                         bool archiving,
                         CellPropertyCollection cellPropertyCollection)
    : Cell(pMutationState,
           pCellCycleModel,
           pSrnModel,
           archiving,
           cellPropertyCollection)
{
}

CellPtr ComplexCell::Divide()
{
    // Check we're allowed to divide
    assert(!IsDead());
    assert(mCanDivide);

    mCanDivide = false;

    // Reset properties of parent cell
    mpCellCycleModel->ResetForDivision();
    mpSrnModel->ResetForDivision();

    // Create copy of cell property collection to modify for daughter cell
    CellPropertyCollection daughter_property_collection = mCellPropertyCollection;

    // Remove the CellId from the daughter cell, as a new one will be assigned in the constructor
    daughter_property_collection.RemoveProperty<CellId>();

    // Copy all cell data (note we create a new object not just copying the pointer)
    assert(daughter_property_collection.HasPropertyType<CellData>());

    // Get the existing copy of the cell data and remove it from the daughter cell
    boost::shared_ptr<CellData> p_cell_data = GetCellData();
    daughter_property_collection.RemoveProperty(p_cell_data);

    // Create a new cell data object using the copy constructor and add this to the daughter cell
    MAKE_PTR_ARGS(CellData, p_daughter_cell_data, (*p_cell_data));
    daughter_property_collection.AddProperty(p_daughter_cell_data);

    // Copy all cell vec data (note we create a new object not just copying the pointer)
    if (daughter_property_collection.HasPropertyType<CellVecData>())
    {
        // Get the existing copy of the cell data and remove it from the daughter cell
        boost::shared_ptr<CellVecData> p_cell_vec_data = GetCellVecData();
        daughter_property_collection.RemoveProperty(p_cell_vec_data);

        // Create a new cell data object using the copy constructor and add this to the daughter cell
        MAKE_PTR_ARGS(CellVecData, p_daughter_cell_vec_data, (*p_cell_vec_data));
        daughter_property_collection.AddProperty(p_daughter_cell_vec_data);
    }

    // Record the cell chemistry for splitting cell data
    AbstractChemistry* p_cell_chemistry = new AbstractChemistry();

    if (static_cast<ChemicalSrnModel*>(mpSrnModel)->SRNType() == "Chemical")
    {
        ChemicalSrnModel* p_srn_model = static_cast<ChemicalSrnModel*>(mpSrnModel);
        p_cell_chemistry->AddChemistry(p_srn_model->GetCellChemistry()); // from SRN
    }

    if (static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel)->CellCycleType()=="Chemical")
    {
        SimpleChemicalThresholdCellCycleModel* p_cc_model = static_cast<SimpleChemicalThresholdCellCycleModel*>(mpCellCycleModel);
        p_cell_chemistry->AddChemistry(p_cc_model->GetThresholdChemistry()); // from cell cycle model
    }

    // Transport property
    if (mCellPropertyCollection.HasProperty<TransportCellProperty>())
    {
        auto transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* p_transport_chemistry = transport_cell_property->GetTransportReactionSystem()->GetCellChemistry();

        // Add the cell chemistry due to transport
        p_cell_chemistry->AddChemistry(p_transport_chemistry);
    }

    // membrane property
    if (mCellPropertyCollection.HasProperty<MembraneCellProperty>())
    {
        auto p_membrane_property = boost::static_pointer_cast<MembraneCellProperty>(mCellPropertyCollection.GetPropertiesType<MembraneCellProperty>().GetProperty());
        AbstractChemistry* p_membrane_chemistry = p_membrane_property->GetMembraneReactionSystem()->GetCellChemistry();

        // Add the cell chemistry due to membrane
        p_cell_chemistry->AddChemistry(p_membrane_chemistry);
    }

    // based on the parent cell data determine the split ratio
    DetermineSplitRatio();

    // share the two cellDatas between the two cells, use p_cell_chemistry
    // split chemical cell data
    double parent_species_concentration = 0.0;
    double new_parent_species_concentration = 0.0;
    double daughter_species_concentration = 0.0;
    unsigned num_chemicals = p_cell_chemistry->GetNumberChemicals();

    for (unsigned i=0; i<num_chemicals; ++i)
    {
        parent_species_concentration = this->GetCellData()->GetItem(p_cell_chemistry->GetChemicalNamesByIndex(i));
        
        if (IsChemicalShared(p_cell_chemistry->GetChemicalNamesByIndex(i)))
        {
            // Chemical concentration is shared upon division
            new_parent_species_concentration = SplitParentCellData(parent_species_concentration);

            daughter_species_concentration = parent_species_concentration - new_parent_species_concentration;
        }
        else
        {
            // Chemical concentration is duplicated upon division
            new_parent_species_concentration = parent_species_concentration;
            daughter_species_concentration = parent_species_concentration;
        }
        
        // Ensure non-negative concentration
        daughter_species_concentration = std::max(daughter_species_concentration, 0.0);
        
        this->GetCellData()->SetItem(p_cell_chemistry->GetChemicalNamesByIndex(i), new_parent_species_concentration);
        p_daughter_cell_data->SetItem(p_cell_chemistry->GetChemicalNamesByIndex(i), daughter_species_concentration);
    }

    // Run through cell properties and create new objects for them

    // Transport property
    if (mCellPropertyCollection.HasProperty<TransportCellProperty>())
    {
        // Parent property
        auto transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* p_transport_chemistry = transport_cell_property->GetTransportReactionSystem()->GetCellChemistry();

        // Create new transport cell property
        daughter_property_collection.RemoveProperty(transport_cell_property);

        // Use copy construtor
        boost::shared_ptr<TransportCellProperty> p_daughter_transport_property(new TransportCellProperty(*transport_cell_property));
        daughter_property_collection.AddProperty(p_daughter_transport_property);

        // Split the properties between the two cells
        transport_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_transport_property->PreparePostDivisionDaughter(*transport_cell_property, mSplitRatio);
    }

    // Membrane property
    if (mCellPropertyCollection.HasProperty<MembraneCellProperty>())
    {
        // Parent property
        auto p_membrane_property = boost::static_pointer_cast<MembraneCellProperty>(mCellPropertyCollection.GetPropertiesType<MembraneCellProperty>().GetProperty());

        // Create new membrane cell property
        daughter_property_collection.RemoveProperty(p_membrane_property);
        boost::shared_ptr<MembraneCellProperty> p_daughter_membrane_property(new MembraneCellProperty(*p_membrane_property));
        daughter_property_collection.AddProperty(p_daughter_membrane_property);

        // Split the properties between the two cells
        p_membrane_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_membrane_property->PreparePostDivisionDaughter(*p_membrane_property, mSplitRatio);
    }

    // Cell analytics property
    if (mCellPropertyCollection.HasProperty<CellAnalyticsProperty>())
    {
        // Parent property
        auto p_analytics_property = boost::static_pointer_cast<CellAnalyticsProperty>(mCellPropertyCollection.GetPropertiesType<CellAnalyticsProperty>().GetProperty());

        // Create new transport cell property
        daughter_property_collection.RemoveProperty(p_analytics_property);

        // Use copy construtor
        boost::shared_ptr<CellAnalyticsProperty> p_daughter_analytics_property(new CellAnalyticsProperty(*p_analytics_property));
        daughter_property_collection.AddProperty(p_daughter_analytics_property);

        // Split the properties between the two cells
        p_analytics_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_analytics_property->PreparePostDivisionDaughter(*p_analytics_property, mSplitRatio);
    }

    // Create daughter cell with modified cell property collection
    CellPtr p_new_cell(new ComplexCell(GetMutationState(), mpCellCycleModel->CreateCellCycleModel(), mpSrnModel->CreateSrnModel(), false, daughter_property_collection));

    // Initialise properties of daughter cell
    p_new_cell->GetCellCycleModel()->InitialiseDaughterCell();
    p_new_cell->GetSrnModel()->InitialiseDaughterCell();

    // Set the daughter cell to inherit the apoptosis time of the parent cell
    p_new_cell->SetApoptosisTime(mApoptosisTime);

    std::vector<double> prime_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);
    std::vector<double> daughter_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);

    AbstractChemistry* p_threshold_chemistry = static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetThresholdChemistry();
    std::string chemical_name;

    for (unsigned i=0; i<num_chemicals; ++i)
    {
        chemical_name = p_cell_chemistry->GetChemicalNamesByIndex(i);

        if (!p_threshold_chemistry->CheckChemical(new AbstractChemical(chemical_name)))
        {
            prime_cell_threshold_species_concentrations[p_threshold_chemistry->GetChemicalIndexByName(chemical_name)] = p_cell_data->GetItem(chemical_name);
        }
    }

    for (unsigned i=0; i<num_chemicals; ++i)
    {
        chemical_name = p_cell_chemistry->GetChemicalNamesByIndex(i);
        if (!p_threshold_chemistry->CheckChemical(new AbstractChemical(chemical_name)))
        {
            daughter_cell_threshold_species_concentrations[p_threshold_chemistry->GetChemicalIndexByName(chemical_name)] = p_daughter_cell_data->GetItem(chemical_name);
        }
    }

    // update chemical cell cycles for each of the daughter cells
    static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->SetSpeciesConcentrations(prime_cell_threshold_species_concentrations);
    static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->SetSpeciesConcentrations(daughter_cell_threshold_species_concentrations);

    // update Ode from cell data for each of the daughter cells
    static_cast<ChemicalSrnModel*>(this->GetSrnModel())->UpdateOdeStatesFromCellData();
    static_cast<ChemicalSrnModel*>(p_new_cell->GetSrnModel())->UpdateOdeStatesFromCellData();
    
    return p_new_cell;
}

void ComplexCell::DetermineSplitRatio()
{
    mSplitRatio = 0.5;
}

double ComplexCell::SplitParentCellData(double current_parent_value)
{
    return current_parent_value*GetSplitRatio();
}

bool ComplexCell::IsChemicalShared(std::string chemical_name)
{
    // Determine whether the chemical concetration is shared between parnt and daughter after division
    bool is_chemical_shared = false;

    for (unsigned i = 0; i < mChemicalNames.size(); ++i)
    {
        // For each cellular chemical
        if (mChemicalNames[i] == chemical_name)
        {
            for (unsigned j = 0; j < mShareKey.size(); ++j)
            {
                // Test if the chemical is to be shared
                if (mChemicalDivsionRules[i] == mShareKey[j])
                {
                    is_chemical_shared = true;
                }
            }
        }
    }

    return is_chemical_shared; 
}

double ComplexCell::GetSplitRatio()
{
    return mSplitRatio;
}

void ComplexCell::SetSplitRatio(double splitRatio)
{
    mSplitRatio = splitRatio;
}

std::vector<std::string> ComplexCell::GetShareKey()
{
    return mShareKey;
}

void ComplexCell::SetShareKey(std::vector<std::string> shareKey)
{
    mShareKey = shareKey;
}

std::vector<std::string> ComplexCell::GetChemicalNames()
{
    return mChemicalNames;
}

void ComplexCell::SetChemicalNames(std::vector<std::string> chemicalNames)
{
    mChemicalNames = chemicalNames;
}

std::vector<std::string> ComplexCell::GetChemicalDivsionRules()
{
    return mChemicalDivsionRules;
}

void ComplexCell::SetChemicalDivsionRules(std::vector<std::string> chemicalRules)
{
    mChemicalDivsionRules = chemicalRules;
}