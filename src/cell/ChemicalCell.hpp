#ifndef CHEMICALCELL_HPP_
#define CHEMICALCELL_HPP_

#include "Cell_virtual.hpp"
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellMutationState.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "CellData.hpp"
#include "CellVecData.hpp"
#include "Cell.hpp"

#include "AbstractCellCycleModel.hpp"
#include "AbstractSrnModel.hpp"
#include "CellPropertyCollection.hpp"

#include "AbstractChemistry.hpp"
#include "MembraneCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "CellAnalyticsProperty.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).
class AbstractSrnModel; // Circular definition (cells need to know about subcellular reaction network models and vice-versa).
class Cell;

// version of cell.hpp in which the cell division properties may be overridden

class ChemicalCell : public Cell
{
protected:
    using Cell::Divide;

    double mSplitRatio = 0.5; // proportion of parent cell volume retained, rest goes to daughter 

public:

    ChemicalCell(boost::shared_ptr<AbstractCellProperty> pMutationState,
         AbstractCellCycleModel* pCellCycleModel,
         AbstractSrnModel* pSrnModel=nullptr,
         bool archiving=false,
         CellPropertyCollection cellPropertyCollection=CellPropertyCollection());

    virtual ~ChemicalCell()
    {
    };

    virtual CellPtr Divide();

    virtual void DetermineSplitRatio();

    virtual double SplitParentCellData(double);

    double GetSplitRatio();

    void SetSplitRatio(double);
};

ChemicalCell::ChemicalCell(
        boost::shared_ptr<AbstractCellProperty> pMutationState,
        AbstractCellCycleModel* pCellCycleModel,
        AbstractSrnModel* pSrnModel,
        bool archiving,
        CellPropertyCollection cellPropertyCollection)

    : Cell(pMutationState, pCellCycleModel, pSrnModel, archiving, cellPropertyCollection)
{
}

CellPtr ChemicalCell::Divide()
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

    // copy cell data

    // Copy all cell data (note we create a new object not just copying the pointer)
    assert(daughter_property_collection.HasPropertyType<CellData>());
    // Get the existing copy of the cell data and remove it from the daughter cell
    boost::shared_ptr<CellData> p_cell_data = GetCellData();
    daughter_property_collection.RemoveProperty(p_cell_data);
    // Create a new cell data object using the copy constructor and add this to the daughter cell
    MAKE_PTR_ARGS(CellData, p_daughter_cell_data, (*p_cell_data));
    daughter_property_collection.AddProperty(p_daughter_cell_data);
    // Copy all cell Vec data (note we create a new object not just copying the pointer)
    if (daughter_property_collection.HasPropertyType<CellVecData>())
    {
        // Get the existing copy of the cell data and remove it from the daughter cell
        boost::shared_ptr<CellVecData> p_cell_vec_data = GetCellVecData();
        daughter_property_collection.RemoveProperty(p_cell_vec_data);
        // Create a new cell data object using the copy constructor and add this to the daughter cell
        MAKE_PTR_ARGS(CellVecData, p_daughter_cell_vec_data, (*p_cell_vec_data));
        daughter_property_collection.AddProperty(p_daughter_cell_vec_data);
    }

    // record the cell chemistry for splitting cell data
    AbstractChemistry* p_cell_chemistry = new AbstractChemistry();

    if (static_cast<ChemicalSrnModel*>(mpSrnModel)->SRNType()=="Chemical")
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
        boost::shared_ptr<TransportCellProperty> p_transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* p_transport_chemistry = p_transport_cell_property->GetTransportReactionSystem()->GetCellChemistry();

        // Add the cell chemistry due to transport
        p_cell_chemistry->AddChemistry(p_transport_chemistry);
    }

    // Membrane property
    if (mCellPropertyCollection.HasProperty<MembraneCellProperty>())
    {
        boost::shared_ptr<MembraneCellProperty> p_membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(mCellPropertyCollection.GetPropertiesType<MembraneCellProperty>().GetProperty());
        AbstractChemistry* membraneChemistry = p_membrane_cell_property->GetMembraneReactionSystem()->GetCellChemistry();

        // Add the cell chemistry due to membrane
        p_cell_chemistry->AddChemistry(membraneChemistry);
    }

    // Based on the parent cell data determine the split ratio
    DetermineSplitRatio();

    // share the two cellDatas between the two cells, use p_cell_chemistry
    // split chemical cell data
    double parent_species_concentration=0.0;
    double new_parent_species_concentration=0.0;
    double daughter_species_concentration=0.0;
    unsigned numChemicals = p_cell_chemistry->GetNumberChemicals();

    for (unsigned i=0; i<numChemicals; i++)
    {
        parent_species_concentration = this->GetCellData()->GetItem(cellChemistry->GetChemicalNamesByIndex(i));
        
        new_parent_species_concentration = SplitParentCellData(parent_species_concentration);

        daughter_species_concentration = parent_species_concentration - new_parent_species_concentration;

        if (daughter_species_concentration<0.0)
        {
            // must be positive concentration
            daughter_species_concentration =0.0;
        }
        
        this->GetCellData()->SetItem(cellChemistry->GetChemicalNamesByIndex(i), new_parent_species_concentration);
        p_daughter_cell_data->SetItem(cellChemistry->GetChemicalNamesByIndex(i), daughter_species_concentration);
    }

    // run through cell properties and create new objects for them 
    
    //daughter_property_collection = mCellPropertyCollection;

    // transport property
    if (mCellPropertyCollection.HasProperty<TransportCellProperty>())
    {
        // parent property
        boost::shared_ptr<TransportCellProperty> p_transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(mCellPropertyCollection.GetPropertiesType<TransportCellProperty>().GetProperty());
        AbstractChemistry* p_transport_chemistry = p_transport_cell_property ->GetTransportReactionSystem()->GetCellChemistry();

        // create new transport cell property
        daughter_property_collection.RemoveProperty(p_transport_cell_property);
        // use copy construtor
        boost::shared_ptr<TransportCellProperty> p_daughter_transport_property(new TransportCellProperty(*p_transport_cell_property));
        daughter_property_collection.AddProperty(p_daughter_transport_property);

        // split the properties betwene the two cells
        p_transport_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_transport_property->PreparePostDivisionDaughter(*p_transport_cell_property, mSplitRatio);
    }

    // Membrane cell property
    if (mCellPropertyCollection.HasProperty<MembraneCellProperty>())
    {
        boost::shared_ptr<MembraneCellProperty> p_membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(mCellPropertyCollection.GetPropertiesType<MembraneCellProperty>().GetProperty());

        // create new membrane cell property
        daughter_property_collection.RemoveProperty(p_membrane_cell_property);
        boost::shared_ptr<MembraneCellProperty> p_daughter_membrane_property(new MembraneCellProperty(*p_membrane_cell_property));
        daughter_property_collection.AddProperty(p_daughter_membrane_property);

        // split the properties betwene the two cells
        p_membrane_cell_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_membrane_property->PreparePostDivisionDaughter(*p_membrane_cell_property, mSplitRatio);
    }

    // Cell analytics property
    if (mCellPropertyCollection.HasProperty<CellAnalyticsProperty>())
    {
        // Parent property
        boost::shared_ptr<CellAnalyticsProperty> p_cell_analytics_property = boost::static_pointer_cast<CellAnalyticsProperty>(mCellPropertyCollection.GetPropertiesType<CellAnalyticsProperty>().GetProperty());

        // Create new transport cell property
        daughter_property_collection.RemoveProperty(p_cell_analytics_property);

        // Use copy construtor
        boost::shared_ptr<CellAnalyticsProperty> p_daughter_cellAnalytics_property(new CellAnalyticsProperty(*p_cell_analytics_property));
        daughter_property_collection.AddProperty(p_daughter_cellAnalytics_property);

        // Split the properties between the two cells
        p_cell_analytics_property->PreparePostDivisionParent(mSplitRatio);
        p_daughter_cellAnalytics_property->PreparePostDivisionDaughter(*p_cell_analytics_property, mSplitRatio);
    }

    // Create new chemical cell

    // Create daughter cell with modified cell property collection
    CellPtr p_new_cell(new ChemicalCell(GetMutationState(), mpCellCycleModel->CreateCellCycleModel(), mpSrnModel->CreateSrnModel(), false, daughter_property_collection));

    // Initialise properties of daughter cell
    p_new_cell->GetCellCycleModel()->InitialiseDaughterCell();
    p_new_cell->GetSrnModel()->InitialiseDaughterCell();

    // Set the daughter cell to inherit the apoptosis time of the parent cell
    p_new_cell->SetApoptosisTime(mApoptosisTime);

    std::vector<double> prime_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);
    std::vector<double> daughter_cell_threshold_species_concentrations(static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->GetNumberThresholdSpecies(),0.0);

    AbstractChemistry* p_threshold_chemistry = static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->GetThresholdChemistry();
    std::string chemical_name;

    for (unsigned i=0; i<numChemicals; i++)
    {
        chemical_name = p_cell_chemistry->GetChemicalNamesByIndex(i);

        if (!p_threshold_chemistry->CheckChemical(new AbstractChemical(chemical_name)))
        {
            prime_cell_threshold_species_concentrations[p_threshold_chemistry->GetChemicalIndexByName(chemical_name)] = p_cell_data->GetItem(chemical_name);
        }
    }

    for (unsigned i=0; i<numChemicals; i++)
    {
        chemical_name = p_cell_chemistry->GetChemicalNamesByIndex(i);
        if (!p_threshold_chemistry->CheckChemical(new AbstractChemical(chemical_name)))
        {
            daughter_cell_threshold_species_concentrations[p_threshold_chemistry->GetChemicalIndexByName(chemical_name)] = p_daughter_cell_data->GetItem(chemical_name);
        }
    }

    // Update chemical cell cycles for each of the daughter cells
    static_cast<SimpleChemicalThresholdCellCycleModel*>(this->GetCellCycleModel())->SetSpeciesConcentrations(prime_cell_threshold_species_concentrations);
    static_cast<SimpleChemicalThresholdCellCycleModel*>(p_new_cell->GetCellCycleModel())->SetSpeciesConcentrations(daughter_cell_threshold_species_concentrations);

    // Update ODE from cell data for each of the daughter cells
    static_cast<ChemicalSrnModel*>(this->GetSrnModel())->UpdateOdeStatesFromCellData();
    static_cast<ChemicalSrnModel*>(p_new_cell->GetSrnModel())->UpdateOdeStatesFromCellData();
    
    return p_new_cell;
}

void ChemicalCell::DetermineSplitRatio()
{
    mSplitRatio = 0.5;
}

double ChemicalCell::SplitParentCellData(double current_parent_value)
{
    return current_parent_value*GetSplitRatio();
}

double ChemicalCell::GetSplitRatio()
{
    return mSplitRatio;
}

void ChemicalCell::SetSplitRatio(double split_ratio)
{
    mSplitRatio = split_ratio;
}

#endif