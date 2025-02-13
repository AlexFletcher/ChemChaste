#ifndef CHEMICALTRACKINGMODIFIER_HPP_
#define CHEMICALTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "ChemicalSrnModel.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellProperty.hpp"
#include "TransportCellProperty.hpp"
#include "MembraneCellProperty.hpp"

/**
 * Class to handle the different methods of altering the chemical concentrations 
 * (state variable) within each cell in a cell based simulation coupled to a 
 * domain based PDE system. The class utlises and update cell properties and 
 * sums all the contributions at the end of a simulation timestep, chekcing the 
 * concentrations against a zeroing threshold.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ChemicalTrackingModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
private:
    double mDeltaError;

    bool mIsCheckConcentration;

public:

    ChemicalTrackingModifier();

    virtual ~ChemicalTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

    void SetupSRNFromCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    void UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    void SetDeltaError(double);

    double GetDeltaError();

    void CheckConcentration(const std::vector<double>&);

    void CheckConcentration(const double&);

    void SetIsCheckConcentration(bool);

    bool GetIsCheckConcentration();
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::ChemicalTrackingModifier()
    : AbstractCellBasedSimulationModifier<SPACE_DIM>()
{
    SetIsCheckConcentration(true);
    SetDeltaError(1e-6);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::~ChemicalTrackingModifier()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    SetupSRNFromCellData(rCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetupSRNFromCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    rCellPopulation.Update();

    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        // update cells ODE states from cellData
        ChemicalSrnModel* p_model = static_cast<ChemicalSrnModel*>(cell_iter->GetSrnModel());
    
        AbstractChemistry* this_cell_chemistry = p_model->GetCellChemistry();
  
        unsigned numChemicals = this_cell_chemistry->GetNumberChemicals();
    
        std::vector<double> this_SRN_concentration_vector(numChemicals,0.0);
        std::string this_chemical_name="";
 
        for (unsigned i=0; i<numChemicals;++i)
        {
            this_chemical_name = this_cell_chemistry->GetChemicalNamesByIndex(i);
            this_SRN_concentration_vector[i] = cell_iter->GetCellData()->GetItem(this_chemical_name);
        }
       
        p_model->GetOdeSystem()->SetStateVariables(this_SRN_concentration_vector);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // Make sure the cell population is updated, i.e cell cycle, cell state
    rCellPopulation.Update();

    // Run through each cell in the population and update the internal chemical concentrations in cellData
    unsigned count = 0;
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        ChemicalSrnModel* p_model = static_cast<ChemicalSrnModel*>(cell_iter->GetSrnModel());
        AbstractChemistry* this_cell_srn_chemistry = p_model->GetCellChemistry();
        unsigned numSrnChemicals = this_cell_srn_chemistry->GetNumberChemicals();
        
        // Retrieve the cell data values for the SRN model 
        double this_concentration = 0;
        std::string this_name = "";
        for (unsigned i=0; i<numSrnChemicals;++i)
        {
            this_name = this_cell_srn_chemistry->GetChemicalNamesByIndex(i);
           
            this_concentration = p_model->GetStateValueByName(this_name);

            CheckConcentration(this_concentration);
            cell_iter->GetCellData()->SetItem(this_name, this_concentration);
        }

        // Update the cell data based on any transport properties
        if (cell_iter->template HasCellProperty<TransportCellProperty>())
        {
            auto transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(cell_iter->template rGetCellPropertyCollection(). template GetPropertiesType<TransportCellProperty>().GetProperty());

            AbstractChemistry* this_cell_transport_chemistry = transport_cell_property->GetTransportReactionSystem()->GetCellChemistry();
            unsigned numTransportChemicals = this_cell_transport_chemistry->GetNumberChemicals();
            
            for (unsigned i=0; i<numTransportChemicals; ++i)
            {
                this_name = this_cell_transport_chemistry->GetChemicalNamesByIndex(i);

                this_concentration = cell_iter->GetCellData()->GetItem(this_name);        
                this_concentration += transport_cell_property->GetReactionCalls() *transport_cell_property->GetChangeInternalCellBoundaryConcentrationByName(this_name);

                CheckConcentration(this_concentration);
                cell_iter->GetCellData()->SetItem(this_name, this_concentration);
            }
        }

        // Update the cell data based on any membrane properties
        if (cell_iter->template HasCellProperty<MembraneCellProperty>())
        {
            auto membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(cell_iter->template rGetCellPropertyCollection(). template GetPropertiesType<MembraneCellProperty>().GetProperty());

            AbstractChemistry* this_cell_membrane_chemistry = membrane_cell_property->GetMembraneReactionSystem()->GetCellChemistry();
            unsigned numMembraneChemicals = this_cell_membrane_chemistry->GetNumberChemicals();
            
            for (unsigned i=0; i<numMembraneChemicals; ++i)
            {
                this_name = this_cell_membrane_chemistry->GetChemicalNamesByIndex(i);
           
                this_concentration = cell_iter->GetCellData()->GetItem(this_name);
                
                this_concentration += membrane_cell_property->GetReactionCalls() *membrane_cell_property->GetChangeInternalCellBoundaryConcentrationByName(this_name);
                CheckConcentration(this_concentration);
                cell_iter->GetCellData()->SetItem(this_name, this_concentration);
            }
        }

        // update the internal cell states for the properties
        if (cell_iter->template HasCellProperty<TransportCellProperty>())
        {
            auto transport_cell_property = boost::static_pointer_cast<TransportCellProperty>(cell_iter->template rGetCellPropertyCollection(). template GetPropertiesType<TransportCellProperty>().GetProperty());

            AbstractChemistry* this_cell_transport_chemistry = transport_cell_property->GetTransportReactionSystem()->GetCellChemistry();
            unsigned numTransportChemicals = this_cell_transport_chemistry->GetNumberChemicals();

            std::vector<double> transportStateVector(numTransportChemicals,0.0);
            
            std::string this_state="";
            for (unsigned i=0; i<numTransportChemicals; ++i)
            {
                this_state = this_cell_transport_chemistry->GetChemicalNamesByIndex(i);
            
                transportStateVector[i] = cell_iter->GetCellData()->GetItem(this_state);
            }
            
            transport_cell_property->UpdateCellConcentrationVector(transportStateVector);
        }

        if (cell_iter->template HasCellProperty<MembraneCellProperty>())
        {
            auto membrane_cell_property = boost::static_pointer_cast<MembraneCellProperty>(cell_iter->template rGetCellPropertyCollection(). template GetPropertiesType<MembraneCellProperty>().GetProperty());

            AbstractChemistry* this_cell_membrane_chemistry = membrane_cell_property->GetMembraneReactionSystem()->GetCellChemistry();
            unsigned numMembraneChemicals = this_cell_membrane_chemistry->GetNumberChemicals();
            std::vector<double> membraneStateVector(numMembraneChemicals,0.0);
            
            std::string this_state="";
            for (unsigned i=0; i<numMembraneChemicals; ++i)
            {
                this_state = this_cell_membrane_chemistry->GetChemicalNamesByIndex(i);
            
                membraneStateVector[i] = cell_iter->GetCellData()->GetItem(this_state);
            }

            membrane_cell_property->UpdateCellConcentrationVector(membraneStateVector);
        }

        // finally update the SRN state variables to the new post transport values
        p_model->UpdateOdeStatesFromCellData();

        count++;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetDeltaError(double deltaError)
{
    mDeltaError = deltaError;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::GetDeltaError()
{
    return mDeltaError;
}

// overload between a vector and a double value
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::CheckConcentration(const std::vector<double>& rY)
{
    // If chemical concentration gets too low then nan can occur, concentration must be +ve
    for (unsigned i=0; i<rY.size(); ++i)
    {
        // due to the discrete nature occasionally rY can evaluate to <0
        // ensure rY >= 0
        if (rY[i]<mDeltaError)
        {
            const_cast<double&>(rY[i]) = 0;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::CheckConcentration(const double& rY)
{
    // if chemical concentration gets too low then nan can occur, concentration must be +ve

    // due to the discrete nature occasionally rY can evaluate to <0
    // ensure rY >= 0
    if (rY<mDeltaError)
    {
        const_cast<double&>(rY) = 0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::SetIsCheckConcentration(bool isCheckConcentration)
{
    mIsCheckConcentration = isCheckConcentration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ChemicalTrackingModifier<ELEMENT_DIM,SPACE_DIM>::GetIsCheckConcentration()
{
    return mIsCheckConcentration;    
}

#endif 