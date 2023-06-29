#ifndef SCHNACKENBERGSRNMODEL_HPP_
#define SCHNACKENBERGSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "SchnackenbergOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

/**
 * A subclass of AbstractOdeSrnModel that includes a Delta-Notch ODE system in 
 * the sub-cellular reaction network.
 *
 * \todo #2752 document this class more thoroughly here
 */
class SchnackenbergSrnModel : public AbstractOdeSrnModel
{
protected:

    /**
     * Protected copy-constructor for use by CreateSrnModel(). The only way for 
     * external code to create a copy of a SRN model is by calling that method, 
     * to ensure that a model of the correct subclass is created. This 
     * copy-constructor helps subclasses to ensure that all member variables are 
     * correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a 
     * daughter cell upon cell division. Note that the parent SRN model will 
     * have had ResetForDivision() called just before CreateSrnModel() is 
     * called, so performing an exact copy of the parent is suitable behaviour. 
     * Any daughter-cell-specific initialisation can be done in 
     * InitialiseDaughterCell().
     *
     * @param rModel the SRN model to copy.
     */
    SchnackenbergSrnModel(const SchnackenbergSrnModel& rModel);

public:

    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver 
     *                   object (allows the use of different ODE solvers)
     */
    SchnackenbergSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return a copy of the current SRN model.
     */
    AbstractSrnModel* CreateSrnModel();

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise();

    /**
     * Overridden SimulateToTime() method for custom behaviour.
     *
     * \todo #2752 say what it does in this class
     */
    void SimulateToCurrentTime();

    /**
     * @return the value of the state variable U
     */
    double GetU();

    /**
     * @return the value of the state variable V
     */
    double GetV();

    /**
     * Output SRN model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSrnModelParameters(out_stream& rParamsFile);
};

#endif /* SCHNACKENBERSRNMODEL_HPP_ */