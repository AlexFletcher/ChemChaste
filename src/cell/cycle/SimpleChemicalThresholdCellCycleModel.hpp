#ifndef SIMPLECHEMICALTHRESHOLDCELLCYCLEMODEL_HPP_
#define SIMPLECHEMICALTHRESHOLDCELLCYCLEMODEL_HPP_

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "AbstractChemistry.hpp"

/**
 * Simple cell cycle model wherein the concentrations of species within the cell 
 * are compared to maximum and minimum thresholds. If a single species is above 
 * the threshols then cell division is initilaed while if below cell death is 
 * triggered.
 */
class SimpleChemicalThresholdCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
{
protected:

    /** \todo Document member variable. */
    double mCurrentStarvationDuration;

    /** \todo Document member variable. */
    double mCurrentStarvationOnsetTime;

    /** \todo Document member variable. */
    double mCriticalStarvationDuration;

    /** \todo Document member variable. */
    AbstractChemistry* mpThresholdChemistry;

    /** \todo Document member variable. */
    std::vector<double> mSpeciesConcentrations;

    /** \todo Document member variable. */
    unsigned mNumberThresholdSpecies;

    /** \todo Document member variable. */
    std::vector<double> mMaxThresholdConcentrationValues;

    /** \todo Document member variable. */
    std::vector<bool> mIsMaximumThresholdSet;

    /** \todo Document member variable. */
    std::vector<double> mMinThresholdConcentrationValues;

    /** \todo Document member variable. */
    std::vector<bool> mIsMinimumThresholdSet;

    /** \todo Implement archiving. */
public:

    /**
     * Default constructor.
     */
    SimpleChemicalThresholdCellCycleModel();

    /**
     * Protected copy-constructor for use by CreateCellCycleModel. The only way 
     * for external code to create a copy of a cell cycle model is by calling 
     * that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member 
     * variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a 
     * daughter cell upon cell division. Note that the parent cell cycle model 
     * will have had ResetForDivision() called just before 
     * CreateCellCycleModel() is called, so performing an exact copy of the 
     * parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    SimpleChemicalThresholdCellCycleModel(const SimpleChemicalThresholdCellCycleModel& rModel);

    /**
     * Destructor.
     */
    virtual ~SimpleChemicalThresholdCellCycleModel();

    /**
     * \todo Document method.
     * 
     * @return
     */ 
    virtual std::string CellCycleType()
    {
        return "Chemical";
    }

    /**
     * \todo Document method.
     * 
     * @param pChemistry
     */
    void SetUp(AbstractChemistry* pChemistry);

    /**
     * \todo Document method.
     */
    void Initialise();

    /**
     * \todo Document method.
     */
    void InitialiseDaughterCell();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    double GetCurrentStarvationDuration() const;

    /**
     * \todo Document method.
     * 
     * @return 
     */
    double GetCurrentStarvationOnsetTime() const;

    /**
     * \todo Document method.
     */
    void UpdateCellCyclePhase();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    bool ReadyToDivide();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * \todo Document method.
     */
    void UpdateStarvationDuration();

    /**
     * \todo Document method.
     * 
     * @param parameters
     * 
     * @return 
     */
    double CellDeathProbability(std::vector<double> parameters);

    /**
     * \todo Document method.
     * 
     * @return 
     */
    bool ConcentrationAboveMaxThreshold();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    bool ConcentrationBelowMinThreshold();

    /**
     * \todo Document method.
     */
    void PrepareForDivision();

    /**
     * \todo Document method.
     * 
     * @param maxThreshold
     */
    void SetMaximumSpeciesThreshold(std::vector<double> maxThreshold);

    /**
     * \todo Document method.
     * 
     * @param minThreshold
     */
    void SetMinimumSpeciesThreshold(std::vector<double> minThreshold);

    /**
     * \todo Document method.
     * 
     * @param threshold
     * @param index
     */
    void SetMaximumThresholdByIndex(double threshold, unsigned index);

    /**
     * \todo Document method.
     * 
     * @param threshold
     * @param index
     */
    void SetMinimumThresholdByIndex(double threshold, unsigned index);

    /**
     * \todo Document method.
     * 
     * @param concentrations
     */
    void SetSpeciesConcentrations(std::vector<double> concentrations);

    /**
     * \todo Document method.
     * 
     * @param threshold
     * @param index
     */
    void SetSpeciesConcentrationsByIndex(double threshold, unsigned index);

    /**
     * \todo Document method.
     * 
     * @param speciesNumber
     */
    void SetNumberThresholdSpecies(unsigned speciesNumber);

    /**
     * \todo Document method.
     * 
     * @param thresholdCheck
     */
    void SetMaximumThresholdCheck(std::vector<bool> thresholdCheck);

    /**
     * \todo Document method.
     * 
     * @param thresholdCheck
     */
    void SetMinimumThresholdCheck(std::vector<bool> thresholdCheck);

    /**
     * \todo Document method.
     * 
     * @param thresholdCheck
     * @param index
     */
    void SetMaximumThresholdCheckByIndex(bool thresholdCheck, unsigned index);

    /**
     * \todo Document method.
     * 
     * @param thresholdCheck
     * @param index
     */
    void SetMinimumThresholdCheckByIndex(bool thresholdCheck, unsigned index);

    /**
     * \todo Document method.
     * 
     * @param pChemistry
     */
    void SetThresholdChemistry(AbstractChemistry* pChemistry);

    /**
     * \todo Document method.
     * 
     * @return 
     */
    AbstractChemistry* GetThresholdChemistry();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    std::vector<double> GetMaximumSpeciesThreshold();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    std::vector<double> GetMinimumSpeciesThreshold();

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    double GetMaximumThresholdByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    double GetMinimumThresholdByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @return 
     */
    std::vector<double> GetSpeciesConcentrations();

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    double GetSpeciesConcentrationsByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @return 
     */
    unsigned GetNumberThresholdSpecies();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    std::vector<bool> GetMaximumThresholdCheck();

    /**
     * \todo Document method.
     * 
     * @return 
     */
    std::vector<bool> GetMinimumThresholdCheck();

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    bool GetMaximumThresholdCheckByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    bool GetMinimumThresholdCheckByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @return 
     */
    double GetCriticalStarvationDuration() const;

    /**
     * \todo Document method.
     * 
     * @param criticalStarvationDuration
     * 
     * @return 
     */
    void SetCriticalStarvationDuration(double criticalStarvationDuration);

    /**
     * \todo Document method.
     * 
     * @param currentStarvationOnsetTime
     * 
     * @return 
     */
    void SetCurrentStarvationOnsetTime(double currentStarvationOnsetTime);
};

#endif /* SIMPLECHEMICALTHRESHOLDCELLCYCLEMODEL_HPP_ */