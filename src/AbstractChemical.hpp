#ifndef ABSTRACTCHEMICAL_HPP_
#define ABSTRACTCHEMICAL_HPP_

#include <string>

/**
 * Abstract class used to contain information about a chemical species in 
 * ChemChaste. Contains the name of the chemical and basic chemcial properties.
 */
class AbstractChemical
{
private:

    /** The name of the chemical. Defaults to empty string. */
    std::string mChemicalName;

    /** The molecular size of the chemical. Defaults to 1.0. */
    double mSize;

    /** The molecular mass of the chemical. Defaults to 1.0. */
    double mMass;

    /** Molecular charge for case of electrical reactions. Defaults to 0. */
    int mValence;

    /**
     * Numerical dimensions for the chemical traits, used in the ODE information 
     * to define the states. Defaults to "non-dim".
     */
    std::string mChemicalDimensions;

    /**
     * Whether the Gibbs formation energy is known. Defaults to false.
     */
    bool mFormationKnown;

    /**
     * The Gibbs formation energy. Defaults to 0.
     */
    double mFormationGibbs;

public:

    /**
     * Default constructor.
     * 
     * @param chemicalName the name of the chemical. Defaults to empty string.
     * @param size The molecular size of the chemical. Defaults to 1.0.
     * @param mass The molecular mass of the chemical. Defaults to 1.0.
     * @param valence Molecular charge for case of electrical reactions. 
     *                Defaults to 0.
     * @param chemicalDimensions Numerical dimensions for the chemical traits. 
     *                           Defaults to 1.0.
     */
    AbstractChemical(std::string chemicalName="", 
                     double size=1.0,
                     double mass=1.0,
                     int valence=0,
                     std::string chemicalDimensions="non-dim");

    /**
     * Destructor.
     */
    virtual ~AbstractChemical();

    /**
     * Set mChemicalName.
     * 
     * @param chemicalName The name of the chemical.
     */
    void SetChemicalName(std::string chemicalName);

    /**
     * Set mSize.
     * 
     * @param size The molecular size of the chemical.
     */
    void SetChemicalSize(double size);

    /**
     * Set mMass.
     * 
     * @param mass The molecular mass of the chemical.
     */
    void SetChemicalMass(double mass);

    /**
     * Set mValence.
     * 
     * @param valence Molecular charge for case of electrical reactions.
     */
    void SetChemicalValence(int valence);
    
    /**
     * Set mChemicalDimensions.
     * 
     * @param chemicalDimensions Numerical dimensions for the chemical traits.
     */
    void SetChemicalDimensions(std::string chemicalDimensions);

    /**
     * Set mFormationKnown.
     * 
     * @param formationKnown Whether the Gibbs formation energy is known.
     */
    void SetChemicalFormationKnown(bool formationKnown);

    /**
     * Set mFormationGibbs.
     * 
     * @param formationGibbs The Gibbs formation energy.
     */
    void SetChemicalFormationGibbs(double formationGibbs);

    /**
     * @return mChemicalName.
     */
    std::string GetChemicalName();

    /**
     * @return mSize.
     */
    double GetChemicalSize();

    /**
     * @return mMass.
     */
    double GetChemicalMass();

    /**
     * @return mCValence.
     */
    int GetChemicalValence();

    /**
     * @return mChemicalDimensions.
     */
    std::string GetChemicalDimensions();

    /**
     * @return mFormationKnown.
     */
    bool IsChemicalFormationKnown();

    /**
     * @return mFormationGibbs.
     */
    double GetChemicalFormationGibbs();

    /**
     * @return "AbstractChemical".
     */
    virtual std::string GetChemicalType();
};

#endif /* ABSTRACTCHEMICAL_HPP_ */
