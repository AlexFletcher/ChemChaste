#ifndef ABSTRACTCHEMISTRY_HPP_
#define ABSTRACTCHEMISTRY_HPP_

#include <iostream>
#include "AbstractChemical.hpp"

/**
 * Abstract property to contain information about the collection of chemical 
 * species in the system container operator for adding two 
 * AbstractReactionSystems to form union for bulk copy contructor for mutating 
 * reaction network.
 * 
 * \todo inherit? : public AbstractChemical
 */
class AbstractChemistry 
{
protected:

    /**
     * Number of chemicals.
     * Initialised to zero in constructor.
     */
    unsigned mNumberChemicals;

    /**
     * Vector of pointers to chemicals. 
     * Initialised to empty vector in constructor.
     */
    std::vector<AbstractChemical*> mChemicalVector;

    /**
     * Vector of names of chemicals. 
     * Initialised to empty vector in constructor.
     */
    std::vector<std::string> mChemicalNames;

    /**
     * Vector of chemical dimensions. 
     * Initialised to empty vector in constructor.
     */
    std::vector<std::string> mChemicalDimensions;
    
public:

    /**
     * Constructor.
     */
    AbstractChemistry();

    /**
     * Destructor.
     */
    virtual ~AbstractChemistry();

    /**
     * Check whether a given chemical is already present, and if it isn't, 
     * add it to mChemicalVector.
     *
     * @param pChemical pointer to a chemical
     */
    virtual void AddChemical(AbstractChemical* pChemical);

    /**
     * @param pChemical pointer to a chemical
     *
     * @return whether pChemical is present.
     */
    virtual bool CheckChemical(AbstractChemical* pChemical);

    /**
     * \todo document method
     *
     * @param pChemical pointer to a chemical
     */
    virtual void UpdateChemicalVectors(AbstractChemical* pChemical);

    /**
     * \todo document method
     *
     * @param chemicalVector vector of pointers to chemicals
     */
    void SetChemicalVector(std::vector<AbstractChemical*> chemicalVector);

    /**
     * \todo document method
     *
     * @param chemicalNames vector of names of chemicals
     */
    void SetChemicalNames(std::vector<std::string> chemicalNames);

    /**
     * \todo document method
     *
     * @param chemicalDimensions vector of chemical dimensions
     */
    void SetChemicalDimensions(std::vector<std::string> chemicalDimensions);

    /**
     * @return mChemicalDimensions
     */
    std::vector<AbstractChemical*> rGetChemicalVector();

    /**
     * @return mChemicalNames
     */
    std::vector<std::string> GetChemicalNames();

    /**
     * \todo throw exception if no chemical has this index
     *
     * @param index \todo document argument
     *
     * @return \todo document return type
     */
    std::string GetChemicalNamesByIndex(unsigned index);

    /**
     * \todo throw exception if no chemical has this name
     *
     * @param name name of a chemical
     *
     * @return pointer to the chemical with this name
     */
    unsigned GetChemicalIndexByName(std::string name);

    /**
     * @return mChemicalDimensions
     */
    std::vector<std::string> GetChemicalDimensions();
    /**
     * \todo throw exception if no chemical has this index
     *
     * @param index \todo document argument
     *
     * @return \todo document return type
     */
    std::string GetChemicalDimensionsByIndex(unsigned index);

    /**
     * \todo document method
     *
     * @return \todo document return type
     */
    virtual std::string GetChemistryType();

    /**
     * \todo document method
     *
     * @param pChemistry pointer to a chemistry
     */
    virtual void AddChemistry(AbstractChemistry* pChemistry);

    /**
     * @return mNumberChemicals
     */
    unsigned GetNumberChemicals();
};

#endif /* ABSTRACTCHEMISTRY_HPP_ */
