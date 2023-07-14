#ifndef CHEMICALCELLPROPERTY_HPP_
#define CHEMICALCELLPROPERTY_HPP_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "AbstractCellProperty.hpp"
#include "StateVariableRegister.hpp"

/**
 * Cell property to handle the basic chemical property required for cell 
 * simulations.
 */
class ChemicalCellProperty : public AbstractCellProperty
{
protected:

    /** Register for the state variables present in the cell. */
    StateVariableRegister* mpStateVariableRegister;

    /** Concentration vector of the cell, would be derived from CellData. */
    std::vector<double> mConcentrationVector;

public:

    /**
     * Constructor.
     */
    ChemicalCellProperty();

    /**
     * Destructor.
     */
    virtual ~ChemicalCellProperty();

    /**
     * \todo document method
     * 
     * @param
     */
    virtual void InitialiseCell(std::vector<std::string> pRegister, 
                                std::vector<double> concentrationVector);

    /**
     * \todo document method
     * 
     * @param
     */
    virtual void InitialiseCell(StateVariableRegister*, std::vector<double>);

    /**
     * \todo document method
     * 
     * @param rConcentrationVector
     */
    virtual void UpdateCellConcentrationVector(
        std::vector<double>& rConcentrationVector);

    /**
     * \todo document method
     * 
     * @param pRegister
     */
    void SetStateVariableRegister(StateVariableRegister* pRegister);

    /**
     * @return mpStateVariableRegister
     */
    StateVariableRegister* GetStateVariableRegister();

    /**
     * @return mConcentrationVector
     */
    std::vector<double> GetCellConcentrationVector();

    /**
     * \todo document method
     * 
     * @param index
     */
    double GetCellConcentrationByIndex(unsigned index);

    /**
     * \todo document method
     * 
     * @param name
     */
    double GetCellConcentrationByName(std::string name);
};

#endif /* CHEMICALCELLPROPERTY_HPP_ */