#ifndef SIMPLECHEMICALTHRESHOLDCELLCYCLEFROMFILE_HPP_
#define SIMPLECHEMICALTHRESHOLDCELLCYCLEFROMFILE_HPP_

#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#include "AbstractChemistry.hpp"
#include "SimpleChemicalThresholdCellCycleModel.hpp"
 
/**
 * \todo Document this class.
 */
class SimpleChemicalThresholdCellCycleFromFile : public SimpleChemicalThresholdCellCycleModel
{
protected:

    using SimpleChemicalThresholdCellCycleModel::SetUp;

    /** \todo Document member variable. */
    std::string mCellCycleFilename;

    /** \todo Document member variable. */
    double mSmallestMaximumThreshold = 1e-6;

public:

    /**
     * Default constructor.
     * 
     * \todo Document method.
     * 
     * @param fileName
    */
    SimpleChemicalThresholdCellCycleFromFile(std::string fileName);

    /**
     * Destructor.
     */
    virtual ~SimpleChemicalThresholdCellCycleFromFile();

    /**
     * \todo Document method.
     */
    virtual void SetUp();

    /**
     * \todo Document method.
     * 
     * @param fileName
     */
    std::vector<std::vector<std::string> > ReadMatrix(std::string fileName);

    /**
     * \todo Document method.
     * 
     * @param line 
     */
    std::vector<std::string> ParseMatrixLineString(std::string line);
};

#endif /* SIMPLECHEMICALTHRESHOLDCELLCYCLEFROMFILE_HPP_ */