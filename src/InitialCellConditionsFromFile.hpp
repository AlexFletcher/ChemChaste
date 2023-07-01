#ifndef INITIALCELLCONDITIONSFROMFILE_HPP_
#define INITIALCELLCONDITIONSFROMFILE_HPP_

#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <stdlib.h> 

#include "RandomNumberGenerator.hpp"
#include "AbstractChemistry.hpp"

/**
 * \todo Document this class.
 */
class InitialCellConditionsFromFile
{
protected:

    std::vector<std::string> mChemicalNamesVector;

    std::vector<double> mConcentrationVector;

    std::vector<bool> mPerturbConcentrationVector;

    std::vector<std::string> mBoolTestDictionary = {"true","True","TRUE","1"};

public:

    InitialCellConditionsFromFile(std::string);

    virtual ~InitialCellConditionsFromFile()
    {
    };

    std::vector<std::vector<std::string>> ReadMatrix(std::string);

    std::vector<std::string> parseMatrixLineString(std::string);

    bool StringToBool(std::string);

    void SetChemicalNamesVector(std::vector<std::string>);

    void SetConcentrationVector(std::vector<double>);

    void SetPerturbConcentrationVector(std::vector<bool>);

    std::vector<std::string> GetChemicalNamesVector();

    std::vector<double> GetConcentrationVector();

    std::vector<bool> GetPerturbConcentrationVector();
};

#endif