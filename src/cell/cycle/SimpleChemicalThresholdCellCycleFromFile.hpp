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
     * @param 
    */
    SimpleChemicalThresholdCellCycleFromFile(std::string);

    /**
     * Destructor.
     */
    virtual ~SimpleChemicalThresholdCellCycleFromFile()
    {
    };

    /**
     * \todo Document method.
     */
    virtual void SetUp();

    /**
     * \todo Document method.
     */
    std::vector<std::vector<std::string> > ReadMatrix(std::string);

    /**
     * \todo Document method.
     */
    std::vector<std::string> parseMatrixLineString(std::string);
};

SimpleChemicalThresholdCellCycleFromFile::SimpleChemicalThresholdCellCycleFromFile(std::string thresholdFilename) 
    : SimpleChemicalThresholdCellCycleModel(),
      mCellCycleFilename(thresholdFilename)
{
}

void SimpleChemicalThresholdCellCycleFromFile::SetUp()
{
    std::vector<std::string> thresholdChemicalVector;
    std::vector<double> maximumThresholdVector;
    std::vector<double> minimumThresholdVector;
    std::vector<bool> maximumThresholdCheck;
    std::vector<bool> minimumThresholdCheck;

    // Read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> thresholdDatabase;
    thresholdDatabase = ReadMatrix(mCellCycleFilename);

    // Parse the matrix
    for (unsigned i=0; i<thresholdDatabase.size(); i++)
    {
        // The first element of the read in database is the identifying names
        thresholdChemicalVector.push_back(thresholdDatabase[i][0]);

        if (thresholdDatabase[i].size() == 3)
        {
            // Data for both the maximum and minimum thresholds

            if (std::stoul(thresholdDatabase[i][1].c_str()) < std::stoul(thresholdDatabase[i][2].c_str()))
            {
                // Then can't possibly divide as will be apoptotic,so don't check
                maximumThresholdVector.push_back(std::stoul(thresholdDatabase[i][1].c_str()));
                minimumThresholdVector.push_back(std::stoul(thresholdDatabase[i][2].c_str()));
                maximumThresholdCheck.push_back(false);
                minimumThresholdCheck.push_back(true);
            }
            else
            {
                maximumThresholdVector.push_back(std::stoul(thresholdDatabase[i][1].c_str()));
                minimumThresholdVector.push_back(std::stoul(thresholdDatabase[i][2].c_str()));
                maximumThresholdCheck.push_back(true);
                minimumThresholdCheck.push_back(true);
            }
        }
        else if (thresholdDatabase[i].size() == 2)
        {
            // Data for only the maximum the thresholds
            maximumThresholdVector.push_back(std::stoul(thresholdDatabase[i][1].c_str()));
            minimumThresholdVector.push_back(0.0);
            maximumThresholdCheck.push_back(true);
            minimumThresholdCheck.push_back(false);
        }
        else
        {
            // Then no data for either the maximum and minimum thresholds
            maximumThresholdVector.push_back(0.0);
            minimumThresholdVector.push_back(0.0);
            maximumThresholdCheck.push_back(false);
            minimumThresholdCheck.push_back(false);
        }
        // Error checking
        if (maximumThresholdVector[i] < mSmallestMaximumThreshold)
        {
            maximumThresholdVector[i] = 0.0;
            maximumThresholdCheck[i] = false;
        }
    }

    // store the information
    AbstractChemistry* pThresholdChemistry = new AbstractChemistry();

    for (unsigned i = 0; i < thresholdChemicalVector.size(); i++)
    {
        pThresholdChemistry->AddChemical(new AbstractChemical(thresholdChemicalVector[i]));
    }
    
    SimpleChemicalThresholdCellCycleModel::SetUp(pThresholdChemistry);
    SetUp(pThresholdChemistry);
    SetMaximumSpeciesThreshold(maximumThresholdVector);
    SetMinimumSpeciesThreshold(minimumThresholdVector);
    SetNumberThresholdSpecies(thresholdChemicalVector.size());
    SetMaximumThresholdCheck(maximumThresholdCheck);
    SetMinimumThresholdCheck(minimumThresholdCheck);
}

std::vector<std::vector<std::string>> SimpleChemicalThresholdCellCycleFromFile::ReadMatrix(std::string filename)
{
    /*
     * Parse a matrix file (.csv) line by line, ignore escape lines, containing 
     * file information that is lines starting with '#' 
     */
    
    std::string line;
    std::ifstream inputFile(filename);

    /*
     * Read all data types as std::string therefore return the matrix of strings 
     * for personalised methods down stream.
     */
    std::vector<std::vector<std::string> > output_matrix = std::vector<std::vector<std::string>>();

    // Check file exists and is openable
    if (inputFile.is_open())
    {
        // Open the matrix file
        while (getline(inputFile, line))
        {
            // While the file still has lines not read, read line left to right, top to bottom
            if (!line.empty())
            {
                if (line.at(0) != '#')
                {
                    output_matrix.push_back(parseMatrixLineString(line));
                }   
            }
        }
        inputFile.close();
        return output_matrix;
    }
    else
    {
        ///\todo replace with EXCEPTION
        std::cout << "Error: Unable to open file: " << filename << std::endl;
        return output_matrix;
    }
}

std::vector<std::string> SimpleChemicalThresholdCellCycleFromFile::parseMatrixLineString(std::string line)
{
    // For a line string in the matrix read, parse into vector data entries based on delimiters ','
    std::vector<std::string> row_vector = std::vector<std::string>();

    // Delimiter, may be modified by further methods
    std::string delim = ",";
    std::string matrixCell;

    // Determine the position of the delimiter
    size_t delimiter_pos = line.find(delim);

    bool IsEndOfLine = false;
    
    while (!IsEndOfLine)
    {
        // While not at the end of the file, sample sub strings from the posiiton of the delimiter
        if (delimiter_pos == std::string::npos)
        {
            IsEndOfLine = true;
        }
        
        // Sample substring from begining of the string to the delimiter positioon, store as data entry
        matrixCell = line.substr(0, delimiter_pos);

        // Eemove the sampled entry from the string
        line = line.substr(delimiter_pos + 1, std::string::npos);

        row_vector.push_back(matrixCell);

        // Update delimiter position
        delimiter_pos = line.find(delim);
    }
    return row_vector;
}

#endif /* SIMPLECHEMICALTHRESHOLDCELLCYCLEFROMFILE_HPP_ */