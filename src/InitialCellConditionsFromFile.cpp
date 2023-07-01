#include "InitialCellConditionsFromFile.hpp"

InitialCellConditionsFromFile::InitialCellConditionsFromFile(std::string conditionFilename) 
{
    std::vector<std::string> chemicalNamesVector;
    std::vector<double> concentrationVector;
    std::vector<bool> perturbConcentrationVector;

    // read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string>> conditionsDatabase;
    conditionsDatabase = ReadMatrix(conditionFilename);

    // parse the matrix
    for (unsigned i=0; i<conditionsDatabase.size(); ++i)
    {
        // the first element of the read in database is the identifying names
        chemicalNamesVector.push_back(conditionsDatabase[i][0]);

        if (conditionsDatabase[0].size()==2)
        {
            // data for species concentration
            concentrationVector.push_back(atof(conditionsDatabase[i][1].c_str()));
        }
        else if (conditionsDatabase[0].size()==3)
        {
            // data for species concentration and whether to perturb the initial value
            
            perturbConcentrationVector.push_back(StringToBool(conditionsDatabase[i][2]));
            if (perturbConcentrationVector[i]==true)
            {
                concentrationVector.push_back(fabs(atof(conditionsDatabase[i][1].c_str()) + RandomNumberGenerator::Instance()->ranf()));
            }
            else
            {
                concentrationVector.push_back(atof(conditionsDatabase[i][1].c_str()));
            }
        }
    }

    // Store the information
    SetChemicalNamesVector(chemicalNamesVector);
    SetConcentrationVector(concentrationVector);
    SetPerturbConcentrationVector(perturbConcentrationVector);
}

std::vector<std::vector<std::string>> InitialCellConditionsFromFile::ReadMatrix(std::string filename)
{
    // parse a matrix file (.csv) line by line, ignore escape line,s containing file information
    // that is lines starting with '#' 
    
    std::string line;
    std::ifstream inputFile(filename);

    // read all data types as std::string therefore return the matrix of strings for personalised
    // methods down stream
    std::vector<std::vector<std::string>> outputMatrix = std::vector<std::vector<std::string>>();

    // check file exists and is openable
    if (inputFile.is_open())
    {
        // open the matrix file
        while (getline(inputFile,line))
        {
            // while the file still has lines not read.
            // read line left to right, top to bottom.
            if (!line.empty())
            {
                if (line.at(0)=='#')
                {
                }
                else
                {
                    outputMatrix.push_back(parseMatrixLineString(line));
                }   
            }
        }
        inputFile.close();

        return outputMatrix;
    }
    else
    {
        ///\todo replace with EXCEPTION
        std::cout << "Error: Unable to open file: "<<filename << std::endl;
        return outputMatrix;
    }
}

std::vector<std::string> InitialCellConditionsFromFile::parseMatrixLineString(std::string line)
{
    // for a line string in the matrix read, parse into vector data entries based on delimiters ','
    std::vector<std::string> rowVector = std::vector<std::string>();

    // delimiter, may be modified by further methods
    std::string delim = ",";
    std::string matrixCell;

    // determine the position of the delimiter
    size_t delimiter_pos=line.find(delim);

    bool IsEndOfLine = false;
    
    while (!IsEndOfLine)
    {
        // while not at the end of the file, sample sub strings from the posiiton of the delimiter
        if (delimiter_pos == std::string::npos)
        {
            IsEndOfLine = true;
        }
        
        // sample substring from begining of the string to the delimiter positioon, store as data entry
        matrixCell = line.substr(0,delimiter_pos);

        // remove the sampled entry from the string
        line = line.substr(delimiter_pos+1,std::string::npos);

        rowVector.push_back(matrixCell);

        // update delimiter position
        delimiter_pos = line.find(delim);
    }
    return rowVector;
}

bool InitialCellConditionsFromFile::StringToBool(std::string test_string)
{
    // take in a test string and determine if a true character string is present

    for (unsigned i=0; i<mBoolTestDictionary.size(); ++i)
    {
        if (test_string == mBoolTestDictionary[i])
        {
            return true;
        }
    }

    return false;
}

void InitialCellConditionsFromFile::SetChemicalNamesVector(std::vector<std::string> namesVector)
{
    mChemicalNamesVector = namesVector;
}

void InitialCellConditionsFromFile::SetConcentrationVector(std::vector<double> concentrationVector)
{
    mConcentrationVector = concentrationVector;
}

void InitialCellConditionsFromFile::SetPerturbConcentrationVector(std::vector<bool> perturbConcentrationVector)
{
    mPerturbConcentrationVector = perturbConcentrationVector;
}

std::vector<std::string> InitialCellConditionsFromFile::GetChemicalNamesVector()
{
    return mChemicalNamesVector;
}

std::vector<double> InitialCellConditionsFromFile::GetConcentrationVector()
{
    return mConcentrationVector;
}

std::vector<bool> InitialCellConditionsFromFile::GetPerturbConcentrationVector()
{
    return mPerturbConcentrationVector;
}