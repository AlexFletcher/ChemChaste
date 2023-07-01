#include "SimpleChemicalThresholdCellCycleFromFile.hpp"

SimpleChemicalThresholdCellCycleFromFile::SimpleChemicalThresholdCellCycleFromFile(std::string thresholdFilename) 
    : SimpleChemicalThresholdCellCycleModel(),
      mCellCycleFilename(thresholdFilename)
{
}

void SimpleChemicalThresholdCellCycleFromFile::SetUp()
{
    std::vector<std::string> threshold_chemical_vector;
    std::vector<double> max_threshold_vector;
    std::vector<double> min_threshold_vector;
    std::vector<bool> max_threshold_check;
    std::vector<bool> min_threshold_check;

    // Read and store the diffusion data as a matrix of strings 
    std::vector<std::vector<std::string> > threshold_database;
    threshold_database = ReadMatrix(mCellCycleFilename);

    // Parse the matrix
    for (unsigned i = 0; i < threshold_database.size(); ++i)
    {
        // The first element of the read in database is the identifying names
        threshold_chemical_vector.push_back(threshold_database[i][0]);

        if (threshold_database[i].size() == 3)
        {
            // Data for both the maximum and minimum thresholds
            if (std::stoul(threshold_database[i][1].c_str()) < std::stoul(threshold_database[i][2].c_str()))
            {
                // Then can't possibly divide as will be apoptotic,so don't check
                max_threshold_vector.push_back(std::stoul(threshold_database[i][1].c_str()));
                min_threshold_vector.push_back(std::stoul(threshold_database[i][2].c_str()));
                max_threshold_check.push_back(false);
                min_threshold_check.push_back(true);
            }
            else
            {
                max_threshold_vector.push_back(std::stoul(threshold_database[i][1].c_str()));
                min_threshold_vector.push_back(std::stoul(threshold_database[i][2].c_str()));
                max_threshold_check.push_back(true);
                min_threshold_check.push_back(true);
            }
        }
        else if (threshold_database[i].size() == 2)
        {
            // Data for only the maximum the thresholds
            max_threshold_vector.push_back(std::stoul(threshold_database[i][1].c_str()));
            min_threshold_vector.push_back(0.0);
            max_threshold_check.push_back(true);
            min_threshold_check.push_back(false);
        }
        else
        {
            // Then no data for either the maximum and minimum thresholds
            max_threshold_vector.push_back(0.0);
            min_threshold_vector.push_back(0.0);
            max_threshold_check.push_back(false);
            min_threshold_check.push_back(false);
        }

        // Error checking
        if (max_threshold_vector[i] < mSmallestMaximumThreshold)
        {
            max_threshold_vector[i] = 0.0;
            max_threshold_check[i] = false;
        }
    }

    // Store the information
    AbstractChemistry* p_threshold_chemistry = new AbstractChemistry();

    for (unsigned i = 0; i < threshold_chemical_vector.size(); ++i)
    {
        p_threshold_chemistry->AddChemical(new AbstractChemical(threshold_chemical_vector[i]));
    }
    
    SimpleChemicalThresholdCellCycleModel::SetUp(p_threshold_chemistry);
    SetUp(p_threshold_chemistry);
    SetMaximumSpeciesThreshold(max_threshold_vector);
    SetMinimumSpeciesThreshold(min_threshold_vector);
    SetNumberThresholdSpecies(threshold_chemical_vector.size());
    SetMaximumThresholdCheck(max_threshold_check);
    SetMinimumThresholdCheck(min_threshold_check);
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
    std::string matrix_cell;

    // Determine the position of the delimiter
    size_t delimiter_pos = line.find(delim);

    bool is_end_of_line = false;
    
    while (!is_end_of_line)
    {
        // While not at the end of the file, sample sub strings from the posiiton of the delimiter
        if (delimiter_pos == std::string::npos)
        {
            is_end_of_line = true;
        }
        
        // Sample substring from begining of the string to the delimiter positioon, store as data entry
        matrix_cell = line.substr(0, delimiter_pos);

        // Eemove the sampled entry from the string
        line = line.substr(delimiter_pos + 1, std::string::npos);

        row_vector.push_back(matrix_cell);

        // Update delimiter position
        delimiter_pos = line.find(delim);
    }
    return row_vector;
}