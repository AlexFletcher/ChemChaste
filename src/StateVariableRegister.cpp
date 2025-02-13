#include "StateVariableRegister.hpp"

StateVariableRegister::StateVariableRegister(std::vector<std::string> VariableNameVector)
    : mStateVariableRegister(VariableNameVector)
{
    mNumberSystemStateVariables = VariableNameVector.size();
}

StateVariableRegister::~StateVariableRegister()
{
}

void StateVariableRegister::UpdateStateVariableRegister(std::vector<std::string> variableRegister)
{
    mStateVariableRegister = variableRegister;
    SetNumStateVariables(mStateVariableRegister.size());
}

std::vector<std::string> StateVariableRegister::GetStateVariableRegisterVector()
{
    return mStateVariableRegister;
}

void StateVariableRegister::AddStateVariableVector(std::vector<std::string> variableRegister)
{
    for (unsigned i=0; i<variableRegister.size(); ++i)
    {
        // check candidate variable is indeed new, prevent duplicates
        if (!IsStateVariablePresent(variableRegister[i]))
        {
            mStateVariableRegister.push_back(variableRegister[i]);
            SetNumStateVariables(mStateVariableRegister.size());
        }
    }
}

void StateVariableRegister::AddStateVariable(std::string new_variable_name)
{
    // check candidate variable is indeed new, prevent duplicates
    if (!IsStateVariablePresent(new_variable_name))
    {
        mStateVariableRegister.push_back(new_variable_name);
        SetNumStateVariables(mStateVariableRegister.size());
    }
}

void StateVariableRegister::RemoveStateVariable(std::string remove_variable_name)
{
    // check candidate variable is indeed present
    if (IsStateVariablePresent(remove_variable_name))
    {
        this->mStateVariableRegister.erase( this->mStateVariableRegister.begin() + RetrieveStateVariableIndex(remove_variable_name));
        SetNumStateVariables(this->mStateVariableRegister.size());
    }
}

bool StateVariableRegister::IsStateVariablePresent(std::string variable_name)
{
    // function to test if a variable is present in the vairable register, used to locate or add new variable
    bool is_present = false;
    for (unsigned variable_index=0; variable_index<this->mStateVariableRegister.size(); variable_index++)
    {
        if (this->mStateVariableRegister[variable_index] == variable_name)
        {
            is_present = true;
            break;
        }
    }
    return is_present;
}

unsigned StateVariableRegister::RetrieveStateVariableIndex(std::string variable_name)
{
    // only use when know variable is present
    unsigned index =0;
    for (unsigned variable_index=0; variable_index<this->mStateVariableRegister.size(); variable_index++)
    {
        if (this->mStateVariableRegister[variable_index] == variable_name)
        {
            index=variable_index;
            break;
        }
    }
    return index;
}

std::string StateVariableRegister::RetrieveStateVariableName(unsigned index)
{
    if (index <= mNumberSystemStateVariables)
    {
        return mStateVariableRegister[index];
    }
    else
    {
        return "Error";   
    }
}

void StateVariableRegister::SetNumStateVariables(unsigned number_of_state_variables)
{
    mNumberSystemStateVariables = number_of_state_variables;
}

unsigned StateVariableRegister::GetNumStateVariables()
{
    return mNumberSystemStateVariables;
}

std::vector<unsigned> StateVariableRegister::FindIndicesInThisRegister(StateVariableRegister* p_new_register)
{
    // the desired variable isn't found in this system, store a number that can't be accessed, i.e larger than the number of variables
    // need to implment check for this
    std::vector<unsigned> matchedIndices(p_new_register->GetNumStateVariables(),this->mNumberSystemStateVariables);

    for (unsigned i=0; i<p_new_register->GetNumStateVariables(); ++i)
    {
        
        if (IsStateVariablePresent(p_new_register->RetrieveStateVariableName(i)))
        {   
            matchedIndices[i] = this->RetrieveStateVariableIndex(p_new_register->RetrieveStateVariableName(i));
        }
    }

    return matchedIndices;
}
    
std::vector<unsigned> StateVariableRegister::FindIndicesInThatRegister(StateVariableRegister* p_new_register)
{
    // the desired variable isn't found in this system, store a number that can't be accessed, i.e larger than the number of variables
    // need to implment check for this
    std::vector<unsigned> matchedIndices(this->mNumberSystemStateVariables,p_new_register->GetNumStateVariables());

    for (unsigned i=0; i<this->mNumberSystemStateVariables; ++i)
    {
        if (p_new_register->IsStateVariablePresent(this->mStateVariableRegister[i]))
        {
            matchedIndices[i] = p_new_register->RetrieveStateVariableIndex(this->mStateVariableRegister[i]);
        }
    }

    return matchedIndices;
}
    
std::vector<std::string> StateVariableRegister::FindCommonNamesInRegisters(StateVariableRegister* p_new_register)
{
    std::vector<std::string> matchedNames;

    for (unsigned i=0; i<this->mNumberSystemStateVariables; ++i)
    {
        if (p_new_register->IsStateVariablePresent(this->mStateVariableRegister[i]))
        {
            matchedNames.push_back( mStateVariableRegister[i]);
        }
    }

    return matchedNames;
}