#include "AbstractChemical.hpp"

AbstractChemical::AbstractChemical(std::string chemicalName, 
    double size,
    double mass,
    int valence,
    std::string chemicalDimensions)
    : mChemicalName(chemicalName),
      mSize(size),
      mMass(mass),
      mValence(valence),
      mChemicalDimensions(chemicalDimensions)
{
    // default to not knowing the formation energy
    SetChemicalFormationGibbs(0.0);
    SetChemicalFormationKnown(false);
}

void AbstractChemical::SetChemicalName(std::string chemicalName)
{
    mChemicalName = chemicalName;
}

void AbstractChemical::SetChemicalSize(double size)
{
    mSize = size;
}

void AbstractChemical::SetChemicalMass(double Mass)
{
    mMass  =Mass;
}

void AbstractChemical::SetChemicalValence(int valence)
{
    mValence = valence;
}

void AbstractChemical::SetChemicalDimensions(std::string chemicalDimensions)
{
    mChemicalDimensions = chemicalDimensions;
}

void AbstractChemical::SetChemicalFormationKnown(bool formationKnown)
{
    mFormationKnown = formationKnown;
}

void AbstractChemical::SetChemicalFormationGibbs(double formationGibbs)
{
    mFormationGibbs = formationGibbs;
}

std::string AbstractChemical::GetChemicalName()
{
    return mChemicalName;
}

double AbstractChemical::GetChemicalSize()
{
    return mSize;
}

double AbstractChemical::GetChemicalMass()
{
    return mMass;
}

int AbstractChemical::GetChemicalValence()
{
    return mValence;
}

std::string AbstractChemical::GetChemicalDimensions()
{
    return mChemicalDimensions;
}

bool AbstractChemical::IsChemicalFormationKnown()
{
    return mFormationKnown;
}

double AbstractChemical::GetChemicalFormationGibbs()
{
    return mFormationKnown;
}

std::string GetChemicalType()
{
    return "AbstractChemical";
}