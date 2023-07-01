#ifndef MASSACTIONREACTION_HPP_
#define MASSACTIONREACTION_HPP_

#include <string>
#include <tuple>
#include <vector>
#include <cmath>

#include "AbstractChemical.hpp"
#include "AbstractReversibleReaction.hpp"

// mass action reactions are generally reversible

/**
 * \todo Document this class.
 */
class MassActionReaction : public AbstractReversibleReaction
{
private:

    using AbstractReversibleReaction::UpdateReaction;
    using AbstractReversibleReaction::UpdateReactionRate;
    using AbstractReversibleReaction::GetReactionType;
    using AbstractReversibleReaction::ParseReactionInformation;

    double mReactionGibbs;
    
    double mRkj = 8.3144598e-3;
    double mTemp = 300; // kelvin

    bool mIsGibbs;

    bool mIsReversible;

    double mGibbsFreeEnergy = 0.0;

    double mForwardReactionRateConstant;

    double mReverseReactionRateConstant;

    std::string mGibbsDelimiter = "deltaG =";

public:

    MassActionReaction( std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(),
                        std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(),
                        std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(),
                        std::vector<unsigned> stoichProducts = std::vector<unsigned>(),
                        bool IsGibbs = false,
                        bool IsReversible = true,
                        double ForwardReactionRateConstant = 1.0,
                        double ReverseReactionRateConstant = 1.0);

    virtual ~MassActionReaction()
    {
    };

    virtual void UpdateReaction();

    virtual void UpdateReactionRate(AbstractChemistry*, const std::vector<double>&);

    virtual std::string GetReactionType();

    virtual void ParseReactionInformation(std::string, bool);

    double GetGibbsFreeEnergy();

    void SetGibbsFreeEnergy(double);

    double CalculateReactionQuotient(AbstractChemistry* , const std::vector<double>& );

    double DeltaGtoKf(double deltaG, double Kr = 1.0);

    double KftodeltaG(double Kf, double Kr = 1.0);

    double CalculateGibbsFromQuotient(double G0, double Q);

    void SetReactionTemperature(double);

    double GetReactionTemperature();

    double GetForwardReactionRateConstant();

    void SetForwardReactionRateConstant(double);
    
    double GetReverseReactionRateConstant();

    void SetReverseReactionRateConstant(double);

    void SetGibbsDelmiter(std::string);

    std::string GetGibbsDelimiter();
};

#endif /* MASSACTIONREACTION_HPP_ */