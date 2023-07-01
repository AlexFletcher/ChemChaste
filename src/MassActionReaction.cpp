#include "MassActionReaction.hpp"

MassActionReaction::MassActionReaction(std::vector<AbstractChemical*> substrates,
                        std::vector<AbstractChemical*> products,
                        std::vector<unsigned> stoichSubstrates,
                        std::vector<unsigned> stoichProducts,
                        bool IsGibbs,
                        bool IsReversible,
                        double ForwardReactionRateConstant, // takes the form of the gibbs energy
                        double ReverseReactionRateConstant)
    :   AbstractReversibleReaction(substrates,
                                   products,
                                   stoichSubstrates,
                                   stoichProducts,
                                   ForwardReactionRateConstant,
                                   ReverseReactionRateConstant),
        mIsGibbs(IsGibbs),
        mIsReversible(IsReversible),
        mForwardReactionRateConstant(ForwardReactionRateConstant),
        mReverseReactionRateConstant(ReverseReactionRateConstant)
{
    if (mIsGibbs)
    {
        mGibbsFreeEnergy = ForwardReactionRateConstant;
    }
    if (!mIsReversible)
    {
        mReverseReactionRateConstant = 0.0;
    }
}

void MassActionReaction::UpdateReaction()
{
    return;
}

void MassActionReaction::UpdateReactionRate(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{

    double kf = mForwardReactionRateConstant;
    double kr = mReverseReactionRateConstant;

    if (mIsGibbs)
    {
        kf = DeltaGtoKf(mGibbsFreeEnergy, kr);
    }

    //double reaction_quotient = CalculateReactionQuotient(systemChemistry, currentSystemConc);

    double forwardFlux=1.0;
    double reverseFlux=1.0;

    std::vector<AbstractChemical*> p_chemical_vector = systemChemistry->rGetChemicalVector();
    unsigned index = 0;
    for (std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
            chem_iter != p_chemical_vector.end();
            ++chem_iter, ++index)
    {
        AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

        for (unsigned j=0; j<mNumSubstrates; j++)
        {
            if (mpSubstrates[j]->GetChemicalName()==p_system_chemical->GetChemicalName())
            {
                forwardFlux *=  std::pow(currentSystemConc[index],mStoichSubstrates[j]);
                break;
            }
        }
        for (unsigned j=0; j<mNumProducts; j++)
        {
            if (mpProducts[j]->GetChemicalName()==p_system_chemical->GetChemicalName())
            {
                reverseFlux *=  std::pow(currentSystemConc[index],mStoichProducts[j]);
                break;
            }
        }
    }  

    SetForwardReactionRate(kf*forwardFlux);
    SetReverseReactionRate(kr*reverseFlux);
}

std::string MassActionReaction::GetReactionType()
{
    return "MassActionReaction";
}


void MassActionReaction::ParseReactionInformation(std::string reaction_information, bool IsReversible=false)
{
    mIsReversible = IsReversible;

    if (!mIsReversible)
    {
        if (reaction_information.find(mIrreversibleRateName) != std::string::npos)
        {

            size_t pos= reaction_information.find(mIrreversibleRateName);
            mForwardReactionRateConstant = atof(reaction_information.substr(pos+mIrreversibleRateName.size()+1,std::string::npos).c_str());
            mReverseReactionRateConstant = 0.0;
        }
    }
    else
    {  
        if (reaction_information.find(mGibbsDelimiter) != std::string::npos)
        {
            size_t pos = reaction_information.find(mGibbsDelimiter);
            std::cout << "Gibbs raw: "<<reaction_information.substr(pos+mGibbsDelimiter.size()+1,std::string::npos).c_str() << std::endl;
            mGibbsFreeEnergy = atof(reaction_information.substr(pos+mGibbsDelimiter.size()+1,std::string::npos).c_str());
            mIsGibbs = true;
            std::cout << "Gibbs translated: " << mGibbsFreeEnergy << std::endl;
        }
        else
        {
            // abstractReversibleReaction set the reaction rates not constants so update
            AbstractReversibleReaction::ParseReactionInformation(reaction_information,mIsReversible);
            SetForwardReactionRateConstant(GetForwardReactionRate());
            SetReverseReactionRateConstant(GetReverseReactionRate());
        }   
    }
}

double MassActionReaction::CalculateReactionQuotient(AbstractChemistry* systemChemistry, const std::vector<double>& currentSystemConc)
{
    double quotient = 1.0;
    
    if (mNumSubstrates ==0 || mNumProducts==0)
    {
        quotient = 0.0;
    }
    else
    {
        double products_concentrations = 1.0;
        double substrates_concentrations = 1.0;
    
        // need to check against the concentration of each chemical in the system
        std::vector<AbstractChemical*> p_chemical_vector = systemChemistry->rGetChemicalVector();
        unsigned index = 0;
        for (std::vector<AbstractChemical*>::iterator chem_iter = p_chemical_vector.begin();
                chem_iter != p_chemical_vector.end();
                ++chem_iter, ++index)
        {
            AbstractChemical *p_system_chemical = dynamic_cast<AbstractChemical*>(*chem_iter);

            for (unsigned j=0; j<mNumSubstrates; j++)
            {
                if (mpSubstrates[j]->GetChemicalName()==p_system_chemical->GetChemicalName())
                {
                    substrates_concentrations *=  std::pow(currentSystemConc[index],mStoichSubstrates[j]);
                    break;
                }
            }
            for (unsigned j=0; j<mNumProducts; j++)
            {
                if (mpProducts[j]->GetChemicalName()==p_system_chemical->GetChemicalName())
                {
                    products_concentrations *=  std::pow(currentSystemConc[index],mStoichProducts[j]);
                    break;
                }
            }
        }    

        quotient = products_concentrations/substrates_concentrations;
    }

    return quotient;
}

double MassActionReaction::GetGibbsFreeEnergy()
{
    return mGibbsFreeEnergy;
}

void MassActionReaction::SetGibbsFreeEnergy(double gibbsEnergy)
{
    mGibbsFreeEnergy = gibbsEnergy;
}

double MassActionReaction::DeltaGtoKf(double deltaG, double Kr)
{
    return Kr*exp(-deltaG/(mRkj*mTemp));
}

double MassActionReaction::KftodeltaG(double Kf, double Kr)
{
    return -mRkj*mTemp*log(Kf/Kr);
}

double MassActionReaction::CalculateGibbsFromQuotient(double G0, double Q)
{
    // Q reaction quotient
    return G0 - 8.3144598e-3*mTemp*log(Q);
}

void MassActionReaction::SetReactionTemperature(double temp)
{
    mTemp = temp;
}

double MassActionReaction::GetReactionTemperature()
{
    return mTemp;
}

double MassActionReaction::GetForwardReactionRateConstant()
{
    return mForwardReactionRateConstant;
}

void MassActionReaction::SetForwardReactionRateConstant(double forwardReactionRateConstant)
{
    mForwardReactionRateConstant = forwardReactionRateConstant;
}

double MassActionReaction::GetReverseReactionRateConstant()
{
    return mReverseReactionRateConstant;
}

void MassActionReaction::SetReverseReactionRateConstant(double reverseReactionRateConstant)
{
    mReverseReactionRateConstant = reverseReactionRateConstant;
}

void MassActionReaction::SetGibbsDelmiter(std::string gibbsDelimiter)
{
    mGibbsDelimiter = gibbsDelimiter;
}

std::string MassActionReaction::GetGibbsDelimiter()
{
    return mGibbsDelimiter;
}