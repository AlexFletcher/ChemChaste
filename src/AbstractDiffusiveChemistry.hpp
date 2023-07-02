#ifndef ABSTRACTDIFFUSIVECHEMISTRY_HPP_
#define ABSTRACTDIFFUSIVECHEMISTRY_HPP_

#include "AbstractChemistry.hpp"
#include "AbstractDiffusiveChemical.hpp"

#include <string>
#include <tuple>
#include <vector>

/**
 * Abstract property to contain information about the collection of chemical 
 * species in the system.
 */
class AbstractDiffusiveChemistry : public AbstractChemistry
{
private:

    /** \todo Document member variable. */
    unsigned mNumberDiffusiveChemicals;
    
    using AbstractChemistry::CheckChemical;
    using AbstractChemistry::AddChemical;
    using AbstractChemistry::UpdateChemicalVectors;
    using AbstractChemistry::AddChemistry;

protected:

    /** \todo Document member variable. */
    std::vector<std::string> mDiffusionChemicalNames;

    /** \todo Document member variable. Maybe this isn't necessary and the data structure of the AbstractDiffusiveChemicals is all that is necessary? */
    std::vector<std::vector<double> > mDiffusivityMatrix;
    
    /** \todo Document member variable. Maybe this isn't necessary and the data structure of the AbstractDiffusiveChemicals is all that is necessary? */
    std::vector<std::vector<std::string> > mDiffusionDomainMatrix;

public:

    /**
     * Constructor.
     */
    AbstractDiffusiveChemistry();

    /**
     * Destructor.
     */
    virtual ~AbstractDiffusiveChemistry();

    /**
     * \todo Document method.
     * 
     * @param pChemical
     */
    virtual void AddChemical(AbstractDiffusiveChemical* pChemical);

    /**
     * \todo Document method.
     * 
     * @param pChemical
     * 
     * @return whether chemical is already in class chemical vector.
     */
    virtual bool CheckChemical(AbstractDiffusiveChemical* pChemical); 

    /**
     * \todo Document method.
     * 
     * @param pChemical
     */
    virtual void UpdateChemicalVectors(AbstractDiffusiveChemical* pChemical);

    /**
     * \todo Document method.
     * 
     * @param pChemical
     * 
     * @return whether pChemical is already in class chemical vector.
     */
    virtual void UpdateDomainVector(AbstractDiffusiveChemical* pChemical);

    /**
     * Function to combine the chemistries in the effort to form an overall 
     * union of chemistries, for use as domain chemistry determines the solver 
     * state varaible ordering. Chemistry classes need not be of the same type, 
     * hence derived classes will need to dertermine the "highest class" define 
     * the highest class to be the parent chemistry.
     */
    virtual void AddChemistry(AbstractDiffusiveChemistry* pChemistry);

    /**
     * \todo Document method.
     * 
     * @return 
     */
    virtual std::string GetChemistryType();

    /**
     * Set mDiffusivityMatrix.
     * 
     * @param diffusivityMatrix the new value for mDiffusivityMatrix
     */
    void SetDiffusivityMatrix(std::vector<std::vector<double> > diffusivityMatrix);

    /**
     * @return mDiffusivityMatrix.
    */
    std::vector<std::vector<double> > GetDiffusivityMatrix();

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    std::vector<double> GetDiffusivityVectorByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @param chemicalName
     * @param domainName
     * 
     * @return 
     */
    double GetDiffusivityValueByChemicalAndDomainName(std::string chemicalName, std::string domainName);

    /**
     * \todo Document method.
     * 
     * @param chemicalName
     * 
     * @return 
     */
    double GetDiffusivityValueByChemicalName(std::string chemicalName);

    /**
     * Set mDiffusionDomainMatrix.
     * 
     * @param diffusionDomains the new value for mDiffusionDomainMatrix
     */
    void SetDiffusionDomainsMatrix(std::vector<std::vector<std::string> > diffusionDomains);

    /**
     * @return mDiffusionDomainMatrix
     */
    std::vector<std::vector<std::string>> GetDiffusionDomainsMatrix();

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    std::vector<std::string> GetDiffusionDomainsVectorByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @return 
     */
    unsigned GetNumberDiffusiveChemicals();

    /**
     * \todo Document method.
     * 
     * @param index
     * 
     * @return 
     */
    std::string GetDiffusiveChemicalNamesByIndex(unsigned index);

    /**
     * \todo Document method.
     * 
     * @param chemicalName
     * 
     * @return 
     */
    unsigned GetDiffusiveChemicalIndexByName(std::string chemicalName);
};

AbstractDiffusiveChemistry::AbstractDiffusiveChemistry()
{
    mDiffusivityMatrix = std::vector<std::vector<double> >();
    mDiffusionDomainMatrix = std::vector<std::vector<std::string> >();
    mNumberDiffusiveChemicals = 0;
}

AbstractDiffusiveChemistry::~AbstractDiffusiveChemistry()
{
}

void  AbstractDiffusiveChemistry::AddChemical(AbstractDiffusiveChemical* pChemical)
{
    /*
     * Add a chemical to the diffusive chemistry, check whether the diffusive 
     * domain of the chemical is new, and if it is new then add it to the 
     * chemical diffusion properties.
     */
    bool is_new_chem = false;
    if (mpChemicalVector.empty())
    {
        is_new_chem = true;
    }
    else
    {
        is_new_chem = CheckChemical(pChemical);
    }
    if (is_new_chem)
    {
        // If candidate is a new chemical then naturally the domain would be new for the chemical
        UpdateChemicalVectors(pChemical);
    }
    else
    {
        // Check whether the existing chemical is in a new domain
        UpdateDomainVector(pChemical);
    }
}

///\todo same structure as the AbstractChemical, remove?
bool AbstractDiffusiveChemistry::CheckChemical(AbstractDiffusiveChemical* pChemical)
{
    bool is_new_chem  = true;
    for (std::vector<AbstractChemical*>::iterator chem_iter = mpChemicalVector.begin();
         chem_iter != mpChemicalVector.end();
         ++chem_iter)
    {
        if (pChemical->GetChemicalName() == dynamic_cast<AbstractDiffusiveChemical*>(*chem_iter)->GetChemicalName())
        {   
            is_new_chem = false;
            break;
        }
    }
    return is_new_chem;
}

void AbstractDiffusiveChemistry::UpdateDomainVector(AbstractDiffusiveChemical* pChemical)
{
    for (std::vector<AbstractChemical*>::iterator chem_iter = mpChemicalVector.begin();
         chem_iter != mpChemicalVector.end();
         ++chem_iter)
    {
        // Assume the new chemical is known in the system
        if (pChemical->GetChemicalName() == dynamic_cast<AbstractDiffusiveChemical*>(*chem_iter)->GetChemicalName())
        {   
            unsigned num_candidate_domains = chemical->GetNumDiffusiveDomains();

            // Chemical whose domain is to be tested has been found in the chemical data structure
            for (unsigned candidate_domain_index = 0; candidate_domain_index < num_candidate_domains; ++candidate_domain_index)
            {
                // Test each of the potential new domains, assume to be new
                bool is_new_domain = true;
                unsigned index_record = 0;
                
                for (unsigned domain_index = 0; domain_index < dynamic_cast<AbstractDiffusiveChemical*>(*chem_iter)->GetChemicalName()->GetNumDiffusiveDomains(); domain_index++)
                {
                    // For each of the current existing domains
                    std::string existing_query_domain = dynamic_cast<AbstractDiffusiveChemical*>(*chem_iter)->GetChemicalName()->GetDiffusiveDomainByIndex(domain_index);
                
                    if (pChemical->GetDiffusiveDomainByIndex(candidate_domain_index) == existing_query_domain)
                    {
                        // New domain has already be recorded
                        is_new_domain = false;
                        index_record = candidate_domain_index;
                        break;
                    }
                }
                if (is_new_domain)
                {
                    // Then domain is new and add to DiffusiveChemical 
                    dynamic_cast<AbstractDiffusiveChemical*>(*chem_iter)->AddDiffusiveDomain(chemical->GetDiffusiveDomainByIndex(index_record),chemical->GetChemicalDiffusivityByIndex(index_record));
                }
            }
            
            break;
        }
    }
}

void AbstractDiffusiveChemistry::UpdateChemicalVectors(AbstractDiffusiveChemical* chemical)
{
    // virtual function to overide with further updates
    // additional diffusion properties
    AbstractChemistry::UpdateChemicalVectors(chemical); // implicit upcasting to AbstractChemcial

    UpdateDomainVector(chemical);

    mDiffusivityMatrix.push_back(chemical->GetChemicalDiffusivityVector());
    mDiffusionDomainMatrix.push_back(chemical->GetDiffusiveDomainVector());
    mDiffusionChemicalNames.push_back(chemical->GetChemicalName());
    mNumberDiffusiveChemicals +=1;
}

void AbstractDiffusiveChemistry::AddChemistry(AbstractDiffusiveChemistry* pChemistry)
{
    // Want to add the chemical vectors while preventing duplicates
    std::vector<AbstractChemical*> chemical_vector = pChemistry->rGetChemicalVector();
    for (std::vector<AbstractChemical*>::iterator chem_iter = chemical_vector.begin(); 
         chem_iter != chemical_vector.end();
         ++chem_iter)
    {
        // Add each of the chemicals from the additional chemistry in turn; this checks for duplicates implicitly
        AddChemical(dynamic_cast<AbstractDiffusiveChemical*>(*chem_iter));
    }
}

std::string AbstractDiffusiveChemistry::GetChemistryType()
{
    // Vvirtual function to be overriden in derived classes, used to identify properties to compy when adding chemistries
    return "AbstractDiffusiveChemistry";
}

void AbstractDiffusiveChemistry::SetDiffusivityMatrix(std::vector<std::vector<double> > diffusivityMatrix)
{
    mDiffusivityMatrix=diffusivityMatrix;
}

std::vector<std::vector<double>> AbstractDiffusiveChemistry::GetDiffusivityMatrix()
{
    return mDiffusivityMatrix;
}

std::vector<double> AbstractDiffusiveChemistry::GetDiffusivityVectorByIndex(unsigned index)
{
    if (index < mNumberDiffusiveChemicals)
    {
        return mDiffusivityMatrix[index];
    }
    else
    {
        std::cout << "Error: AbstractDiffusiveChemistry::GetDiffusivityByIndex(unsigned index), index out of bounds" << std::endl;
        std::vector<double> returnVec(1, 0.0);
        return returnVec;
    } 
}

double AbstractDiffusiveChemistry::GetDiffusivityValueByChemicalAndDomainName(std::string chemicalName, std::string domainName)
{
    // Return the diffusivity value for a given species with a given domain name
    unsigned index = GetDiffusiveChemicalIndexByName(chemicalName);

    std::vector<std::string> domainsVector = GetDiffusionDomainsVectorByIndex(index);
    std::vector<double> diffusivityVector = GetDiffusivityVectorByIndex(index);

    bool is_domain_found = false;
    for (unsigned domain_index = 0; domain_index < domainsVector.size(); domain_index++)
    {
        if (domainName == domainsVector[domain_index])
        {
            is_domain_found = true;
            return diffusivityVector[domain_index];
        }
    }
    if (!is_domain_found)
    {
        // The domain is unaccounted for in the AbstractDiffusiveChemical data structure
        return 0.0; // i.e not diffusive
    }
    return 0.0;
}

double AbstractDiffusiveChemistry::GetDiffusivityValueByChemicalName(std::string chemicalName)
{
    unsigned index = GetDiffusiveChemicalIndexByName(chemicalName);
    return GetDiffusivityVectorByIndex(index)[0];
}

void AbstractDiffusiveChemistry::SetDiffusionDomainsMatrix(std::vector<std::vector<std::string> > diffusionDomains)
{
    mDiffusionDomainMatrix = diffusionDomains;
}

std::vector<std::vector<std::string>> AbstractDiffusiveChemistry::GetDiffusionDomainsMatrix()
{
    return mDiffusionDomainMatrix;
}

std::vector<std::string> AbstractDiffusiveChemistry::GetDiffusionDomainsVectorByIndex(unsigned index)
{
    if (index < mNumberDiffusiveChemicals)
    {
        return mDiffusionDomainMatrix[index];
    }
    else
    {
        std::cout << "Error: AbstractDiffusiveChemistry::GetDiffusionDomainsByIndex(unsigned index), index out of bounds" << std::endl;
        std::vector<std::string> returnVec(1, "Null");
        return returnVec;
    } 
}

unsigned AbstractDiffusiveChemistry::GetNumberDiffusiveChemicals()
{
    return mNumberDiffusiveChemicals;
}

std::string AbstractDiffusiveChemistry::GetDiffusiveChemicalNamesByIndex(unsigned index)
{
    if (index < mNumberDiffusiveChemicals)
    {
        return mDiffusionChemicalNames[index];
    }
    else
    {
        std::cout << "Error: AbstractDiffusiveChemistry::GetDiffusiveChemicalNamesByIndex(unsigned index), index out of bounds" << std::endl;
        return "Error";
    }
}

unsigned AbstractDiffusiveChemistry::GetDiffusiveChemicalIndexByName(std::string chemicalName)
{
    bool is_found = false;
    for (unsigned index = 0; index < mNumberDiffusiveChemicals; index++)
    {
        if (chemicalName == mDiffusionChemicalNames[index])
        {
            is_found = true;
            return index;
        }

    }
    if (!is_found)
    {
        return mNumberDiffusiveChemicals; // should act out of bounds for future methods, yielding non-diffusive chemical
    }
    return mNumberDiffusiveChemicals;
}

#endif /* ABSTRACTDIFFUSIVECHEMISTRY_HPP_ */