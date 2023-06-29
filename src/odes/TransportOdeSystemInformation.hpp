#ifndef TRANSPORTODESYSTEMINFORMATION_HPP_
#define TRANSPORTODESYSTEMINFORMATION_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractTransportReactionSystem.hpp"

/**
 * \todo Document class.
 */
template<class ODE_SYSTEM>
class TransportOdeSystemInformation : public AbstractOdeSystemInformation
{
private:

    /** \todo Document member variable. */
    static boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > mpInstance;

    /** \todo Document member variable. */
    AbstractTransportReactionSystem* mpReactionSystem;

protected:

    /** \todo not defined anywhere? */
    TransportOdeSystemInformation(const TransportOdeSystemInformation<ODE_SYSTEM>&);

    /** \todo not defined anywhere? */
    TransportOdeSystemInformation& operator= (const TransportOdeSystemInformation<ODE_SYSTEM>&);

    /** \todo Document method. */
    void Initialise();

    /** \todo Document method. */
    AbstractTransportReactionSystem* GetReactionSystem();

    /** \todo Document method. */
    void SetReactionSystem(AbstractTransportReactionSystem* pReactionSystem);

public:
  
    /** \todo Document method. */
    static boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > Instance();

    /**
     * Default constructor; calls Initialise.
     *
     * Designed to be used as follows by ODE system classes in their 
     * constructors:
     * mpSystemInfo.reset(new TransportOdeSystemInformation<CLASS>);
     */
    TransportOdeSystemInformation(AbstractTransportReactionSystem* pReactionSystem);

    ///\todo Add destructor?
};

template<class ODE_SYSTEM>
TransportOdeSystemInformation<ODE_SYSTEM>::TransportOdeSystemInformation(AbstractTransportReactionSystem* pReactionSystem)
    : mpReactionSystem(pReactionSystem)
{
    assert(mpInstance == nullptr);
    TransportOdeSystemInformation<ODE_SYSTEM>::Initialise();
}

template<class ODE_SYSTEM>
void TransportOdeSystemInformation<ODE_SYSTEM>::Initialise()
{
}

template<class ODE_SYSTEM>
boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > TransportOdeSystemInformation<ODE_SYSTEM>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new TransportOdeSystemInformation<ODE_SYSTEM>);
        mpInstance->Initialise();
    }
    return mpInstance;
}

template<class ODE_SYSTEM>
AbstractTransportReactionSystem* TransportOdeSystemInformation<ODE_SYSTEM>::GetReactionSystem()
{
    return mpReactionSystem;
}

template<class ODE_SYSTEM>
void TransportOdeSystemInformation<ODE_SYSTEM>::SetReactionSystem(AbstractTransportReactionSystem* pReactionSystem)
{
    mpReactionSystem  = pReactionSystem;
}

template<class ODE_SYSTEM>
boost::shared_ptr<TransportOdeSystemInformation<ODE_SYSTEM> > TransportOdeSystemInformation<ODE_SYSTEM>::mpInstance;

#endif /* TRANSPORTODESYSTEMINFORMATION_HPP_ */