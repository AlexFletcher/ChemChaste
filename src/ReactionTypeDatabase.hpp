#ifndef REACTIONTYPEDATABASE_HPP_
#define REACTIONTYPEDATABASE_HPP_

#include <string>
#include <vector>

//======================================================//
//              CHEMICAL REACTION HEADERS               //   
//======================================================//

// chemical reactions (ReactionTablet)
#include "AbstractReaction.hpp"
#include "AbstractReversibleReaction.hpp"
#include "MassActionReaction.hpp"
#include "SpectatorDependentReaction.hpp"
#include "MichaelisMentenReaction.hpp"

//======================================================//
//             TRANSPORT REACTION HEADERS               //   
//======================================================//

// transport reactions (TransportTablet)
#include "AbstractTransportReaction.hpp"
#include "AbstractTransportOutReaction.hpp"
#include "AbstractReversibleTransportReaction.hpp"
#include "MassActionTransportReaction.hpp"

//======================================================//
//              MEMBRANE REACTION HEADERS               //   
//======================================================//

#include "AbstractMembraneReaction.hpp"
#include "AbstractReversibleMembraneReaction.hpp"
#include "MassActionCoupledMembraneReaction.hpp"

//======================================================//
//                  CHEMICAL REACTIONS                  //   
//======================================================//

void ReactionTablet(AbstractReaction*& p_reaction,  std::string reactionType = "", std::vector<AbstractChemical*> substrates = std::vector<AbstractChemical*>(), std::vector<AbstractChemical*> products = std::vector<AbstractChemical*>(), std::vector<unsigned> stoichSubstrates = std::vector<unsigned>(), std::vector<unsigned> stoichProducts = std::vector<unsigned>(), std::string reactionInformation = "", bool isReversible = false, AbstractChemistry* p_systemChemistry = new AbstractChemistry())
{
    if (reactionType == "ZerothOrderReaction")
    {
        delete p_reaction;
        p_reaction = new AbstractReaction(substrates, products, stoichSubstrates, stoichProducts);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "ZerothOrderReversibleReaction")
    {
        delete p_reaction;
        p_reaction = new AbstractReversibleReaction(substrates, products, stoichSubstrates, stoichProducts);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "MassActionReaction")
    {        
        delete p_reaction;
        p_reaction = new MassActionReaction(substrates, products, stoichSubstrates, stoichProducts);

        // need to have it in the dynamic cast form to access function of downcast form, object natively of upcast form
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "SpectatorDependentReaction")
    {
        delete p_reaction;
        p_reaction = new SpectatorDependentReaction(substrates, products, stoichSubstrates, stoichProducts, p_systemChemistry);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "MichaelisMentenReaction")
    {
        delete p_reaction;
        p_reaction = new MichaelisMentenReaction(substrates, products, stoichSubstrates, stoichProducts, p_systemChemistry);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    // add more else if (reactionType == "...")
    else
    {
        std::cout << "ReactionTablet Error: "<<reactionType<<" is unknown to the ReactionTypeDataBase." << std::endl;
    }
};

//======================================================//
//                  TRANSPORT REACTIONS                 //   
//======================================================//

void TransportTablet(AbstractTransportReaction*& p_reaction,  std::string reactionType = "", std::vector<AbstractChemical*> bulkReactionSpecies = std::vector<AbstractChemical*>(), std::vector<AbstractChemical*> cellReactionSpecies = std::vector<AbstractChemical*>(), std::vector<unsigned> stoichBulk = std::vector<unsigned>(), std::vector<unsigned> stoichCell = std::vector<unsigned>(), std::string reactionInformation = "", bool isReversible = false, AbstractChemistry* p_systemChemistry = new AbstractChemistry())
{
    if (reactionType == "ZerothOrderTransportIntoCell")
    {
        delete p_reaction;
        p_reaction = new AbstractTransportReaction(bulkReactionSpecies, cellReactionSpecies, stoichBulk, stoichCell);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "ZerothOrderTransportOutOfCell")
    {
        delete p_reaction;
        p_reaction = new AbstractTransportOutReaction(bulkReactionSpecies, cellReactionSpecies, stoichBulk, stoichCell);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "ZerothOrderReversibleTransport")
    {
        delete p_reaction;
        p_reaction = new AbstractReversibleTransportReaction(bulkReactionSpecies, cellReactionSpecies, stoichBulk, stoichCell);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "MassActionTransportReaction")
    {
        delete p_reaction;
        p_reaction = new MassActionTransportReaction(bulkReactionSpecies, cellReactionSpecies, stoichBulk, stoichCell);

        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else
    {
        std::cout << "TransportTablet Error: "<<reactionType<<" is unknown to the ReactionTypeDataBase." << std::endl;
    }
};

//======================================================//
//                  MEMBRANE REACTIONS                  //   
//======================================================//

void MembraneTablet(AbstractMembraneReaction*& p_reaction,  std::string reactionType = "", std::vector<AbstractChemical*> bulkSubstrates = std::vector<AbstractChemical*>(), std::vector<AbstractChemical*> bulkProducts = std::vector<AbstractChemical*>(), std::vector<AbstractChemical*> cellSubstrates = std::vector<AbstractChemical*>(), std::vector<AbstractChemical*> cellProducts = std::vector<AbstractChemical*>(), std::vector<unsigned> stoichBulkSubstrates = std::vector<unsigned>(), std::vector<unsigned> stoichBulkProducts = std::vector<unsigned>(), std::vector<unsigned> stoichCellSubstrates = std::vector<unsigned>(), std::vector<unsigned> stoichCellProducts = std::vector<unsigned>(), std::string reactionInformation = "", bool isReversible = false, AbstractChemistry* p_bulkChemistry = new AbstractChemistry(), AbstractChemistry* p_cellChemistry = new AbstractChemistry())
{
    if (reactionType == "ZerothOrderCoupledMembrane")
    {
        delete p_reaction;
        p_reaction = new AbstractMembraneReaction(bulkSubstrates, bulkProducts, cellSubstrates, cellProducts, stoichBulkSubstrates, stoichBulkProducts,stoichCellSubstrates,stoichCellProducts);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);
    }
    else if (reactionType == "ZerothOrderReversibleMembrane")
    {
        delete p_reaction;
        p_reaction = new AbstractReversibleMembraneReaction(bulkSubstrates, bulkProducts, cellSubstrates, cellProducts, stoichBulkSubstrates, stoichBulkProducts,stoichCellSubstrates,stoichCellProducts);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);

    }
    else if (reactionType == "MassActionCoupledMembraneReaction")
    {
        delete p_reaction;
        p_reaction = new MassActionCoupledMembraneReaction(bulkSubstrates, bulkProducts, cellSubstrates, cellProducts, stoichBulkSubstrates, stoichBulkProducts,stoichCellSubstrates,stoichCellProducts);
        p_reaction->ParseReactionInformation(reactionInformation, isReversible);

    }
    else
    {
        std::cout << "MembraneTablet Error: "<<reactionType<<" is unknown to the ReactionTypeDataBase." << std::endl;
    }
}

#endif