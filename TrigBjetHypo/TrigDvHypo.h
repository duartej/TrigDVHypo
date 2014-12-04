// ************************************************
//
// NAME:     TrigDvHypo.h
// PACKAGE:  Trigger/TrigHypothesis/TrigDvHypo
//
// AUTHOR:   Andrea Coccaro
// EMAIL:    Andrea.Coccaro@ge.infn.it
// 
// ************************************************

#ifndef TRIGDVHYPO_H
#define TRIGDVHYPO_H

#include "TrigInterfaces/HypoAlgo.h"


class TrigEFBjetContainer;
class TrigTimerSvc;
class TriggerElement;

/**
 * @brief Hypo class for HLT Displaced vertex multi-track selection.
 *
 * @author Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
 *         (based in the TrigBjetHypo class)
 *
 * This is the base hypothesis class for the displaced vertex selection 
 * Selection can be performed using as handles from the secondary vertices
 * found in the TrigDvFex class.
 *
 * The event is selected using:
 *   - 
 *   -
 *   -
 *   -
 *  
 */

#include<string>

class TrigEFBjetContainer;  // Necesito esto?

class TrigDvHypo : public HLT::HypoAlgo 
{
    public:
      	/** @brief Constructor. */
      	TrigDvHypo (const std::string&, ISvcLocator*);
      	/** @brief Destructor. */
      	~TrigDvHypo ();
	
      	HLT::ErrorCode hltInitialize();
      	HLT::ErrorCode hltExecute(const HLT::TriggerElement*, bool&);
      	HLT::ErrorCode hltFinalize(); 

    private:
       	/** @brief Total execution time of TrigDvFex class. */
      	TrigTimer *m_totTimer;
	
      	/** @brief Pointer to TrigEFBjet collection. */
      	TrigEFBjetContainer* m_trigEFBjetColl;
      
	/** @brief The jet collection to be used ("EFJet" or "" for default config), "SplitJet" for new config */
	std::string m_jetKey;
	
       	/** @brief DeclareProperty: Displaced-vertex track multiplicity cut */ 
      	int m_ntrackDV;
       	
     	/** @brief DeclareProperty: Displaced-vertex mass cut (using the pion hypothesis.. */ 
      	float m_massDV;

	/** @brief DeclareProperty: Displaced-vertex minimum distance to any primary vertex*/
	float m_rDV;

       	/** @brief DeclareProperty: if acceptAll flag is set to true, every event is taken. */ 
      	bool m_acceptAll;

      	/** @brief DeclareProperty: string corresponding to the trigger level in which the algorithm is running. */
      	std::string m_instance;
	
      	/** @brief to check the beam spot flag status. */
      	bool m_useBeamSpotFlag;
	
      	/** @brief Cut counter. */
      	//unsigned short int m_cutCounter; --> try to use bitmask to extract the cuts used
	int m_cutCounter;
};

#endif

