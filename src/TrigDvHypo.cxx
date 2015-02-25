// ************************************************
//
// NAME:     TrigDvHypo.cxx
// PACKAGE:  Trigger/TrigHypothesis/TrigDvHypo
//
// AUTHOR:   Jordi Duarte-Campderros (based in TrigBjetHypo)
// EMAIL:    jorge.duarte.campderros@cern.ch
// 
// ************************************************

#include "TrigBjetHypo/TrigDvHypo.h"

#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "TrigSteeringEvent/TrigPassBits.h"

#include "InDetBeamSpotService/IBeamCondSvc.h"

#include "TrigParticle/TrigEFBjetContainer.h"

//** ----------------------------------------------------------------------------------------------------------------- **//

TrigDvHypo::TrigDvHypo(const std::string& name, ISvcLocator* pSvcLocator) :
    HLT::HypoAlgo(name, pSvcLocator),
    m_totTimer(0),
    m_trigEFBjetColl(0),
    m_jetKey(""),
    m_ntrackDV(-1),
    m_massDV(-1),
    m_rDV(-1),
    m_acceptAll(false),
    m_instance(""),
    m_useBeamSpotFlag(false),
    m_cutCounter(0)
{
    declareProperty("JetKey",     m_jetKey   = ""     );
    
    declareProperty("NTrackDV",   m_ntrackDV  = 2      );
    declareProperty("MassDV",     m_massDV   = 6.*CLHEP::GeV );
    declareProperty("DistanceDV", m_rDV      = 4.*CLHEP::mm );

    declareProperty ("AcceptAll", m_acceptAll         );
    declareProperty ("Instance",  m_instance          );
    
    declareProperty ("UseBeamSpotFlag", m_useBeamSpotFlag = false);
    
    declareMonitoredVariable("CutCounter",    m_cutCounter,      AutoClear);
    declareMonitoredVariable("nTracksperRoI", m_tracksPerRoI, AutoClear);
}


//** ----------------------------------------------------------------------------------------------------------------- **//


TrigDvHypo::~TrigDvHypo() {}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvHypo::hltInitialize() 
{
    //* Get message service *//
    ATH_MSG_INFO("Initializing TrigDvHypo");

    ATH_MSG_INFO("declareProperty review:\n" 
	    << " NTrackDV        = " << m_ntrackDV   << "\n"
	    << " MassDV          = " << m_massDV    << "\n"
	    << " DistanceDV      = " << m_rDV       << "\n"
	    << " AcceptAll       = " << m_acceptAll << "\n"
	    << " Instance        = " << m_instance  << "\n"
	    << " UseBeamSpotFlag = " << m_useBeamSpotFlag);
  
    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvHypo::hltExecute(const HLT::TriggerElement* outputTE, bool& pass) 
{
    ATH_MSG_DEBUG("Executing TrigDvHypo");
    //* AcceptAll declare property setting *//
    if (m_acceptAll)
    {
        ATH_MSG_DEBUG("REGTEST: AcceptAll property is set: taking all "
                << "events and applying the selection only for saving the TrigPassBits");
    }
    else 
    {
        ATH_MSG_DEBUG("REGTEST: AcceptAll property not set: applying the "
                << "selection and saving the TrigPassBits");
    }
    //* initialise monitoring variables *//
    m_cutCounter = -1;

    //* Get RoI descriptor *//
    const TrigRoiDescriptor* roiDescriptor = 0;
    HLT::ErrorCode stat = getFeature(outputTE, roiDescriptor, m_jetKey);
    if (stat == HLT::OK) 
    {
        ATH_MSG_DEBUG("Using outputTE: " << "RoI id " 
                << roiDescriptor->roiId()<< ", Phi = " 
                <<  roiDescriptor->phi() << ", Eta = " 
                << roiDescriptor->eta());
    }
    else
    {
        ATH_MSG_WARNING("No RoI for this Trigger Element ");    
        return HLT::NAV_ERROR;
    }
  
    //* Define TrigPassBits for b-jets *//
    TrigPassBits *bitsEF=0;
    //* Retrieve xAOD b-jet object *//
    const TrigEFBjetContainer* trigEFBjetContainer=0;
    if(getFeature(outputTE, trigEFBjetContainer, "EFBjetDvFex") != HLT::OK) 
    {
        ATH_MSG_WARNING("Failed to get TrigEFBjetContainer");
        pass = false;
        return HLT::OK;
    }

    if(!trigEFBjetContainer)
    {
        ATH_MSG_DEBUG("Empty TrigEFBjetContainer");
        pass = false;
        return HLT::OK;
    }

    ATH_MSG_DEBUG("Got EFBjetContainer with " << trigEFBjetContainer->size() << " EFBjet (DV modified) object");
 
    if(trigEFBjetContainer->size() > 1) 
    {
        ATH_MSG_ERROR("More than one EFBjet-SV object to analyse: this should never happen");
    	return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::NAV_ERROR);
    }
    else if(trigEFBjetContainer->size() == 0) 
    {
        ATH_MSG_ERROR("No EFBjet-SV object to analyse: this should never happen");
        return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::NAV_ERROR);
    }
    
    //* Retrieve beamspot information *//
    if(m_useBeamSpotFlag) 
    {
      	IBeamCondSvc* m_iBeamCondSvc; 
    	StatusCode sc = service("BeamCondSvc", m_iBeamCondSvc);
      	if(sc.isFailure() || m_iBeamCondSvc == 0) 
        {
            ATH_MSG_WARNING("Could not retrieve Beam Conditions Service. ");
        }
        else 
        {
           int beamSpotStatus = 0;        
           int beamSpotBitMap = m_iBeamCondSvc->beamStatus();    
	    
           // To be promoted to a function BjetHelper..
           beamSpotStatus = ((beamSpotBitMap & 0x4) == 0x4);  
           if(beamSpotStatus) 	    
           {
               beamSpotStatus = ((beamSpotBitMap & 0x3) == 0x3);
           }
           else
           {
               m_cutCounter=0;
               pass = false;      
               return HLT::OK;
           }
        }
    }
  
   
    //* Add TrigPassBits for EF SV *//
    bitsEF = HLT::makeTrigPassBits(trigEFBjetContainer);

    //* to separate bad input TE and true behaviour *//
    //m_cutCounter=1;
  
    bool result = true;
    //* Loop over EFBjets and perform cut */
    for(auto &  trigEFBjet : *trigEFBjetContainer)
    {
        m_tracksPerRoI= trigEFBjet->xNVtx();
        
        result = true;
        // Apply the sequential cuts:
	    // number of tracks
        if(trigEFBjet->xNVtx() < m_ntrackDV)
        {
            m_cutCounter = 1;
            result *= false;
        }
        
        // Invariant mass of the SV
        if(trigEFBjet->xMVtx() < m_massDV)
        {
            m_cutCounter = 2;
            result *= false;
        }
        
        // distance to the PV
	    if(trigEFBjet->xSV() < m_rDV)
        {
            m_cutCounter = 3; 
            result *= false;
        }
    
        if(result)
        {
            HLT::markPassing(bitsEF, trigEFBjet,trigEFBjetContainer);
            pass = true;
        }
    }

    if(pass)
    {
        m_cutCounter = 4;
    }
    //* Print trigger decision *//
    if(m_acceptAll) 
    {
        ATH_MSG_DEBUG("REGTEST: Trigger decision is 1");
    } 
    else 
    {
        ATH_MSG_DEBUG("REGTEST: Trigger decision is " << pass);
    }

    //* Print TrigPassBits to outputTE *//
    // Nota que estic cambiant EFBjets -->  DvEFBjets
    if(attachBits(outputTE, bitsEF, "DvEFBjets") != HLT::OK) 
    {
    	ATH_MSG_ERROR("Problem attaching TrigPassBits for DV");
    }
    
    if(m_acceptAll) 
    {
    	m_cutCounter = 4;
    	pass = true;
    }
    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvHypo::hltFinalize() 
{
    ATH_MSG_INFO("Finalizing TrigDvHypo");

    return HLT::OK;
}

