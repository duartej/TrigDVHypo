// ************************************************
//
// NAME:     TrigDvFex.cxx
// PACKAGE:  Trigger/TrigHypothesis/TrigBjetHypo
// 
// ************************************************

#include "TrigBjetHypo/TrigDvFex.h"
#include "TaggerHelper.h"

// Eigen library (Amg::error)
#include "EventPrimitives/EventPrimitivesHelpers.h"

#include "Particle/TrackParticleContainer.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

// --> Changing TrigEFBjetContainer --> TrigInDetTrackCollection ??
#include "TrigParticle/TrigEFBjetContainer.h"

#include "TrigSteeringEvent/TrigOperationalInfo.h"


#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODBase/IParticle.h"

//** ---------------------------------------------------------------------- **//
const double TOLERANCE=1e-30;

TrigDvFex::TrigDvFex(const std::string& name, ISvcLocator* pSvcLocator) :
  HLT::FexAlgo(name, pSvcLocator),
  m_trigEFBjetColl(0),
  m_trigBjetJetInfo(0),
  m_totTracks(0),
  m_totSelTracks(0)
{
  declareProperty ("AlgoId",             m_algo);
  declareProperty ("Instance",           m_instance);
  declareProperty ("JetKey",             m_jetKey     = ""); //"" needed for default config, SplitJet for new config
  declareProperty ("PriVtxKey",          m_priVtxKey  = "EFHistoPrmVtx");

  declareProperty ("onlinemon",          m_mon_online = true);
  declareProperty ("validation",         m_mon_validation = true);
  
  declareProperty ("HistoPrmVtxAtEF",    m_histoPrmVtxAtEF    = true);
  declareProperty ("UseEtaPhiTrackSel",  m_useEtaPhiTrackSel  = false);
  
  declareProperty ("UseJetDirection",    m_useJetDirection);

  declareProperty ("TrkSel_Chi2",        m_trkSelChi2         = 0.0);
  declareProperty ("TrkSel_BLayer",      m_trkSelBLayer       = 1);
  declareProperty ("TrkSel_PixHits",     m_trkSelPixHits      = 2);
  declareProperty ("TrkSel_SiHits",      m_trkSelSiHits       = 4);
  declareProperty ("TrkSel_D0",          m_trkSelD0           = 300.0*CLHEP::mm);
  declareProperty ("TrkSel_Z0",          m_trkSelZ0           = 300.0*CLHEP::mm);
  declareProperty ("TrkSel_Pt",          m_trkSelPt           = 4.0*CLHEP::GeV);

  declareMonitoredStdContainer("trk_a0",            m_mon_trk_a0,        AutoClear);
  declareMonitoredStdContainer("trk_a0_sel",        m_mon_trk_a0_sel,    AutoClear);
  declareMonitoredStdContainer("trk_S(a0)_sel",     m_mon_trk_Sa0_sel,   AutoClear);
  declareMonitoredStdContainer("trk_z0",            m_mon_trk_z0,        AutoClear);
  declareMonitoredStdContainer("trk_z0_sel",        m_mon_trk_z0_sel,    AutoClear);
  declareMonitoredStdContainer("trk_z0_sel_PV",     m_mon_trk_z0_sel_PV, AutoClear);
  declareMonitoredStdContainer("trk_S(z0)_sel",     m_mon_trk_Sz0_sel,   AutoClear);
  declareMonitoredStdContainer("trk_prob",          m_mon_trk_prob,      AutoClear);

  declareMonitoredVariable    ("roi_nTracks",       m_totTracks);
  declareMonitoredVariable    ("roi_nTracks_sel",   m_totSelTracks);
  declareMonitoredStdContainer("roi_stepsToSelect", m_listCutApplied, AutoClear);
  declareMonitoredObject      ("roi_selectedTracks", *this, &TrigDvFex::totSelectedTracks);

  declareMonitoredVariable    ("roi_deltaEtaJet",       m_deltaEtaJet,       AutoClear);
  declareMonitoredVariable    ("roi_deltaPhiJet",       m_deltaPhiJet,       AutoClear);
  declareMonitoredVariable    ("roi_deltaEtaTrkJet",    m_deltaEtaTrkJet,    AutoClear);
  declareMonitoredVariable    ("roi_deltaPhiTrkJet",    m_deltaPhiTrkJet,    AutoClear);
  declareMonitoredVariable    ("roi_deltaEtaJetTrkJet", m_deltaEtaJetTrkJet, AutoClear);
  declareMonitoredVariable    ("roi_deltaPhiJetTrkJet", m_deltaPhiJetTrkJet, AutoClear);

  m_taggerHelper = new TaggerHelper(msg(), msgLvl());
}


//** ----------------------------------------------------------------------------------------------------------------- **//


TrigDvFex::~TrigDvFex() 
{
  if(m_taggerHelper)            delete m_taggerHelper;
  if(m_trigBjetJetInfo)         delete m_trigBjetJetInfo;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::hltInitialize() 
{
    // Get message service
    ATH_MSG_INFO("Initializing TrigDvFex, version " << PACKAGE_VERSION);
    // declareProperty overview
    ATH_MSG_DEBUG("declareProperty review:");

    ATH_MSG_DEBUG(" JetKey = "               << m_jetKey );
    ATH_MSG_DEBUG(" PriVtxKey = "            << m_priVtxKey );

    ATH_MSG_DEBUG(" AlgoId = "              << m_algo ); 
    ATH_MSG_DEBUG(" Instance = "            << m_instance );
 
    ATH_MSG_DEBUG(" HistoPrmVtxAtEF = "     << m_histoPrmVtxAtEF );
    ATH_MSG_DEBUG(" UseEtaPhiTrackSel = "   << m_useEtaPhiTrackSel );
    ATH_MSG_DEBUG(" UseJetDirection = "     << m_useJetDirection );

    ATH_MSG_DEBUG(" TrkSel_Chi2 = "     << m_trkSelChi2 ); 
    ATH_MSG_DEBUG(" TrkSel_BLayer = "   << m_trkSelBLayer ); 
    ATH_MSG_DEBUG(" TrkSel_SiHits = "   << m_trkSelSiHits ); 
    ATH_MSG_DEBUG(" TrkSel_D0 = "       << m_trkSelD0 ); 
    ATH_MSG_DEBUG(" TrkSel_Z0 = "       << m_trkSelZ0 ); 
    ATH_MSG_DEBUG(" TrkSel_Pt = "       << m_trkSelPt ); 

    m_trigBjetJetInfo    = new TrigBjetJetInfo();

    return HLT::OK;
}



//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::getTrackCollection(const xAOD::TrackParticleContainer*& pointerToEFTrackCollections,
                 const HLT::TriggerElement* whateverTE)
{
    pointerToEFTrackCollections = 0;

    std::vector<const xAOD::TrackParticleContainer*> vectorOfEFTrackCollections;
    HLT::ErrorCode status = getFeatures(whateverTE, vectorOfEFTrackCollections, "");//"bJetTracks"); 
    if(status != HLT::OK) 
    {
       ATH_MSG_ERROR("Failed to get TrackParticleContainer from the trigger element");
       return status;
    } 
    ATH_MSG_DEBUG("Got " << vectorOfEFTrackCollections.size() << " TrackParticleContainer");
    
    if( vectorOfEFTrackCollections.size() > 1 )
    {
        ATH_MSG_ERROR("The vector of  xAOD::TrackParticleContainer have more than 1 element!");
        return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::NAV_ERROR);
    }
    else if( vectorOfEFTrackCollections.size() < 1 )
    {
        ATH_MSG_DEBUG("The vector of xAOD:TrackParticleContainer have none element!");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }

    pointerToEFTrackCollections = vectorOfEFTrackCollections[0];
    
    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//
bool TrigDvFex::efTrackSel(const xAOD::TrackParticle*& track, unsigned int i) 
{
    uint8_t nBlayerHits = 0;
    uint8_t nPixHits    = 0;  
    uint8_t nSCTHits    = 0; 
  
    track->summaryValue(nBlayerHits, xAOD::numberOfBLayerHits);
    track->summaryValue(nPixHits,    xAOD::numberOfPixelHits);
    track->summaryValue(nSCTHits,    xAOD::numberOfSCTHits);
    
    int   nSiHits = nPixHits + nSCTHits;
    float theta   = track->theta();
    float qOverPt = track->qOverP()/TMath::Sin(theta); 
    float pT      = (1.0/qOverPt);
    float d0      = track->d0();
    float z0      = track->z0();


    ATH_MSG_VERBOSE( "efTrackSel method\n" <<
	    "  Track number "    << i+1  << " to be selected must be:\n" <<
       	    "    Pt "            << fabs(pT)                      << " >= " << m_trkSelPt << "\n" <<
       	    "    d0 "            << fabs(d0)                      << " >= " << m_trkSelD0 << "\n" <<
	   // "    z0*sin(theta) " << fabs(z0-zv)*TMath::Sin(theta) << " <= " << m_trkSelZ0 << "\n" <<
       	    "    bLayer "        << (int)nBlayerHits              << " >= " << m_trkSelBLayer << "\n" <<
       	    "    pixelHit "      << (int)nPixHits                 << " >= " << m_trkSelPixHits<< "\n" <<
            "    SiHit "         << (int)nSiHits                  << " >= " << m_trkSelSiHits << "\n" <<
       	    "    Prob(chi2) "    << TMath::Prob(track->chiSquared(), (int)nSiHits*3-5) << " > " << m_trkSelChi2);
  
    if(m_useEtaPhiTrackSel) 
    {
    	if(fabs(track->eta() - m_trigBjetJetInfo->etaRoI()) > 0.2) 
        {
            ATH_MSG_DEBUG("  track " << i+1 << " is not selected (eta matching)");
            m_listCutApplied.push_back(CutListMonitor::RoIEtaMatching); 
            return false;
        }
        
        if(fabs(m_taggerHelper->phiCorr(m_taggerHelper->phiCorr(track->phi()) - m_trigBjetJetInfo->phiRoI())) > 0.2) 
        {
            ATH_MSG_DEBUG( "  track " << i+1 << " is not selected (phi matching)");
            m_listCutApplied.push_back(CutListMonitor::RoIPhiMatching); 
            return false;
        }
    }

    if(fabs(pT) < m_trkSelPt) 
    {
        ATH_MSG_DEBUG("  track " << i+1 << " not selected (pT cut)");
        m_listCutApplied.push_back(CutListMonitor::PtCut); 
        return false;
    }
  
    if(fabs(d0) < m_trkSelD0) 
    {
        ATH_MSG_DEBUG("  track " << i+1 << " not selected (d0 cut)");
        m_listCutApplied.push_back(CutListMonitor::D0Cut); 
        return false;
    }

    /*if(fabs(z0-zv)*TMath::Sin(theta) > m_trkSelZ0) 
    {
       ATH_MSG_DEBUG("  track " << i+1 << " not selected (z0 cut)");
       m_listCutApplied.push_back(CutListMonitor::Z0Cut); 
       return false;
    }*/

    if(nBlayerHits < m_trkSelBLayer)
    {
       ATH_MSG_DEBUG("  track " << i+1 << " not selected (missing b-layer hit)");
       m_listCutApplied.push_back(CutListMonitor::BLayerHitsCut); 
       return false;
    }

    if(nPixHits < m_trkSelPixHits) 
    { 
    	ATH_MSG_DEBUG("  track " << i+1 << " not selected (too few pixel hits)");
    	m_listCutApplied.push_back(CutListMonitor::PixelHitsCut); 
	    return false;
    }
   
    if(nSiHits < m_trkSelSiHits) 
    {
	    ATH_MSG_DEBUG("  track " << i+1 << " not selected (too few silicon hits)");
	    m_listCutApplied.push_back(CutListMonitor::SiHitsCut); 
	    return false;
    }
  
    if(TMath::Prob(track->chiSquared(), (int)nSiHits*3-5) <= m_trkSelChi2) 
    {
        ATH_MSG_DEBUG("  track " << i+1 << " not selected (chi2 cut)");
        m_listCutApplied.push_back(CutListMonitor::Chi2Cut); 
        return false;
    }

    ATH_MSG_DEBUG("    track " << i+1 << " is selected");
   
    m_listCutApplied.push_back(CutListMonitor::SELECTED);
    // Monitoring stuff -------------------------------------------------
    if(m_mon_validation || m_mon_online)
    {
        m_mon_trk_a0_sel.push_back(track->d0());
        m_mon_trk_z0_sel.push_back(track->z0());
        m_mon_trk_z0_sel_PV.push_back(track->z0());//-m_trigBjetPrmVtxInfo->zPrmVtx());
        
        // ez0
	    const float errIP1D = Amg::error(track->definingParametersCovMatrix(),1);
        // ed0
	    const float errIP2D = Amg::error(track->definingParametersCovMatrix(),0);
    	const float z0 = track->z0();//;-m_trigBjetPrmVtxInfo->zPrmVtx();

        float phiJetObject = 0.0;
        float etaJetObject = 0.0; 
        
        if(m_useJetDirection == 1) 
        {
            // Using HLT Jet 
            phiJetObject = m_trigBjetJetInfo->phiJet();
            etaJetObject = m_trigBjetJetInfo->etaJet();
        }
        else if(m_useJetDirection == 2) 
        {
            // Using HLT track-based jet
            phiJetObject = m_trigBjetJetInfo->phiTrkJet();
            etaJetObject = m_trigBjetJetInfo->etaTrkJet();
        } 
        else if(m_useJetDirection == 3) 
        {
            // Using the LvL1 jet RoI 
	        phiJetObject = m_trigBjetJetInfo->phiRoI();
            etaJetObject = m_trigBjetJetInfo->etaRoI();
        }
        const float d0Sign = m_taggerHelper->signedD0(track->d0(), track->phi(), phiJetObject);
        const float z0Sign = m_taggerHelper->signedZ0(z0, track->eta(), etaJetObject);
        //const float sigmaBeamSpot = (m_trigBjetPrmVtxInfo->xBeamSpotWidth()+
        //        m_trigBjetPrmVtxInfo->yBeamSpotWidth())/2.0;
        float IP1D = 0.0;
        float IP2D = 0.0;
        if(fabs(errIP1D) > TOLERANCE)
        {
            IP1D = z0Sign/sqrt(errIP1D*errIP1D);
        }
        if( (fabs(errIP2D) > TOLERANCE) )//|| sigmaBeamSpot)
        {
            IP2D = d0Sign/sqrt(IP2D*errIP2D);//+sigmaBeamSpot*sigmaBeamSpot);
        }
        m_mon_trk_Sz0_sel.push_back(IP1D);
        m_mon_trk_Sa0_sel.push_back(IP2D);
    }
    
    return true;
}

//** ----------------------------------------------------------------------------------------------------------------- **//
HLT::ErrorCode TrigDvFex::getJet(const xAOD::Jet * & jetreturn, const HLT::TriggerElement* inputTE)
{
    // Initialize
    jetreturn = 0;

    //const xAOD::JetContainer* jets(0);
    const xAOD::JetContainer* jets = 0;
    HLT::ErrorCode ec = getFeature(inputTE, jets, m_jetKey);
    if(ec!=HLT::OK) 
    {
        ATH_MSG_WARNING("Failed to get JetCollection");
        return HLT::OK;
    } 
    
    ATH_MSG_DEBUG("Obtained JetContainer");
  
    if(jets == 0)
    {
        ATH_MSG_WARNING("Jet collection pointer is 0");
        return HLT::OK;
        //return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::UNKNOWN);;
    }
 
    std::vector<const xAOD::Jet*> theJets(jets->begin(), jets->end());
    std::size_t njets = theJets.size();
    ATH_MSG_DEBUG("  * with size: " << njets);
  
    if( njets == 0 )
    {
        ATH_MSG_DEBUG("so, JetCollection is empty");
        return HLT::OK;
        //return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }
    else if(njets > 1)
    {
        ATH_MSG_DEBUG("Something is wrong, it should not be more than one jet");
	    return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::NAV_ERROR);
    }
    // The jet passed back by reference 
    jetreturn=theJets[0];
  
    if(msgLvl() <= MSG::DEBUG)
    {
        for(const xAOD::Jet* aJet : theJets) 
        { 
            static int i=0;   
            ATH_MSG_DEBUG("  ["<<i<<"-Jet]: Et=" << aJet->p4().Et()/CLHEP::GeV << " (GeV),"
                   << "  Eta=" << aJet->p4().Eta() << " Phi= " << aJet->p4().Phi() 
                << " Phi (corrected with helper)= " << m_taggerHelper->phiCorr(aJet->p4().Phi())) ;
        }
    }
    return HLT::OK;
}
 

// CAVEAT!! Be careful with the use of inputTE which maybye doesn't have all the features which is suppose
// to have. The ouputTE is going to contain, in some cases, references to objects related to just processed
// RoI (in the previous sequences, which are actually executing this sequence...). So, if you find problems
// just use outputTE, which is safer (although inputTE is slightly faster, marginaly in fact)
HLT::ErrorCode TrigDvFex::hltExecute(const HLT::TriggerElement* inputTE, HLT::TriggerElement* outputTE) 
{
    ATH_MSG_DEBUG("Executing TrigDvFex");

    // Clear and initialize data members
    m_totSelTracks = 0;
    m_totTracks    = 0;
    m_setPVInfo    = false;
  
    m_trigBjetJetInfo->clear();

    // -----------------------------------
    // Get RoI descriptor
    // -----------------------------------
    const TrigRoiDescriptor* roiDescriptor = 0;
    if(getFeature(inputTE, roiDescriptor, m_jetKey) != HLT::OK) 
    {
        ATH_MSG_DEBUG("No feature for this Trigger Element");    	
    	return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::NAV_ERROR);
    }
    ATH_MSG_DEBUG("Using TE: " << "RoI id " << roiDescriptor->roiId()
	    << ", Phi = " <<  roiDescriptor->phi() << ", Eta = " << roiDescriptor->eta());
    
    // -----------------------------------
    // Get Track collection 
    // -----------------------------------
    const xAOD::TrackParticleContainer *  pointerToEFTrackCollections = 0;
    HLT::ErrorCode tpstatus = getTrackCollection(pointerToEFTrackCollections,outputTE);
    if(tpstatus != HLT::OK)
    {
        ATH_MSG_DEBUG("No track collection retrieved... "); 
    } 
    ATH_MSG_DEBUG("Track collection retrieved");
    // Track info
    std::vector<TrigBjetTrackInfo> trigBjetTrackInfoVector;
    m_trigBjetTrackInfoVector = &trigBjetTrackInfoVector;
    // Get number of reconstructed tracks in this RoI
    if(pointerToEFTrackCollections)
    {
        m_totTracks = m_taggerHelper->getTrackNumber(pointerToEFTrackCollections);
        for(unsigned int j = 0; j < m_totTracks; ++j) 
        {
            const xAOD::TrackParticle* track = (*pointerToEFTrackCollections)[j];
            if(m_mon_validation)
            {
                m_mon_trk_a0.push_back(track->d0());
                m_mon_trk_z0.push_back(track->z0());
            }
        
            if(efTrackSel(track, j)) 
            {
                m_totSelTracks++;
                TrigBjetTrackInfo trigBjetTrackInfo(track);
            
                float d0Corr=0, z0Corr=0;
                d0Corr=track->d0(); 
                z0Corr=track->z0();
                trigBjetTrackInfo.setIPCorr(d0Corr, z0Corr);
                ATH_MSG_DEBUG(trigBjetTrackInfo);
                trigBjetTrackInfoVector.push_back(trigBjetTrackInfo);
      	    }
        }
    }
    else
    {
        m_totTracks = 0;
    }  

    
    // --- JetInfo
    // Jets: remember in some previous step, each jet has been associated to
    // a RoI, so the RoI and the jet information is the same
    m_trigBjetJetInfo->setEtaPhiRoI(roiDescriptor->eta(), m_taggerHelper->phiCorr(roiDescriptor->phi()));
    const xAOD::Jet * thejet = 0;
    HLT::ErrorCode xaodc = getJet(thejet,inputTE);
    if( xaodc != HLT::OK )
    {
        return xaodc;
    }	
    
    if(thejet)
    {
        m_trigBjetJetInfo->setEtaPhiJet(thejet->p4().Eta(),m_taggerHelper->phiCorr(thejet->p4().Phi()));
        m_trigBjetJetInfo->setEtJet(thejet->p4().Et());
    }

  
    // New container to
    m_trigEFBjetColl = new TrigEFBjetContainer();
    // -----------------------------------
    // Create TrigEFBjet and attach feature
    // -----------------------------------
    // Note that the meaning of some characteristics are tuned to 
    // DV case (decay length significance is actually the decay length)
    // Make sense to put more than one TrigEFBjet object (for instance,
    // one for each SV, if there are more than one...
    TrigEFBjet* trigEFBjet = new TrigEFBjet(roiDescriptor->roiId(), 
	    m_trigBjetJetInfo->etaJet(), m_trigBjetJetInfo->phiJet(),
	    0, 0, 0, 0.0,  m_trigBjetJetInfo->etJet(),
	    -1,-1,-1, -1,-1,  // Note, i can use this to fill other stuff if I need 
	    10.0*CLHEP::mm, 10.0*CLHEP::GeV,  //TRACK-
        0.0, m_totSelTracks);   //TRACK-

    trigEFBjet->validate(true);
    m_trigEFBjetColl->push_back(trigEFBjet);
    
    if(!m_trigEFBjetColl) 
    {
    	ATH_MSG_ERROR("Feature TrigEFBjetContainer not found");
       	return HLT::ErrorCode(HLT::Action::ABORT_JOB, HLT::Reason::BAD_JOB_SETUP);
    }
    
    HLT::ErrorCode stat = attachFeature(outputTE, m_trigEFBjetColl, "EFBjetDvFex");
    if(stat != HLT::OK) 
    {
	    ATH_MSG_DEBUG("Failed to attach TrigEFBjetContainer to navigation");
	    return stat;
    }
    
    // Summary info to print if it is in the proper message level    
    if(msgLvl() <= MSG::DEBUG) 
    {
      	const EventInfo* pEventInfo = 0;
     	if( !store() || store()->retrieve(pEventInfo).isFailure() ) 
        {
            ATH_MSG_DEBUG("Failed to get EventInfo ");
        } 
        else 
        {
            ATH_MSG_DEBUG("TrigDvFex::hltExecute END");
            ATH_MSG_DEBUG("DV summary (Run " << pEventInfo->event_ID()->run_number()    
                    << "; Event " << pEventInfo->event_ID()->event_number() << ")");
            ATH_MSG_DEBUG("REGTEST:  RoI " << roiDescriptor->roiId() << ", Phi = "   
                    << roiDescriptor->phi() << ", Eta = "   << roiDescriptor->eta());
            ATH_MSG_DEBUG("REGTEST:  Tracks: " << m_totTracks << " reconstructed and " 
                    << m_totSelTracks <<" selected");
            if(pointerToEFTrackCollections)
            {
                for(const xAOD::TrackParticle * track: *pointerToEFTrackCollections)
                {
                    ATH_MSG_DEBUG("REGTEST:     d0: " << track->d0()/CLHEP::mm << " [mm], eta:"
                           << track->eta() << ", phi:" << track->phi() << ", origin: (" 
                           << track->vx()/CLHEP::mm << "," << track->vy()/CLHEP::mm << ","
                           << track->vz()/CLHEP::mm << ") [mm]" );
                }
            }
        }
    }

    return HLT::OK;
}

//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::hltFinalize() 
{
    ATH_MSG_INFO("Finalizing TrigDvFex");
    return HLT::OK;
}



