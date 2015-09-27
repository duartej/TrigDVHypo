// ************************************************
//
// NAME:     TrigDvFex.cxx
// PACKAGE:  Trigger/TrigHypothesis/TrigBjetHypo
// 
// ************************************************

#include "TrigBjetHypo/TrigDvFex.h"

// Eigen library (Amg::error)
#include "EventPrimitives/EventPrimitivesHelpers.h"

#include "Particle/TrackParticleContainer.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

// --> Changing TrigEFBjetContainer --> TrigInDetTrackCollection ??
#include "TrigParticle/TrigEFBjetContainer.h"

//JDC ---< TESTING!!
#include "TrkRIO_OnTrack/RIO_OnTrack.h"
//#include "TrkPrepRawData/PrepRawDataContainer.h"
#include "IRegionSelector/IRegSelSvc.h"
#include "IRegionSelector/IRoiDescriptor.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"

// FIXME: TO BE DEPRECATED
#include "TrigSteeringEvent/TrigOperationalInfo.h"

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODBase/IParticle.h"

//** ---------------------------------------------------------------------- **//
const double TOLERANCE=1e-30;
const double HALF_PI = M_PI/2.0;


TrigDvFex::TrigDvFex(const std::string& name, ISvcLocator* pSvcLocator) :
  HLT::FexAlgo(name, pSvcLocator),
  m_regionSelector("RegSelSvc",this->name()),
  m_pixcontainer_v(0),
  m_sctcontainer_v(0),
  m_trtcontainer_v(0),
  m_is_pix_filled(false),
  m_is_sct_filled(false),
  m_is_trt_filled(false),
  m_are_prds_filled(false),
  m_trigEFBjetColl(0),
  m_trigBjetJetInfo(0),
  m_totTracks(0),
  m_totSelTracks(0),
  m_unusedhits(0),
  m_unusedhits_fraction(0),
  m_deltaEtaTrk(0),
  m_deltaPhiTrk(0)
{
  declareProperty ("AlgoId",             m_algo);

  // Services
  declareProperty ("JetKey",             m_jetKey     = ""); //"" needed for default config, SplitJet for new config

  declareProperty ("onlinemon",          m_mon_online = true);
  declareProperty ("validation",         m_mon_validation = true);
  
  declareProperty ("UseEtaPhiTrackSel",  m_useEtaPhiTrackSel  = false);
  
  declareProperty ("TrkSel_Chi2",        m_trkSelChi2         = 0.0);
  declareProperty ("TrkSel_BLayer",      m_trkSelBLayer       = 1);
  declareProperty ("TrkSel_PixHits",     m_trkSelPixHits      = 2);
  declareProperty ("TrkSel_SiHits",      m_trkSelSiHits       = 4);
  declareProperty ("TrkSel_D0",          m_trkSelD0           = 300.0*CLHEP::mm);
  declareProperty ("TrkSel_Z0",          m_trkSelZ0           = 300.0*CLHEP::mm);
  declareProperty ("TrkSel_Pt",          m_trkSelPt           = 4.0*CLHEP::GeV);

  declareMonitoredStdContainer("trk_a0",            m_mon_trk_a0,        AutoClear);
  declareMonitoredStdContainer("trk_a0_sel",        m_mon_trk_a0_sel,    AutoClear);
  declareMonitoredStdContainer("trk_Sa0_sel",       m_mon_trk_Sa0_sel,   AutoClear);
  declareMonitoredStdContainer("trk_z0",            m_mon_trk_z0,        AutoClear);
  declareMonitoredStdContainer("trk_z0_sel",        m_mon_trk_z0_sel,    AutoClear);
  declareMonitoredStdContainer("trk_Sz0_sel",       m_mon_trk_Sz0_sel,   AutoClear);

  declareMonitoredStdContainer("trk_theta",         m_mon_trk_theta,     AutoClear);
  declareMonitoredStdContainer("trk_pt",            m_mon_trk_pt,       AutoClear);
  declareMonitoredStdContainer("trk_qOverPt",       m_mon_trk_qOverPt,   AutoClear);
  declareMonitoredStdContainer("trk_eta",           m_mon_trk_eta,       AutoClear);

  declareMonitoredStdContainer("trk_BlayerHits",    m_mon_trk_BlayerHits,AutoClear);
  declareMonitoredStdContainer("trk_PixHits",       m_mon_trk_PixHits,   AutoClear);
  declareMonitoredStdContainer("trk_SCTHits",       m_mon_trk_SCTHits,   AutoClear);
  declareMonitoredStdContainer("trk_TRTHits",       m_mon_trk_TRTHits,   AutoClear);
  
  declareMonitoredStdContainer("trk_PixHoles",      m_mon_trk_PixHoles,  AutoClear);
  declareMonitoredStdContainer("trk_SCTHoles",      m_mon_trk_SCTHoles,  AutoClear);
  declareMonitoredStdContainer("trk_TRTHoles",      m_mon_trk_TRTHoles,  AutoClear);
  
  declareMonitoredStdContainer("trk_PixSharedHits", m_mon_trk_PixSharedHits,   AutoClear);
  declareMonitoredStdContainer("trk_SCTSharedHits", m_mon_trk_SCTSharedHits,   AutoClear);
  //declareMonitoredStdContainer("trk_TRTSharedHits", m_mon_trk_TRTSharedHits,   AutoClear);
  
  declareMonitoredStdContainer("trk_theta_sel",      m_mon_trkSel_theta,     AutoClear);
  declareMonitoredStdContainer("trk_pt_sel",         m_mon_trkSel_pt,       AutoClear);
  declareMonitoredStdContainer("trk_qOverPt_sel",    m_mon_trkSel_qOverPt,   AutoClear);
  declareMonitoredStdContainer("trk_eta_sel",        m_mon_trkSel_eta,       AutoClear);

  declareMonitoredStdContainer("trk_BlayerHits_sel", m_mon_trkSel_BlayerHits,AutoClear);
  declareMonitoredStdContainer("trk_PixHits_sel",    m_mon_trkSel_PixHits,   AutoClear);
  declareMonitoredStdContainer("trk_SCTHits_sel",    m_mon_trkSel_SCTHits,   AutoClear);
  declareMonitoredStdContainer("trk_TRTHits_sel",    m_mon_trkSel_TRTHits,   AutoClear);
  
  declareMonitoredStdContainer("trk_PixHoles_sel",   m_mon_trkSel_PixHoles,  AutoClear);
  declareMonitoredStdContainer("trk_SCTHoles_sel",   m_mon_trkSel_SCTHoles,  AutoClear);
  declareMonitoredStdContainer("trk_TRTHoles_sel",   m_mon_trkSel_TRTHoles,  AutoClear);
  
  declareMonitoredStdContainer("trk_PixSharedHits_sel", m_mon_trkSel_PixSharedHits,   AutoClear);
  declareMonitoredStdContainer("trk_SCTSharedHits_sel", m_mon_trkSel_SCTSharedHits,   AutoClear);
  //declareMonitoredStdContainer("trk_TRTSharedHits_sel", m_mon_trkSel_TRTSharedHits,   AutoClear);

  declareMonitoredVariable    ("roi_nTracks",       m_totTracks);
  declareMonitoredVariable    ("roi_nTracks_sel",   m_totSelTracks);
  declareMonitoredStdContainer("roi_stepsToSelect", m_listCutApplied, AutoClear);
  declareMonitoredObject      ("roi_selectedTracks", *this, &TrigDvFex::totSelectedTracks);

  declareMonitoredVariable    ("roi_unusedhits" ,         m_unusedhits,         AutoClear);
  declareMonitoredVariable    ("roi_unusedhits_fraction", m_unusedhits_fraction,AutoClear);
  
  declareMonitoredVariable    ("roi_deltaEtaTrk",    m_deltaEtaTrk,    AutoClear);
  declareMonitoredVariable    ("roi_deltaPhiTrk",    m_deltaPhiTrk,    AutoClear);

}


//** ----------------------------------------------------------------------------------------------------------------- **//


TrigDvFex::~TrigDvFex() 
{
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

    ATH_MSG_DEBUG(" AlgoId = "              << m_algo ); 
 
    ATH_MSG_DEBUG(" TrkSel_Chi2 = "     << m_trkSelChi2 ); 
    ATH_MSG_DEBUG(" TrkSel_BLayer = "   << m_trkSelBLayer ); 
    ATH_MSG_DEBUG(" TrkSel_SiHits = "   << m_trkSelSiHits ); 
    ATH_MSG_DEBUG(" TrkSel_D0 = "       << m_trkSelD0 ); 
    ATH_MSG_DEBUG(" TrkSel_Z0 = "       << m_trkSelZ0 ); 
    ATH_MSG_DEBUG(" TrkSel_Pt = "       << m_trkSelPt ); 

    m_trigBjetJetInfo    = new TrigBjetJetInfo();

    // Retrieving Region Selector Tool
    if ( m_regionSelector.retrieve().isFailure() ) 
    {
        ATH_MSG_FATAL("Unable to retrieve RegionSelector tool "
                << m_regionSelector.type());
        return HLT::ErrorCode(HLT::Action::ABORT_JOB, HLT::Reason::BAD_JOB_SETUP);
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

    // Clear and initialize data members ---> put it in a fucntion and be sure is DONE for all variables!
    // FIXME
    m_etaRoI = -999;
    m_phiRoI = -999;
    m_totSelTracks = 0;
    m_totTracks    = 0;
    m_unusedhits   = 0;
    m_unusedhits_fraction = 0.0;
    m_setPVInfo    = false;
  
    m_trigBjetJetInfo->clear();

    // -----------------------------------
    // Get RoI descriptor
    // -----------------------------------
    const TrigRoiDescriptor* roiDescriptor = 0;
    if(getFeature(inputTE, roiDescriptor, m_jetKey) != HLT::OK) 
    {
        ATH_MSG_ERROR("No feature for this Trigger Element");    	
    	return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::NAV_ERROR);
    }
    m_etaRoI = roiDescriptor->eta();
    m_phiRoI = phiCorr(roiDescriptor->phi());
    ATH_MSG_DEBUG("Using TE: " << "RoI id " << roiDescriptor->roiId()
	    << ", Phi = " <<  m_phiRoI << ", Eta = " << m_etaRoI);

    // -----------------------------------
    // Get Hit related collection
    // -----------------------------------
    // update the list of PRDs if needed ---> FIXME? Put it inside getPRDs?
    HLT::ErrorCode ec = updatePRDs(TrigDvFexDet::PIXEL);
    if( ec != HLT::OK )
    {
        return ec;
    }
    ec = updatePRDs(TrigDvFexDet::SCT);
    if( ec != HLT::OK )
    {
        return ec;
    }
    ec = updatePRDs(TrigDvFexDet::TRT);
    if( ec != HLT::OK )
    {
        return ec;
    }
    // Get the PRDs available in the RoI 
    const PrdsAtRoI prds_at_roi = getPRDsfromID(roiDescriptor);
    
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
    std::unordered_set<const Trk::PrepRawData*> prds_at_tracks;
    if(pointerToEFTrackCollections)
    {
        m_totTracks = pointerToEFTrackCollections->size();
        for(unsigned int j = 0; j < m_totTracks; ++j) 
        {
            const xAOD::TrackParticle* track = (*pointerToEFTrackCollections)[j];
            // updating the prds used by the tracks
            updateTrackPRDs(track,prds_at_tracks);

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
    ATH_MSG_DEBUG("Found " << m_totTracks << " tracks in the RoI " 
           << " using a total of " << prds_at_tracks.size() << " PRDs");

    // Obtain the number of unused hits in the RoI
    std::map<DETID,int> unusedhits = { {PIXEL,0}, {SCT,0},{TRT,0} } ;
    int unusedhits_blayer = 0;
    
    unsigned int _usedHits = 0;
    // Was used in the prds' tracks -->> change code..?..
    // Pixel
    for(auto & currentPRD: prds_at_roi.pixel)
    {
        if( prds_at_tracks.find(currentPRD) == prds_at_tracks.end() )
        {
            ++(unusedhits[PIXEL]);
            // Is in the B-layer?
            if( prds_at_roi.blayer.find(currentPRD) != prds_at_roi.blayer.end() )
            {
                ++unusedhits_blayer;
            }
        }
        else
        {
            ++_usedHits;
        }
    }
    // SCT        
    for(auto & currentPRD: prds_at_roi.sct)
    {
        if( prds_at_tracks.find(currentPRD) == prds_at_tracks.end() )
        {
            ++(unusedhits[SCT]);
        }
        else
        {
            ++_usedHits;
        }
    }
    // TRT
    for(auto & currentPRD: prds_at_roi.trt)
    {
        if( prds_at_tracks.find(currentPRD) == prds_at_tracks.end() )
        {
            ++(unusedhits[TRT]);
        }
        else
        {
            ++_usedHits;
        }
    }

    m_unusedhits = (unusedhits[PIXEL]+unusedhits[SCT]+unusedhits[TRT]);

    
    if( prds_at_roi.size() > 0)
    {
        m_unusedhits_fraction = m_unusedhits/((float)(prds_at_roi.size()));
    }
    else
    {
        m_unusedhits_fraction = 0.0;
    }

    // blayer
    float unusedhits_blayer_frac = 0.0;
    if( prds_at_roi.blayer.size() > 0 )
    {
        unusedhits_blayer_frac = ((float)unusedhits_blayer)/((float)prds_at_roi.blayer.size());
    }
    
    // debug information
    if(msgLvl() <= MSG::DEBUG) 
    {
        const unsigned int prds_pix_total = prds_at_roi.pixel.size();
        const unsigned int prds_sct_total = prds_at_roi.sct.size();
        const unsigned int prds_trt_total = prds_at_roi.trt.size();
        
        ATH_MSG_DEBUG("Number of unused hits: " << m_unusedhits 
                << " (from a total of " << prds_at_roi.size() << " in the RoI, i.e. " 
                << (m_unusedhits_fraction*100.0) << "%)");
        ATH_MSG_DEBUG("Number of used hits: " << _usedHits);
        ATH_MSG_DEBUG("  * Pixel: " << (unusedhits[PIXEL]) << " - "
                << ((prds_pix_total != 0) ? ((float)(unusedhits[PIXEL]))/((float)prds_pix_total)*100. : 0) << "%");
        ATH_MSG_DEBUG("   * SCT  : " << (unusedhits[SCT])  << " - "
                << ((prds_sct_total != 0) ? ((float)(unusedhits[SCT]))/((float)prds_sct_total)*100. : 0) << "%");
        ATH_MSG_DEBUG("   * TRT  : " << (unusedhits[TRT])  << " - "
                << ((prds_trt_total != 0) ? ((float)(unusedhits[TRT]))/((float)prds_trt_total)*100. : 0) << "%");
    }
    // --- JetInfo
    // Jets: remember in some previous step, each jet has been associated to
    // a RoI, so the RoI and the jet information is the same
    m_trigBjetJetInfo->setEtaPhiRoI(roiDescriptor->eta(), phiCorr(roiDescriptor->phi()));
    const xAOD::Jet * thejet = 0;
    HLT::ErrorCode xaodc = getJet(thejet,inputTE);
    if( xaodc != HLT::OK )
    {
        return xaodc;
    }	
    
    if(thejet)
    {
        m_trigBjetJetInfo->setEtaPhiJet(thejet->p4().Eta(),thejet->p4().Phi());
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
    TrigEFBjet* trigEFBjet = new TrigEFBjet(
        roiDescriptor->roiId(),  // roi
	    m_trigBjetJetInfo->etaJet(), // eta
        m_trigBjetJetInfo->phiJet(), // phi
	    0,                        //  TrackCollection
        0,                        //  PriVtxCollection
        0,                        //  SecVtxCollection
        prds_at_roi.pixel.size(), // prmVtx
        m_trigBjetJetInfo->etJet(), // ptJet
	    unusedhits[PIXEL],     // xComb 
        prds_at_roi.sct.size(),// xIP1d
        unusedhits[SCT],       // xIP2d
        prds_at_roi.trt.size(),// xIP3d
        unusedhits[TRT],       // xChi2 
	    10.0*CLHEP::mm,        // xSv 
        10.0*CLHEP::GeV,       // xMvtx
        unusedhits_blayer_frac,// xEvtx
        m_totSelTracks);       // xNvtx

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
HLT::ErrorCode TrigDvFex::hltEndEvent() 
{
    ATH_MSG_INFO("hltEndEvent: End event resetting");
    // resetting the pointers
    m_pixcontainer_v = 0;
    m_is_pix_filled  = false;

    m_sctcontainer_v = 0;
    m_is_sct_filled  = false;

    m_trtcontainer_v = 0;
    m_is_trt_filled  = false;

    m_are_prds_filled= false;

    return HLT::OK;
}

//** ----------------------------------------------------------------------------------------------------------------- **//
HLT::ErrorCode TrigDvFex::hltFinalize() 
{
    ATH_MSG_INFO("Finalizing TrigDvFex");
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
    
    // New
    uint8_t nPixHoles     = 0;
    uint8_t nPixSharedHits= 0;
    track->summaryValue(nPixHoles,      xAOD::numberOfPixelHoles);
    track->summaryValue(nPixSharedHits, xAOD::numberOfPixelSharedHits);

    uint8_t nSCTHoles     = 0;
    uint8_t nSCTSharedHits= 0;
    track->summaryValue(nSCTHoles,      xAOD::numberOfSCTHoles);
    track->summaryValue(nSCTSharedHits, xAOD::numberOfSCTSharedHits);

    uint8_t nTRTHits      = 0;
    uint8_t nTRTHoles     = 0;
    //uint8_t nTRTSharedHits= 0;
    track->summaryValue(nTRTHits,       xAOD::numberOfTRTHits);
    track->summaryValue(nTRTHoles,      xAOD::numberOfTRTHoles);
    //track->summaryValue(nTRTSharedHits, xAOD::numberOfTRTSharedHits);

    
    
    const int   nSiHits  = nPixHits + nSCTHits;
    // const int   totalHits= nSiHits + nTRTHits; --< Not used by the moment
    const float theta    = track->theta();
    const float phi      = phiCorr(track->phi());
    const float eta      = track->eta();
    const float qOverPt  = track->qOverP()/TMath::Sin(theta); 
    //const float pT      = (1.0/qOverPt);
    const float pT       = track->pt();
    const float d0       = track->d0();
    const float z0       = track->z0();


    ATH_MSG_VERBOSE( "efTrackSel method\n" <<
	    "  Track number "    << i+1  << " to be selected must be:\n" <<
       	    "    Pt [GeV]:"  << fabs(pT)/CLHEP::GeV           << " >= " << m_trkSelPt/CLHEP::GeV << "\n" <<
       	    "    d0 [mm]:"   << fabs(d0)/CLHEP::mm            << " >= " << m_trkSelD0/CLHEP::mm << "\n" <<
	   // "    z0*sin(theta) " << fabs(z0-zv)*TMath::Sin(theta) << " <= " << m_trkSelZ0 << "\n" <<
       	    "    bLayer "        << (int)nBlayerHits              << " >= " << m_trkSelBLayer << "\n" <<
       	    "    pixelHit "      << (int)nPixHits                 << " >= " << m_trkSelPixHits<< "\n" <<
            "    SiHit "         << (int)nSiHits                  << " >= " << m_trkSelSiHits << "\n" <<
       	    "    Prob(chi2) "    << TMath::Prob(track->chiSquared(), (int)nSiHits*3-5) << " > " << m_trkSelChi2);
    
    // Monitoring stuff: Validation-----------------------------------------------
    if(m_mon_validation)
    {
        m_mon_trk_a0.push_back(d0/CLHEP::mm);
        m_mon_trk_z0.push_back(z0/CLHEP::mm);
        m_mon_trk_theta.push_back(theta);
        m_mon_trk_pt.push_back(pT/CLHEP::GeV);
        m_mon_trk_qOverPt.push_back(qOverPt);
        m_mon_trk_eta.push_back(eta);

        // Quality
        m_mon_trk_BlayerHits.push_back((int)nBlayerHits);
        m_mon_trk_PixHits.push_back((int)nPixHits);
        m_mon_trk_SCTHits.push_back((int)nSCTHits);
        m_mon_trk_TRTHits.push_back((int)nTRTHits);
        m_mon_trk_PixHoles.push_back((int)nPixHoles);
        m_mon_trk_SCTHoles.push_back((int)nSCTHoles);
        m_mon_trk_TRTHoles.push_back((int)nTRTHoles);
        m_mon_trk_PixSharedHits.push_back((int)nPixSharedHits);
        m_mon_trk_SCTSharedHits.push_back((int)nSCTSharedHits);
        //m_mon_trk_TRTSharedHits.push_back((int)nTRTSharedHits);
    }
        
    if(m_useEtaPhiTrackSel) 
    {
    	if(std::abs(eta - m_etaRoI) > 0.2) 
        {
            ATH_MSG_DEBUG("  track " << i+1 << " is not selected (eta matching)");
            m_listCutApplied.push_back(CutListMonitor::RoIEtaMatching); 
            return false;
        }
        
        if(std::abs(phiCorr(phi- m_phiRoI)) > 0.2) 
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
    if(m_mon_validation )
    {
        m_mon_trkSel_theta.push_back(theta);
        m_mon_trkSel_pt.push_back(pT/CLHEP::GeV);
        m_mon_trkSel_qOverPt.push_back(qOverPt);
        m_mon_trkSel_eta.push_back(eta);

        // Quality
        m_mon_trkSel_BlayerHits.push_back((int)nBlayerHits);
        m_mon_trkSel_PixHits.push_back((int)nPixHits);
        m_mon_trkSel_SCTHits.push_back((int)nSCTHits);
        m_mon_trkSel_TRTHits.push_back((int)nTRTHits);
        m_mon_trkSel_PixHoles.push_back((int)nPixHoles);
        m_mon_trkSel_SCTHoles.push_back((int)nSCTHoles);
        m_mon_trkSel_TRTHoles.push_back((int)nTRTHoles);
        m_mon_trkSel_PixSharedHits.push_back((int)nPixSharedHits);
        m_mon_trkSel_SCTSharedHits.push_back((int)nSCTSharedHits);
        //m_mon_trkSel_TRTSharedHits.push_back((int)nTRTSharedHits);
        
        // RoI related
        m_deltaEtaTrk = m_etaRoI-eta;
        m_deltaPhiTrk = phiCorr(m_phiRoI-phi);
    }   

    if(m_mon_validation || m_mon_online)
    {
        m_mon_trk_a0_sel.push_back(d0/CLHEP::mm);
        m_mon_trk_z0_sel.push_back(z0/CLHEP::mm);
        
        // ez0
	    const float errIP1D = Amg::error(track->definingParametersCovMatrix(),1);
        // ed0
	    const float errIP2D = Amg::error(track->definingParametersCovMatrix(),0);
        
        // Significance d0, z0: Note that the beam spot is not included meaning
        // that the significance is going to grow (because of the inherent width of
        // the beam)/
        float IP1D = 0.0;
        float IP2D = 0.0;
        if(fabs(errIP1D) > TOLERANCE)
        {
            IP1D = z0/sqrt(errIP1D*errIP1D);
        }
        if( (fabs(errIP2D) > TOLERANCE) )
        {
            IP2D = d0/sqrt(IP2D*errIP2D);
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
                   << "  Eta=" << aJet->p4().Eta() << " Phi= " << aJet->p4().Phi());
        }
    }
    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//
// Correct the phi angle to be defined between [-PI,PI]
float TrigDvFex::phiCorr(const float & phi) const
{
    float phicorr = phi;
    if(phicorr < -M_PI)
    { 
        phicorr += 2.0*M_PI;
    }
    else if(phi > M_PI) 
    {
        phicorr -= 2.0*M_PI;
    }
    
    return phicorr;
}

// Update the full list of PRDs per Event. The full filling job is done only
// once per Event, because the PRDs collection belong to an Event (the same
// collections for the same Event's RoIs)
HLT::ErrorCode TrigDvFex::updatePRDs(const TrigDvFexDet & detsystem)
{
    // check if it was already filled ( with the first RoI of the event)
    if( m_are_prds_filled )
    {
        return HLT::OK;
    }

    std::string colname;
    StatusCode sc(0);
    // Getting the prepared raw data
    if( detsystem == TrigDvFexDet::PIXEL ) 
    {
        // Just do it at the beginning of the event
        if( m_is_pix_filled )
        {
            return HLT::OK;
        }
        colname = "PixelTrigClusters"; // m_pixhitsname;
        m_is_pix_filled = true;
        sc = evtStore()->retrieve(m_pixcontainer_v,colname);
    }
    else if( detsystem == TrigDvFexDet::SCT )
    {
        // Just do it at the beginning of the event
        if( m_is_sct_filled )
        {
            return HLT::OK;
        }
        colname = "SCT_TrigClusters"; // m_scthitsname;
        m_is_sct_filled = true;
        sc = evtStore()->retrieve(m_sctcontainer_v,colname);
    }
    else if( detsystem == TrigDvFexDet::TRT )
    {
        // Just do it at the beginning of the event
        if( m_is_trt_filled )
        {
            return HLT::OK;
        }
        colname = "TRT_TrigDriftCircles"; // m_trthitsname;
        sc = evtStore()->retrieve(m_trtcontainer_v,colname);
        m_is_trt_filled = true;
    }
    else
    {
        ATH_MSG_ERROR("Inconsistent key for the detector enum of TrigDvFex algorithm!");
        return HLT::ErrorCode(HLT::Action::ABORT_JOB, HLT::Reason::BAD_JOB_SETUP);
    }

    if(sc.isFailure())
    {
        ATH_MSG_ERROR("Not found the '" << colname << "' collection '");
        return HLT::ErrorCode(HLT::Action::ABORT_JOB, HLT::Reason::BAD_JOB_SETUP);
    }
    
    m_are_prds_filled = m_is_pix_filled && m_is_sct_filled && m_is_trt_filled;

    // Just in DEBUG mode show some info
    if(msgLvl() <= MSG::DEBUG && m_are_prds_filled )  
    {
        ATH_MSG_INFO("Total number of PrepRawData extracted from the Event");
        // Pixel
        unsigned int nPRDsPix = 0;
        for(auto & pixhitcollection: * m_pixcontainer_v)
        {
            if( pixhitcollection == 0)
            {
                continue;
            }
            nPRDsPix += pixhitcollection->size();
        }
        ATH_MSG_INFO("Got " << nPRDsPix << " PRDs from Pixel"); 
        // SCT
        unsigned int nPRDsSCT = 0;
        for(auto & scthitcollection: *m_sctcontainer_v)
        {
            if( scthitcollection == 0 )
            {
                continue;
            }
            nPRDsSCT += scthitcollection->size();
        }
        ATH_MSG_INFO("Got " << nPRDsSCT << " PRDs from SCT" );
        // TRT
        unsigned int nPRDsTRT= 0;
        for(auto & trthitcollection: *m_trtcontainer_v)
        {
            if(trthitcollection == 0)
            {
                continue;
            }
            nPRDsTRT += trthitcollection->size();
        }
        ATH_MSG_INFO("Got " << nPRDsTRT << " PRDs from TRT");
    }

    return HLT::OK;
}


// Extract the PRDs which belongs to detector elements placed at 
// the current RoI 
PrdsAtRoI TrigDvFex::getPRDsfromID(const TrigRoiDescriptor * roi)
{
    // =====================================================
    // Obtaining the hits inside this RoI 
    
    // Initialize vector of PRDs
    //std::vector<const Trk::PrepRawData*> prds_at_roi; 
    PrdsAtRoI prds_at_roi;
    
    // Pixel hash id's
    std::vector<IdentifierHash> listOfPixIds;
    m_regionSelector->DetHashIDList(PIXEL,*roi,listOfPixIds);
    ATH_MSG_DEBUG("Found " << listOfPixIds.size() << " PIXEL det. elements");
    // Pixel hits
    InDet::SiClusterContainer::const_iterator endPixContainer = m_pixcontainer_v->end();
    for(auto & id_pix: listOfPixIds)
    {
        InDet::SiClusterContainer::const_iterator currentPixCol = m_pixcontainer_v->indexFind(id_pix);
        {
            // This pixel prds element is not in the PRDs collection
            // or is empty
            if(currentPixCol == endPixContainer || (*currentPixCol) == 0)
            {
                continue;
            }

            const InDet::SiClusterCollection * pixcollection = *currentPixCol;
            for( auto & pixhit: * pixcollection )
            {
                prds_at_roi.pixel.insert(pixhit);
                if( pixhit->detectorElement()->isBlayer() )
                {
                    prds_at_roi.blayer.insert(pixhit);
                }
            }
        }
    }
    ATH_MSG_DEBUG("Number of Pixel PRDs at this RoI: " << prds_at_roi.pixel.size()
            << "(" << prds_at_roi.blayer.size() << " B-layers)"); 
    
    // SCT has id's
    std::vector<IdentifierHash> listOfSCTIds;
    m_regionSelector->DetHashIDList(SCT,*roi,listOfSCTIds);
    ATH_MSG_DEBUG("Found " << listOfSCTIds.size() << " SCT det. elements");
    // SCT hits
    InDet::SiClusterContainer::const_iterator endSCTContainer = m_sctcontainer_v->end();
    for(auto & id_sct: listOfSCTIds)
    {
        InDet::SiClusterContainer::const_iterator currentSCTCol = m_sctcontainer_v->indexFind(id_sct);
        {
            // This SCT prds element is not in the PRDs collection
            // or is empty
            if(currentSCTCol == endSCTContainer || (*currentSCTCol) == 0)
            {
                continue;
            }

            const InDet::SiClusterCollection * sctcollection = *currentSCTCol;
            
            for( auto & scthit: * sctcollection )
            {
                prds_at_roi.sct.insert(scthit);
            }
        }
    }
    ATH_MSG_DEBUG("Number of SCT PRDs at this RoI: " << prds_at_roi.sct.size());
    
    // TRT hash id's
    std::vector<IdentifierHash> listOfTRTIds;
    m_regionSelector->DetHashIDList(TRT,*roi,listOfTRTIds);
    ATH_MSG_DEBUG("Found " << listOfTRTIds.size() << " TRT det. elements");
    // TRT hits
    InDet::TRT_DriftCircleContainer::const_iterator endTRTContainer = m_trtcontainer_v->end();
    for(auto & id_trt: listOfTRTIds)
    {
        InDet::TRT_DriftCircleContainer::const_iterator currentTRTCol = m_trtcontainer_v->indexFind(id_trt);
        {
            // This TRT prds element is not in the PRDs collection
            // or is empty
            if(currentTRTCol == endTRTContainer || (*currentTRTCol) == 0)
            {
                continue;
            }

            const InDet::TRT_DriftCircleCollection * trtcollection = *currentTRTCol;
            
            for( auto & trthit: * trtcollection )
            {
                //prds_at_roi.push_back(trthit);
                prds_at_roi.trt.insert(trthit);
            }
        }
    }
    ATH_MSG_DEBUG("Number of TRT PRDs at this RoI: " << prds_at_roi.trt.size());
    ATH_MSG_INFO("Total Number of PRDs at this RoI: " << prds_at_roi.size());
    
    return prds_at_roi;
}


// Update the set containing the PRDs used by the tracks, with the i-track
void TrigDvFex::updateTrackPRDs(const xAOD::TrackParticle * track, std::unordered_set<const Trk::PrepRawData*> & prds)
{ 
    // Retrieve the Track link
    const Trk::Track * trk_track = track->track();
    if( track != 0 )
    {
        for(DataVector<const Trk::MeasurementBase>::const_iterator it = trk_track->measurementsOnTrack()->begin();
                it != trk_track->measurementsOnTrack()->end(); ++it)
        {
            const Trk::RIO_OnTrack * rot = dynamic_cast<const Trk::RIO_OnTrack*>(*it);
            if( rot != 0 )
            {
                prds.insert(rot->prepRawData());
            }
        }
    }
    else
    { 
        ATH_MSG_WARNING("No link to Trk::Track!! Cannot fill the 'Unused hits' variable"); 
    }
    ATH_MSG_INFO("Number of PRDs on this track: " << trk_track->measurementsOnTrack()->size());
}

