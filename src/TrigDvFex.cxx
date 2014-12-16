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

//#include "TrigInDetEvent/TrigInDetTrackCollection.h"
#include "Particle/TrackParticleContainer.h"
//#include "TrigInDetEvent/TrigVertexCollection.h"
#include "VxSecVertex/VxSecVertexInfo.h"
#include "VxSecVertex/VxSecVKalVertexInfo.h"
//#include "TrigCaloEvent/TrigT2Jet.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"
//#include "EventPrimitives/EventPrimitives.h"
//#include "EventPrimitives/EventPrimitivesHelpers.h"

#include "TrigParticle/TrigEFBjetContainer.h"

#include "TrigSteeringEvent/TrigOperationalInfo.h"

#include "InDetBeamSpotService/IBeamCondSvc.h"

//#include "TrigInDetEvent/TrigVertex.h"

//#include "xAODBTagging/BTaggingContainer.h"
//#include "xAODBTagging/BTagging.h"
//#include "xAODBTagging/BTaggingAuxContainer.h"

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/VertexAuxContainer.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODBase/IParticle.h"

//** ---------------------------------------------------------------------- **//
const double TOLERANCE=1e-30;

TrigDvFex::TrigDvFex(const std::string& name, ISvcLocator* pSvcLocator) :
  HLT::FexAlgo(name, pSvcLocator),
  m_trackJetFinderTool("TrigTrackJetFinderTool",this),      	
  m_trigEFBjetColl(0),
  m_trigBjetPrmVtxInfo(0),
  m_trigBjetSecVtxInfo(0),
  m_trigBjetJetInfo(0),
  m_totTracks(0),
  m_totSelTracks(0)
{
  declareProperty ("AlgoId",             m_algo);
  declareProperty ("Instance",           m_instance);
  declareProperty ("JetKey",             m_jetKey     = ""); //"" needed for default config, SplitJet for new config
  declareProperty ("PriVtxKey",          m_priVtxKey  = "EFHistoPrmVtx"); //Does this still work with default config?

  declareProperty ("onlinemon",          m_mon_online = true);
  declareProperty ("validation",         m_mon_validation = true);
  
  declareProperty ("UseBeamSpotFlag",    m_useBeamSpotFlag    = false);
  declareProperty ("SetBeamSpotWidth",   m_setBeamSpotWidth   = 0.05);

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
  if(m_trigBjetPrmVtxInfo)      delete m_trigBjetPrmVtxInfo;
  if(m_trigBjetSecVtxInfo)      delete m_trigBjetSecVtxInfo;
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
 
    ATH_MSG_DEBUG(" UseBeamSpotFlag = "     << m_useBeamSpotFlag ); 
    ATH_MSG_DEBUG(" SetBeamSpotWidth = "    << m_setBeamSpotWidth );

    ATH_MSG_DEBUG(" HistoPrmVtxAtEF = "     << m_histoPrmVtxAtEF );
    ATH_MSG_DEBUG(" UseEtaPhiTrackSel = "   << m_useEtaPhiTrackSel );
    ATH_MSG_DEBUG(" UseJetDirection = "     << m_useJetDirection );

    ATH_MSG_DEBUG(" TrkSel_Chi2 = "     << m_trkSelChi2 ); 
    ATH_MSG_DEBUG(" TrkSel_BLayer = "   << m_trkSelBLayer ); 
    ATH_MSG_DEBUG(" TrkSel_SiHits = "   << m_trkSelSiHits ); 
    ATH_MSG_DEBUG(" TrkSel_D0 = "       << m_trkSelD0 ); 
    ATH_MSG_DEBUG(" TrkSel_Z0 = "       << m_trkSelZ0 ); 
    ATH_MSG_DEBUG(" TrkSel_Pt = "       << m_trkSelPt ); 

    // Retrieve TrigTrackJetFinder tool
    StatusCode sc = m_trackJetFinderTool.retrieve();
    if(sc.isFailure()) 
    {
      ATH_MSG_FATAL("Failed to locate tool " << m_trackJetFinderTool);
      return HLT::BAD_JOB_SETUP;
    }
    ATH_MSG_INFO("Retrieved tool " << m_trackJetFinderTool);
  
    m_trigBjetPrmVtxInfo = new TrigBjetPrmVtxInfo();
    m_trigBjetSecVtxInfo = new TrigBjetSecVtxInfo();
    m_trigBjetJetInfo    = new TrigBjetJetInfo();

    return HLT::OK;
}



//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::getTrackCollection(const xAOD::TrackParticleContainer*& pointerToEFTrackCollections,
                 const HLT::TriggerElement* whateverTE)
{
    std::vector<const xAOD::TrackParticleContainer*> vectorOfEFTrackCollections;
    HLT::ErrorCode status = getFeatures(whateverTE, vectorOfEFTrackCollections, ""); 
    if(status != HLT::OK) 
    {
       ATH_MSG_ERROR("Failed to get TrackParticleContainer from the trigger element");
    } 
    else if(msgLvl() <= MSG::DEBUG) 
    { 
       ATH_MSG_DEBUG("Got " << vectorOfEFTrackCollections.size() << " TrackParticleContainer");
    }
    
    if(vectorOfEFTrackCollections.begin() == vectorOfEFTrackCollections.end())
    {
       pointerToEFTrackCollections = 0;
       return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::MISSING_FEATURE);
    }
    
    pointerToEFTrackCollections = *(vectorOfEFTrackCollections.begin());
    return HLT::OK;
}

//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::getPrmVtxCollection(const xAOD::VertexContainer*& pointerToEFPrmVtxCollections, 
	const HLT::TriggerElement* whateverTE, const std::string & vtxkey) 
{
    std::vector<const xAOD::VertexContainer*> vectorOfEFPrmVtxCollections;
    HLT::ErrorCode status = getFeatures(whateverTE, vectorOfEFPrmVtxCollections, vtxkey);
    if(status != HLT::OK) 
    {
        ATH_MSG_ERROR("Failed to get xAOD::VertexContainer from the trigger element");
    } 
    else
    {
        ATH_MSG_DEBUG("Got " << vectorOfEFPrmVtxCollections.size() << " xAOD::VertexContainer");
    }
    
    std::vector<const xAOD::VertexContainer*>::iterator pPrmVtxColl    = vectorOfEFPrmVtxCollections.begin();
    std::vector<const xAOD::VertexContainer*>::iterator lastPrmVtxColl = vectorOfEFPrmVtxCollections.end();

    for(; pPrmVtxColl != lastPrmVtxColl; ++pPrmVtxColl)
    {
        ATH_MSG_DEBUG("Size of pPrmVtxColl = " << (*pPrmVtxColl)->size());
        if((*pPrmVtxColl)->size() == 0)
        {
            continue;
        }
        ATH_MSG_VERBOSE("xAOD::VertexContainer with label " << (*pPrmVtxColl)->front()->vertexType());    
        if( (*pPrmVtxColl)->front()->vertexType() != xAOD::VxType::PriVtx )
        {
            continue;
        }
        ATH_MSG_DEBUG("Selected collection with Primary Vertex label ");
        ATH_MSG_DEBUG("First PV has z-position = " <<  (*pPrmVtxColl)->front()->z());
        
        pointerToEFPrmVtxCollections = *pPrmVtxColl;
        return HLT::OK;
    }
    // Didn't found any PV 
    pointerToEFPrmVtxCollections = 0;
    ATH_MSG_DEBUG("Not found any Primary Vertex valid collection!");
    // return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::MISSING_FEATURE);
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::getSecVtxCollection(const Trk::VxSecVertexInfoContainer*& pointerToSecVtxCollections,
       	const HLT::TriggerElement* whateverTE) 
{
    std::vector<const Trk::VxSecVertexInfoContainer*> vectorOfSecVtxCollections;
    // FIXME:: Harcoded name, change for a property
    HLT::ErrorCode status = getFeatures(whateverTE, vectorOfSecVtxCollections, "SecVtxInfo");
    if(status != HLT::OK) 
    {
       ATH_MSG_ERROR("Failed to get VxSecVertexInfoContainer from the trigger element");
       return status;
    }
    ATH_MSG_DEBUG("Got " << vectorOfSecVtxCollections.size() << " VxSecVertexInfoContainer");
  
    std::vector<const Trk::VxSecVertexInfoContainer*>::iterator pSecVtxColl    = vectorOfSecVtxCollections.begin();
    std::vector<const Trk::VxSecVertexInfoContainer*>::iterator lastSecVtxColl = vectorOfSecVtxCollections.end();
  
    // JDC:: What's the logic? A vector of vectors of sec. vtx. 
    for( ; pSecVtxColl != lastSecVtxColl; ++pSecVtxColl)
    {
	    if((*pSecVtxColl)->size() == 0)
        {   
            continue;
        }
        const std::vector<xAOD::Vertex*> vectorOfVxCandidates = (*pSecVtxColl)->front()->vertices();
    
        if(vectorOfVxCandidates.size() == 0)
        {
            continue;
        }
        
        ATH_MSG_VERBOSE("VxSecVertexInfoContainer with label type " << 
                (*vectorOfVxCandidates.begin())->vertexType());
        
        if((*vectorOfVxCandidates.begin())->vertexType() == 2) 
        {
            ATH_MSG_DEBUG("Selected collection with Secondary Vertex label"); 
            ATH_MSG_DEBUG("First SV has z-position = " << (*vectorOfVxCandidates.begin())->z());
            // Just need the first one? I'm assuming the previous alg, put all the SV in the same set
            pointerToSecVtxCollections = *pSecVtxColl;
            return HLT::OK;
        }
    }
    // Didn't found any SV
    pointerToSecVtxCollections = 0;
    ATH_MSG_DEBUG("Not found any Secondary Vertex valid collection!");
    return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    //return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::MISSING_FEATURE);
}


//** ----------------------------------------------------------------------------------------------------------------- **//

// This method should use a new data-member called m_svcollection, in order to be sure we properly
// set this collection in the getSecVtx or maybe this method should be absorbed in the getSecVtxCollection
// method?? This will avoid a lot of cross-checsk and ifs...
HLT::ErrorCode TrigDvFex::setSecVtxInfo(const Trk::VxSecVertexInfoContainer*& pointerToEFSecVtxCollections, 
					  const xAOD::Vertex* & pvselected) 
{
    if(!pvselected)
    {
        ATH_MSG_DEBUG("No primary vertex collection sent when extracting sec vtx info");
        return HLT::OK;
    }

    if(!pointerToEFSecVtxCollections) 
    {
        ATH_MSG_DEBUG("No secondary vertex collection sent when extracting sec vtx info");
        return HLT::OK; // FIXME I need this, so you should stop the FEX
    }
  
    const Trk::VxSecVertexInfo* m_secVertexInfo = (*pointerToEFSecVtxCollections)[0];
    if(!m_secVertexInfo) 
    {   
       ATH_MSG_DEBUG("No secondary vertex when extracting sec vtx info");
       return HLT::OK; // FIXME I need this so you should stop the FEX
    }
  
    const Trk::VxSecVKalVertexInfo * myVKalSecVertexInfo = dynamic_cast<const Trk::VxSecVKalVertexInfo*>(m_secVertexInfo);
 
    if(!myVKalSecVertexInfo) 
    {
        ATH_MSG_DEBUG("The cast to VKal secondary vertex went wrong, the pointer is zero");
        return HLT::OK;
    }
  
    const std::vector<xAOD::Vertex*> & myVertices = myVKalSecVertexInfo->vertices();
   
    if(myVertices.size() == 0)
    {
        return HLT::OK;
    }
 
    m_trigBjetSecVtxInfo->setVtxMass(myVKalSecVertexInfo->mass());
    m_trigBjetSecVtxInfo->setEnergyFraction(myVKalSecVertexInfo->energyFraction());
    m_trigBjetSecVtxInfo->setN2TrkVtx(myVKalSecVertexInfo->n2trackvertices());
    
    int NTracksInSV=0;
    
    if(myVertices.size()>1) 
    {
        ATH_MSG_WARNING("Secondary vertex from InDetVKalVxInJet has more than one vertex, is this ok?");
    }
    
    for(xAOD::Vertex * vertexIt : myVertices) 
    {
        if(!(vertexIt)) 
        {
    	    ATH_MSG_DEBUG("Secondary vertex from InDetVKalVxInJet has zero pointer. Skipping this vtx..");
    	    continue;
        }
        ATH_MSG_DEBUG("Secondary Vertex Candidate at (" 
    	    << vertexIt->position().x() << "," 
    	    << vertexIt->position().y() << "," 
         	    << vertexIt->position().z());
    
        // Check we have vxTrackAtVertex
        std::vector<Trk::VxTrackAtVertex> * myTracks = 0;
        if( vertexIt->vxTrackAtVertexAvailable() )
        {
            myTracks = &(vertexIt->vxTrackAtVertex());
            NTracksInSV+=myTracks->size();
        }
        else
        {
            ATH_MSG_WARNING("TrackAtVertex object not available! Setting 0 tracks attached to this vertex");
        }
    }
    m_trigBjetSecVtxInfo->setNTrksInVtx(NTracksInSV);
    
    //Calculate decay length and significance here
    ///Use the same utilities as in InDetVKalVxInJet
    //Note that this never happend because the pvselected should be filled at least at 000
    if(!pvselected) 
    {
        ATH_MSG_ERROR("Not filled properly the Primary Vertex (needed at least one default at (0,0,0)");
        return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::UNKNOWN);
    }
  
    ATH_MSG_DEBUG("Primary vertex for decay length (" 
    	<< pvselected->position().x() << "," 
    	<< pvselected->position().y() << "," 
    	<< pvselected->position().z() << ") and error (" 
    	<< pvselected->covariancePosition()(0,0) << "," 
    	<< pvselected->covariancePosition()(1,1) << "," 
    	<< pvselected->covariancePosition()(2,2) << ")");
    
    //Toggle
    m_trigBjetSecVtxInfo->isValid(true);

    // Getting the first vertex??
    const Amg::Vector3D SecVrt = myVertices[0]->position();

    double distance = (SecVrt - pvselected->position()).mag();
    // Note, storing the distance 3D between SV and PV, watch out with the name 
    // of the method!!
    m_trigBjetSecVtxInfo->setDecayLengthSignificance(distance);
    if(msgLvl() <= MSG::DEBUG) 
    {
        double dist2D = (SecVrt - pvselected->position()).perp();            
	    ATH_MSG_DEBUG("Calculating secondary vertex decay length with primary vertex at (" 
                << pvselected->position().x()/CLHEP::mm << "," << pvselected->position().y()/CLHEP::mm
                << "," << pvselected->position().z()/CLHEP::mm << ") and sec vtx at ("
                << SecVrt.x()/CLHEP::mm << "," << SecVrt.y()/CLHEP::mm << "," << SecVrt.z()/CLHEP::mm 
                <<  ") which gives 3D decay length " << distance/CLHEP::mm 
                << " and 2D(R/phi) decay length " 
                << dist2D/CLHEP::mm << " [all in mm]");
    }
    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//
bool TrigDvFex::efTrackSel(const xAOD::TrackParticle*& track, unsigned int i) 
{
    float zv = m_trigBjetPrmVtxInfo->zPrmVtx();
    
    uint8_t nBlayerHits = 0;
    uint8_t nPixHits    = 0;  
    uint8_t nSCTHits    = 0; 
  
    track->summaryValue(nBlayerHits, xAOD::numberOfBLayerHits);
    track->summaryValue(nPixHits,    xAOD::numberOfPixelHits);
    track->summaryValue(nSCTHits,    xAOD::numberOfSCTHits);
    
    int   nSiHits = nPixHits + nSCTHits; //summary->get(Trk::numberOfPixelHits)+summary->get(Trk::numberOfSCTHits);
    float theta   = track->theta();
    float qOverPt = track->qOverP()/TMath::Sin(theta); 
    float pT      = (1.0/qOverPt);
    float d0      = track->d0();
    float z0      = track->z0();

    // FIX FOR REDEFINED IP REFERENCE (ATR-9051)
    // m_taggerHelper->IPCorr(track->measuredPerigee()->parameters()[Trk::d0], 
    //               track->measuredPerigee()->parameters()[Trk::z0], 
    //                    d0,z0,track->phi(), track->eta(), pT, 
    //                    m_trigBjetPrmVtxInfo->xBeamSpot(), m_trigBjetPrmVtxInfo->yBeamSpot()); 
    //d0=track->measuredPerigee()->parameters()[Trk::d0]; // THIS LINE WAS DOING NOTHING!
    //z0=track->measuredPerigee()->parameters()[Trk::z0]; 
    //             // THIS LINE WAS DOING NOTHING! //+m_trigBjetPrmVtxInfo->zBeamSpot(); 
    // END FIX 

    ATH_MSG_VERBOSE( "efTrackSel method\n" <<
	    "  Track number "    << i+1  << " to be selected must be:\n" <<
       	    "    Pt "            << fabs(pT)                      << " >= " << m_trkSelPt << "\n" <<
       	    "    d0 "            << fabs(d0)                      << " <= " << m_trkSelD0 << "\n" <<
	    "    z0*sin(theta) " << fabs(z0-zv)*TMath::Sin(theta) << " <= " << m_trkSelZ0 << "\n" <<
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
  
    if(fabs(d0) > m_trkSelD0) 
    {
        ATH_MSG_DEBUG("  track " << i+1 << " not selected (d0 cut)");
        m_listCutApplied.push_back(CutListMonitor::D0Cut); 
        return false;
    }

    if(fabs(z0-zv)*TMath::Sin(theta) > m_trkSelZ0) 
    {
       ATH_MSG_DEBUG("  track " << i+1 << " not selected (z0 cut)");
       m_listCutApplied.push_back(CutListMonitor::Z0Cut); 
       return false;
    }

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
        m_mon_trk_z0_sel_PV.push_back(track->z0()-m_trigBjetPrmVtxInfo->zPrmVtx());
        
        // ez0
	    const float errIP1D = Amg::error(track->definingParametersCovMatrix(),1);
        // ed0
	    const float errIP2D = Amg::error(track->definingParametersCovMatrix(),0);
    	const float z0 = track->z0()-m_trigBjetPrmVtxInfo->zPrmVtx();

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
        const float sigmaBeamSpot = (m_trigBjetPrmVtxInfo->xBeamSpotWidth()+
                m_trigBjetPrmVtxInfo->yBeamSpotWidth())/2.0;
        float IP1D = 0.0;
        float IP2D = 0.0;
        if(fabs(errIP1D) > TOLERANCE)
        {
            IP1D = z0Sign/sqrt(errIP1D*errIP1D);
        }
        if((fabs(errIP2D) > TOLERANCE) || sigmaBeamSpot)
        {
            IP2D = d0Sign/sqrt(IP2D*errIP2D+sigmaBeamSpot*sigmaBeamSpot);
        }
        m_mon_trk_Sz0_sel.push_back(IP1D);
        m_mon_trk_Sa0_sel.push_back(IP2D);
    }
    
    return true;
}

//** ----------------------------------------------------------------------------------------------------------------- **//
HLT::ErrorCode TrigDvFex::checkxAODJets(const HLT::TriggerElement* inputTE)
{
    //const xAOD::JetContainer* jets(0);
    const xAOD::JetContainer* jets = 0;
    HLT::ErrorCode ec = getFeature(inputTE, jets, m_jetKey);
    if(ec!=HLT::OK) 
    {
        ATH_MSG_WARNING("Failed to get JetCollection");
        return ec;
    } 
    
    ATH_MSG_DEBUG("Obtained JetContainer");
    ATH_MSG_DEBUG("pass 2 " << &jets);
  
    if(jets == 0)
    {
        ATH_MSG_WARNING("Jet collection pointer is 0");
        return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::UNKNOWN);;
    }
 
    std::vector<const xAOD::Jet*> theJets(jets->begin(), jets->end());
    std::size_t njets = theJets.size();
    ATH_MSG_DEBUG("pass 2| jet size: " << njets);
  
    if( njets == 0 )
    {
        ATH_MSG_DEBUG("JetCollection is empty");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }
    ATH_MSG_DEBUG("JetCollection contains " << njets <<"jets");
  
    if(njets > 1)
    {
        ATH_MSG_DEBUG("Something is wrong, it should not be more than one jet");
	    return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::NAV_ERROR);
    }
  
    if(msgLvl() <= MSG::DEBUG)
    {
        for(const xAOD::Jet* aJet : theJets) 
        { 
            static int i=0;   
            ATH_MSG_DEBUG("  ["<<i<<"-Jet]: Et=" << aJet->p4().Et() << "(GeV), Eta=" << aJet->p4().Eta());
        }
    }
    return HLT::OK;
}
 


//** ----------------------------------------------------------------------------------------------------------------- **//

// CAVEAT!! Be careful with the use of inputTE which maybye doesn't have all the features which is suppose
// to have. The ouputTE is going to contain, in some cases, references to objects related to just processed
// RoI (in the previous sequences, which are actually executing this sequence...). So, if you find problems
// just use outputTE, which is safer (although inputTE is slightly faster, marginaly in fact)
HLT::ErrorCode TrigDvFex::hltExecute(const HLT::TriggerElement* inputTE, HLT::TriggerElement* outputTE) 
{
    // --> FIXME
    // WORKFLOW: 
    //          1. Extract ROI
    //          2. Extract Secondary Vertices from algorithm
    //               |--- InDet::TrigVxSecondaryCombo/TrigVxSecondaryCombo_Bjet_EF  (Package: InDetTrigVxSecondary) 
    //             being the inputs TE 'HLT_j55_eta_jsplit_EFID' and 'HLT_super_EFID_prmVtxCombo' and the output
    //             [HLT_j55_eta_jsplit_EFID__superVtx]   
    //             Need to know the Secondary vertices input, etc.. (I know there is also a secondary vertex info)
    //          3. NO NEED From PV, either Tracks, either anything else!!!
    //
    ATH_MSG_DEBUG("Executing TrigDvFex");

    // Clear and initialize data members
    m_totSelTracks = 0;
    m_totTracks    = 0;
  
    m_trigBjetPrmVtxInfo->clear();
    m_trigBjetSecVtxInfo->clear();
    m_trigBjetJetInfo->clear();
    // Track info
    std::vector<TrigBjetTrackInfo> trigBjetTrackInfoVector;
    m_trigBjetTrackInfoVector = &trigBjetTrackInfoVector;

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
    // Get EF jets (Do I need them?, yes when we work with the trigger DV+object)
    // But, let's change the order, firts tracks: this should be a function depending
    // of the object muon, electron, MET, JET
    // -----------------------------------
    float m_et_EFjet = 0;
    if(m_instance == "EF") 
    {
        std::vector<const TrigOperationalInfo*> m_vectorOperationalInfo;
        if(getFeatures(inputTE, m_vectorOperationalInfo, "EFJetInfo") != HLT::OK) 
        {
            ATH_MSG_WARNING("Failed to get TrigOperationalInfo");
            return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
        }
        else 
        {
            ATH_MSG_DEBUG("Number of TrigOperationalInfo objects: " << m_vectorOperationalInfo.size());
        }

        // -----------------------------------
        // Get operational info
        // -----------------------------------
        std::vector<const TrigOperationalInfo*>::const_iterator m_operationalInfo;
        for(m_operationalInfo=m_vectorOperationalInfo.begin(); 
    	    m_operationalInfo!=m_vectorOperationalInfo.end(); ++m_operationalInfo) 
        {
            if( (*m_operationalInfo)->defined("EFJetEt")==1 ) 
            {
                unsigned int m_etSize = (*m_operationalInfo)->infos().first.size();
                if(m_etSize!=1) 
                {
                    ATH_MSG_WARNING("More than one Et threshold associated to the same EF jet");
                    return HLT::ErrorCode(HLT::Action::ABORT_CHAIN, HLT::Reason::NAV_ERROR);
                }
                m_et_EFjet = (*m_operationalInfo)->get("EFJetEt");
            }
        }
    }

    // Set properties of the EF-jet 
    m_trigBjetJetInfo->setEtaPhiJet(roiDescriptor->eta(), m_taggerHelper->phiCorr(roiDescriptor->phi()));
    m_trigBjetJetInfo->setEtaPhiRoI(roiDescriptor->eta(), m_taggerHelper->phiCorr(roiDescriptor->phi()));
    m_trigBjetJetInfo->setEtJet(m_et_EFjet);
  
    //xAOD jets from TE
    ATH_MSG_DEBUG( "pass 1 " << m_et_EFjet);
    HLT::ErrorCode xaodc = checkxAODJets(inputTE);
    if( xaodc != HLT::OK )
    {
        return xaodc;
    }	
  
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Ready to get the tracks & vertices
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // -----------------------------------
    // Retrieve beamspot information
    // -----------------------------------
    IBeamCondSvc* m_iBeamCondSvc = 0;
    StatusCode sc = service("BeamCondSvc", m_iBeamCondSvc);
  
    if(sc.isFailure() || m_iBeamCondSvc == 0) 
    {
    	m_iBeamCondSvc = 0;
        ATH_MSG_WARNING("Could not retrieve Beam Conditions Service. ");
    } 
    else 
    {
     	Amg::Vector3D m_beamSpot = m_iBeamCondSvc->beamPos();
     	int m_beamSpotBitMap = m_iBeamCondSvc->beamStatus();
     	// Check if beam spot is from online algorithms
        int m_beamSpotStatus = ((m_beamSpotBitMap & 0x4) == 0x4);
     	// Check if beam spot fit converged
	    if(m_beamSpotStatus)
        {
            // and update
            m_beamSpotStatus = ((m_beamSpotBitMap & 0x3) == 0x3);
        }
        ATH_MSG_DEBUG("Beam spot from service: x=" << m_beamSpot.x() 
                << ", y=" << m_beamSpot.y() << ", z=" << m_beamSpot.z() 
                << ", tiltXZ=" << m_iBeamCondSvc->beamTilt(0) << ", tiltYZ=" 
                << m_iBeamCondSvc->beamTilt(1) << ", sigmaX=" 
                << m_iBeamCondSvc->beamSigma(0) << ", sigmaY=" 
                << m_iBeamCondSvc->beamSigma(1) << ", sigmaZ=" 
                << m_iBeamCondSvc->beamSigma(2) << ", status=" << m_beamSpotStatus);
	    
        // Update Primary vertex info class with beam spot position, tilt and width
        m_trigBjetPrmVtxInfo->setBeamSpot(m_beamSpot.x(), m_beamSpot.y(), m_beamSpot.z());
        m_trigBjetPrmVtxInfo->setBeamSpotTilt(m_iBeamCondSvc->beamTilt(0), m_iBeamCondSvc->beamTilt(1));    
        m_trigBjetPrmVtxInfo->setBeamSpotWidth(m_iBeamCondSvc->beamSigma(0), 
                m_iBeamCondSvc->beamSigma(1), m_iBeamCondSvc->beamSigma(2));
        m_trigBjetPrmVtxInfo->setBeamSpotStatus(m_beamSpotStatus);

        ATH_MSG_DEBUG(*m_trigBjetPrmVtxInfo);
    }
   
    // -----------------------------------
    // Create collections
    // -----------------------------------
    m_trigEFBjetColl = new TrigEFBjetContainer();
   
    // Create pointers to collections
    const xAOD::TrackParticleContainer * pointerToEFTrackCollections = 0;

    // Create pointers to TrigVertex collection
    const xAOD::VertexContainer *          pointerToEFPrmVtxCollections = 0;
    const Trk::VxSecVertexInfoContainer *  pointerToEFSecVtxCollections = 0;

    // -----------------------------------
    // Get EF track collection 
    // -----------------------------------
    HLT::ErrorCode status = getFeature(outputTE, pointerToEFTrackCollections);
    if(status != HLT::OK) 
    {
     	ATH_MSG_DEBUG("No HLT track collection retrieved");
        return status; // Should I return?
    } 
    ATH_MSG_DEBUG("HLT track collection retrieved");
    
    // -----------------------------------
    // Get secondary vertex collection 
    // -----------------------------------
    HLT::ErrorCode retSVstatus = getSecVtxCollection(pointerToEFSecVtxCollections, inputTE);
    if(retSVstatus != HLT::OK)
    {
        ATH_MSG_DEBUG("No secondary vertex collection retrieved");
        return retSVstatus;
    } 
    ATH_MSG_DEBUG("Secondary vertex collection retrieved");
  
    // -----------------------------------
    // Get primary vertex collection
    // -----------------------------------
    float m_xPrmVtx=0, m_yPrmVtx=0, m_zPrmVtx=0;
    std::string vtxlabel;
    if(m_histoPrmVtxAtEF)     // PV from TrigT2HistoPrmVtx
    {
        vtxlabel=m_priVtxKey;
    }
    else                      // PV from ID tracking
    {
        vtxlabel="";
    }
    // retrieve the vtx collection
    HLT::ErrorCode retPVstatus = getPrmVtxCollection(pointerToEFPrmVtxCollections, inputTE,vtxlabel);
    if(retPVstatus != HLT::OK)
    {
        ATH_MSG_DEBUG("No primary vertex collection retrieved");
	    //return retPVstatus;
    }
    else
    {
        ATH_MSG_DEBUG("Primary vertex collection retrieved");
    }

    // Aux var: PV with higher sum_{tracks} pt^2
    const xAOD::Vertex *pvselected = 0;
    // Protect against null pointers
    if(pointerToEFPrmVtxCollections) 
    {
        // Protect against empty vectors
	    if(pointerToEFPrmVtxCollections->size()==1) 
	    {
            pvselected = ((*pointerToEFPrmVtxCollections)[0]);
        }
	    else if(pointerToEFPrmVtxCollections->size() > 1)
        {
            // Aux map to order
	        std::map<float,const xAOD::Vertex *> auxPVOrder;
	        // Get as PV, the one with higher sum pt^2 of the tracks
	        for(const xAOD::Vertex * prmVertex : *pointerToEFPrmVtxCollections)
            {
                const size_t ntracks = prmVertex->nTrackParticles(); 
                float sumpt2 = 0.0;
		        for( size_t i = 0; i < ntracks; ++i)
                {
                    const double trackpt= prmVertex->trackParticle(i)->pt();
                    sumpt2 += (trackpt*trackpt);
                }
                auxPVOrder[sumpt2] = prmVertex;
            }
            // Higher sum_{tracks} pt^2
            pvselected = auxPVOrder.rbegin()->second;
        }
        else
        {
            ATH_MSG_DEBUG("Empty primary vertex collection retrieved");       
        }
    }
    else   
    {
        ATH_MSG_DEBUG("No primary vertex collection retrieved");
    }
    
    if(pvselected)
    {
        m_zPrmVtx = pvselected->z();
        if(!m_histoPrmVtxAtEF)
        {
            m_xPrmVtx = (float)(pvselected->x());
            m_yPrmVtx = (float)(pvselected->y());
        }
    }
    
    // -----------------------------------
    // Apply beam spot correction for tilt
    // ----------------------------------- 
    float m_xBeamSpot = m_trigBjetPrmVtxInfo->xBeamSpot() + 
        tan(m_trigBjetPrmVtxInfo->xBeamSpotTilt()) * (m_zPrmVtx-m_trigBjetPrmVtxInfo->zBeamSpot());
    
    float m_yBeamSpot = m_trigBjetPrmVtxInfo->yBeamSpot() + 
        tan(m_trigBjetPrmVtxInfo->yBeamSpotTilt()) * (m_zPrmVtx-m_trigBjetPrmVtxInfo->zBeamSpot());
    
    float m_zBeamSpot = m_trigBjetPrmVtxInfo->zBeamSpot();
    
    m_trigBjetPrmVtxInfo->setBeamSpot(m_xBeamSpot, m_yBeamSpot, m_zBeamSpot);  // Redundant?
    m_trigBjetPrmVtxInfo->setPrmVtx(m_xPrmVtx, m_yPrmVtx, m_zPrmVtx);
    
    ATH_MSG_DEBUG(*m_trigBjetPrmVtxInfo);
    
    m_trackJetFinderTool->clear();
    ATH_MSG_DEBUG("Set input  z-vtx to trackjet tool " << m_trigBjetPrmVtxInfo->zPrmVtx());
    m_trackJetFinderTool->inputPrimaryVertexZ(m_trigBjetPrmVtxInfo->zPrmVtx());
    ATH_MSG_DEBUG("Done set input  z-vtx to trackjet tool " << m_trigBjetPrmVtxInfo->zPrmVtx());
    
    // Get number of reconstructed tracks in this RoI
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
            ATH_MSG_DEBUG("  " << trigBjetTrackInfo);
            trigBjetTrackInfoVector.push_back(trigBjetTrackInfo);
            //m_trackJetFinderTool->addTrack(track, j);
      	}
    }
    
    std::vector<int> tracksTrackJet;
    float etaTrackJet, phiTrackJet;
    
    m_trackJetFinderTool->findJet(tracksTrackJet, etaTrackJet, phiTrackJet);
    if(etaTrackJet != -99 && phiTrackJet != -99) 
    {
     	m_trigBjetJetInfo->setEtaPhiTrkJet(etaTrackJet, phiTrackJet);
    } 
    else 
    {
     	m_trigBjetJetInfo->setEtaPhiTrkJet(m_trigBjetJetInfo->etaRoI(), m_trigBjetJetInfo->phiRoI());
        ATH_MSG_DEBUG("eta Jet = eta RoI");
    }
    ATH_MSG_DEBUG(*m_trigBjetJetInfo);
  
    // -----------------------------------
    // For monitoring
    // -----------------------------------
    m_deltaEtaJet       = m_trigBjetJetInfo->etaRoI()-m_trigBjetJetInfo->etaJet();
    m_deltaPhiJet       = m_trigBjetJetInfo->phiRoI()-m_trigBjetJetInfo->phiJet();
    m_deltaEtaTrkJet    = m_trigBjetJetInfo->etaRoI()-m_trigBjetJetInfo->etaTrkJet();
    m_deltaPhiTrkJet    = m_trigBjetJetInfo->phiRoI()-m_trigBjetJetInfo->phiTrkJet();
    m_deltaEtaJetTrkJet = m_trigBjetJetInfo->etaJet()-m_trigBjetJetInfo->etaTrkJet();
    m_deltaPhiJetTrkJet = m_trigBjetJetInfo->phiJet()-m_trigBjetJetInfo->phiTrkJet();
    
    // Get secondary vertex information at EF --> Do this first!!
    HLT::ErrorCode statSVInfo = setSecVtxInfo(pointerToEFSecVtxCollections,pvselected);

    // Some info to print if it is in the proper message level
    if(msgLvl() <= MSG::DEBUG) 
    {
      	const EventInfo* pEventInfo = 0;
     	if( !store() || store()->retrieve(pEventInfo).isFailure() ) 
        {
            ATH_MSG_DEBUG("Failed to get EventInfo ");
        } 
        else 
        {
            ATH_MSG_DEBUG("DV summary (Run " << pEventInfo->event_ID()->run_number() 
                    << "; Event " << pEventInfo->event_ID()->event_number() << ")");
            ATH_MSG_DEBUG("REGTEST:  RoI " << roiDescriptor->roiId() << ", Phi = "   
                    << roiDescriptor->phi() << ", Eta = "   << roiDescriptor->eta());
            ATH_MSG_DEBUG("REGTEST:  Tracks: " << m_totTracks << " reconstructed and " 
                    << m_totSelTracks <<" selected");

            if(pointerToEFPrmVtxCollections) 
            {
                ATH_MSG_DEBUG("REGTEST:  Primary vertex: " 
                        << pointerToEFPrmVtxCollections->size() << " reconstructed"
                        << ", (x,y,z) = (" << m_trigBjetPrmVtxInfo->xPrmVtx() << "," 
                        << m_trigBjetPrmVtxInfo->yPrmVtx() << "," 
                        << m_trigBjetPrmVtxInfo->zPrmVtx() << ")");
            }
            if(pointerToEFSecVtxCollections)
            {
                ATH_MSG_DEBUG("REGTEST:  Secondary vertex: " 
                        << pointerToEFSecVtxCollections->size() << " reconstructed");
            }
            else	    
            {
                ATH_MSG_DEBUG("REGTEST:  Secondary vertex: 0 reconstructed");
            }
            ATH_MSG_DEBUG("REGTEST:  SV Decay Length: " 
                    <<  m_trigBjetSecVtxInfo->decayLengthSignificance()/CLHEP::mm << " mm, SV mass: " 
                    <<  m_trigBjetSecVtxInfo->vtxMass()/CLHEP::GeV << " GeV, SV efrac: " 
                    <<  m_trigBjetSecVtxInfo->energyFraction()/CLHEP::GeV
                    << " GeV, SV 2-track vertex multiplicity " << m_trigBjetSecVtxInfo->n2TrkVtx());
            ATH_MSG_DEBUG("REGTEST: List weights stored probability and likelihood objects:");
        }
    }
    
    // Attach a Primary Vertex and Secondary Vertex ... that'all we need
    // -----------------------------------
    // -----------------------------------
    // Create TrigEFBjet and attach feature
    // -----------------------------------
    // Note that the meaning of some characteristics are tuned to 
    // DV case: * Put a PV (a coll with the PV selected)
    TrigEFBjet* trigEFBjet = new TrigEFBjet(roiDescriptor->roiId(), 
	    m_trigBjetJetInfo->etaJet(), m_trigBjetJetInfo->phiJet(),
	    0, 0, 0, m_trigBjetPrmVtxInfo->zPrmVtx(), m_trigBjetJetInfo->etJet(),
	    -1,-1,-1, -1,-1,  // Note, i can use this to fill other stuff if I need 
	    m_trigBjetSecVtxInfo->decayLengthSignificance(), m_trigBjetSecVtxInfo->vtxMass(), 
	    m_trigBjetSecVtxInfo->energyFraction(), m_trigBjetSecVtxInfo->n2TrkVtx()); 
    
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

    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::hltFinalize() 
{
    ATH_MSG_INFO("Finalizing TrigDvFex");
    return HLT::OK;
}



