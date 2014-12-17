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
  // Note that "T2PrimaryVertex" belongs to "EFHistoPrmVtx" but has been used 
  // by the TrigVxSecondayrCombo algorithm to do all the calculation related the SV.
  // Also, is chosen a PV with better fit quality if there is more than on PV, altough
  // the criteria should be the highest sum_{track} pt^2...
  declareProperty ("PriVtxKey",          m_priVtxKey  = "T2PrimaryVertex");

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
	const HLT::TriggerElement* whateverTE, const std::string & vtxkey, const bool & ispvfromsvalg) 
{
    std::vector<const xAOD::VertexContainer*> vectorOfEFPrmVtxCollections;
    HLT::ErrorCode status = getFeatures(whateverTE, vectorOfEFPrmVtxCollections, vtxkey);
    if(status != HLT::OK) 
    {
        ATH_MSG_ERROR("Failed to get xAOD::VertexContainer from the trigger element");
        return status;
    } 
    ATH_MSG_DEBUG("Got " << vectorOfEFPrmVtxCollections.size() << " xAOD::VertexContainer");
    
    if( vectorOfEFPrmVtxCollections.size() > 1 )
    {
        ATH_MSG_ERROR("The vector of  xAOD::VertexContainer have more than 1 element!");
        return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::NAV_ERROR);
    }

    const xAOD::VertexContainer * vectorOfPv = vectorOfEFPrmVtxCollections[0];
    // Check if we have a unique valid collection (when extracted from the SV-builder algorithm
    if( ispvfromsvalg && vectorOfPv->size() != 1)
    {
        ATH_MSG_ERROR("The PV '" << vtxkey << "' coming from the SV-builder algorithm,"
               " have more than 1 element. But by construction, shouldn't!");
        return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::NAV_ERROR);
    }
       
    for( auto & pv : *vectorOfPv)
    {
        // This check is probably redundant, but nevertheless...
        if( pv->vertexType() != xAOD::VxType::VertexType::PriVtx )
        {
            continue;
        }
        ATH_MSG_DEBUG("Selected collection with Primary Vertex label ");
        ATH_MSG_DEBUG("First PV has z-position = " <<  pv->z());
        
        pointerToEFPrmVtxCollections = vectorOfPv;
        return HLT::OK;        
    }    
    // Didn't found any PV collection
    pointerToEFPrmVtxCollections = 0;
    ATH_MSG_DEBUG("Not found any Primary Vertex valid collection!");
    return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
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
    
    if( vectorOfSecVtxCollections.size() > 1 )
    {
        ATH_MSG_ERROR("The vector of VxSecVertexInfoContainer have more than 1 element!");
        return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::NAV_ERROR);
    }
    else if( vectorOfSecVtxCollections.size() < 1 )
    {
        ATH_MSG_ERROR("The vector of VxSecVertexInfoContainer have none element!");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }

    const Trk::VxSecVertexInfoContainer * vectorOfSvInfo = vectorOfSecVtxCollections[0];
    // Check if each VxSecVertexInfo contains a valid set of vertices, note that
    // Trk::VxSecVertexInfo is a base class with no useful methods except of the 
    // vertices one
    for( auto & svinfo : *vectorOfSvInfo)
    {
        const std::vector<xAOD::Vertex*> vectorOfSv = svinfo->vertices();
        if( vectorOfSv.size() == 0 )
        {
            continue;
        }
        // This check is probably redundant, but nevertheless...
        if( vectorOfSv[0]->vertexType() == xAOD::VxType::VertexType::SecVtx)
        {
            //So, the vectorOfSvInfo contains at least one collection 
            // of Secondary Vertices... storing and return
            ATH_MSG_DEBUG("Selected vertex-collection with Secondary Vertex label"); 
            ATH_MSG_DEBUG("The first SV has z-position = " << vectorOfSv[0]->z());
            pointerToSecVtxCollections = vectorOfSvInfo;
            return HLT::OK;
        }
    }
    // Didn't found any valid secondary vertex collection, the chain should abort?
    pointerToSecVtxCollections = 0;
    ATH_MSG_DEBUG("Not found any Secondary Vertex valid collection!");
    return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
}


//** ----------------------------------------------------------------------------------------------------------------- **//

// This method should use a new data-member called m_svcollection, in order to be sure we properly
// set this collection in the getSecVtx or maybe this method should be absorbed in the getSecVtxCollection
// method?? This will avoid a lot of cross-checsk and ifs...
HLT::ErrorCode TrigDvFex::setSecVtxInfo(const Trk::VxSecVertexInfoContainer*& pointerToEFSecVtxCollections, 
					  const xAOD::Vertex* & pvselected) 
{
    // Redundant... but 
    if(!pvselected)
    {
        ATH_MSG_DEBUG("No primary vertex collection sent when extracting sec vtx info");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }

    // Also redundant..
    if(!pointerToEFSecVtxCollections) 
    {
        ATH_MSG_DEBUG("No secondary vertex collection sent when extracting sec vtx info");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }
  
    const Trk::VxSecVertexInfo* m_secVertexInfo = (*pointerToEFSecVtxCollections)[0];
    // Also redundant..
    if(!m_secVertexInfo) 
    { 
        ATH_MSG_DEBUG("No secondary vertex when extracting sec vtx info");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }
  
    // Downcast from the base Trk::VxSecVertexInfo 
    const Trk::VxSecVKalVertexInfo * myVKalSecVertexInfo = dynamic_cast<const Trk::VxSecVKalVertexInfo*>(m_secVertexInfo);
 
    if(!myVKalSecVertexInfo) 
    {
        ATH_MSG_DEBUG("The cast to VKal secondary vertex went wrong, the pointer is zero");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }
  
    const std::vector<xAOD::Vertex*> & myVertices = myVKalSecVertexInfo->vertices();
   
    if(myVertices.size() == 0)
    {
        ATH_MSG_DEBUG("The vertices size associated to the myVKalSecVertexInfo is 0!" 
                << " Content of the myVKalSecVertexInfo " << std::endl 
                << "   - mass: " << myVKalSecVertexInfo->mass()/CLHEP::GeV << " [GeV] " << std::endl
                << "   - energyFraction: " << myVKalSecVertexInfo->energyFraction() << " " << std::endl
                << "   - energyTrkInJet: " << myVKalSecVertexInfo->energyTrkInJet() << " " << std::endl
                << "   - number of 2Tracks-vertices: " << myVKalSecVertexInfo->n2trackvertices() );
        return HLT::OK;
    }
 
    m_trigBjetSecVtxInfo->setVtxMass(myVKalSecVertexInfo->mass());
    m_trigBjetSecVtxInfo->setEnergyFraction(myVKalSecVertexInfo->energyFraction());
    m_trigBjetSecVtxInfo->setN2TrkVtx(myVKalSecVertexInfo->n2trackvertices());
    
    int NTracksInSV=0;
    
    if(myVertices.size()>1) 
    {
        ATH_MSG_WARNING("Secondary vertex from InDetVKalVxInJet has more than one vertex,"
                << " they have been merged into just one (same mass, energies, etc..");
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
    ATH_MSG_DEBUG("Primary vertex for decay length (" 
    	<< pvselected->position().x() << "," 
    	<< pvselected->position().y() << "," 
    	<< pvselected->position().z() << ") and error (" 
    	<< pvselected->covariancePosition()(0,0) << "," 
    	<< pvselected->covariancePosition()(1,1) << "," 
    	<< pvselected->covariancePosition()(2,2) << ")");
    //Toggle
    m_trigBjetSecVtxInfo->isValid(true);

    // Store the SV with higher mass in the collection
    // Getting the first vertex?? --> or a mean/weigthed mean? maybe
    /*float x = 0.;
    float y = 0.;
    float z = 0.;
    const unsigned int nv = myVertices.size();
    for(auto & sv: myVertices)
    {
        x += sv->position->x();
        y += sv->position->y();
        z += sv->position->z();
    }
    const Amg::Vector3D SecVrt(x/float(nv),y/float(nv),z/float(nv));*/

    const Amg::Vector3D SecVrt = myVertices[0]->position();

    double distance = (SecVrt - pvselected->position()).mag();
    // Note, storing the distance 3D between SV and PV, watch out with the name 
    // of the method!!
    m_trigBjetSecVtxInfo->setDecayLengthSignificance(distance);
    if(msgLvl() <= MSG::DEBUG) 
    {
        double dist2D = (SecVrt - pvselected->position()).perp();            
	    ATH_MSG_DEBUG("Calculated secondary vertex decay length with primary vertex at (" 
                << pvselected->position().x()/CLHEP::mm << "," << pvselected->position().y()/CLHEP::mm
                << "," << pvselected->position().z()/CLHEP::mm << ") and sec. vtx at ("
                << SecVrt.x()/CLHEP::mm << "," << SecVrt.y()/CLHEP::mm << "," << SecVrt.z()/CLHEP::mm 
                <<  ")");
        ATH_MSG_DEBUG("    * 3D decay length: " << distance/CLHEP::mm << " [mm]"); 
        ATH_MSG_DEBUG("    * 2D(R/phi) decay length: " << dist2D/CLHEP::mm << " [mm]");
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
 

//** ----------------------------------------------------------------------------------------------------------------- **//
HLT::ErrorCode TrigDvFex::setBeamSpotRelated()
{
    // This function should be called after filling the Primary vertex, so check it
    if(!m_setPVInfo)
    {
        ATH_MSG_ERROR("This function cannot be called before setting the" 
                << " Primary Vertex info (TrigPrmVtxInfo class) and informed"
                << " by the private 'm_setPVInfo'. This is a logical error which"
                << " should be fixed. Check the TrigDvFex::hltExecute method.");
	    return HLT::ErrorCode(HLT::Action::ABORT_JOB,HLT::Reason::UNKNOWN);
    }

    IBeamCondSvc* iBeamCondSvc = 0;
    StatusCode sc = service("BeamCondSvc", iBeamCondSvc);
  
    if(sc.isFailure() || iBeamCondSvc == 0) 
    {
    	iBeamCondSvc = 0;
        ATH_MSG_WARNING("Could not retrieve Beam Conditions Service. ");
        return HLT::OK;
    } 
     	
    Amg::Vector3D beamSpot = iBeamCondSvc->beamPos();
    int beamSpotBitMap = iBeamCondSvc->beamStatus();
    // Check if beam spot is from online algorithms
    int beamSpotStatus = ((beamSpotBitMap & 0x4) == 0x4);
    // Check if beam spot fit converged
    if(beamSpotStatus)
    {
        // and update
        beamSpotStatus = ((beamSpotBitMap & 0x3) == 0x3);
    }
    ATH_MSG_DEBUG("Beam spot from service: x=" << beamSpot.x()
            << ", y=" << beamSpot.y() << ", z=" << beamSpot.z() 
            << ", tiltXZ=" << iBeamCondSvc->beamTilt(0) << ", tiltYZ=" 
            << iBeamCondSvc->beamTilt(1) << ", sigmaX=" 
            << iBeamCondSvc->beamSigma(0) << ", sigmaY=" 
            << iBeamCondSvc->beamSigma(1) << ", sigmaZ=" 
            << iBeamCondSvc->beamSigma(2) << ", status=" << beamSpotStatus);

    // Update Primary vertex info class with beam spot position, tilt and width
    m_trigBjetPrmVtxInfo->setBeamSpotTilt(iBeamCondSvc->beamTilt(0), iBeamCondSvc->beamTilt(1));    
    m_trigBjetPrmVtxInfo->setBeamSpotWidth(iBeamCondSvc->beamSigma(0), 
            iBeamCondSvc->beamSigma(1), iBeamCondSvc->beamSigma(2));
    m_trigBjetPrmVtxInfo->setBeamSpotStatus(beamSpotStatus);
    
    ATH_MSG_DEBUG(*m_trigBjetPrmVtxInfo);
    
    // -----------------------------------
    // Apply beam spot correction for tilt
    // ----------------------------------- 
    const float xBeamSpot_corr = beamSpot.x() +
        tan(m_trigBjetPrmVtxInfo->xBeamSpotTilt())*(m_trigBjetPrmVtxInfo->zPrmVtx()-beamSpot.z());
    
    const float yBeamSpot_corr = beamSpot.y() +
        tan(m_trigBjetPrmVtxInfo->yBeamSpotTilt())*(m_trigBjetPrmVtxInfo->zPrmVtx()-beamSpot.z());
    
    const float zBeamSpot_corr = beamSpot.z();
   
    m_trigBjetPrmVtxInfo->setBeamSpot(xBeamSpot_corr, yBeamSpot_corr, zBeamSpot_corr);
    // -- Primary vertex info class filled!
    
    ATH_MSG_DEBUG(*m_trigBjetPrmVtxInfo);
    return HLT::OK;
    
    /*m_trackJetFinderTool->clear();
    ATH_MSG_DEBUG("Set input  z-vtx to trackjet tool " << m_trigBjetPrmVtxInfo->zPrmVtx());
    m_trackJetFinderTool->inputPrimaryVertexZ(m_trigBjetPrmVtxInfo->zPrmVtx());
    ATH_MSG_DEBUG("Done set input  z-vtx to trackjet tool " << m_trigBjetPrmVtxInfo->zPrmVtx());
    
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
                ATH_MSG_DEBUG("  " << trigBjetTrackInfo);
                trigBjetTrackInfoVector.push_back(trigBjetTrackInfo);
      	    }
        }
    }
    else
    {
        m_totTracks = 0;
    }*/
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
    m_setPVInfo    = false;
  
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
    // Get secondary vertex collection 
    // -----------------------------------
    const Trk::VxSecVertexInfoContainer *  pointerToEFSecVtxCollections = 0;

    HLT::ErrorCode retsvstatus = getSecVtxCollection(pointerToEFSecVtxCollections, inputTE);
    if(retsvstatus != HLT::OK)
    {
        ATH_MSG_DEBUG("No sec:ondary vertex collection retrieved... Stopping execution");
        return retsvstatus;
    } 
    ATH_MSG_DEBUG("Secondary vertex collection retrieved");

    // -----------------------------------
    // Get primary vertex collection
    // -----------------------------------
    const xAOD::VertexContainer *   pointerToEFPrmVtxCollections = 0;
    
    std::string vtxlabel;
    bool ispvfromsvalgo=true;
    if(m_histoPrmVtxAtEF)     // PV from TrigT2HistoPrmVtx
    {
        vtxlabel=m_priVtxKey;
    }
    else                      // PV from ID tracking
    {
        vtxlabel="";
        ispvfromsvalgo=false;
    }
    // retrieve the vtx collection
    HLT::ErrorCode retpvstatus = getPrmVtxCollection(pointerToEFPrmVtxCollections, inputTE,vtxlabel,ispvfromsvalgo);
    if(retpvstatus != HLT::OK)
    {
        ATH_MSG_DEBUG("No primary vertex collection retrieved, this is an indication" 
                << " that there is no Secondary vertex collection either. Exiting...");
	    return retpvstatus;
    }
    ATH_MSG_DEBUG("Primary vertex collection retrieved");

    // Get the PV. If the PV collection is coming from the secondary vertex builder
    // algorithm, there is just one selected PV. Otherwise, it could be that there
    // is more than one, so the used criteria is taking the PV with 
    // highest sum_{tracks} pt^2
    const xAOD::Vertex *pvselected = 0;
    if(pointerToEFPrmVtxCollections)
    {
        const size_t npvc = pointerToEFPrmVtxCollections->size();
        if(npvc == 1)
        {
            pvselected = (*pointerToEFPrmVtxCollections)[0];
        }
        else if(npvc > 1)
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
            // Higher sum_{tracks} pt^2 (remember, the highest is the last one
            // in a std::map 
            pvselected = auxPVOrder.rbegin()->second;
        }
        else
        {
            ATH_MSG_DEBUG("Empty primary vertex collection retrieved. Exiting execution...");
            return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
        }
    }
    else   
    {
        ATH_MSG_DEBUG("No primary vertex collection retrieved. Exiting execution...");
        return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE);
    }

    // Filling the relevant info of the Secondary vertex
    retpvstatus = setSecVtxInfo(pointerToEFSecVtxCollections,pvselected);
    if(retpvstatus != HLT::OK)
    {
        ATH_MSG_DEBUG("Problems filling the Secondary Vertex related info"); 
	    return retpvstatus;
    }
    
    // Filling the Auxiliar (relevant) helper classes
    // --- PrmVtxInfo
    m_trigBjetPrmVtxInfo->setPrmVtx(pvselected->x(),pvselected->y(),pvselected->z());
    m_setPVInfo=true;
    // Updated the pv info with beamspot
    const HLT::ErrorCode statusbs = setBeamSpotRelated();
    if( statusbs != HLT::OK )
    {
        return statusbs;
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
                    << " GeV, SV 2-track vertex multiplicity: " 
                    << m_trigBjetSecVtxInfo->n2TrkVtx() 
                    << " , SV tracks multiplicity: " << m_trigBjetSecVtxInfo->nTrksInVtx());
        }
    }

    return HLT::OK;
    /*
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Ready to get the tracks & vertices
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // -----------------------------------
    // Retrieve beamspot information
    // -----------------------------------
    // NOT NEEDED!!! Right?
   
    // Create pointers to collections
    const xAOD::TrackParticleContainer * pointerToEFTrackCollections = 0;

    // -----------------------------------
    // Get EF track collection 
    // -----------------------------------
    HLT::ErrorCode status = getFeature(outputTE, pointerToEFTrackCollections);
    if(status != HLT::OK) 
    {
     	ATH_MSG_DEBUG("No HLT track collection retrieved");
        //return HLT::ErrorCode(HLT::Action::CONTINUE,HLT::Reason::MISSING_FEATURE); // Should I return?
    } 
    else
    {
        ATH_MSG_DEBUG("HLT track collection retrieved");
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
                ATH_MSG_DEBUG("  " << trigBjetTrackInfo);
                trigBjetTrackInfoVector.push_back(trigBjetTrackInfo);
      	    }
        }
    }
    else
    {
        m_totTracks = 0;
    }
    
    //std::vector<int> tracksTrackJet;
    //float etaTrackJet, phiTrackJet;
    
    // Not using it...
    //m_trackJetFinderTool->findJet(tracksTrackJet, etaTrackJet, phiTrackJet);
    //if(etaTrackJet != -99 && phiTrackJet != -99) 
    //{
    // 	m_trigBjetJetInfo->setEtaPhiTrkJet(etaTrackJet, phiTrackJet);
    //} 
    //else 
    //{
    m_trigBjetJetInfo->setEtaPhiTrkJet(m_trigBjetJetInfo->etaRoI(), m_trigBjetJetInfo->phiRoI());
    ATH_MSG_DEBUG("eta Jet = eta RoI");
    //}
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

    return HLT::OK;*/
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigDvFex::hltFinalize() 
{
    ATH_MSG_INFO("Finalizing TrigDvFex");
    return HLT::OK;
}



