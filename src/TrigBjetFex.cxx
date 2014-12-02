// ************************************************
//
// NAME:     TrigBjetFex.cxx
// PACKAGE:  Trigger/TrigHypothesis/TrigBjetHypo
// 
// ************************************************

#include "TrigBjetHypo/TrigBjetFex.h"
#include "TrigBjetHypo/TrigBjetTagger.h"
#include "TrigBjetHypo/TuningLikelihood.h"

#include "TaggerHelper.h"

#include "TrigInDetEvent/TrigInDetTrackCollection.h"
#include "Particle/TrackParticleContainer.h"
#include "TrigInDetEvent/TrigVertexCollection.h"
#include "VxVertex/VxContainer.h"
#include "VxSecVertex/VxSecVertexInfo.h"
#include "VxSecVertex/VxSecVKalVertexInfo.h"
#include "TrigCaloEvent/TrigT2Jet.h"
#include "TrigSteeringEvent/TrigRoiDescriptor.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"
#include "EventPrimitives/EventPrimitives.h"
#include "EventPrimitives/EventPrimitivesHelpers.h"

#include "TrigParticle/TrigEFBjetContainer.h"

#include "TrigSteeringEvent/TrigOperationalInfo.h"

#include "InDetBeamSpotService/IBeamCondSvc.h"

#include "TrigInDetEvent/TrigVertex.h"

#include "xAODBTagging/BTaggingContainer.h"
#include "xAODBTagging/BTagging.h"
#include "xAODBTagging/BTaggingAuxContainer.h"

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/VertexAuxContainer.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODBase/IParticle.h"

//** ----------------------------------------------------------------------------------------------------------------- **//


TrigBjetFex::TrigBjetFex(const std::string& name, ISvcLocator* pSvcLocator) :
  HLT::FexAlgo(name, pSvcLocator),
  m_trackJetFinderTool("TrigTrackJetFinderTool",this),      	
  m_trigEFBjetColl(0),
  m_trigBjetTagger(0),
  m_constTrigBjetTagger(0),
  m_trigBjetPrmVtxInfo(0),
  m_trigBjetSecVtxInfo(0),
  m_trigBjetJetInfo(0),
  m_totTracks(0),
  m_totSelTracks(0),
  m_tuningLikelihoodIP1D(0),
  m_tuningLikelihoodIP2D(0),
  m_tuningLikelihoodIP3D(0),
  m_tuningLikelihoodIP1D_lowSiHits(0),
  m_tuningLikelihoodIP2D_lowSiHits(0),
  m_tuningLikelihoodIP3D_lowSiHits(0),
  m_tuningLikelihoodMVtx(0),
  m_tuningLikelihoodEVtx(0),
  m_tuningLikelihoodNVtx(0),
  m_tuningLikelihoodSV(0)
{
  declareProperty ("AlgoId",             m_algo);
  declareProperty ("Instance",           m_instance);
  declareProperty ("Taggers",            m_taggers);
  declareProperty ("JetKey",             m_jetKey     = ""); //"" needed for default config, SplitJet for new config
  declareProperty ("PriVtxKey",          m_priVtxKey  = "EFHistoPrmVtx"); //Does this still work with default config?

  declareProperty ("par_0_MC",           m_par_0_MC);
  declareProperty ("par_1_MC",           m_par_1_MC);
  declareProperty ("par_0_DT",           m_par_0_DT);
  declareProperty ("par_1_DT",           m_par_1_DT);

  declareProperty ("SizeIP1D",           m_sizeIP1D);
  declareProperty ("bIP1D",              m_bIP1D);
  declareProperty ("uIP1D",              m_uIP1D);
  declareProperty ("SizeIP2D",           m_sizeIP2D);
  declareProperty ("bIP2D",              m_bIP2D);
  declareProperty ("uIP2D",              m_uIP2D);
  declareProperty ("SizeIP3D",           m_sizeIP3D);
  declareProperty ("bIP3D",              m_bIP3D);
  declareProperty ("uIP3D",              m_uIP3D);

  declareProperty ("useLowSiHits",                 m_useLowSiHits = false);
  declareProperty ("SizeIP1D_lowSiHits",           m_sizeIP1D_lowSiHits);
  declareProperty ("bIP1D_lowSiHits",              m_bIP1D_lowSiHits);
  declareProperty ("uIP1D_lowSiHits",              m_uIP1D_lowSiHits);
  declareProperty ("SizeIP2D_lowSiHits",           m_sizeIP2D_lowSiHits);
  declareProperty ("bIP2D_lowSiHits",              m_bIP2D_lowSiHits);
  declareProperty ("uIP2D_lowSiHits",              m_uIP2D_lowSiHits);
  declareProperty ("SizeIP3D_lowSiHits",           m_sizeIP3D_lowSiHits);
  declareProperty ("bIP3D_lowSiHits",              m_bIP3D_lowSiHits);
  declareProperty ("uIP3D_lowSiHits",              m_uIP3D_lowSiHits);

  declareProperty ("SizeMVtx",           m_sizeMVtx);
  declareProperty ("bMVtx",              m_bMVtx);
  declareProperty ("uMVtx",              m_uMVtx);
  declareProperty ("SizeEVtx",           m_sizeEVtx);
  declareProperty ("bEVtx",              m_bEVtx);
  declareProperty ("uEVtx",              m_uEVtx);
  declareProperty ("SizeNVtx",           m_sizeNVtx);
  declareProperty ("bNVtx",              m_bNVtx);
  declareProperty ("uNVtx",              m_uNVtx);
  declareProperty ("SizeSV",             m_sizeSV);
  declareProperty ("bSV",                m_bSV);
  declareProperty ("uSV",                m_uSV);

  declareProperty ("UseBeamSpotFlag",    m_useBeamSpotFlag    = false);
  declareProperty ("SetBeamSpotWidth",   m_setBeamSpotWidth   = 0.05);

  declareProperty ("UseParamFromData",   m_useParamFromData   = false);

  declareProperty ("UseErrIPParam",      m_useErrIPParam      = false);
  declareProperty ("HistoPrmVtxAtEF",    m_histoPrmVtxAtEF    = true);
  declareProperty ("UseEtaPhiTrackSel",  m_useEtaPhiTrackSel  = false);

  declareProperty ("UseJetDirection",    m_useJetDirection);
  declareProperty ("RetrieveHLTJets",    m_retrieveHLTJets    = true);
  declareProperty ("TagHLTJets",         m_tagHLTJets         = 0);

  declareProperty ("TrkSel_Chi2",        m_trkSelChi2         = 0.001);
  declareProperty ("TrkSel_BLayer",      m_trkSelBLayer       = 1);
  declareProperty ("TrkSel_PixHits",     m_trkSelPixHits      = 2);
  declareProperty ("TrkSel_SiHits",      m_trkSelSiHits       = 4);
  declareProperty ("TrkSel_D0",          m_trkSelD0           = 1*CLHEP::mm);
  declareProperty ("TrkSel_Z0",          m_trkSelZ0           = 2*CLHEP::mm);
  declareProperty ("TrkSel_Pt",          m_trkSelPt           = 1*CLHEP::GeV);

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
  declareMonitoredObject      ("roi_selectedTracks", *this, &TrigBjetFex::totSelectedTracks);

  declareMonitoredVariable    ("roi_deltaEtaJet",       m_deltaEtaJet,       AutoClear);
  declareMonitoredVariable    ("roi_deltaPhiJet",       m_deltaPhiJet,       AutoClear);
  declareMonitoredVariable    ("roi_deltaEtaTrkJet",    m_deltaEtaTrkJet,    AutoClear);
  declareMonitoredVariable    ("roi_deltaPhiTrkJet",    m_deltaPhiTrkJet,    AutoClear);
  declareMonitoredVariable    ("roi_deltaEtaJetTrkJet", m_deltaEtaJetTrkJet, AutoClear);
  declareMonitoredVariable    ("roi_deltaPhiJetTrkJet", m_deltaPhiJetTrkJet, AutoClear);

  declareMonitoredObject("X(IP1D)", m_constTrigBjetTagger, &TrigBjetTagger::getXIP1D);
  declareMonitoredObject("X(IP2D)", m_constTrigBjetTagger, &TrigBjetTagger::getXIP2D);
  declareMonitoredObject("X(IP3D)", m_constTrigBjetTagger, &TrigBjetTagger::getXIP3D);
  declareMonitoredObject("X(SVTX)", m_constTrigBjetTagger, &TrigBjetTagger::getXSVTX);
  declareMonitoredObject("X(COMB)", m_constTrigBjetTagger, &TrigBjetTagger::getXCOMB);
  declareMonitoredObject("X(CHI2)", m_constTrigBjetTagger, &TrigBjetTagger::getXCHI2);

  m_taggerHelper = new TaggerHelper(msg(), msgLvl());
}


//** ----------------------------------------------------------------------------------------------------------------- **//


TrigBjetFex::~TrigBjetFex() {

  if (m_taggerHelper)            delete m_taggerHelper;
  if (m_trigBjetTagger)          delete m_trigBjetTagger;
  if (m_trigBjetPrmVtxInfo)      delete m_trigBjetPrmVtxInfo;
  if (m_trigBjetSecVtxInfo)      delete m_trigBjetSecVtxInfo;
  if (m_trigBjetJetInfo)         delete m_trigBjetJetInfo;
  if (m_tuningLikelihoodIP1D)    delete m_tuningLikelihoodIP1D;
  if (m_tuningLikelihoodIP2D)    delete m_tuningLikelihoodIP2D;
  if (m_tuningLikelihoodIP3D)    delete m_tuningLikelihoodIP3D;
  if (m_tuningLikelihoodIP1D_lowSiHits)    delete m_tuningLikelihoodIP1D_lowSiHits;
  if (m_tuningLikelihoodIP2D_lowSiHits)    delete m_tuningLikelihoodIP2D_lowSiHits;
  if (m_tuningLikelihoodIP3D_lowSiHits)    delete m_tuningLikelihoodIP3D_lowSiHits;
  if (m_tuningLikelihoodMVtx)    delete m_tuningLikelihoodMVtx;
  if (m_tuningLikelihoodEVtx)    delete m_tuningLikelihoodEVtx;
  if (m_tuningLikelihoodNVtx)    delete m_tuningLikelihoodNVtx;
  if (m_tuningLikelihoodSV)      delete m_tuningLikelihoodSV;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigBjetFex::hltInitialize() {
  
  // Get message service
      if (msgLvl() <= MSG::INFO) 
      msg() << MSG::INFO << "Initializing TrigBjetFex, version " << PACKAGE_VERSION << endreq;

    // declareProperty overview
        if (msgLvl() <= MSG::DEBUG) {
          msg() << MSG::DEBUG << "declareProperty review:" << endreq;

	  msg() << MSG::DEBUG << "JetKey = "               << m_jetKey << endreq;
	  msg() << MSG::DEBUG << "PriVtxKey = "            << m_priVtxKey << endreq;

          msg() << MSG::DEBUG << " AlgoId = "              << m_algo << endreq; 
          msg() << MSG::DEBUG << " Instance = "            << m_instance << endreq;
 
          msg() << MSG::DEBUG << " UseBeamSpotFlag = "     << m_useBeamSpotFlag << endreq; 
          msg() << MSG::DEBUG << " SetBeamSpotWidth = "    << m_setBeamSpotWidth << endreq;

          msg() << MSG::DEBUG << " UseParamFromData = "    << m_useParamFromData << endreq; 

          msg() << MSG::DEBUG << " Taggers = "             << m_taggers << endreq; 
          msg() << MSG::DEBUG << " UseErrIPParam = "       << m_useErrIPParam << endreq; 
          msg() << MSG::DEBUG << " UseJetDirection = "     << m_useJetDirection << endreq; 
          msg() << MSG::DEBUG << " RetrieveHLTJets = "     << m_retrieveHLTJets << endreq; 
          msg() << MSG::DEBUG << " TagHLTJets = "          << m_tagHLTJets << endreq;
          msg() << MSG::DEBUG << " HistoPrmVtxAtEF = "     << m_histoPrmVtxAtEF << endreq;
          msg() << MSG::DEBUG << " UseEtaPhiTrackSel = "   << m_useEtaPhiTrackSel << endreq;

          msg() << MSG::DEBUG << " JetProb 0 MC = "      << m_par_0_MC << endreq; 
          msg() << MSG::DEBUG << " JetProb 1 MC = "      << m_par_1_MC << endreq; 
          msg() << MSG::DEBUG << " JetProb 0 DT = "      << m_par_0_DT << endreq; 
          msg() << MSG::DEBUG << " JetProb 1 DT = "      << m_par_1_DT << endreq; 

          msg() << MSG::DEBUG << " SizeIP1D = "          << m_sizeIP1D << endreq; 
          msg() << MSG::DEBUG << " bIP1D = "             << m_bIP1D << endreq; 
          msg() << MSG::DEBUG << " uIP1D = "             << m_uIP1D << endreq; 
          msg() << MSG::DEBUG << " SizeIP2D = "          << m_sizeIP2D << endreq;
          msg() << MSG::DEBUG << " bIP2D = "             << m_bIP2D << endreq; 
          msg() << MSG::DEBUG << " uIP2D = "             << m_uIP2D << endreq;  
          msg() << MSG::DEBUG << " SizeIP3D = "          << m_sizeIP3D << endreq; 
          msg() << MSG::DEBUG << " bIP3D = "             << m_bIP3D << endreq; 
          msg() << MSG::DEBUG << " uIP3D = "             << m_uIP3D << endreq; 

          msg() << MSG::DEBUG << " SizeIP1D_lowSiHits = "  << m_sizeIP1D_lowSiHits << endreq; 
          msg() << MSG::DEBUG << " bIP1D_lowSiHits = "     << m_bIP1D_lowSiHits << endreq; 
          msg() << MSG::DEBUG << " uIP1D_lowSiHits = "     << m_uIP1D_lowSiHits << endreq; 
          msg() << MSG::DEBUG << " SizeIP2D_lowSiHits = "  << m_sizeIP2D_lowSiHits << endreq;
          msg() << MSG::DEBUG << " bIP2D_lowSiHits = "     << m_bIP2D_lowSiHits << endreq; 
          msg() << MSG::DEBUG << " uIP2D_lowSiHits = "     << m_uIP2D_lowSiHits << endreq;  
          msg() << MSG::DEBUG << " SizeIP3D_lowSiHits = "  << m_sizeIP3D_lowSiHits << endreq; 
          msg() << MSG::DEBUG << " bIP3D_lowSiHits = "     << m_bIP3D_lowSiHits << endreq; 
          msg() << MSG::DEBUG << " uIP3D_lowSiHits = "     << m_uIP3D_lowSiHits << endreq; 

          msg() << MSG::DEBUG << " SizeIP1D = "          << m_sizeIP1D << endreq; 
          msg() << MSG::DEBUG << " bIP1D = "             << m_bIP1D << endreq; 
          msg() << MSG::DEBUG << " uIP1D = "             << m_uIP1D << endreq; 

          msg() << MSG::DEBUG << " TrkSel_Chi2 = "     << m_trkSelChi2 << endreq; 
          msg() << MSG::DEBUG << " TrkSel_BLayer = "   << m_trkSelBLayer << endreq; 
          msg() << MSG::DEBUG << " TrkSel_SiHits = "   << m_trkSelSiHits << endreq; 
          msg() << MSG::DEBUG << " TrkSel_D0 = "       << m_trkSelD0 << endreq; 
          msg() << MSG::DEBUG << " TrkSel_Z0 = "       << m_trkSelZ0 << endreq; 
          msg() << MSG::DEBUG << " TrkSel_Pt = "       << m_trkSelPt << endreq; 

          msg() << MSG::DEBUG << " SizeMVtx = "       << m_sizeMVtx << endreq; 
          msg() << MSG::DEBUG << " bMVtx = "          << m_bMVtx << endreq; 
          msg() << MSG::DEBUG << " uMVtx = "          << m_uMVtx << endreq;
          msg() << MSG::DEBUG << " SizeEVtx = "       << m_sizeEVtx << endreq; 
          msg() << MSG::DEBUG << " bEVtx = "          << m_bEVtx << endreq; 
          msg() << MSG::DEBUG << " uEVtx = "          << m_uEVtx << endreq;  
          msg() << MSG::DEBUG << " SizeNVtx = "       << m_sizeNVtx << endreq; 
          msg() << MSG::DEBUG << " bNVtx = "          << m_bNVtx << endreq; 
          msg() << MSG::DEBUG << " uNVtx = "          << m_uNVtx << endreq;
          msg() << MSG::DEBUG << " SizeSV = "         << m_sizeSV << endreq; 
          msg() << MSG::DEBUG << " bSV = "            << m_bSV << endreq; 
          msg() << MSG::DEBUG << " uSV = "            << m_uSV << endreq;  
        }

    m_trigBjetTagger = new TrigBjetTagger(this, msg(), msgLvl());
    m_constTrigBjetTagger  = const_cast<const TrigBjetTagger*>(m_trigBjetTagger);

    if (msgLvl() <= MSG::DEBUG) 
      msg() << MSG::DEBUG << "Retrieving tuning likelihoods." << endreq;

    m_tuningLikelihoodIP1D = new TuningLikelihood(&m_sizeIP1D[0], &m_bIP1D[0], &m_uIP1D[0], m_sizeIP1D.size());
    m_trigBjetTagger->fillLikelihoodMap("IP1D", m_tuningLikelihoodIP1D);
    
    m_tuningLikelihoodIP2D = new TuningLikelihood(&m_sizeIP2D[0], &m_bIP2D[0], &m_uIP2D[0], m_sizeIP2D.size());   
    m_trigBjetTagger->fillLikelihoodMap("IP2D", m_tuningLikelihoodIP2D);

    m_tuningLikelihoodIP3D = new TuningLikelihood(&m_sizeIP3D[0], &m_bIP3D[0], &m_uIP3D[0], m_sizeIP3D.size()); 
    m_trigBjetTagger->fillLikelihoodMap("IP3D", m_tuningLikelihoodIP3D);

    if (m_useLowSiHits) {
      m_tuningLikelihoodIP1D_lowSiHits = new TuningLikelihood(&m_sizeIP1D_lowSiHits[0], &m_bIP1D_lowSiHits[0], &m_uIP1D_lowSiHits[0], m_sizeIP1D_lowSiHits.size());
      m_trigBjetTagger->fillLikelihoodMap("IP1D_lowSiHits", m_tuningLikelihoodIP1D_lowSiHits);

      m_tuningLikelihoodIP2D_lowSiHits = new TuningLikelihood(&m_sizeIP2D_lowSiHits[0], &m_bIP2D_lowSiHits[0], &m_uIP2D_lowSiHits[0], m_sizeIP2D_lowSiHits.size());
      m_trigBjetTagger->fillLikelihoodMap("IP2D_lowSiHits", m_tuningLikelihoodIP2D_lowSiHits);

      m_tuningLikelihoodIP3D_lowSiHits = new TuningLikelihood(&m_sizeIP3D_lowSiHits[0], &m_bIP3D_lowSiHits[0], &m_uIP3D_lowSiHits[0], m_sizeIP3D_lowSiHits.size());
      m_trigBjetTagger->fillLikelihoodMap("IP3D_lowSiHits", m_tuningLikelihoodIP3D_lowSiHits);
    }

    m_tuningLikelihoodMVtx = new TuningLikelihood(&m_sizeMVtx[0], &m_bMVtx[0], &m_uMVtx[0], m_sizeMVtx.size()); 
    m_trigBjetTagger->fillLikelihoodMap("MVTX", m_tuningLikelihoodMVtx);
    
    m_tuningLikelihoodEVtx = new TuningLikelihood(&m_sizeEVtx[0], &m_bEVtx[0], &m_uEVtx[0], m_sizeEVtx.size()); 
    m_trigBjetTagger->fillLikelihoodMap("EVTX", m_tuningLikelihoodEVtx);

    m_tuningLikelihoodNVtx = new TuningLikelihood(&m_sizeNVtx[0], &m_bNVtx[0], &m_uNVtx[0], m_sizeNVtx.size()); 
    m_trigBjetTagger->fillLikelihoodMap("NVTX", m_tuningLikelihoodNVtx);

    m_tuningLikelihoodSV = new TuningLikelihood(&m_sizeSV[0], &m_bSV[0], &m_uSV[0], m_sizeSV.size()); 
    m_trigBjetTagger->fillLikelihoodMap("SVTX", m_tuningLikelihoodSV);

    // Retrieve TrigTrackJetFinder tool
        StatusCode sc = m_trackJetFinderTool.retrieve();
    if(sc.isFailure()) {
      msg() << MSG::FATAL << "Failed to locate tool " << m_trackJetFinderTool << endreq;
      return HLT::BAD_JOB_SETUP;
    } else
      msg() << MSG::INFO << "Retrieved tool " << m_trackJetFinderTool << endreq;
  
    m_trigBjetPrmVtxInfo = new TrigBjetPrmVtxInfo();
    m_trigBjetSecVtxInfo = new TrigBjetSecVtxInfo();
    m_trigBjetJetInfo    = new TrigBjetJetInfo();

    return HLT::OK;
}



//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigBjetFex::getCollection(const xAOD::TrackParticleContainer*& pointerToEFTrackCollections, const HLT::TriggerElement* outputTE) {

  std::vector<const xAOD::TrackParticleContainer*> vectorOfEFTrackCollections;

  HLT::ErrorCode status = getFeatures(outputTE, vectorOfEFTrackCollections, ""); 

  if (status != HLT::OK) {
    msg() << MSG::ERROR << "Failed to get TrackParticleContainer from the trigger element" << endreq;
  } else if (msgLvl() <= MSG::DEBUG) 
    msg() << MSG::DEBUG << "Got " << vectorOfEFTrackCollections.size() << " TrackParticleContainer" << endreq;
  
  std::vector<const xAOD::TrackParticleContainer*>::iterator pTrackColl    = vectorOfEFTrackCollections.begin();
  std::vector<const xAOD::TrackParticleContainer*>::iterator lastTrackColl = vectorOfEFTrackCollections.end();
  
  if (pTrackColl == lastTrackColl) {
    pointerToEFTrackCollections = 0;
    return HLT::ERROR;
  } else {
    pointerToEFTrackCollections = *pTrackColl;
    return HLT::OK;
  }

}

//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigBjetFex::getPrmVtxCollection(const xAOD::VertexContainer*& pointerToEFPrmVtxCollections, 
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
    return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::MISSING_FEATURE);
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigBjetFex::getSecVtxCollection(const Trk::VxSecVertexInfoContainer*& pointerToSecVtxCollections,
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
    return HLT::ErrorCode(HLT::Action::ABORT_CHAIN,HLT::Reason::MISSING_FEATURE);
}


//** ----------------------------------------------------------------------------------------------------------------- **//


  //HLT::ErrorCode TrigBjetFex::getSecVtxInfo(const Trk::VxSecVertexInfoContainer*& pointerToEFSecVtxCollections , const xAOD::VertexContainer*& pointerToEFPrmVtxCollections, const TrigVertexCollection*& pointerToPrmVtxCollections) {

HLT::ErrorCode TrigBjetFex::getSecVtxInfo(const Trk::VxSecVertexInfoContainer*& pointerToEFSecVtxCollections, 
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
    	return HLT::OK;
    }
  
    const Trk::VxSecVertexInfo* m_secVertexInfo = (*pointerToEFSecVtxCollections)[0];
    if(!m_secVertexInfo) 
    {   
       ATH_MSG_DEBUG("No secondary vertex when extracting sec vtx info");
       return HLT::OK;
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
    std::vector<xAOD::Vertex*>::const_iterator verticesIt=myVertices.begin();
    std::vector<xAOD::Vertex*>::const_iterator verticesEnd=myVertices.end();
    
    if(myVertices.size()>1) 
    {
        ATH_MSG_WARNING("Secondary vertex from InDetVKalVxInJet has more than one vertex, is this ok?");
    }
    
    for(xAOD::Vertex * vertexIt : myVertices) 
    {
        if(!(*vertexIt)) 
        {
    	ATH_MSG_DEBUG("Secondary vertex from InDetVKalVxInJet has zero pointer. Skipping this vtx..");
    	continue;
        }
        ATH_MSG_DEBUG("VxCandidate at (" 
    	    << (*vertexIt)->position().x() << "," 
    	    << (*vertexIt)->position().y() << "," 
         	    << (*vertexIt)->position().z());
    
        // Check we have vxTrackAtVertex
        std::vector<Trk::VxTrackAtVertex> * myTracks = 0;
        if( (*vertex)->vxTrackAtVertexAvailable() )
        {	
    	myTracks = &((*vertexIt)->vxTrackAtVertex());
    	NTracksInSV+=myTracks->size();
        }
        else
        {
    	ATH_MSG_WARNING("TrackAtVertex object not available! Setting 0 tracks attached to this vertex");
        }
    }
    m_trigBjetSecVtxInfo->setNTrksInVtx(NTracksInSV);
  
    //Calculate decay length and significance here
    //Use the same utilities as in InDetVKalVxInJet
    //Some gymnastics needed...ugh...
    if(!pPrmVrt) 
    {
        return HLT::OK;
    }
  
    ATH_MSG_DEBUG("Primary vertex for decay length (" 
    	<< pvselected->position().x() << "," 
    	<< pvselected->position().y() << "," 
    	<< pvselected->position().z() << ") and error (" 
    	<< pvselected->covariancePosition()(0,0) << "," 
    	<< pvselected->covariancePosition()(1,1) << "," 
    	<< pvselected->covariancePosition()(2,2) << ")");
    
    /* Not needed
    // Needs some Eigen migration magic here
    CLHEP::HepLorentzVector jetDirection;
  
    if (m_useJetDirection == 1) 
    {
        jetDirection.setX(cos(m_trigBjetJetInfo->phiJet()));
        jetDirection.setY(sin(m_trigBjetJetInfo->phiJet()));
        jetDirection.setZ(sinh(m_trigBjetJetInfo->etaJet()));
    } 
    else if(m_useJetDirection == 2) 
    {
        jetDirection.setX(cos(m_trigBjetJetInfo->phiTrkJet()));
        jetDirection.setY(sin(m_trigBjetJetInfo->phiTrkJet()));
        jetDirection.setZ(sinh(m_trigBjetJetInfo->etaTrkJet()));
    } 
    else if(m_useJetDirection == 3) 
    {
        jetDirection.setX(cos(m_trigBjetJetInfo->phiRoI()));
        jetDirection.setY(sin(m_trigBjetJetInfo->phiRoI()));
        jetDirection.setZ(sinh(m_trigBjetJetInfo->etaRoI()));
    }
    
    //Gymnastics to get sec vtx position and error matrix in correct format for later
    // Getting the first vertex??
    const Amg::Vector3D SecVrt = myVertices[0]->position();

    std::vector<double> SecVrtErr(6,0);
    SecVrtErr[0] = myVertices[0]->covariancePosition()(0,0);
    SecVrtErr[1] = myVertices[0]->covariancePosition()(0,1);
    SecVrtErr[2] = myVertices[0]->covariancePosition()(1,1); 
    SecVrtErr[3] = myVertices[0]->covariancePosition()(0,2);
    SecVrtErr[4] = myVertices[0]->covariancePosition()(1,2);
    SecVrtErr[5] = myVertices[0]->covariancePosition()(2,2); 
    
    //Calculate the jet-vrt direction for the signed decay length calculation
    double JetVrtDir = jetDirection.px()*(SecVrt.x()-pvselected->position().x())
	+ jetDirection.py()*(SecVrt.y()-pvselected->position().y())
	+ jetDirection.pz()*(SecVrt.z()-pvselected->position().z());

    //Decay length
    float sign3D = 0.0;
    double dist3D = m_taggerHelper->VrtVrtDist(*pvselected,SecVrt,SecVrtErr,sign3D);

    if(JetVrtDir < 0)
    {
	sign3D = -sign3D;
    }

    //m_trigBjetSecVtxInfo->setDecayLengthSignificance(sign3D);  */
    //Toggle
    m_trigBjetSecVtxInfo->isValid(true);

    double distance = (SecVrt - pvseleted->position()).mag();
    // Note, storing the distance 3D between SV and PV, watch out with the name 
    // of the method!!
    m_trigBjetSecVtxInfo->setDecayLengthSignificance(distance);
    if(msgLvl() <= MSG::DEBUG) 
    {
	/*if(fabs(distance-dist3D)>0.0001) 
	{
	    ATH_MSG_DEBUG("decay length distance do not match among tools: tool " 
		    << dist3D << " manual " << distance);
	}*/
        double dist2D = (SecVrt - pvselected->position()).perp();            
	ATH_MSG_DEBUG("Calculating secondary vertex decay length with primary vertex at (" 
  		<< pvselected->position().x()/CLHEP::mm << "," << pvselected->position().y()/CLHEP::mm
		<< "," << pvselected->position().z()/CLHEP::mm << ") and sec vtx at ("
   		<< SecVrt.x()/CLHEP::mm << "," << SecVrt.y()/CLHEP::mm << "," << SecVrt.z()/CLHEP::mm 
	//	<< ") and jet direction (px,py,pz) = ("
	//       	<< jetDirection.px() << "," << jetDirection.py() << "," << jetDirection.pz() 
		<<  ") which gives 3D decay length " << distance/CLHEP::mm 
	//	<< " (VrtVrtDist tool " 
	//	<< dist3D << ") and 3D significance " << sign3D  
		<< " and 2D(R/phi) decay length " 
		<< dist2D/CLHEP::mm << " [all in mm]");
    }
    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


bool TrigBjetFex::efTrackSel(const xAOD::TrackParticle*& track, unsigned int i) 
{
    float zv = m_trigBjetPrmVtxInfo->zPrmVtx();
    //const Trk::TrackSummary *summary = track->trackSummary();
    //const Trk::FitQuality *quality   = track->fitQuality();
    
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
    return true;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigBjetFex::hltExecute(const HLT::TriggerElement* inputTE, HLT::TriggerElement* outputTE) 
{
    ATH_MSG_DEBUG("Executing TrigBjetFex");

    // Clear and initialize data members
    m_totSelTracks = 0;
    m_totTracks    = 0;
  
    m_trigBjetPrmVtxInfo->clear();
    m_trigBjetSecVtxInfo->clear();
    m_trigBjetJetInfo->clear();
    // JDC:: Why are you asigning this address?
    std::vector<TrigBjetTrackInfo> trigBjetTrackInfoVector;
    m_trigBjetTrackInfoVector = &trigBjetTrackInfoVector;

    // This is really horrible... 
    ATH_MSG_VERBOSE("Printing out inputTE will get rid of the compilation warning: " << inputTE);
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
    if(m_retrieveHLTJets && m_instance == "EF") 
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
    		// JDC:: I think by construction, this case it cannot be... not sure
		//     	 look at TrigOperationalInfo implementation
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
  
    const xAOD::JetContainer* jets(0); // JDC:: Check this initialization, shoudn't be "const ... * jet = 0;"?
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
   
    xAOD::BTaggingAuxContainer trigBjetAuxContainer;
    m_trigBTaggingContainer = new xAOD::BTaggingContainer();
    m_trigBTaggingContainer->setStore(&trigBjetAuxContainer);
    
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
    xAOD::Vertex *pvselected = 0;
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

    m_trigBjetPrmVtxInfo->setBeamSpot(m_xBeamSpot, m_yBeamSpot, m_zBeamSpot);
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
     	m_mon_trk_a0.push_back(track->d0());
    	m_mon_trk_z0.push_back(track->z0());
     	
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
    
    //std::vector<std::string>::iterator pTagger    = m_taggers.begin();
    //std::vector<std::string>::iterator lastTagger = m_taggers.end();
    
    /*if(m_useBeamSpotFlag && !m_trigBjetPrmVtxInfo->beamSpotStatus()) 
    {
	ATH_MSG_DEBUG("Beam spot status flag set to " << m_trigBjetPrmVtxInfo->beamSpotStatus() 
		<< ". Discriminant weights are not computed.");
   	ATH_MSG_DEBUG("Beam spot flag set to " << m_useBeamSpotFlag << 
		". Discriminant weights are not computed.");
     	m_listCutApplied.push_back(CutListMonitor::BeamSpotWrongStatus);
    	//m_trigBjetTagger->getWeights();
    } 
    else if(m_trigBjetPrmVtxInfo->xBeamSpotWidth()>m_setBeamSpotWidth || 
	    m_trigBjetPrmVtxInfo->yBeamSpotWidth()>m_setBeamSpotWidth) 
    {
     	ATH_MSG_DEBUG("Beam spot width is more than " << m_setBeamSpotWidth 
		<< "um. Discriminant weights are not computed.");
    	m_listCutApplied.push_back(CutListMonitor::BeamSpotTooWide);
    	//m_trigBjetTagger->getWeights();
    }
    else 
    {
  	ATH_MSG_DEBUG("Computing discriminant weights using taggers: " << m_taggers <<
	       	" and using calibration from " << (m_useParamFromData==0 ? "MC" : "data") 
		<< " for CHI2");
      	//Look for a sec vertex?
	bool retrieveSV = false;
    	for( ; pTagger != lastTagger; ++pTagger) 
	{
	    if((*pTagger).find("VTX") != std::string::npos) 
	    {
	      	retrieveSV = true;
	       	break;
	    }
    	}
     	if(retrieveSV) 
	{
	    // Get secondary vertex information at EF from TrigBjetFex::getSecVtxInfo
	    if(getSecVtxInfo(pointerToEFSecVtxCollections, pointerToEFPrmVtxCollections) != HLT::OK)
	    { 
		ATH_MSG_DEBUG("No EF SV information retrieved from TrigBjetFex::getSecVtxInfo");
	    }
	    ATH_MSG_DEBUG(*m_trigBjetSecVtxInfo);
     	}
	// JDC:: Not need the weights, just the values
    	//m_trigBjetTagger->getWeights(m_trigBjetTrackInfoVector, 
	//	m_trigBjetPrmVtxInfo, m_trigBjetSecVtxInfo, m_trigBjetJetInfo);
    }*/
    
    // Get secondary vertex information at EF --> Do this first!!
    HLT::ErrorCode statSVInfo = getSecVtxInfo(pointerToEFSecVtxCollections,pvselected);

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
	    ATH_MSG_DEBUG("Bjet slice summary (Run " << pEventInfo->event_ID()->run_number() 
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
       
	    /*pTagger = m_taggers.begin();
	    for ( ; pTagger != lastTagger; ++pTagger)
	    {
		ATH_MSG_DEBUG("REGTEST:  X(" << (*pTagger) << ") = " 
			<< m_trigBjetTagger->taggersXMap((*pTagger)));
	    }*/
	}
    }

    // Attach a Primary Vertex and Secondary Vertex ... that'all we need
    // -----------------------------------
    // -----------------------------------
    // Create TrigEFBjet and attach feature
    // -----------------------------------
    // Note that the meaning of some characteristics are tuned to 
    // DV case
    TrigEFBjet* trigEFBjet = new TrigEFBjet(roiDescriptor->roiId(), 
	    m_trigBjetJetInfo->etaJet(), m_trigBjetJetInfo->phiJet(),
	    0, 0, 0, m_trigBjetPrmVtxIne fo->zPrmVtx(), m_trigBjetJetInfo->etJet(),
	    -1,-1,-1, -1,-1,  // Note, i can use this to fill 
	    m_trigBjetSecVtxInfo->decayLengthSignificance(), m_trigBjetSecVtxInfo->vtxMass(), 
	    m_trigBjetSecVtxInfo->energyFraction(), m_trigBjetSecVtxInfo->n2TrkVtx()); 
    
    trigEFBjet->validate(true);
    m_trigEFBjetColl->push_back(trigEFBjet);
    
    if(!m_trigEFBjetColl) 
    {
	ATH_MSG_ERROR("Feature TrigEFBjetContainer not found");
       	return HLT::ErrorCode(HLT::Action::ABORT_JOB, HLT::Reason::BAD_JOB_SETUP);
    }
    
    HLT::ErrorCode stat = attachFeature(outputTE, m_trigEFBjetColl, "EFBjetFex");
    if(stat != HLT::OK) 
    {
	ATH_MSG_DEBUG("Failed to attach TrigEFBjetContainer to navigation");
	return stat;
    }

    // -----------------------------------
    // Create xAOD::BTagging and attach feature
    // -----------------------------------
    xAOD::BTagging * newBTag = new xAOD::BTagging();
    m_trigBTaggingContainer->push_back(newBTag);
    newBTag->setSV1_pu(m_trigBjetTagger->taggersPuMap("MVTX")*
	    m_trigBjetTagger->taggersPuMap("NVTX")*m_trigBjetTagger->taggersPuMap("EVTX"));
    newBTag->setSV1_pb(m_trigBjetTagger->taggersPbMap("MVTX")*
	    m_trigBjetTagger->taggersPbMap("NVTX")*m_trigBjetTagger->taggersPbMap("EVTX"));
    newBTag->setIP2D_pu(m_trigBjetTagger->taggersPuMap("IP2D"));
    newBTag->setIP2D_pb(m_trigBjetTagger->taggersPbMap("IP2D"));
    
    newBTag->setIP3D_pu(m_trigBjetTagger->taggersPuMap("IP3D"));
    newBTag->setIP3D_pb(m_trigBjetTagger->taggersPbMap("IP3D"));

    ATH_MSG_DEBUG("IP2D u/b: " << m_trigBjetTagger->taggersPuMap("IP2D") << "/" 
	    << m_trigBjetTagger->taggersPbMap("IP2D") << "   IP3D u/b: " 
	    << m_trigBjetTagger->taggersPuMap("IP3D") << "/" << m_trigBjetTagger->taggersPbMap("IP3D")
    	    << "   SV1 u/b: " << m_trigBjetTagger->taggersPuMap("MVTX")*
	            m_trigBjetTagger->taggersPuMap("NVTX")*m_trigBjetTagger->taggersPuMap("EVTX") 

	    << "/" <<  m_trigBjetTagger->taggersPbMap("MVTX")*
	    m_trigBjetTagger->taggersPbMap("NVTX")*m_trigBjetTagger->taggersPbMap("EVTX"));
    if(!m_trigBTaggingContainer) 
    {
    	ATH_MSG_ERROR("Feature BTaggingContainer not found");
       	return HLT::ErrorCode(HLT::Action::ABORT_JOB, HLT::Reason::BAD_JOB_SETUP);
    }
  
    stat = attachFeature(outputTE, m_trigBTaggingContainer, "HLTBjetFex");
    if(stat != HLT::OK) 
    {
	ATH_MSG_DEBUG("Failed to attach BTaggingContainer to navigation");
    	return stat;
    }*/
    return HLT::OK;
}


//** ----------------------------------------------------------------------------------------------------------------- **//


HLT::ErrorCode TrigBjetFex::hltFinalize() {

  if (msgLvl() <= MSG::INFO) 
    msg() << MSG::INFO << "Finalizing TrigBjetFex" << endreq;

  return HLT::OK;
}



