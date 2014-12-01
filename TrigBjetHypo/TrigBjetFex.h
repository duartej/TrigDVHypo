// ************************************************
//
// NAME:     TrigBjetFex.h
// PACKAGE:  Trigger/TrigHypothesis/TrigBjetHypo
//
// AUTHOR:   Jordi Duarte-Campderros (Based in TrigBjetFex
//           class)
// EMAIL:    jorge.duarte.campderros@cern.ch
// 
// ************************************************

#ifndef TRIGBJETFEX_H
#define TRIGBJETFEX_H

#include "GaudiKernel/ToolHandle.h"

#include "TrigInterfaces/FexAlgo.h"

#include "TrigBjetHypo/TrigBjetTagger.h"
#include "TrigBjetHypo/TuningLikelihood.h"
#include "TrigBjetHypo/TrigBjetDataHelper.h"

#include "TrigTrackJetFinderTool/ITrigTrackJetFinderTool.h"
#include "VxSecVertex/VxSecVertexInfo.h"



class MsgStream;
class TrigEFBjetContainer;
class TriggerElement;

class TrigInDetTrackCollection;
class TrackParticleContainer;
class TrigVertexCollection;
class VxContainer;
class TrigBjetTagger;
class TrigRoiDescriptor;

class TaggerHelper;
class TrigBjetTrackInfo;
class TrigBjetPrmVtxInfo;
class TrigBjetSecVtxInfo;
class TrigBjetJetInfo;

namespace Rec {
  class TrackParticleContainer;
  class TrackParticle;
}

/**
 * @brief FEX class for the b-jet weight construction to be used by subsequent hypothesis algorithm.
 *
 * @author Andrea Coccaro <Andrea.Coccaro@ge.infn.it>
 *
 * This is the base feature extraction class for the HLT b-jet slice. 
 * It computes relevant quantities based on track impact parameters and secondary vertices, creates Bjet object and attaches it to TE.
 */


class TrigBjetFex: public HLT::FexAlgo {
  
 public:

  /** @brief For monitoring purposes the TrigBjetTagger class has been declared as friend. */
  friend class TrigBjetTagger;

  /** @brief Constructor. */
  TrigBjetFex(const std::string&, ISvcLocator*);
  /** @brief Destructor. */
  ~TrigBjetFex();

  HLT::ErrorCode hltInitialize();
  HLT::ErrorCode hltExecute(const HLT::TriggerElement*, HLT::TriggerElement*);
  HLT::ErrorCode hltFinalize(); 

 private:

  std::string m_jetKey;
  std::string m_priVtxKey;

  /** @brief To retrieve track collections reconstructed at EF and stored in TrackParticleContainer. */
  HLT::ErrorCode getCollection(const xAOD::TrackParticleContainer*&, const HLT::TriggerElement*);
  /** @brief To retrieve primary vertex collections reconstructed and stored in TrigVertexCollection. */
  //  HLT::ErrorCode getPrmVtxCollection(const xAOD::VertexContainer*&, const HLT::TriggerElement*);
  /** @brief To retrieve primary vertex collections reconstructed and stored in VxContainer. */
  HLT::ErrorCode getPrmVtxCollection(const xAOD::VertexContainer*&, const HLT::TriggerElement*);
  /** @brief To retrieve secondary vertex collections reconstructed and stored in Trk::VxSecVertexInfoContainer. */
  HLT::ErrorCode getSecVtxCollection(const Trk::VxSecVertexInfoContainer*&, const HLT::TriggerElement*);
  /** @brief To retrieve the EF secondary vertex information. */
  //HLT::ErrorCode getSecVtxInfo(const Trk::VxSecVertexInfoContainer*&, const xAOD::VertexContainer*&, const TrigVertexCollection*&);
  HLT::ErrorCode getSecVtxInfo(const Trk::VxSecVertexInfoContainer*&, const xAOD::VertexContainer*&);

  /** @brief To select EF tracks. */
  bool efTrackSel(const xAOD::TrackParticle*&, unsigned int);

  /** @brief TrigTrackJetFinder tool. */
  ToolHandle<ITrigTrackJetFinderTool> m_trackJetFinderTool;
  /** @brief Pointer to TrigEFBjet collection. */
  TrigEFBjetContainer* m_trigEFBjetColl;
  /** @brief Pointer to xAOD::BTaggingContainer. */
  xAOD::BTaggingContainer* m_trigBTaggingContainer;

  /** @brief Pointer to TaggerHelper class. */ 
  TaggerHelper* m_taggerHelper;

  /** @brief Pointer to TrigBjetTagger class. */ 
  TrigBjetTagger* m_trigBjetTagger;
  /** @brief Const pointer to TrigBjetTagger class (for monitoring). */ 
  const TrigBjetTagger* m_constTrigBjetTagger;

  /** @brief Internal EDM class for primary vertex and beam-spot information, defined in TrigBjetDataHelper. */ 
  TrigBjetPrmVtxInfo* m_trigBjetPrmVtxInfo;
  /** @brief Internal EDM class for secondary vertex information, defined in TrigBjetDataHelper. */ 
  TrigBjetSecVtxInfo* m_trigBjetSecVtxInfo;
  /** @brief Internal EDM class for jet direction information, defined in TrigBjetDataHelper. */ 
  TrigBjetJetInfo*    m_trigBjetJetInfo;
  /** @brief Internal EDM class for track parameter information, defined in TrigBjetDataHelper. */ 
  std::vector<TrigBjetTrackInfo>* m_trigBjetTrackInfoVector;

  /** @brief track reconstruction algorithm. */
  int m_algo;
  /** @brief list of computed taggers. */
  std::vector<std::string> m_taggers;
  /** @brief switch to use impact parameter errors parametrization. */
  bool m_useErrIPParam;
  /** @brief switch to use primary vertex computed with TrigT2HistoPrmVtx or InDetTrigPriVxFinder algorithm. */
  bool m_histoPrmVtxAtEF;
  /** @brief switch to perform track-RoI eta/phi matching when selecting reconstructed tracks. */
  bool m_useEtaPhiTrackSel;
  /** @brief switch to estimate the track impact parameter sign using the HLT jet direction (1) or the HLT track-jet direction (2) or the LVL1 jet RoI direction (3). */
  unsigned int m_useJetDirection;
  /** @brief switch to retrieve the HLT jets. True (false) for the 2012 (2011) setup. */
  bool m_retrieveHLTJets;
  /** @brief switch to tag the HLT jets. */
  unsigned int m_tagHLTJets;

  /** @brief string corresponding to the trigger level in which the algorithm is running. */
  std::string m_instance;

  /** @brief to check the beam spot flag status and start track selection accordingly. */
  bool m_useBeamSpotFlag;
  /** @brief to set the maximum beam spot width to enable the selection. */
  float m_setBeamSpotWidth;

  /** @brief to use JetProb parametrization derived from data. */
  bool m_useParamFromData;

  /** @brief lower bound of the chi square of the reconstructed track (to perform track selection). */
  float m_trkSelChi2;
  /** @brief lower bound of the number of hits on b-layer of reconstructed track (to perform track selection). */
  int   m_trkSelBLayer;
  /** @brief lower bound of the number of hits in pixel detector of reconstructed track (to perform track selection). */
  int   m_trkSelPixHits;
  /** @brief lower bound of the number of hits in silicon detectors of reconstructed track (to perform track selection). */
  int   m_trkSelSiHits;
  /** @brief upper bound of transverse impact parameter of reconstructed track (to perform track selection). */
  float m_trkSelD0;
  /** @brief upper bound of longitudinal impact parameter of reconstructed track (to perform track selection). */
  float m_trkSelZ0;
  /** @brief lower bound of pT of the reconstructed track (to perform track selection). */
  float m_trkSelPt;

  //////////////////////
  //* for monitoring *//
  //////////////////////

  /** @brief to monitor track selection. */
  std::vector<float> m_listCutApplied;

  /** @brief to monitor track quantities and effect ot the selection. */
  std::vector<float> m_mon_trk_a0, m_mon_trk_a0_sel, m_mon_trk_Sa0_sel;
  std::vector<float> m_mon_trk_z0, m_mon_trk_z0_sel, m_mon_trk_z0_sel_PV, m_mon_trk_Sz0_sel;
  std::vector<float> m_mon_trk_prob;

  /** @brief To retrieve selected tracks in percentage. */ 
  float totSelectedTracks() const {
    if (!m_totTracks) return -0.1;
    else return (float)m_totSelTracks/(float)m_totTracks;
  };

  /** @brief Number of reconstructed tracks per RoI. */
  float m_totTracks;
  /** @brief Number of selected tracks per RoI. */
  float m_totSelTracks;

  /** @brief Delta eta between the LVL1 jet RoI and the EF jet. */
  float m_deltaEtaJet;
  /** @brief Delta phi between the LVL1 jet RoI and the EF jet. */
  float m_deltaPhiJet;
  /** @brief Delta eta between the LVL1 jet RoI and the EF track-jet. */
  float m_deltaEtaTrkJet;
  /** @brief Delta phi between the LVL1 jet RoI and the EF track-jet. */
  float m_deltaPhiTrkJet;
  /** @brief Delta eta between the EF jet and the EF track-jet. */
  float m_deltaEtaJetTrkJet;
  /** @brief Delta phi between the EF jet and the EF track-jet. */
  float m_deltaPhiJetTrkJet;

  //////////////////////
  //* for calibration *//
  //////////////////////

  /** @brief Pointer to TuningLikelihood class for IP1D likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodIP1D;
  /** @brief Pointer to TuningLikelihood class for IP2D likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodIP2D;
  /** @brief Pointer to TuningLikelihood class for IP3D likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodIP3D;
  /** @brief Pointer to TuningLikelihood class for IP1D likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodIP1D_lowSiHits;
  /** @brief Pointer to TuningLikelihood class for IP2D likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodIP2D_lowSiHits;
  /** @brief Pointer to TuningLikelihood class for IP3D likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodIP3D_lowSiHits;
  /** @brief Pointer to TuningLikelihood class for MVtx likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodMVtx;
  /** @brief Pointer to TuningLikelihood class for EVtx likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodEVtx;
  /** @brief Pointer to TuningLikelihood class for NVtx likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodNVtx;
  /** @brief Pointer to TuningLikelihood class for SV likelihood tagger. */ 
  TuningLikelihood* m_tuningLikelihoodSV;

  /** @brief necessary for calibration constants for one-dimensional likelihood taggers based on impact parameters. */
  std::vector<float> m_sizeIP1D, m_bIP1D, m_uIP1D;  
  /** @brief necessary for calibration constants for two-dimensional likelihood taggers based on impact parameters. */
  std::vector<float> m_sizeIP2D, m_bIP2D, m_uIP2D;
  /** @brief necessary for calibration constants for three-dimensional likelihood taggers based on impact parameters. */
  std::vector<float> m_sizeIP3D, m_bIP3D, m_uIP3D;

  /** @brief necessary for calibration constants for one-dimensional likelihood taggers based on impact parameters. */
  std::vector<float> m_sizeIP1D_lowSiHits, m_bIP1D_lowSiHits, m_uIP1D_lowSiHits;
  /** @brief necessary for calibration constants for two-dimensional likelihood taggers based on impact parameters. */
  std::vector<float> m_sizeIP2D_lowSiHits, m_bIP2D_lowSiHits, m_uIP2D_lowSiHits;
  /** @brief necessary for calibration constants for three-dimensional likelihood taggers based on impact parameters. */
  std::vector<float> m_sizeIP3D_lowSiHits, m_bIP3D_lowSiHits, m_uIP3D_lowSiHits;

  bool m_useLowSiHits;

  /** @brief necessary for calibration constants for one-dimensional likelihood taggers based on mass of secondary vertex. */
  std::vector<float> m_sizeMVtx, m_bMVtx, m_uMVtx;
  /** @brief necessary for calibration constants for one-dimensional likelihood taggers based on energy secondary vertex. */
  std::vector<float> m_sizeEVtx, m_bEVtx, m_uEVtx;
  /** @brief necessary for calibration constants for one-dimensional likelihood taggers based on number of tracks linked to secondary vertex. */
  std::vector<float> m_sizeNVtx, m_bNVtx, m_uNVtx;
  /** @brief necessary for calibration constants for three-dimensional likelihood taggers based on secondary vertex. */
  std::vector<float> m_sizeSV, m_bSV, m_uSV;

  /** @brief necessary for calibration constants for BjetProbability */
  std::vector<float> m_par_0_MC, m_par_1_MC, m_par_0_DT, m_par_1_DT;
};

#endif
