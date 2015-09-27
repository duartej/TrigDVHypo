// ************************************************
//
// NAME:     TrigDvFex.h
// PACKAGE:  Trigger/TrigHypothesis/TrigBjetHypo
//
// AUTHOR:   Jordi Duarte-Campderros 
//           (Based in TrigBjetFex class)
// EMAIL:    jorge.duarte.campderros@cern.ch
// 
// ************************************************

#ifndef TRIGDVFEX_H
#define TRIGDVFEX_H

#include "GaudiKernel/ToolHandle.h"
// JDC TESTING
#include "GaudiKernel/ServiceHandle.h"

#include "TrigInterfaces/FexAlgo.h"
//#include "TrigSteeringEvent/Enums.h"

#include "TrigBjetHypo/TrigBjetDataHelper.h"

// JDC TESTING
#include "InDetPrepRawData/TRT_DriftCircleContainer.h"
#include "InDetPrepRawData/SiClusterContainer.h"

// system
#include<vector>
#include<unordered_set>
#include<string>


//class HLT::TriggerElement;  Enable by default

// Not real classes, just typedefs
//class xAOD::TrackParticle;
//class xAOD::TrackParticleContainer;
//class xAOD::VertexContainer;
//
//class Rec::TrackParticle; --> enabled by default
//class Rec::TrackParticleContainer; --> enabled by default

//class Trk::VxSecVertexInfoContainer;  --> typedef in VxSecVertexInfo.h

class TrigEFBjetContainer;

class TrigBjetTrackInfo;
class TrigBjetJetInfo;

// JDC TESTING
class IRegSelSvc;
class TrigRoiDescriptor;
namespace Trk
{
    class PrepRawData;
}
//namespace InDet
//{
//    class SiClusterContainer;
//    class TRT_DriftCircleContainer;
//}

enum class TrigDvFexDet 
{
    PIXEL = 1,
    SCT   = 2,
    TRT   = 3
};

struct PrdsAtRoI
{
    std::unordered_set<const Trk::PrepRawData*> pixel;
    std::unordered_set<const Trk::PrepRawData*> blayer;
    std::unordered_set<const Trk::PrepRawData*> sct;
    std::unordered_set<const Trk::PrepRawData*> trt;

    inline void clear() 
    {
        pixel.clear();
        blayer.clear();
        sct.clear();
        trt.clear();
    };

    inline unsigned int size() const
    {
        return (pixel.size()+sct.size()+trt.size());
    };
};

/**
 * @brief FEX class for the displaced vertex searches (r_{DV} < 300 mm) to be used 
 *        by subsequent hypothesis algorithm.
 *
 * @author Jordi Duarte-Campderros <jorge.duarte.campderros@cern.ch>
 *
 * This is the base feature extraction class for the HLT DV-searches triggers. 
 * The trigger strategy is based in the RoI-level reconstruction of the tracks.
 * Each RoI (L1 jet) defines a eta+-DeltaEta and phi+-DeltaPhi which is used to look 
 * after seeds at the Si (Pixel and SCT) and/or TRT detectors. The seeds obtained in
 * that RoI are used to reconstruct tracks compatibles with the RoI eta and phi (i.e. 
 * the tracks pass through this eta+-deltaEta,phi+-deltaPhi).
 *
 * --- TO BE DELETED: NEEDS TO DECIDE WHICH QUANTITES TO STORED ---
 * The Fex computes relevant quantities based on track impact pa and secondary vertices, creates a 
 * modified Bjet object and attaches it to TE. Note that the Bjet created stores the decay length instead
 * of the significance of the decay length. Also, the Primary Vertex (with higher sum_{tracks}pT^2) and
 * the secondary vertices found in the event are filled in the Bjet object
 */


class TrigDvFex: public HLT::FexAlgo 
{
    public:
        /** @brief Constructor. */
        TrigDvFex(const std::string& name, ISvcLocator * isvclocator);
        /** @brief Destructor. */
        ~TrigDvFex();
    
        HLT::ErrorCode hltInitialize();
        HLT::ErrorCode hltExecute(const HLT::TriggerElement* input, HLT::TriggerElement* output);
        HLT::ErrorCode hltFinalize(); 
        // JDC 
        HLT::ErrorCode hltEndEvent();
    
    private:
        std::string m_jetKey;
 
        /** @brief Enable monitoring histograms (online) */
        bool m_mon_online;
        /** @brief Enable monitoring histograms (validation) */
        bool m_mon_validation;

        // JDC TESTING --> 
        ServiceHandle<IRegSelSvc>  m_regionSelector;
        //! hits-related info
        const InDet::SiClusterContainer * m_pixcontainer_v;
        const InDet::SiClusterContainer * m_sctcontainer_v;
        const InDet::TRT_DriftCircleContainer * m_trtcontainer_v;

        bool m_is_pix_filled;
        bool m_is_sct_filled;
        bool m_is_trt_filled;
        bool m_are_prds_filled;
    

        /** @brief To retrieve track collections reconstructed at EF and stored in TrackParticleContainer. */
        HLT::ErrorCode getTrackCollection(const xAOD::TrackParticleContainer*&, const HLT::TriggerElement*);
        /** @brief Retrieved the xAOD::Jet collection attached to the TE. */
        HLT::ErrorCode getJet(const xAOD::Jet * &,const HLT::TriggerElement*);

        /** @brief Helper function to define the azimuthal angle between [-PI,PI] */
        float phiCorr(const float & phi) const;
    
        /** @brief To select EF tracks. */
        bool efTrackSel(const xAOD::TrackParticle*&, unsigned int);

        /** @brief Function retrieving the Hits (PrepRawData) found in the inner detector at a given RoI 
         * and also return the element index of the first element of the SCT and TRT, PIXEL is always 0*/
        /*std::vector<const Trk::PrepRawData*> getPRDsfromID(const TrigRoiDescriptor * roi,
                std::pair<unsigned int,unsigned int> & sct_trt_indices);*/
        PrdsAtRoI getPRDsfromID(const TrigRoiDescriptor * roi);
        /** @brief update the Prepared raw data collection for the given event */
        HLT::ErrorCode updatePRDs(const TrigDvFexDet & detectortype);
        
        /** @brief update the Prepared raw data collection for the given track */
        void updateTrackPRDs(const xAOD::TrackParticle * track, std::unordered_set<const Trk::PrepRawData*> & prds);

        /** @brief Pointer to TrigEFBjet collection. */
        TrigEFBjetContainer* m_trigEFBjetColl;
    
        
        /** @brief Internal EDM class for jet direction information, defined in TrigBjetDataHelper. */ 
        TrigBjetJetInfo*    m_trigBjetJetInfo;
        /** @brief Internal EDM class for track parameter information, defined in TrigBjetDataHelper. */ 
        std::vector<TrigBjetTrackInfo>* m_trigBjetTrackInfoVector;
   
        /** @brief check variable to control PV retrieval */
        bool m_setPVInfo;

        /** @brief track reconstruction algorithm. */
        int m_algo;
        /** @brief switch to perform track-RoI eta/phi matching when selecting reconstructed tracks. */
        bool m_useEtaPhiTrackSel;

        /** @brief lower bound of the chi square of the reconstructed track (to perform track selection). */
        float m_trkSelChi2;
        /** @brief lower bound of the number of hits on b-layer of reconstructed track (to perform track selection). */
        int   m_trkSelBLayer;
        /** @brief lower bound of the number of hits in pixel detector of reconstructed track 
         * (to perform track selection). */
        int   m_trkSelPixHits;
        /** @brief lower bound of the number of hits in silicon detectors of reconstructed track 
         * (to perform track selection). */
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
        std::vector<float> m_mon_trk_z0, m_mon_trk_z0_sel, m_mon_trk_Sz0_sel;
        std::vector<float> m_mon_trk_theta, m_mon_trk_pt, m_mon_trk_qOverPt, m_mon_trk_eta;
        std::vector<float> m_mon_trkSel_theta, m_mon_trkSel_pt, m_mon_trkSel_qOverPt, m_mon_trkSel_eta;
        /** @brief quality track */
        std::vector<int> m_mon_trk_BlayerHits, m_mon_trk_PixHits, m_mon_trk_SCTHits, m_mon_trk_TRTHits;
        std::vector<int> m_mon_trk_PixHoles, m_mon_trk_SCTHoles, m_mon_trk_TRTHoles;
        std::vector<int> m_mon_trk_PixSharedHits, m_mon_trk_SCTSharedHits;//, m_mon_trk_TRTSharedHits;

        /** @brief quality track passing d0 cut */
        std::vector<int> m_mon_trkSel_BlayerHits, m_mon_trkSel_PixHits, m_mon_trkSel_SCTHits, m_mon_trkSel_TRTHits;
        std::vector<int> m_mon_trkSel_PixHoles, m_mon_trkSel_SCTHoles, m_mon_trkSel_TRTHoles;
        std::vector<int> m_mon_trkSel_PixSharedHits, m_mon_trkSel_SCTSharedHits;//, m_mon_trkSel_TRTSharedHits;

    
        /** @brief To retrieve selected tracks in percentage. */ 
        inline float totSelectedTracks() const 
        {
            if(!m_totTracks)
            {
                return -0.1;
            }
            else
            {
                return (float)m_totSelTracks/(float)m_totTracks;
            }
        };

        /** @brief Eta of the RoI */
        float m_etaRoI;
        /** @brief RoI phi */
        float m_phiRoI;

        /** @brief Number of reconstructed tracks per RoI. */
        float m_totTracks;
        /** @brief Number of selected tracks per RoI. */
        float m_totSelTracks;

        /** @brief number of total unused hits **/
        float m_unusedhits;
        /** @brief fraction of unused hits from the total hits **/
        float m_unusedhits_fraction;

        /** @brief Delta eta between the LVL1 jet RoI and the EF track-jet. */
        float m_deltaEtaTrk;
        /** @brief Delta phi between the LVL1 jet RoI and the EF track-jet. */
        float m_deltaPhiTrk;
};

#endif
