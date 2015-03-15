from TrigMonitorBase.TrigGenericMonitoringToolConfig import defineHistogram, TrigGenericMonitoringToolConfig

class TrigDvFexMonitoring(TrigGenericMonitoringToolConfig):
    def __init__ (self, name="TrigDvFexMonitoring"):
        super(TrigDvFexMonitoring, self).__init__(name)

        self.Histograms += [ defineHistogram('trk_z0_sel', type='TH1F', 
                title="longitudinal IP for selected tracks;z_{0} [mm]",xbins=301, xmin=-300, xmax=300) ]
        self.Histograms += [ defineHistogram('trk_Sz0_sel', type='TH1F', 
                title="longitudinal IP significance for selected tracks",xbins=301, xmin=-300.0, xmax=300.0) ]
        
        self.Histograms += [ defineHistogram('trk_a0_sel', type='TH1F', 
            title="transverse IP for selected tracks;d_{0} [mm]",xbins=101, xmin=-10, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_Sa0_sel', type='TH1F', 
                title="transverse IP significance for selected tracks",xbins=401, xmin=-1500.0, xmax=1500.0) ]
        
        self.Histograms += [ defineHistogram('roi_nTracks', type='TH1F', 
                title="Number of tracks per RoI before track selection",xbins=50, xmin=-0.5, xmax=49.5) ]
        self.Histograms += [ defineHistogram('roi_nTracks_sel', type='TH1F', 
                title="Number of tracks per RoI after track selection",xbins=50, xmin=-0.5, xmax=49.5) ]
        self.Histograms += [ defineHistogram('roi_selectedTracks', type='TH1F', 
            title="Selected tracks per RoI (in percentage). Note that -1 means no "\
                    "tracks found in the RoI",xbins=20, xmin=-0.2, xmax=1.1) ]
    
        self.Histograms += [ defineHistogram('roi_deltaEtaTrk', type='TH1F',
            title="Delta eta between the RoI and the HLT track-jet",xbins=100, xmin=-1.0, xmax=1.0) ]
        self.Histograms += [ defineHistogram('roi_deltaPhiTrk', type='TH1F',
            title="Delta phi between the RoI and the HLT track-jet",xbins=100, xmin=-1.0, xmax=1.0) ]
        
class TrigEFDvFexValidationMonitoring(TrigDvFexMonitoring):
    def __init__ (self, name="TrigEFDvFexValidationMonitoring"):
        super(TrigEFDvFexValidationMonitoring, self).__init__(name)

        self.defineTarget("Validation")
        
        # Tracks
        self.Histograms += [ defineHistogram('trk_a0', type='TH1F',
            title="transverse IP",xbins=101, xmin=-10, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_z0', type='TH1F',
            title="longitudinal IP",xbins=301, xmin=-300, xmax=300) ]
        
        
        self.Histograms += [ defineHistogram('trk_theta', type='TH1F', 
                title="#theta of the tracks at the perigee;#theta",xbins=100, xmin=0.0, xmax=3.1415) ]
        self.Histograms += [ defineHistogram('trk_pt', type='TH1F', 
                title="p_{T} of the tracks at the perigee;p_{T} [GeV]",xbins=201, xmin=0.0, xmax=200.0) ]
        self.Histograms += [ defineHistogram('trk_qOverPt', type='TH1F', 
                title="signed curvature of the tracks at the perigee;q/p_{T} [1.0/MeV]",xbins=200, xmin=-0.004, xmax=0.004) ]
        self.Histograms += [ defineHistogram('trk_eta', type='TH1F', 
                title="#eta of the tracks at the perigee;#eta",xbins=100, xmin=-2.6, xmax=2.6) ]
        
        self.Histograms += [ defineHistogram('trk_BlayerHits', type='TH1I', 
                title="Tracks hits at the B-layer (first pixel layer)",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_PixHits', type='TH1I', 
                title="Tracks hits at the pixel)",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_SCTHits', type='TH1I', 
                title="Tracks hits at the SCT",xbins=16, xmin=-0.5, xmax=15) ]
        self.Histograms += [ defineHistogram('trk_TRTHits', type='TH1I', 
                title="Tracks hits at the TRT",xbins=51, xmin=-0.5, xmax=50) ]
        
        self.Histograms += [ defineHistogram('trk_PixHoles', type='TH1I', 
                title="Number of pixel layers with absence of hits (pixel))",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_SCTHoles', type='TH1I', 
                title="Number of SCT layers with absence of hits",xbins=16, xmin=-0.5, xmax=15) ]
        self.Histograms += [ defineHistogram('trk_TRTHoles', type='TH1I', 
                title="Number of TRT holes",xbins=51, xmin=-0.5, xmax=50) ]
        
        self.Histograms += [ defineHistogram('trk_PixSharedHits', type='TH1I', 
                title="Number of pixel hits used by more than one track",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_SCTSharedHits', type='TH1I', 
                title="Number of SCT hits used by more than one track",xbins=16, xmin=-0.5, xmax=15) ]
        #self.Histograms += [ defineHistogram('trk_TRTSharedHits', type='TH1I', 
        #        title="Number of TRT hits used by more than one track",xbins=51, xmin=-0.5, xmax=50) ]
        
        self.Histograms += [ defineHistogram('trk_theta_sel', type='TH1F', 
                title="#theta of the tracks at the perigee;#theta",xbins=100, xmin=0.0, xmax=3.1415) ]
        self.Histograms += [ defineHistogram('trk_pt_sel', type='TH1F', 
                title="p_{T} of the tracks at the perigee;p_{T} [GeV]",xbins=201, xmin=0.0, xmax=200.0) ]
        self.Histograms += [ defineHistogram('trk_qOverPt_sel', type='TH1F', 
                title="signed curvature of the tracks at the perigee;q/p_{T} [1.0/MeV]",xbins=200, xmin=-0.004, xmax=0.004) ]
        self.Histograms += [ defineHistogram('trk_eta_sel', type='TH1F', 
                title="#eta of the tracks at the perigee;#eta",xbins=100, xmin=-4.5, xmax=4.5) ]
        
        self.Histograms += [ defineHistogram('trk_BlayerHits_sel', type='TH1I', 
                title="Tracks hits at the B-layer (first pixel layer)",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_PixHits_sel', type='TH1I', 
                title="Tracks hits at the pixel)",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_SCTHits_sel', type='TH1I', 
                title="Tracks hits at the SCT",xbins=16, xmin=-0.5, xmax=15) ]
        self.Histograms += [ defineHistogram('trk_TRTHits_sel', type='TH1I', 
                title="Tracks hits at the TRT",xbins=51, xmin=-0.5, xmax=50) ]
        
        self.Histograms += [ defineHistogram('trk_PixHoles_sel', type='TH1I', 
                title="Number of pixel layers with absence of hits (pixel))",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_SCTHoles_sel', type='TH1I', 
                title="Number of SCT layers with absence of hits",xbins=16, xmin=-0.5, xmax=15) ]
        self.Histograms += [ defineHistogram('trk_TRTHoles_sel', type='TH1I', 
                title="Number of TRT holes",xbins=51, xmin=-0.5, xmax=50) ]
        
        self.Histograms += [ defineHistogram('trk_PixSharedHits_sel', type='TH1I', 
                title="Number of pixel hits used by more than one track",xbins=11, xmin=-0.5, xmax=10) ]
        self.Histograms += [ defineHistogram('trk_SCTSharedHits_sel', type='TH1I', 
                title="Number of SCT hits used by more than one track",xbins=16, xmin=-0.5, xmax=15) ]
        #self.Histograms += [ defineHistogram('trk_TRTSharedHits_sel', type='TH1I', 
        #        title="Number of TRT hits used by more than one track",xbins=51, xmin=-0.5, xmax=50) ]

        self.Histograms += [ defineHistogram('roi_stepsToSelect',
            type='TH1F', title="Steps to select tracks",xbins=12, xmin=0.0, xmax=12,
            labels='BS flag status:BS width:eta matching:phi matching:pT cut:d0 cut:z0'+
            ' cut:b-layer hit cut:pixel hit cut:silicon hit cut:chi2 cut:selected') ]
        
        self.Histograms += [ defineHistogram('roi_deltaEtaTrk,roi_deltaPhiTrk',
                                              type='TH2F',
                                              title="#Delta#eta vs. #Delta#phi between the RoI and the HLT-track",
                                              xbins=100,xmin=-1.0,xmax=1.0,
                                              ybins=100,ymin=-1.0,ymax=1.0)
                                              ]




class TrigEFDvFexOnlineMonitoring(TrigDvFexMonitoring):
    def __init__ (self, name="TrigEFDvFexOnlineMonitoring"):
        super(TrigEFDvFexOnlineMonitoring, self).__init__(name)

        self.defineTarget("Online")
        
        self.Histograms += [ defineHistogram('roi_stepsToSelect', type='TH1F',
            title="Steps to select tracks",xbins=12, xmin=0.0, xmax=12,
            labels='BS flag status:BS width:eta matching:phi matching:pT cut:d0 cut:z0'+
            ' cut:b-layer hit cut:pixel hit cut:silicon hit cut:chi2 cut:selected') ]


