from TrigMonitorBase.TrigGenericMonitoringToolConfig import defineHistogram, TrigGenericMonitoringToolConfig

class TrigDvFexMonitoring(TrigGenericMonitoringToolConfig):
    def __init__ (self, name="TrigDvFexMonitoring"):
        super(TrigDvFexMonitoring, self).__init__(name)

        self.Histograms += [ defineHistogram('trk_a0_sel', type='TH1F', 
                title="transverse IP for selected tracks",xbins=200, xmin=-300, xmax=300) ]
        self.Histograms += [ defineHistogram('trk_S(a0)_sel', type='TH1F', 
                title="transverse IP significance for selected tracks",xbins=400, xmin=-20.0, xmax=20.0) ]

        self.Histograms += [ defineHistogram('trk_z0_sel', type='TH1F', 
                title="longitudinal IP for selected tracks",xbins=200, xmin=-400, xmax=400) ]
        self.Histograms += [ defineHistogram('trk_z0_sel_PV', type='TH1F', 
                title="longitudinal IP w.r.t. primary vertex for selected tracks",xbins=200, xmin=-300, xmax=300) ]
        self.Histograms += [ defineHistogram('trk_S(z0)_sel', type='TH1F', 
                title="longitudinal IP significance for selected tracks",xbins=400, xmin=-20.0, xmax=20.0) ]
        self.Histograms += [ defineHistogram('trk_prob', type='TH1F', 
                title="track probability estimated by JetProb",xbins=40, xmin=0.0, xmax=1.0) ]
        self.Histograms += [ defineHistogram('roi_nTracks', type='TH1F', 
                title="Number of tracks per RoI before track selection",xbins=50, xmin=-0.5, xmax=49.5) ]
        self.Histograms += [ defineHistogram('roi_nTracks_sel', type='TH1F', 
                title="Number of tracks per RoI after track selection",xbins=50, xmin=-0.5, xmax=49.5) ]
        self.Histograms += [ defineHistogram('roi_selectedTracks', type='TH1F', 
            title="Selected tracks for CHI2 tagger (in percentage)",xbins=20, xmin=-0.2, xmax=1.1) ]
    
        self.Histograms += [ defineHistogram('roi_deltaEtaJet', type='TH1F', 
                title="Delta eta between the LVL1 jet RoI and the HLT jet",xbins=40, xmin=-1.0, xmax=1.0) ]
        self.Histograms += [ defineHistogram('roi_deltaPhiJet', type='TH1F', 
                title="Delta phi between the LVL1 jet RoI and the HLT jet",xbins=40, xmin=-1.0, xmax=1.0) ]
        self.Histograms += [ defineHistogram('roi_deltaEtaTrkJet', type='TH1F', 
        title="Delta eta between the LVL1 jet RoI and the HLT track-jet",xbins=40, xmin=-1.0, xmax=1.0) ]
        self.Histograms += [ defineHistogram('roi_deltaPhiTrkJet', type='TH1F',
            title="Delta phi between the LVL1 jet RoI and the HLT track-jet",xbins=40, xmin=-1.0, xmax=1.0) ]
        self.Histograms += [ defineHistogram('roi_deltaEtaJetTrkJet', type='TH1F',
            title="Delta eta between the HLT jet and the HLT track-jet",xbins=40, xmin=-1.0, xmax=1.0) ]
        self.Histograms += [ defineHistogram('roi_deltaPhiJetTrkJet', type='TH1F',
            title="Delta phi between the HLT jet and the HLT track-jet",xbins=40, xmin=-1.0, xmax=1.0) ]



class TrigEFDvFexValidationMonitoring(TrigDvFexMonitoring):
    def __init__ (self, name="TrigEFDvFexValidationMonitoring"):
        super(TrigEFDvFexValidationMonitoring, self).__init__(name)

        self.defineTarget("Validation")
        
        self.Histograms += [ defineHistogram('trk_a0', type='TH1F',
            title="transverse IP",xbins=200, xmin=-300, xmax=300) ]
        self.Histograms += [ defineHistogram('trk_z0', type='TH1F',
            title="longitudinal IP",xbins=200, xmin=-400, xmax=400) ]

        self.Histograms += [ defineHistogram('roi_stepsToSelect',
            type='TH1F', title="Steps to select tracks",xbins=12, xmin=0.0, xmax=12,
            labels='BS flag status:BS width:eta matching:phi matching:pT cut:d0 cut:z0'+
            ' cut:b-layer hit cut:pixel hit cut:silicon hit cut:chi2 cut:selected') ]



class TrigEFDvFexOnlineMonitoring(TrigDvFexMonitoring):
    def __init__ (self, name="TrigEFDvFexOnlineMonitoring"):
        super(TrigEFDvFexOnlineMonitoring, self).__init__(name)

        self.defineTarget("Online")
        
        self.Histograms += [ defineHistogram('roi_stepsToSelect', type='TH1F',
            title="Steps to select tracks for CHI2 tagger",xbins=12, xmin=0.0, xmax=12,
            labels='BS flag status:BS width:eta matching:phi matching:pT cut:d0 cut:z0'+
            ' cut:b-layer hit cut:pixel hit cut:silicon hit cut:chi2 cut:selected') ]


