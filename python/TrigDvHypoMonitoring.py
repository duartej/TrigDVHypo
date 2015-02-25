from TrigMonitorBase.TrigGenericMonitoringToolConfig import defineHistogram, TrigGenericMonitoringToolConfig

class TrigDvHypoMonitoring(TrigGenericMonitoringToolConfig):
    def __init__ (self, name="TrigDvHypoMonitoring"):
        super(TrigDvHypoMonitoring, self).__init__(name)
        # ----> remove this line and change labels and ranges  
        self.Histograms += [ defineHistogram('CutCounter', type='TH1F', title="Bjet Hypo cut counter",
            xbins=6, xmin=-1.5, xmax=4.5,labels='No SV- obj:No usable beam spot:Rejected (nTrackVtx):'+\
                    'Rejected (Mass):Rejected (decay length):Accepted') ]


class TrigEFDvHypoValidationMonitoring(TrigDvHypoMonitoring):
    def __init__ (self, name="TrigDvEFBjetHypoValidationMonitoring"):
        super(TrigEFDvHypoValidationMonitoring, self).__init__(name)

        self.defineTarget("Validation")
        self.Histograms += [ defineHistogram('nTracksperRoI', 
                                            type='TH1F', 
                                            title="Number of selected tracks per RoI",
                                            xbins=16, xmin=-0.5, xmax=15) 
                                            ]



class TrigEFDvHypoOnlineMonitoring(TrigDvHypoMonitoring):
    def __init__ (self, name="TrigDvEFBjetHypoOnlineMonitoring"):
        super(TrigEFDvHypoOnlineMonitoring, self).__init__(name)

        self.defineTarget("Online")


        


