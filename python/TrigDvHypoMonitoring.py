from TrigMonitorBase.TrigGenericMonitoringToolConfig import defineHistogram, TrigGenericMonitoringToolConfig

class TrigDvHypoMonitoring(TrigGenericMonitoringToolConfig):
    def __init__ (self, name="TrigDvHypoMonitoring"):
        super(TrigDvHypoMonitoring, self).__init__(name)
        # ----> remove this line and change labels and ranges  
        self.Histograms = [ defineHistogram('CutCounter', type='TH1F', title="Bjet Hypo cut counter",
            xbins=4, xmin=-1.5, xmax=2.5,labels='No Bjet obj:No usable beam spot:Rejected:Accepted') ]


class TrigEFDvHypoValidationMonitoring(TrigDvHypoMonitoring):
    def __init__ (self, name="TrigDvEFBjetHypoValidationMonitoring"):
        super(TrigEFDvHypoValidationMonitoring, self).__init__(name)

        self.defineTarget("Validation")



class TrigEFDvHypoOnlineMonitoring(TrigDvHypoMonitoring):
    def __init__ (self, name="TrigDvEFBjetHypoOnlineMonitoring"):
        super(TrigEFDvHypoOnlineMonitoring, self).__init__(name)

        self.defineTarget("Online")


        


