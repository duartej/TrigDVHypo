from TrigBjetHypo.TrigBjetHypoConf import TrigDvFex

from AthenaCommon.Logging import logging
from AthenaCommon.SystemOfUnits import mm, GeV

def getDvFexInstance( instance, algo):
    return DvFex( instance, algo=algo, name="EFDvFex_"+algo )

class DvFex(TrigDvFex):
    __slots__ = []
    
    def __init__(self, instance, algo, name):
        super( DvFex, self ).__init__( name )
        
        mlog = logging.getLogger('BjetHypoConfig.py')  ## SURE??
        
        AllowedInstances = ["EF"]
        AllowedAlgos     = ["EFID"]
        
	# Pre-configuring
        if instance not in AllowedInstances :
            mlog.error("Instance "+instance+" is not supported!")
            return None
        
        self.JetKey = ""
        self.PriVtxKey = "EFHistoPrmVtx"
        
        self.AlgoId = None
        
	if instance=="EF" and algo=="EFID":
	    self.AlgoId = 1
        
        if not self.AlgoId:
	    mlog.error("AlgoId is wrongly set!")
            return None
        
        if instance=="EF" :
            self.Instance = "EF"
        
        self.UseBeamSpotFlag    = False
        self.SetBeamSpotWidth   = 1.0*mm
        
        self.UseEtaPhiTrackSel  = False
        
        if algo=="EFID" :
            self.TrkSel_Chi2    = 0.0
            self.TrkSel_BLayer  = 1
            self.TrkSel_PixHits = 2
            self.TrkSel_SiHits  = 4
            self.TrkSel_D0      = 300.0*mm
            self.TrkSel_Z0      = 300.0*mm
            self.TrkSel_Pt      = 4.0*GeV

        from TrigTimeMonitor.TrigTimeHistToolConfig import TrigTimeHistToolConfig
        time = TrigTimeHistToolConfig("TimeHistogramForTrigBjetHypo")
        time.TimerHistLimits = [0,2]
        
        self.AthenaMonTools = [ time ]
        
	if instance=="EF" :
            from TrigBjetHypo.TrigDvFexMonitoring import TrigEFDvFexValidationMonitoring, TrigEFDvFexOnlineMonitoring
            validation = TrigEFDvFexValidationMonitoring()
            online     = TrigEFDvFexOnlineMonitoring()    
	    self.AthenaMonTools+=[validation],online]



