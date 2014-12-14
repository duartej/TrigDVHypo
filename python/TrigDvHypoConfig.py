from TrigBjetHypo.TrigBjetHypoConf import TrigDvHypo

from AthenaCommon.Logging import logging
from AthenaCommon.SystemOfUnits import GeV

def getDvHypoInstance( instance ):
    return DvHypo( instance=instance, name=instance+"DvHypo")

def getDvHypoNoCutInstance( instance):
    return DvHypoNoCut( instance=instance, name=instance+"DvHypoNoCut" )

def getDvHypoSplitInstance( instance ):
    return DvHypoSplit( instance=instance, name=instance+"DvHypoSplit")

def getDvHypoSplitNoCutInstance( instance):
    return DvHypoSplitNoCut( instance=instance, name=instance+"DvHypoSplitNoCut" )

# FIXME::: Include the cuts in here... or not

class DvHypo (TrigDvHypo):
    __slots__ = []
    
    def __init__(self, instance, name):
        super( DvHypo, self ).__init__( name )
        
        mlog = logging.getLogger('DvHypoConfig.py')
        
        AllowedInstances = ["EF"]
        
        if instance not in AllowedInstances :
            mlog.error("Instance "+instance+" is not supported!")
            return None
        

        if instance=="EF" :
            self.AcceptAll = False
            self.Instance  = "EF"
            self.UseBeamSpotFlag = False

        self.JetKey = "SplitJet"  # FIXME: EVERYTHING SPLIT
        
        from TrigTimeMonitor.TrigTimeHistToolConfig import TrigTimeHistToolConfig
        time = TrigTimeHistToolConfig("TimeHistogramForTrigDvHypo")
        time.TimerHistLimits = [0,0.4]
        
        self.AthenaMonTools = [ time ]
        
        if instance=="EF" :
            from TrigBjetHypo.TrigDvHypoMonitoring import TrigEFDvHypoValidationMonitoring, TrigEFDvHypoOnlineMonitoring
            validation = TrigEFDvHypoValidationMonitoring()
            online     = TrigEFDvHypoOnlineMonitoring()
            self.AthenaMonTools += [validation,online]
            

class DvHypoNoCut (TrigDvHypo):
    __slots__ = []
    
    def __init__(self, instance, name):
        super( DvHypoNoCut, self ).__init__( name )
        
        mlog = logging.getLogger('DvHypoConfig.py')
                
        AllowedInstances = ["EF"]

        self.JetKey = "SplitJet"  # FIXME EVERY THING SPLIT!!
        
        if instance in AllowedInstances :
            from TrigTimeMonitor.TrigTimeHistToolConfig import TrigTimeHistToolConfig
            time = TrigTimeHistToolConfig("TimeHistogramForTrigDvHypo")
            time.TimerHistLimits = [0,0.4]

            self.AthenaMonTools = [ time ]
            
            if instance=="EF" :
                self.AcceptAll = True
                self.Instance  = "EF"
                self.UseBeamSpotFlag = False
                from TrigBjetHypo.TrigDvHypoMonitoring import TrigEFDvHypoValidationMonitoring, TrigEFDvHypoOnlineMonitoring
                validation = TrigEFDvHypoValidationMonitoring()
                online     = TrigEFDvHypoOnlineMonitoring()
                self.AthenaMonTools += [ validation, online ]
        
        else :
            mlog.error("Instance "+instance+" is not supported!")
            return None


### Split instances


class DvHypoSplit (TrigDvHypo):
    __slots__ = []
    
    def __init__(self, instance, name):
        super( DvHypoSplit, self ).__init__( name )
        
        mlog = logging.getLogger('DvHypoConfig.py')
        
        AllowedInstances = ["EF"]
        
        if instance not in AllowedInstances :
            mlog.error("Instance "+instance+" is not supported!")
            return None
        
        if instance=="EF" :
            self.AcceptAll = False
            self.Instance  = "EF"
            self.UseBeamSpotFlag = False

        self.JetKey = "SplitJet"
        
        from TrigTimeMonitor.TrigTimeHistToolConfig import TrigTimeHistToolConfig
        time = TrigTimeHistToolConfig("TimeHistogramForTrigDvHypo")
        time.TimerHistLimits = [0,0.4]
        
        self.AthenaMonTools = [ time ]

        if instance=="EF" :
            from TrigBjetHypo.TrigDvHypoMonitoring import TrigEFDvHypoValidationMonitoring, TrigEFDvHypoOnlineMonitoring
            validation = TrigEFDvHypoValidationMonitoring()
            online     = TrigEFDvHypoOnlineMonitoring()
            self.AthenaMonTools += [ validation, online ]
        
            


class DvHypoSplitNoCut (TrigDvHypo):
    __slots__ = []
    
    def __init__(self, instance, name):
        super( DvHypoSplitNoCut, self ).__init__( name )
        
        mlog = logging.getLogger('DvHypoConfig.py')
                
        AllowedInstances = ["EF"]

        self.JetKey = "SplitJet"
        
        if instance in AllowedInstances :
            from TrigTimeMonitor.TrigTimeHistToolConfig import TrigTimeHistToolConfig
            time = TrigTimeHistToolConfig("TimeHistogramForTrigDvHypo")
            time.TimerHistLimits = [0,0.4]

            self.AthenaMonTools = [ time ]
            
            if instance=="EF" :
                self.AcceptAll = True
                self.Instance  = "EF"
                self.UseBeamSpotFlag = False
                from TrigBjetHypo.TrigDvHypoMonitoring import TrigEFDvHypoValidationMonitoring, TrigEFDvHypoOnlineMonitoring
                validation = TrigEFDvHypoValidationMonitoring()
                online     = TrigEFDvHypoOnlineMonitoring()
                self.AthenaMonTools += [validation,online]        

        else :
            mlog.error("Instance "+instance+" is not supported!")
            return None


        

