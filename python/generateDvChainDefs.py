
__author__  = 'J.Duarte-Campderros (based in code from M.Backes, C.Bernius)'
__doc__="Definition of displaced-vertex chains" 

from AthenaCommon.Logging import logging
logging.getLogger().info("Importing %s",__name__)
logDvDef = logging.getLogger("TriggerMenu.bjet.generateDvChainDefs.py")

from TriggerMenu.menu.MenuUtils import *

import pprint
pp = pprint.PrettyPrinter(indent=4, depth=8)


def __rechaindict__(c):
    """To be deprecated, just for the testing procedure
    """
    from TriggerMenu.menu.DictFromChainName import DictFromChainName
    dfcn = DictFromChainName()

    pl1 = []
    for pch in c['chainParts']:
        pl1.append(pch['L1item'])

    newname = c['chainName'].replace('dv_','').replace('TestChain','j')
    nchlist = [ newname ,c['chainCounter'],c['L1item'],pl1,c['stream'],
            c['groups'],c['EBstep'] ]
    
    return dfcn.getChainDict(nchlist)


###########################################################################
###########################################################################
# Note that the chainDict can be generated from the ChainName list (see in 
# TriggerCommon/TriggerMenu/phython/menu/MC_pp_v5.py you'll see a lot of lists
# You can use the class DictFromChainName in 
# TriggerCommon/TriggerMenu/phython/menu/DictFromChainName.py to obtain the
# chain dict--> you should modify the analysisShortName method or hardcoded 
# the chainDict in order to accept a new chain name... By the moment I'm useing
# the jX_bperf_split as handlers to our triggers!
def generateChainDefs(__chainDict):
    # TESTING
    logDvDef.warning("TESTING DV-Triggers: using Test Slice!")
    logDvDef.warning("TESTING DV-Triggers: Re-building the chainDict from '%s'" 
            % __chainDict['chainName'])
    logDvDef.warning("TESTING DV-Triggers: This is a provisional behaviour while testing")
    chainDict = __rechaindict__(__chainDict)
    # TESTING

    from copy import deepcopy
    chainDict_orig = deepcopy(chainDict)
    logDvDef.debug("Initial Displaced-vertex chainDict: \n %s" % (pp.pformat(chainDict)))

    #----------------------------------------------------------------------------
    # --- preparing the dictionary for jet chain generation ---
    n_chainParts = len(chainDict['chainParts'])
   
    if(n_chainParts >1):
        jetchainDict = _prepareJetChainDict(chainDict)
    elif (n_chainParts == 1):
        logDvDef.debug('Only one chain part dictionary for this chain.')
        jetchainDict = chainDict
    #if(n_chainParts > 0):
    #    jetchainDict = chainDict
    else:
        raise RuntimeError('No chain parts found for chain %s' % (chainDict['chainName']))
    #----------------------------------------------------------------------------

    #----------------------------------------------------------------------------
    # --- obtain jet chain def used as inputTE for first   sequence ---
    # ----> This can be substitute depending the input TE RoI (muon,egamma)
    from TriggerMenu.jet.generateJetChainDefs import generateChainDefs as genJetChainDefs
    thejThrChainDef = get_j35_ChainDef()
    theBJChainDef =  genJetChainDefs(thejThrChainDef)
    logDvDef.debug("ChainDef for dv-jet chain: \n %s" % (str(theBJChainDef)))
    ##----------------------------------------------------------------------------
    # ----> the proper chain
    primaryobjectchainname = jetchainDict['chainName'].replace('dv_','').replace('TestChain','j')
    # -- new dictionary for pure jet (muon,egamma) object trigger
    pochaindict = deepcopy(chainDict)
    pochaindict['chainName'] = primaryobjectchainname
    theAllJetChainDef =  genJetChainDefs(pochaindict)
    logDvDef.debug("Jet Chaindef for dv-jet chain: \n %s" % (str(theAllJetChainDef)))
    ##----------------------------------------------------------------------------

    ##----------------------------------------------------------------------------
    ## --- now merge chain defs so that last TE is from jet serving as input to SV seq ---
    ## DO NOT CHANGE THE MERGING ORDER PLEASE, OTHERWISE WRONG INPUTTE IS SELECTED
    theJetChainDef = mergeChainDefs([theAllJetChainDef, theBJChainDef], 'serial', -1, False)
    logDvDef.debug("Merged chainDef for dv-jet chain: \n %s" % (str(theJetChainDef)))
    ##----------------------------------------------------------------------------

    ##----------------------------------------------------------------------------
    ## --- filtering b-tagged jets chainParts ---
    #listofChainDicts = splitChainDict(chainDict_orig)
    #logDvDef.debug("Split dv-jet chainDict: \n %s" % (pp.pformat(listofChainDicts)))
    
    
    #theJetChainDef =  genJetChainDefs(pochaindict)
    #logDvDef.debug("ChainDef for dv-jet chain: \n %s" % (str(theJetChainDef)))

    # --------------------------------------------------------------------------
    # --- adding to the jet (muon/egamma) sequence, the dv-relateda
    listofChainDicts = splitChainDict(chainDict_orig)
    nchaindict = len(listofChainDicts)
    theListOfChainDefs = []
    for subChainDict in listofChainDicts:
        theJetChainDef1 = deepcopy(theJetChainDef)
        theDvChainDef = buildDvChains(theJetChainDef1, subChainDict, True, nchaindict)    
        theListOfChainDefs += [theDvChainDef] 

    logDvDef.debug("----------------- Beginning of final individual chainDefs for dv-jet chains printout -----------------")
    for chainDef in theListOfChainDefs:
        logDvDef.debug(str(chainDef))        
    logDvDef.debug("----------------- End of final individual chainDefs for dv-jet chains printout -----------------")
    
    if len(theListOfChainDefs)>1:
        theFinalChainDef = mergeChainDefs(theListOfChainDefs,strategy="parallel",\
                offset=-1,preserveL2EFOrder=True,removeDuplicateTEs=True)
    else:
        theFinalChainDef = theListOfChainDefs[0]
        
    logDvDef.debug("----------------- Beginning of final merged chainDefs for b-jet chains printout -----------------")
    logDvDef.debug(str(theFinalChainDef))
    logDvDef.debug("----------------- End of final merged chainDefs for b-jet chains printout -----------------")

    return theFinalChainDef
    
    #----------------------------------------------------------------------------
    # --- obtain jet chain def used as inputTE for first btag sequence ---
    # ----> This can be substitute depending the input TE RoI (muon,egamma)
    #thejThrChainDef = get_j35_ChainDef()
    #if ('j0_bperf' in chainDict['chainName']):
    #    thejThrChainDef['L1item'] = 'L1_J10'
    #    thejThrChainDef['chainName'] = 'j0'
    #    for part in thejThrChainDef['chainParts']:
    #        part['chainPartName'] = 'j0'
    #        part['threshold'] = '0'
    #    if ('L1MU10' in chainDict['chainName']):
    #        thejThrChainDef['L1item'] = 'L1_MU10'
    #    if ('L1RD0_EMPTY' in chainDict['chainName']):
    #        thejThrChainDef['L1item'] = 'L1_RD0_EMPTY'

    #from TriggerMenu.jet.generateJetChainDefs import generateChainDefs as genJetChainDefs
    #theBJChainDef =  genJetChainDefs(thejThrChainDef)
    #logDvDef.debug("ChainDef for dv-jet chain: \n %s" % (str(theBJChainDef)))
    ##----------------------------------------------------------------------------


    ##----------------------------------------------------------------------------
    ## --- build the jet chain, then pass JetChainDef and bjetchainDictionaries to build bjet chains ---
    #theAllJetChainDef =  genJetChainDefs(jetchainDict)
    #logDvDef.debug("Jet Chaindef for dv-jet chain: \n %s" % (str(theAllJetChainDef)))
    ##----------------------------------------------------------------------------

    ##----------------------------------------------------------------------------
    ## --- now merge chain defs so that last TE is from jet serving as input to bjet seq ---
    ## DO NOT CHANGE THE MERGING ORDER PLEASE, OTHERWISE WRONG INPUTTE IS SELECTED
    #theJetChainDef = mergeChainDefs([theAllJetChainDef, theBJChainDef], 'serial', -1, False)
    #logDvDef.debug("Merged chainDef for dv-jet chain: \n %s" % (str(theJetChainDef)))
    ##----------------------------------------------------------------------------

    ##----------------------------------------------------------------------------
    ## --- filtering b-tagged jets chainParts ---
    #listofChainDicts = splitChainDict(chainDict_orig)
    #logDvDef.debug("Split dv-jet chainDict: \n %s" % (pp.pformat(listofChainDicts)))

#    bjetchainDicts = [cdict for cdict in listofChainDicts if cdict['chainParts']['bTag']] 
#    logDvDef.debug("Final dv-jet chainDict: \n %s" % (pp.pformat(bjetchainDicts)))
#
#    theListOfChainDefs = []
#    for subChainDict in bjetchainDicts:
#        theJetChainDef1 = deepcopy(theJetChainDef)
#        theDvChainDef = buildDvChains(theJetChainDef1, subChainDict, True, len(bjetchainDicts))    
#        theListOfChainDefs += [theDvChainDef] 
#
#
#    logDvDef.debug("----------------- Beginning of final individual chainDefs for dv-jet chains printout -----------------")
#    for chainDef in theListOfChainDefs:
#        logDvDef.debug(str(chainDef))        
#
#    logDvDef.debug("----------------- End of final individual chainDefs for dv-jet chains printout -----------------")
#        
#
#    if len(theListOfChainDefs)>1:
#        theFinalChainDef = mergeChainDefs(theListOfChainDefs,strategy="parallel",offset=-1,preserveL2EFOrder=True,removeDuplicateTEs=True)
#    else:
#        theFinalChainDef = theListOfChainDefs[0]
#
#    logDvDef.debug("----------------- Beginning of final merged chainDefs for b-jet chains printout -----------------")
#    logDvDef.debug(str(theFinalChainDef))
#
#
#    logDvDef.debug("----------------- End of final merged chainDefs for b-jet chains printout -----------------")
#        
#
#    return theFinalChainDef


###########################################################################
###########################################################################
def buildDvChains(jchaindef,dvjetdict,doAtL2AndEF=True,numberOfSubChainDicts=1):
    inputTEsEF = jchaindef.signatureList[-1]['listOfTriggerElements'][0]
    print 'MEOW inputTEsEF', inputTEsEF
    print 'MEOW all inpTEs', jchaindef.signatureList

    L2ChainName = "L2_" + dvjetdict['chainName']
    EFChainName = "EF_" + dvjetdict['chainName']
    HLTChainName = "HLT_" + dvjetdict['chainName']   
    topoAlgs = dvjetdict["topo"]

    dvjetparts = dvjetdict['chainParts']

    # Always use the split version
    #if ('split' in dvjetparts['bConfig']):
    theDvChainDef = myDvConfig_split(jchaindef, dvjetdict, inputTEsEF,numberOfSubChainDicts) 
    theDvChainDef.chain_name = 'HLT_'+dvjetdict['chainName']
    #else:
    #    theDvChainDef = myDvConfig1(jchaindef, dvjetdict, inputTEsEF,numberOfSubChainDicts) 
    #    theDvChainDef.chain_name = 'HLT_'+dvjetdict['chainName']

    return theDvChainDef


###################################################################################
###################################################################################

def myDvConfig_split(theChainDef, chainDict, inputTEsEF,numberOfSubChainDicts=1):
    from AthenaCommon.SystemOfUnits import GeV,mm

    L2ChainName = "L2_" + chainDict['chainName']
    EFChainName = "EF_" + chainDict['chainName']
    HLTChainName = "HLT_" + chainDict['chainName']   

    chainParts = chainDict['chainParts']
    btagthresh = chainParts['threshold']
    btagmult = chainParts['multiplicity']
    #btagcut = chainParts['bTag']
    #btagcut = btagcut[1:]
    #btracking = chainParts['bTracking']
    btracking="EFID"
    chainParts['bTracking']=btracking


    #-----------------------------------------------------------------------------------
    # Import of algs
    #-----------------------------------------------------------------------------------

    # jet hypo (can re-use their code)
    from TrigJetHypo.TrigJetHypoConfig import EFJetHypo
    theJetHypo = EFJetHypo(btagthresh, 0.0, 2.5)

    ################

    #jet splitting: 
    # Split the jet collection (stored as "TrigJetRec") into separate TE if:
    #   * the pt_jet is > thershold (default: 15 GeV)
    # Attachs per each splitted-jet:
    #   * TrigRoiDescriptor (called "SplitJet") per default, can be changed with JetOutputKey
    #   * xAOD::JetContainer  (called "SplitJet", can be changed with the same property)
    #   * TrigOperationalInfo (called "EFJetInfo", hardcoded) with just a "EFJetEt" info
    # 
    from TrigBjetHypo.TrigJetSplitterAllTEConfig import getJetSplitterAllTEInstance
    theJetSplit=getJetSplitterAllTEInstance()
    #theJetSplit.JetOuputKey="SplitJetDv" #es necesario? quizas si para asegurarme que utilizo l
    #                       # los  mios? Es posible compartir objetos entre diferentes secuencias??

    #--------------------

    # super ROI building: 
    #  Creates a super RoI from a Jet collection:
    #  Note that is the opposite to the Split FeX (above), and accepts and attach the same
    #  elements, but it takes all the RoI of each jet and it builds only one (superRoI), so
    #  there is one RoI for all the jet collection... IS this pretty equivalent to all the
    #  info in the Hadronic CAL? 
    # NEW-TRACKS from TrigBjetHypo.TrigSuperRoiBuilderAllTEConfig import getSuperRoiBuilderAllTEInstance
    # NEW-TRACKS theSuperRoi=getSuperRoiBuilderAllTEInstance()

    #--------------------

    # tracking:: I don't need tracking... or yes.. 
    from InDetTrigRecExample.EFInDetConfig import TrigEFIDSequence
    if 'EFID' in btracking:
    	# Trigger Track reconstruction:
	# Note the constructor:: (seqName,slicename, seqType)
        #         algorithmName_seqName_EFID will be the name of the algorithms called
	#         by this sequence (bjet,InsideOut): [PixelClustering,SCTClustering,...
	#         Probably here I can change the per default values
        theBjet_tracks = TrigEFIDSequence("Bjet","bjet","InsideOut").getSequence()   
        #theBjet_tracks = TrigEFIDSequence("Bjet","bjet","OutsideIn").getSequence()   --> OLD
        ## -- NEW TRYING TO INVERT SOME CUTS--

    else:
        from TrigFastTrackFinder.TrigFastTrackFinder_Config import TrigFastTrackFinder_Jet
        theTrigFastTrackFinder_Jet = [TrigFastTrackFinder_Jet()]
        from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
        theFastTrackFinderxAOD = TrigInDetSequence("Bjet","bjet","FastxAOD").getSequence()
        theBjet_tracks = TrigEFIDSequence("Bjet","bjet","DataPrep").getSequence()            
        theBjet_tracks += theTrigFastTrackFinder_Jet
        theBjet_tracks += theFastTrackFinderxAOD  # I think this is not needed, see below
        theBjet_tracks += TrigEFIDSequence("Bjet","bjet","InsideOutMerged").getSequence()            

    #--------------------
    # vertex tracking
    # The fast part of the new tracking (no precision stuff)
    # Note that there is available: TrigFastTrackFinder_X, where 
    #  X = Muon/eGamma/muonIso/Tau/TauCore/TauIso//Bphysics/FullScan/BeamSpot/Tile/
    #     FullScan_ZF_Only/Cosmic
    # Configuration: Trigger/TrigAlgorithms/TrigFastTrackFinder/python/TrigFastTrackFinder_Config.py
    # NEW TRACKS-from TrigFastTrackFinder.TrigFastTrackFinder_Config import TrigFastTrackFinder_Jet
    # NEW TRACKS-theTrigFastTrackFinder_Jet = [TrigFastTrackFinder_Jet()]
    # NEW TRACKS-from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
    # NEW TRACKS-theFastTrackFinderxAOD = TrigInDetSequence("BjetPrmVtx","bjet","FastxAOD").getSequence()
    # NEW TRACKS-theVertex_tracks = TrigEFIDSequence("BjetPrmVtx","bjet","DataPrep").getSequence()            
    # NEW TRACKS-theVertex_tracks += theTrigFastTrackFinder_Jet
    # NEW TRACKS-theVertex_tracks += theFastTrackFinderxAOD # does this convert to xAOD?

    # NEW TRACKS-#theVertex_tracks = TrigEFIDSequence("Bjet","bjet","InsideOut").getSequence()        

    #--------------------

    # primary vertexing
    # NEW TRACKS-from TrigT2HistoPrmVtx.TrigT2HistoPrmVtxAllTEConfig import EFHistoPrmVtxAllTE_Jet
    # NEW TRACKS-from TrigT2HistoPrmVtx.TrigT2HistoPrmVtxComboConfig import EFHistoPrmVtxCombo_Jet

    # NEW TRACKS-#--------------------

    # NEW TRACKS-# secondary vertexing
    # NEW TRACKS-from InDetTrigVxSecondary.InDetTrigVxSecondary_LoadTools import TrigVxSecondaryCombo_EF
    # NEW TRACKS-theVxSecondary = TrigVxSecondaryCombo_EF()
    # ---> TESTING FIXME
    #--- Redoing the Vx
    # NEW TRACKS-from AthenaCommon.AppMgr import ToolSvc
    # NEW TRACKS-from InDetVKalVxInJetTool.InDetVKalVxInJetToolConf import InDet__InDetVKalVxInJetTool
    # NEW TRACKS-theVxSecondary.SecVtxFinderList = []
    # NEW TRACKS-InDetVKalVxInJetTool_v2 = InDet__InDetVKalVxInJetTool(name = "InDetEFVKalVxInJetTool_v2",
    # NEW TRACKS-        getNegativeTail = False,
    # NEW TRACKS-        CutChi2         = 999999,
    # NEW TRACKS-        CutSctHits      = 0, 
    # NEW TRACKS-        CutPixelHits    = 0,
    # NEW TRACKS-        CutSiHits       = 0, #was 4
    # NEW TRACKS-        CutBLayHits     = 0,
    # NEW TRACKS-        CutSharedHits   = 0,#2,
    # NEW TRACKS-        FillHist        = True,
    # NEW TRACKS-        CutPt           = 0.5*GeV,
    # NEW TRACKS-        CutA0           = 100*mm,
    # NEW TRACKS-        CutZVrt         = 1500*mm,
    # NEW TRACKS-        MultiVertex     = True,
    # NEW TRACKS-        MultiWithPrimary= True,
    # NEW TRACKS-        #MultiWithOneTrkVrt=True,
    # NEW TRACKS-        VertexMergeCut  = 4.0*mm,
    # NEW TRACKS-        #OutputLevel      = DEBUG
    # NEW TRACKS-        )
    # NEW TRACKS-ToolSvc += InDetVKalVxInJetTool_v2
    # NEW TRACKS-theVxSecondary.SecVtxFinderList += [ InDetVKalVxInJetTool_v2 ]
    #theVxSecondary.SecVtxFinderList[0].FillHist=True
    #theVxSecondary.SecVtxFinderList[0].FillHist=True
    #theVxSecondary.SecVtxFinderList[0].CutPixelHits=0
    #theVxSecondary.SecVtxFinderList[0].CutPt=0.5*GeV
    #theVxSecondary.SecVtxFinderList[0].CutA0=100*mm
    #theVxSecondary.SecVtxFinderList[0].CutZVrt=1500*mm
    #print theVxSecondary.SecVtxFinderList[0]
    #from InDetTrigVxSecondary.InDetTrigVxSecondaryMonitoring import InDetTrigVxSecondaryValidationMonitoring

    #--------------------

    # bjet fex
    #if ('boffperf' in chainParts['bTag']):
    #    from TrigBjetHypo.TrigBtagFexConfig import getBtagFexInstance
    #    theBjetFex = getBtagFexInstance("EF","2012","EFID")
    #else:
    #    from TrigBjetHypo.TrigBjetFexConfig  import getBjetFexSplitInstance
    #    theBjetFex = getBjetFexSplitInstance("EF","2012","EFID")

    #if ('bperf' in chainParts['bTag'] or 'boffperf' in chainParts['bTag']):
    #    from TrigBjetHypo.TrigBjetHypoConfig import getBjetHypoSplitNoCutInstance
    #    theBtagReq = getBjetHypoSplitNoCutInstance("EF")
    #else:
    #    from TrigBjetHypo.TrigBjetHypoConfig import getBjetHypoSplitInstance
    #    theBtagReq = getBjetHypoSplitInstance("EF","2012", btagcut)

    # dv fex
    from TrigBjetHypo.TrigDvFexConfig import getDvFexInstance
    theDvFex = getDvFexInstance("EF","EFID")
    # dv hypo
    from TrigBjetHypo.TrigDvHypoConfig import getDvHypoInstance
    theDvHypo = getDvHypoInstance("EF")

    #-----------------------------------------------------------------------------------
    # TE naming
    #-----------------------------------------------------------------------------------
    
    if 'EFID' in btracking: 
        tracking = 'EFID'
    else: 
        tracking = "IDTrig"

    #jetHypoTE = "HLT_j"+btagthresh+"_eta"    #  Jets eta/pt cut
    jetHypoTE = "HLT_j35_eta"                #  Changed!!
    jetSplitTE = jetHypoTE+"_jsplit"         #  Jets splitted in several TE
    jetTrackTE = jetSplitTE+"_"+tracking     #  Jets reconstructed at HLT
    superTE = "HLT_super"                    #  One SuperRoI only
    superTrackingTE = superTE+"_"+tracking   #  Tracking in the SuperRoI
    prmVertexTE = superTrackingTE+"_prmVtx"  #  PrmVtx obtained with the SuperRoI tracking
    comboPrmVtxTE = prmVertexTE+"Combo"      #  PrmVtx as above, but using other algo
    secVtxTE = jetTrackTE+"__"+"superVtx"    #  Secondary Vertices
    lastTEout = "HLT_"+chainParts['chainPartName'] if numberOfSubChainDicts>1 else EFChainName

    #-----------------------------------------------------------------------------------
    # sequence assembling
    #-----------------------------------------------------------------------------------

    # Remember-> algorithm instance, inputTE, outputTE
    #theChainDef.addSequence(theJetHypo,     inputTEsEF, jetHypoTE) --> Change the first sequence: j35
    #theChainDef.addSequence(theJetSplit,    jetHypoTE,  jetSplitTE) --> Change the first sequence: j35
    theChainDef.addSequence(theJetSplit,inputTEsEF,jetSplitTE)
    theChainDef.addSequence(theBjet_tracks, jetSplitTE, jetTrackTE)

    #theChainDef.addSequence(theSuperRoi,      inputTEsEF, superTE)
    #theChainDef.addSequence(theVertex_tracks, superTE,    superTrackingTE)

    #theChainDef.addSequence([EFHistoPrmVtxAllTE_Jet()], superTrackingTE, prmVertexTE)
    #theChainDef.addSequence([EFHistoPrmVtxCombo_Jet()], [superTrackingTE,prmVertexTE], comboPrmVtxTE)    
    #theChainDef.addSequence(theVxSecondary, [jetTrackTE, comboPrmVtxTE], secVtxTE)
    #theChainDef.addSequence([theDvFex, theDvHypo], secVtxTE, lastTEout)#--> HASTA AQUI --->
    theChainDef.addSequence([theDvFex, theDvHypo], jetTrackTE, lastTEout)
    theChainDef.addSignature(theChainDef.signatureList[-1]['signature_counter']+1, [lastTEout]*int(btagmult))

    return theChainDef

###################################################################################

def myDvConfig1(theChainDef, chainDict, inputTEsEF,numberOfSubChainDicts=1):
    print "*******WARNING********"
    print "This Configuration should not be used!!"
    print "Use instead the '_split' version of the chain"
    print "*****END-WARNING******"
    L2ChainName = "L2_" + chainDict['chainName']
    EFChainName = "EF_" + chainDict['chainName']
    HLTChainName = "HLT_" + chainDict['chainName']   

    chainParts = chainDict['chainParts']
    btagthresh = chainParts['threshold']
    btagmult = chainParts['multiplicity']
    btagcut = chainParts['bTag']
    btagcut = btagcut[1:]
    btracking = chainParts['bTracking']

    #import fexes/hypos
    # -- Using al the TE separately
    from TrigBjetHypo.TrigEFBjetSequenceAllTEConfig import getEFBjetAllTEInstance
    ef_bjetSequence=getEFBjetAllTEInstance()

    from TrigBjetHypo.TrigBjetEtHypoConfig          import getBjetEtHypoInstance
    ef_ethypo_startseq = getBjetEtHypoInstance("EF","StartSequence","35GeV")

    from InDetTrigRecExample.EFInDetConfig import TrigEFIDSequence
    if 'EFID' in btracking:
        ef_bjet_tracks = TrigEFIDSequence("Bjet","bjet","InsideOut").getSequence()        
    else:
        from TrigFastTrackFinder.TrigFastTrackFinder_Config import TrigFastTrackFinder_Jet
        theTrigFastTrackFinder_Jet = [TrigFastTrackFinder_Jet()]
        from TrigInDetConf.TrigInDetSequence import TrigInDetSequence
        theFastTrackFinderxAOD = TrigInDetSequence("Bjet","bjet","FastxAOD").getSequence()
        ef_bjet_tracks = TrigEFIDSequence("Bjet","bjet","DataPrep").getSequence()            
        ef_bjet_tracks += theTrigFastTrackFinder_Jet
        ef_bjet_tracks += theFastTrackFinderxAOD
        ef_bjet_tracks += TrigEFIDSequence("Bjet","bjet","InsideOutMerged").getSequence()            


    from TrigT2HistoPrmVtx.TrigT2HistoPrmVtxAllTEConfig import EFHistoPrmVtxAllTE_Jet
    from TrigT2HistoPrmVtx.TrigT2HistoPrmVtxComboConfig import EFHistoPrmVtxCombo_Jet

    from InDetTrigVxSecondary.InDetTrigVxSecondary_LoadTools import TrigVxSecondary_EF
    ef_VxSecondary_EF = TrigVxSecondary_EF()

    from TrigBjetHypo.TrigBjetEtHypoConfig          import getBjetEtHypoInstance
    ef_EtHypo_Btagging = getBjetEtHypoInstance("EF","Btagging", btagthresh+"GeV")

    # DV FEX
    from TrigBjetHypo.TrigDvFexConfig import getDvFexInstance
    ef_dvfex = getDvFexInstance("EF","EFID")
    # dv hypo
    from TrigBjetHypo.TrigDvHypoConfig import getDvHypoInstance
    ef_hypo = getDvHypoInstance("EF")
    #if ('boffperf' in chainParts['bTag']):
    #    from TrigBjetHypo.TrigBtagFexConfig import getBtagFexInstance
    #    ef_bjet = getBtagFexInstance("EF","2012","EFID")
    #else:
    #    from TrigBjetHypo.TrigBjetFexConfig  import getBjetFexInstance
    #    ef_bjet = getBjetFexInstance("EF","2012","EFID")
    #
    #if ('bperf' in chainParts['bTag'] or 'boffperf' in chainParts['bTag']):
    #    from TrigBjetHypo.TrigBjetHypoConfig import getBjetHypoNoCutInstance
    #    ef_hypo = getBjetHypoNoCutInstance("EF")
    #else:
    #    from TrigBjetHypo.TrigBjetHypoConfig import getBjetHypoInstance
    #    ef_hypo = getBjetHypoInstance("EF","2012", btagcut)


    #------- 2012 EF Sequences based on j35 intput TE-------
    # TE naming
    ef2 ='HLT_BjetSeed'
    ef3 ='HLT_BjetSeed_EtCut%sGeV' % btagthresh
    if ('EFID' in chainParts['bTracking']):
        ef4 ='HLT_BjetSeed_EtCut%sGeV_EFID'  % btagthresh
        ef5 ='HLT_BjetSeed_EtCut%sGeV_AllTEPrmVtx_EFID'  % btagthresh
        ef6 ='HLT_BjetSeed_EtCut%sGeV_ComboPrmVtx_EFID'  % btagthresh
    else:
        ef4 ='HLT_BjetSeed_EtCut%sGeV_IDTrig'  % btagthresh
        ef5 ='HLT_BjetSeed_EtCut%sGeV_AllTEPrmVtx_IDTrig'  % btagthresh
        ef6 ='HLT_BjetSeed_EtCut%sGeV_ComboPrmVtx_IDTrig'  % btagthresh
    if (btagmult == '1'):
        ef7 = 'EF_b%s_%s_%s_VxSecondaryAndBTagHypo' % (btagthresh, btagcut, chainParts['chainPartName'].replace("_"+chainParts['bTracking'],""), )
    else:
        ef7 = 'EF_%sb%s_%s_%s_VxSecondaryAndBTagHypo' % (btagmult, btagthresh, btagcut, chainParts['chainPartName'].replace("_"+chainParts['bTracking'],""))

    theChainDef.addSequence([ef_bjetSequence], inputTEsEF, ef2)
    theChainDef.addSequence(ef_ethypo_startseq, ef2, ef3)
    theChainDef.addSequence(ef_bjet_tracks, ef3 ,ef4)
    theChainDef.addSequence([EFHistoPrmVtxAllTE_Jet()], ef4, ef5)
    theChainDef.addSequence([EFHistoPrmVtxCombo_Jet()], [ef4, ef5], ef6)
    #theChainDef.addSequence([ef_EtHypo_Btagging], ef6, ef7) 
    theChainDef.addSequence([ef_VxSecondary_EF,ef_EtHypo_Btagging], ef6, ef7) 
    lastTEout = "EF_"+chainParts['chainPartName'] if numberOfSubChainDicts>1 else EFChainName
    theChainDef.addSequence([ef_dvfex, ef_hypo], ef7, lastTEout)

    theChainDef.addSignature(theChainDef.signatureList[-1]['signature_counter']+1, [lastTEout]*int(btagmult))

    

    return theChainDef




###########################################################################
###########################################################################
# !!!!!!PLEASE DO NOT CHANGE THIS UNLESS NECESSARY!!!!
# This funciton checks the chainDict of the chain and merges chain parts 
# if only different in multiplicity and b-tag
# to then send off to the jet slice
##########################################################################
def _prepareJetChainDict(cdict):
    
    # -- collecting thresholds in a list to check if there are duplicates -- 
    thresholds = [part['threshold'] for part in cdict['chainParts']]
    from collections import Counter
    counts = Counter(thresholds)
    duplicates = [val for val, count in counts.items() if count > 1]

    if duplicates:
        logDvDef.info('Duplicated thresholds in the jet and b-jet part')
        
        # -- splitting chainparts dictioanries according to bjet and jet -- 
        bjetchainParts =[] 
        jetchainParts  = []
        for part in cdict['chainParts']:
            if not part['bTag']:
                jetchainParts.append(part)
            elif part['bTag']:
                bjetchainParts.append(part)
            else: 
                raise RuntimeError('I do not know what to do with this chain part %s' %(part['chainPartName']))
        
        # -- Loop over jet and bjet chain parts
        # -- pick out the ones with common thresholds
        # -- compare the jet building parts and if nothing differnet, adapt the multiplicity and chainPartName
        mergedbjetpart = []
        for jpart in jetchainParts:
            for bjindex, bjpart in enumerate(bjetchainParts):
                if (int(jpart['threshold']) == int(bjpart['threshold'])):                    
                    # -- check that all jet properties are the same as bjet 
                    # -- except for multiplicity, chainPartName and bTag
                    # -- this is bad hardcoding better clean up ...not
                    allowedDiff1 = set(['multiplicity', 'chainPartName', 'bTag','bTracking'])
                    allowedDiff2 = set(['multiplicity', 'chainPartName', 'bTag'])
                    allowedDiff3 = set(['multiplicity', 'chainPartName', 'bTag','bConfig'])
                    allowedDiff1_noMult = set(['chainPartName', 'bTag','bTracking'])
                    allowedDiff2_noMult = set(['chainPartName', 'bTag'])
                    allowedDiff3_noMult = set(['chainPartName', 'bTag','bConfig'])
                    
                    s1 = set(jpart)
                    s2 = set(bjpart)
                    diff = set(k for k in s1 & s2 if jpart[k] != bjpart[k])

                    if ((diff == allowedDiff1) or (diff == allowedDiff1_noMult) or (diff == allowedDiff2) or (diff == allowedDiff2_noMult) or (diff == allowedDiff3) or (diff == allowedDiff3_noMult)):
                        jpart['multiplicity'] = str(int(jpart['multiplicity'])+int(bjpart['multiplicity']))
                        jpart['chainPartName'] = jpart['multiplicity']+jpart['trigType']+jpart['threshold']
                        mergedbjetpart.append(bjindex)
                    else:
                        logDvDef.info("Jet and underlying jet chain from bjet are not the same despite the same thresholds. Ignore and keep separate dictionaries. The difference is %s " % str(diff))
                        continue



        # -- remove those bjet chainParts that have been merged with the jets
        for index in mergedbjetpart:
            bjetchainParts.pop(index)
        
        # -- modify the chainParts of the dictionary by creating a new one with 
        # -- the updated jetchainParts and remaining bjetchainParts
        cdict['chainParts'] = []
        for jpart in jetchainParts:
            cdict['chainParts'].append(jpart) 
        for bjpart in bjetchainParts:
            cdict['chainParts'].append(bjpart)


    else: 
        logDvDef.info('No duplicated thresholds in the jet/b-jet chain')
    
    logDvDef.debug("Prepared b-jet chainDict: \n %s" % (pp.pformat(cdict)))

    return cdict



###########################################################################
def get_j35_ChainDef():
    # HACK TO GET j35 chains!!!
    j35_chainDict = { 'EBstep': 1, 
                      'L1item': 'L1_J20', 
                      'chainCounter': 21, 
                      'chainName': 'j35', 
                      'chainParts': [ { 'L1item': '', 
                                        'addInfo': [], 
                                        'bTag': 'bmedium', 
                                        'calib': 'had', 
                                        'jetCalib': 'subjes', 
                                        'chainPartName': 'j35', 
                                        'dataType': 'tc', 
                                        'etaRange': '0eta320', 
                                        'extra': '', 
                                        'multiplicity': '1', 
                                        'recoAlg': 'a4', 
                                        'scan': 'FS', 
                                        'signature': 'Jet', 
                                        'threshold': '35', 
                                        'topo': [], 
                                        'trigType': 'j'}, ], 
                      'groups': ['RATE:MultiJet', 'BW:Jets'], 
                      'signature':'Jet', 
                      'signatures': '', 
                      'stream': ['Jet'], 
                      'topo': []} 

    return j35_chainDict

def get_lastTE_j35(ChainDef):
    inputTEsEF = ChainDef.signatureList[-1]['listOfTriggerElements']
    return inputTEsEF
###########################################################################


