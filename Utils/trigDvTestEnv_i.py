#!/usr/bin/env python
# Should be launched with athena.py

#--- Put a check here!!

def sequenceinfo(dictchain):
    """ Just print the different steps of a trigger chain,
    specifying the algorithms involved and the input and 
    output names of the TE.
    """
    k = 0
    for s in dictchain.sequenceList:
        print "Step::%i  [\033[1;34m%s\033[1;m]%4s [\033[1;32m%s\033[1;m]   " % (k,s['input'],"-->",s['output'])
        #alg = "    ::--- Algorithm(s): \n"
        alg = ""
        if type(s['algorithm']) == list:
            for j in s['algorithm']:
                alg += "     |--- %s  (Package: %s)\n" % (j.getFullName(),j.getDlls())
            alg = alg[:-1]
        else:
            j = s['algorithm']
            alg = "     |--- %s  (Package: %s) " % (j.getFullName(),j.getDlls())
        print alg
        k+=1


include("RecExCond/RecExCommon_flags.py")

from TriggerMenu.test.TestSliceFlags import * 
#from TriggerMenu.bjet.BjetSliceFlags import * 
#from TriggerJobOpts.TriggerFlags import TriggerFlags  <-- already imported in the line above

# Test the current chains in the Bjet slice taking the MC_pp_v5 menu 
# --> Too slow
#TriggerFlags.triggerMenuSetup='MC_pp_v5'
#TriggerFlags.readHLTconfigFromXML=False
#TriggerFlags.readLVL1configFromXML=False
#TriggerFlags.BjetSlice.setAll()
#
#
#def dvonly():
#    TriggerFlags.Slices_all_setOff()
#    TriggerFlags.BjetSlice.setAll()
#
#
#from TriggerMenu import useNewTriggerMenu
#useNewTM = useNewTriggerMenu()
#from TriggerMenu.menu.GenerateMenu import GenerateMenu
#m=GenerateMenu.overwriteSignaturesWith(dvonly)

# --- Or including our own chain names...
PhysicsStream='Main'
TriggerFlags.TestSlice.signatures = [
#        ['j55_boffperf',                  8,    'L1_J20',[],  [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j55_bperf',                    10,    'L1_J20',[],  [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j45_bperf_3j45',               11,    'L1_3J15',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j45_bperf_3j45_L13J150ETA24',  12,    'L1_3J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j45_bperf_3j45_L13J20',        13,    'L1_3J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j175_bmedium',                 15,    'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j75_bmedium_3j75',             16,    'L1_4J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['2j55_bmedium_2j55',            17,    'L1_4J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['2j45_bmedium_3j45',            18,    'L1_5J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j175_bmedium_j60_bmedium',     19,    'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j300_bloose',                   9,    'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#
#        # split configuration: --> using bperf to associate the 'bTag' key needed in GenerateMenu
        ['dv_TestChain45',               901,    'L1_3J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet','BW:Jets'],-1],
        ['dv_TestChain55',               902,    'L1_3J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet','BW:Jets'],-1],
        ['dv_3TestChain45_L13J150ETA24', 903,    'L1_3J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet','BW:Jets'],-1],
        ['dv_3TestChain45_L13J20',       904,    'L1_3J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
        ['dv_1TestChain75',              905,    'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
        ['dv_3TestChain75',              906,    'L1_4J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
        ['dv_2TestChain55',              907,    'L1_4J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
        ['dv_3TestChain45',              908,    'L1_5J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
        ['dv_2TestChain60_1TestChain175',909,   'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
        ['dv_TestChain300',              910,   'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j45_bperf_split_3j45',               33,    'L1_3J15',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j45_bperf_split_3j45_L13J150ETA24',  34,    'L1_3J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j45_bperf_split_3j45_L13J20',        35,    'L1_3J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j175_bmedium_split',                 36,    'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j75_bmedium_split_3j75',             37,    'L1_4J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['2j55_bmedium_split_2j55',            38,    'L1_4J20',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['2j45_bmedium_split_3j45',            39,    'L1_5J15.0ETA24',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j175_bmedium_split_j60_bmedium_split',40,   'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],
#        ['j300_bloose_split',                   41,   'L1_J100',[], [PhysicsStream], ['RATE:MultiJet', 'BW:Jets'], -1],

        ]

from TriggerMenu.menu import DictFromChainName as chext
gdch = chext.DictFromChainName()
# Getting the chains dictds
#chains ={}
#for i in TriggerFlags.BjetSlice.signatures():
#    chains[i[0]] = gdch.getChainDict(i)
chains = dict(map(lambda x: (x[0],gdch.getChainDict(x)),TriggerFlags.TestSlice.signatures()))

from TriggerMenu import useNewTriggerMenu
useNewTM = useNewTriggerMenu()
#log.info("Using new TriggerMenu: %r" % useNewTM)

from AthenaCommon.Include import include
from TriggerMenu.test import generateDvChainDefs
chDefs = {}
for trname,signature in chains.iteritems():
    chDefs[trname]=generateDvChainDefs.generateChainDefs(signature)

# --- Final test
# --- Very final test
#def onlyTest():
#    TriggerFlags.Slices_all_setOff()
#from TriggerMenu.menu import GenerateMenu
#gm = GenerateMenu.GenerateMenu()
##chDefsFM= {}
##for trname,signature in chains.iteritems():
##    chDefsFM[trname]=gm.getChainDef(signature)
#GenerateMenu.GenerateMenu.overwriteSignaturesWith(onlyTest)
#gm.generate()

print "Ready to check things!!! Remember chain definitions availables"
print "at 'chDefs' object (chainDefs). Use 'sequenceinfo' function to print info about them"
print "ALSO: 'chains' (chainDict) object available"
