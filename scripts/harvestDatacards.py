import CombineHarvester.CombineTools.ch as ch
from argparse import ArgumentParser
import yaml
from CombineHarvester.Combine_HtautauCP.helpers import *
from CombineHarvester.Combine_HtautauCP.systematics import AddSMRun3Systematics

# HI
description = '''This script makes datacards with CombineHarvester for performing tau ID SF measurments.'''
parser = ArgumentParser(prog="harvesterDatacards",description=description,epilog="Success!")
parser.add_argument('-c', '--config', dest='config', type=str, default='configs/harvestDatacards.yml', action='store', help="set config file")
args = parser.parse_args()

with open(args.config, 'r') as file:
   setup = yaml.safe_load(file)

chans = setup['channels']
if chans == 'all': chans = ['tt'] # only using tt channel for now but can add mt and et later
else: chans = chans.split(',')


output_folder = setup['output_folder']
input_folder = setup['input_folder']
mergeSymBins = setup['mergeSymBins'] # use this option to specify if we want to flatten and/or symmetrise distributions
# TODO: implement this in this script based on the extracted shapes rather than using the additional pre-processing step as we did for Run-2

Run2 = False
if Run2:
    # define background processes
    bkg_procs = ['ZTT','ZL','TTT','VVT','jetFakes']
    
    # define signal processes, which are the same for every channel
    sig_procs = {}
    sig_procs['ggH'] = ['ggH_sm_htt','ggH_ps_htt','ggH_mm_htt']
    sig_procs['qqH'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt','WH_sm_htt','WH_ps_htt','WH_mm_htt','ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']
    
    # define categories which can depend on the channel
    cats = {}
    cats['tt'] = [
            (1, 'tt_2018_zttEmbed'),
            (2, 'tt_2018_jetFakes'),
            (3, 'tt_2018_higgs_Rho_Rho'),
            (4, 'tt_2018_higgs_0A1_Rho_and_0A1_0A1'),
            (5, 'tt_2018_higgs_A1_Rho'),
            (6, 'tt_2018_higgs_A1_A1_PolVec'),
            (7, 'tt_2018_higgs_Pi_Rho_Mixed'),
            (8, 'tt_2018_higgs_Pi_Pi'),
            (9, 'tt_2018_higgs_Pi_A1_Mixed'),
            (10,'tt_2018_higgs_Pi_0A1_Mixed'),
            (11,'tt_2018_higgs_A1_0A1'),
            ]

else: 
    # define background processes
    bkg_procs = ['ZTT','ZL','TTT','VVT','QCD','ZJ','TTJ','VVJ','W']

    # define signal processes, which are the same for every channel
    sig_procs = {}
    sig_procs['ggH'] = ['ggH_sm_prod_sm_htt','ggH_ps_prod_sm_htt','ggH_mm_prod_sm_htt']
    sig_procs['qqH'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt']

    # define categories which can depend on the channel
    cats = {}
    cats['tt'] = [
            (1, 'tt_mva_tau'),
            (2, 'tt_mva_fake'),
            (3, 'tt_higgs_rhorho'),
            (4, 'tt_higgs_rhoa11pr'),
            (5, 'tt_higgs_rhoa1'),
            (6, 'tt_higgs_a1a1'),
            (7, 'tt_higgs_pirho'),
            (8, 'tt_higgs_pipi'),
            (9, 'tt_higgs_pia1'),
            (10,'tt_higgs_pia11pr'),
            (11,'tt_higgs_a11pra1'),
            ]



# Create an empty CombineHarvester instance
cb = ch.CombineHarvester()

# Add processes and observations
for chn in chans:
    # Adding Data,Signal Processes and Background processes to the harvester instance
    cb.AddObservations(['*'], ['htt'], ['13p6TeV'], [chn], cats[chn])
    cb.AddProcesses(['*'], ['htt'], ['13p6TeV'], [chn], bkg_procs, cats[chn], False)
    cb.AddProcesses(['125'], ['htt'], ['13p6TeV'], [chn], sig_procs['ggH'], cats[chn], True)
    cb.AddProcesses(['125'], ['htt'], ['13p6TeV'], [chn], sig_procs['qqH'], cats[chn], True)

# TODO: systematics to be added here
cb = AddSMRun3Systematics(cb)

# Populating Observation, Process and Systematic entries in the harvester instance
for chn in chans:
    if Run2: filename = '%s/htt_%s.inputs-sm-13TeV.root' % (input_folder,chn)
    else: filename = '%s/added_histo_Run2Bins.root' % (input_folder)
    print (">>>   file %s" % (filename))
    cb.cp().channel([chn]).process(bkg_procs).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
    for sig_proc in sig_procs.values(): 
        cb.cp().channel([chn]).process(sig_proc).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC")

ch.SetStandardBinNames(cb)

#TODO: setup bbb's here

if not mergeSymBins:
    # If not flattening/symmetrising then add bbb uncerts using autoMC stats
    cb.SetAutoMCStats(cb, 0., 1, 1)
#else:
    # if not then use old method for BBBs to allow correlations to be taken into account
    # TODO: We could ask combine experts if it is possible to modify autoMCStats to allow for correlations which should be faster 
    # TODO: implement old method for BBBs


# Implement fixes for negative bins and yields

# Zero negetive bins
print(green(">>> Zeroing negative bins"))
cb.ForEachProc(NegativeBins)

print(green(">>> Zeroing negative yields"))
cb.ForEachProc(NegativeYields)

# Write datacards
print(green(">>> Writing datacards..."))
datacardtxt  = "%s/$TAG/$BIN.txt" % (output_folder)
datacardroot = "%s/$TAG/common/$BIN_input.root" % (output_folder)
writer = ch.CardWriter(datacardtxt,datacardroot)
writer.SetVerbosity(1)
writer.SetWildcardMasses([ ])
writer.WriteCards("cmb", cb)
# Cards per category
writer.WriteCards("rhorho",   cb.cp().channel({"tt"}).bin_id({1,2,3}))
writer.WriteCards("rhoa11pr", cb.cp().channel({"tt"}).bin_id({1,2,4}))
writer.WriteCards("rhoa1",    cb.cp().channel({"tt"}).bin_id({1,2,5}))
writer.WriteCards("a1a1",     cb.cp().channel({"tt"}).bin_id({1,2,6}))
writer.WriteCards("pirho",    cb.cp().channel({"tt"}).bin_id({1,2,7}))
writer.WriteCards("pipi",     cb.cp().channel({"tt"}).bin_id({1,2,8}))
writer.WriteCards("pia1",     cb.cp().channel({"tt"}).bin_id({1,2,9}))
writer.WriteCards("pia11pr",  cb.cp().channel({"tt"}).bin_id({1,2,10}))
writer.WriteCards("a11pra1",  cb.cp().channel({"tt"}).bin_id({1,2,11}))

