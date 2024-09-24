import CombineHarvester.CombineTools.ch as ch
from argparse import ArgumentParser
import yaml
from python.helpers import *

# HI
description = '''This script makes datacards with CombineHarvester for performing tau ID SF measurments.'''
parser = ArgumentParser(prog="harvesterDatacards",description=description,epilog="Success!")
parser.add_argument('-c', '--config', dest='config', type=str, default='config/harvestDatacards.yml', action='store', help="set config file")
args = parser.parse_args()

with open(args.config, 'r') as file:
   setup = yaml.safe_load(file)

chans = setup['channels']
if chans == 'all': chans = ['tt'] # only using tt channel for now but can add mt and et later
else: chans = chans.split(',')

output_folder = setup['output_folder']
input_folder = setup['input_folder']

# define background processes
bkg_procs = ['ZTT','ZL','TTT','VVT','jetFakes']

# define signal processes, which are the same for every channel
sig_procs = {}
sig_procs['ggH'] = ['ggH_sm_htt','ggH_ps_htt','ggH_mm_htt']
sig_procs['qqH'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt','WH_sm_htt','WH_ps_htt','WH_mm_htt','ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']

# define categories which can depend on the channel
cats = {}
cats['tt'] = [
        (1, 'zttEmbed'),
        (2, 'jetFakes'),
        (3, 'higgs_Rho_Rho'),
        (4, 'higgs_0A1_Rho_and_0A1_0A1'),
        (5, 'higgs_A1_Rho'),
        (6, 'higgs_A1_A1_PolVec'),
        (7, 'higgs_Pi_Rho_Mixed'),
        (8, 'higgs_Pi_Pi'),
        (9, 'higgs_Pi_A1_Mixed'),
        (10,'tt_2016_higgs_Pi_0A1_Mixed'),
        (11,'tt_2016_higgs_A1_0A1'),
        ]


# Create an empty CombineHarvester instance
cb = CombineHarvester()

# Add processes and observations
for chn in channs:
    # Adding Data,Signal Processes and Background processes to the harvester instance
    cb.AddObservations(['*'], ['htt'], ['13p6TeV'], [chn], cats[chn])
    cb.AddProcesses(['*'], ['htt'], ['13p6TeV'], [chn], bkg_procs, cats[chn], False)
    cb.AddProcesses(['125'], ['htt'], ['13p6TeV'], [chn], sig_procs['ggH'], cats[chn], True)
    cb.AddProcesses(['125'], ['htt'], ['13p6TeV'], [chn], sig_procs['qqH'], cats[chn], True)

# TODO: systematics to be added here

# Populating Observation, Process and Systematic entries in the harvester instance
for chn in chans:
    filename = '%s/htt_%s.inputs-sm-13TeV.root' % (input_folder,chan)
    print (">>>   file %s" % (filename))
    cb.cp().channel([chn]).process(bkg_procs).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
    cb.cp().channel([chn]).process(sig_procs).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")

#SetStandardBinNames(cb) # needed??

#TODO: setup bbb's here


# Write datacards
print(green(">>> writing datacards..."))
datacardtxt  = "%s/cmb/$BIN.txt" % (output_folder)
datacardroot = "%s/cmb/common/$BIN_input.root" % (output_folder)
writer = CardWriter(datacardtxt,datacardroot)
writer.SetVerbosity(1)
writer.SetWildcardMasses([ ])
writer.WriteCards("cmb", cb)
