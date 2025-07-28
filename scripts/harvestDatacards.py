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
merge_mode = setup['merge_mode'] # use this option to specify if we want to flatten and/or symmetrise distributions
# 0: no merging, 1: merge symmetrised bins, 2: Run-2 style merging
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

if merge_mode == 2:
    flat_cats = ['tt_higgs_rhorho', 'tt_higgs_rhoa11pr', 'tt_higgs_rhoa1', 'tt_higgs_pirho', 'tt_higgs_pia11pr', 'tt_higgs_a11pra1']
    sym_cats = ['tt_higgs_a1a1', 'tt_higgs_pipi', 'tt_higgs_pia1']
elif merge_mode == 1: 
    flat_cats = []
    sym_cats = ['tt_higgs_rhorho', 'tt_higgs_rhoa11pr', 'tt_higgs_rhoa1', 'tt_higgs_pirho', 'tt_higgs_pia11pr', 'tt_higgs_a11pra1', 'tt_higgs_a1a1', 'tt_higgs_pipi', 'tt_higgs_pia1']
else:
    flat_cats = []
    sym_cats = []

## Populating Observation, Process and Systematic entries in the harvester instance

for chn in chans:
    if Run2: filename = '%s/htt_%s.inputs-sm-13TeV.root' % (input_folder,chn)
    else: filename = '%s/added_histo_Run2Bins-mergeXbins.root' % (input_folder)
    print (">>>   file %s" % (filename))
    cb.cp().channel([chn]).backgrounds().process(['ZL']).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC") # add data shapes
    if merge_mode == 0: 
        cb.cp().channel([chn]).process(bkg_procs).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
        for sig_proc in sig_procs.values(): cb.cp().channel([chn]).process(sig_proc).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC")
    else:
        for cat in cats[chn]:
            if cat[1] in flat_cats:
                cb.cp().channel([chn]).bin_id([cat[0]]).process(bkg_procs).process(['QCD'],False).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS_flat", "$BIN/$PROCESS_$SYSTEMATIC_flat")
                # jetFakes and signal are symmetrised rather than flattened
                cb.cp().channel([chn]).bin_id([cat[0]]).process(bkg_procs).process(['QCD']).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS_sym", "$BIN/$PROCESS_$SYSTEMATIC_sym")
                for sig_proc in sig_procs.values(): cb.cp().channel([chn]).process(sig_proc).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS$MASS_sym", "$BIN/$PROCESS$MASS_$SYSTEMATIC_sym")
            elif cat[1] in sym_cats:
                cb.cp().channel([chn]).bin_id([cat[0]]).process(bkg_procs).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS_sym", "$BIN/$PROCESS_$SYSTEMATIC_sym")
                for sig_proc in sig_procs.values(): cb.cp().channel([chn]).process(sig_proc).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS$MASS_sym", "$BIN/$PROCESS$MASS_$SYSTEMATIC_sym")
            else:
                cb.cp().channel([chn]).bin_id([cat[0]]).process(bkg_procs).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
                for sig_proc in sig_procs.values(): cb.cp().channel([chn]).bin_id([cat[0]]).process(sig_proc).era(['13p6TeV']).ExtractShapes(filename, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC")

ch.SetStandardBinNames(cb)

#TODO: setup bbb's here

def MatchingProcess(first, second):
    return (
        first.bin()      == second.bin() and
        first.process()  == second.process() and
        first.signal()   == second.signal() and
        first.analysis() == second.analysis() and
        first.era()      == second.era() and
        first.channel()  == second.channel() and
        first.bin_id()   == second.bin_id() and
        first.mass()     == second.mass()
    )

if merge_mode == 0:
    # If not flattening/symmetrising then add bbb uncerts using autoMC stats
    cb.SetAutoMCStats(cb, 0., 1, 1)
else:
    cb.cp().bin_id([1,2]).SetAutoMCStats(cb, 0., 1, 1) # use ausoMCstats for background categories since these don't merge bins
    # if not then use old method for BBBs to allow correlations to be taken into account
    # TODO: We could ask combine experts if it is possible to modify autoMCStats to allow for correlations which should be faster 
    # TODO: implement old method for BBBs

    #bbb_fakes = ch.BinByBinFactory()
    #bbb_fakes.SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_bbb_bin_$#") # this needs to have "_bbb_bin_" in the pattern for the mergeXbbb option to work
    #bbb_fakes.SetAddThreshold(0.)
    #bbb_fakes.SetMergeThreshold(0.5)
    #bbb_fakes.SetFixNorm(False)
    #bbb_fakes.MergeBinErrors(cb.cp().bin_id([1,2], False).backgrounds().process(['QCD']))
    #bbb_fakes.AddBinByBin(cb.cp().bin_id([1,2], False).backgrounds().process(['QCD']), cb)

    bbb_real = ch.BinByBinFactory()
    bbb_real.SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_bbb_bin_$#") # this needs to have "_bbb_bin_" in the pattern for the mergeXbbb option to work
    bbb_real.SetAddThreshold(0.)
    #bbb_real.SetMergeThreshold(0.5)
    bbb_real.SetMergeThreshold(1.0)
    bbb_real.SetFixNorm(False)
    #bbb_real.MergeBinErrors(cb.cp().bin_id([1,2], False).backgrounds().process(['ZTT']).process(['QCD'],False))
    #bbb_real.AddBinByBin(cb.cp().bin_id([1,2], False).backgrounds().process(['ZTT']).process(['QCD'],False), cb)
    bbb_real.MergeBinErrors(cb.cp().backgrounds().process(['QCD'],False))
    bbb_real.AddBinByBin(cb.cp().backgrounds().process(['QCD'],False), cb)


    # As we merged the x-axis bins then we need to rename the bbb uncertainties so that they are correlated properly
    # First we will deal with the catogiries with flat background when all phi_CP bins are merged into 1
    # we need to hardcode the bin number for the xbins
    # Each vector element i corresponds to the number of xbins for bin i+1
    # If these numbers aren't set correctly the method won't work so be careful!
    # Note that the merging is now only performed for the templates that have a flat distribution

    tt_nxbins = [1, 1, 10, 4, 4, 1, 10, 1, 1, 4, 4]  # tt channel binning, set to 1 if no flattening is applied
    mt_nxbins = [1, 1, 10, 1, 4, 4]                 # mt (and et) channel binning, set to 1 if no flattening is applied

    for chan in chans:
        bins = tt_nxbins if ch == "tt" else mt_nxbins
        for i, nxbins in enumerate(bins):
            if nxbins <= 1:
                continue
            print(f"Merging bbb uncertainties for {chan} channel for category {i+1}, nxbins set to {nxbins}")

            def process_callback(proc):
                nominal = proc.ClonedShape().Clone()

                def syst_callback(syst):
                    old_name = syst.name()
                    match_proc = MatchingProcess(proc, syst)

                    if match_proc and "_bbb_bin_" in old_name:
                        bin_num = int(old_name.split('_bbb_bin_')[-1])

                        if (bin_num-1) % nxbins == 0:
                            print('bin_num again:', bin_num)
                            nonum_name = old_name.replace(f"_bbb_bin_{bin_num}", '_bbb_bin_')
                            shape_u_new = syst.ClonedShapeU().Clone()
                            shape_d_new = syst.ClonedShapeD().Clone()
                            shape_u_new.Scale(syst.value_u()) ###
                            shape_d_new.Scale(syst.value_d()) ###
                            shape_u_new.Add(nominal, -1)
                            shape_d_new.Add(nominal, -1)
                            names = []
                            for j in range(bin_num + 1, bin_num + nxbins):
                                names.append(f"{nonum_name}{j}")

                            def merge_syst_shapes(s):
                                shape_u_temp = s.ClonedShapeU().Clone()
                                shape_d_temp = s.ClonedShapeD().Clone()
                                shape_u_temp.Scale(s.value_u()) ###
                                shape_d_temp.Scale(s.value_d()) ###
                                shape_u_temp.Add(nominal, -1)
                                shape_d_temp.Add(nominal, -1)
                                shape_u_new.Add(shape_u_temp)
                                shape_d_new.Add(shape_d_temp)

                            cb.cp().syst_name(names).ForEachSyst(merge_syst_shapes)

                            shape_u_new.Add(nominal)
                            shape_d_new.Add(nominal)

                            syst.set_shapes(shape_u_new, shape_d_new, nominal)
                            
                            for n in names:
                                cb.FilterSysts(lambda s: s.name() == n)
                            
                cb.cp().ForEachSyst(syst_callback)

            #cb.cp().channel([chan]).bin_id([i+1]).backgrounds().process(["QCD"], False).ForEachProc(process_callback)
            #TODO: need to see why the process where everything gets merged into appears to be random


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

