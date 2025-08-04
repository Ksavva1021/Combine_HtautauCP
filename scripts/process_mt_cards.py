import ROOT
import os
from hadd_cp_datacards import hadd_root_files

datacard_name = 'mt_2022_2023preBPix.root'


def roll_histograms(hists):
    # get total number of bins of all histograms
    total_bins = sum(hist.GetNbinsX() for hist in hists)
    # create a new histogram with the total number of bins
    rolled_hist = ROOT.TH1F(hists[0].GetName().split('_cat')[0], '', total_bins, 0, total_bins)

    # Fill the rolled histogram
    bin_offset = 0
    for hist in hists:
        for bin in range(1, hist.GetNbinsX() + 1):
            rolled_hist.SetBinContent(bin + bin_offset, hist.GetBinContent(bin))
        bin_offset += hist.GetNbinsX()

    return rolled_hist

# define a function that combines histograms from ..._cati directories into a single 1D histogram
def merge_seperate_cats(input_file, output_file):
    input_file = ROOT.TFile(input_file, 'READ')
    output_file = ROOT.TFile(output_file, 'RECREATE')

    # Create a dictionary to hold histograms for each directory
    combined_dirs_histograms = {}

    # Loop through all keys in the input file
    for key in input_file.GetListOfKeys():
        obj = key.ReadObj()
        print(key.GetName())
        if isinstance(obj, ROOT.TDirectoryFile):
            dir_name = key.GetName()
            #hist_name = obj.GetName()

            if 'cat' in dir_name:
                dir_name_nocat = '_'.join(dir_name.split('_')[2:])
                cat_name = '_'+dir_name.split('_')[1]
            else: 
                dir_name_nocat = dir_name
                cat_name = ''

            if dir_name_nocat not in combined_dirs_histograms:
                combined_dirs_histograms[dir_name_nocat] = {}
            for hist_key in obj.GetListOfKeys():
                hist_name = hist_key.GetName()
                hist = hist_key.ReadObj()
                # Check if it's a histogram
                if isinstance(hist, ROOT.TH1):
                    # If it's the first time we encounter this histogram, clone it
                    if hist_name not in combined_dirs_histograms[dir_name_nocat]:
                        combined_dirs_histograms[dir_name_nocat][hist_name] = []
                    hist_clone = hist.Clone()
                    hist_clone.SetDirectory(0)
                    hist_clone.SetName(f'{hist_name}{cat_name}')
                    combined_dirs_histograms[dir_name_nocat][hist_name].append(hist_clone)

    # Write the combined histograms to the output file
    for dir_name, histograms in combined_dirs_histograms.items():
        output_file.mkdir(dir_name)
        output_file.cd(dir_name)
        print(f'Writing directory {dir_name}')
        print(histograms)
        for hist_name, hists in histograms.items():
            print(f'Writing histogram {hist_name} with {len(hists)} histograms')
            print(f'Histograms: {hists}')
            print(f'Histogram names: {[hist.GetName() for hist in hists]}')
            if len(hists) > 1:
                # sort histograms based on _catX at the end of their names
                hists.sort(key=lambda x: int(x.GetName().split('_')[-1].replace('cat', '')))
                print('!!!!')
                print(f'Histogram names: {[hist.GetName() for hist in hists]}')
                rolled_hist = roll_histograms(hists)
                rolled_hist.Write(hist_name)
            else:
                hists[0].Write(hist_name)
        #print(f'Writing directory {dir_name} with {len(histograms)} histograms')
        #print(histograms)

    output_file.Close()
    input_file.Close()

merge_seperate_cats(datacard_name, datacard_name.replace('.root', '_merged.root'))