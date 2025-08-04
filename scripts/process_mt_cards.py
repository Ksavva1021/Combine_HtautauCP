import ROOT
import os
from hadd_cp_datacards import hadd_root_files

dirname = 'mt_datacards_desy'



eras = ['preEE', 'postEE','preBPix']

# first hadd together all the datacards inside the era directories
for era in eras:
    era_dir = os.path.join(dirname, era)
    if not os.path.exists(era_dir):
        print(f"Directory {era_dir} does not exist.")
        continue
    
    os.system(f'hadd -f {era_dir}/mt_datacard_{era}.root {era_dir}/shapes__cat_*.root')

# now combine eras into a single file
input_files = [os.path.join(dirname, era, f'mt_datacard_{era}.root') for era in eras] 
output_file = os.path.join(dirname, 'mt_datacard_erascombined.root')

dir_combinations = {
    'mt_mva_fake' : ['fake'],
    'mt_mva_tau' : ['tautau'],
}

for i in [0, 1, 2]:
    dir_combinations[f'mt_mupi_cat{i}'] = [f'_cat{i}_tau2pi']
    dir_combinations[f'mt_murho_cat{i}'] = [f'_cat{i}_tau2rho']
    dir_combinations[f'mt_mua11pr_cat{i}'] = [f'_cat{i}_tau2a1']
    dir_combinations[f'mt_mua1_cat{i}'] = [f'_cat{i}_tau2a1_3pr']

hadd_root_files(input_files, output_file, dir_combinations)

def roll_histograms(hists):
    # get total number of bins of all histograms
    total_bins = sum(hist.GetNbinsX() for hist in hists)
    # create a new histogram with the total number of bins
    rolled_hist = ROOT.TH1F(hists[0].GetName(), '', total_bins, 0, total_bins)

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

            if dir_name[:-1].endswith('cat'):
                dir_name_nonum = dir_name[:-1]
                cat_name = f'cat{dir_name[-1]}'

                if dir_name_nonum not in combined_dirs_histograms:
                    combined_dirs_histograms[dir_name_nonum] = {}

                for hist_key in obj.GetListOfKeys():
                    hist_name = hist_key.GetName()
                    hist = hist_key.ReadObj()

                    # Check if it's a histogram
                    if isinstance(hist, ROOT.TH1):
                        # If it's the first time we encounter this histogram, clone it
                        if hist_name not in combined_dirs_histograms[dir_name_nonum]:
                            combined_dirs_histograms[dir_name_nonum][hist_name] = []
                        hist_clone = hist.Clone()
                        hist_clone.SetDirectory(0)
                        hist_clone.SetName(f'{hist_name}_{cat_name}')
                        combined_dirs_histograms[dir_name_nonum][hist_name].append(hist_clone)

            #else: # write the directory to the output file as is
            #    output_file.cd()
            #    obj.Write()



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
                rolled_hist = roll_histograms(hists)
                rolled_hist.Write(hist_name)
            else:
                hists[0].Write(hist_name)
        #print(f'Writing directory {dir_name} with {len(histograms)} histograms')
        #print(histograms)

    output_file.Close()
    input_file.Close()

merge_seperate_cats(output_file, output_file.replace('.root', '_merged.root'))