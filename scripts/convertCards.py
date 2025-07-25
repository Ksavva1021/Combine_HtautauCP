import ROOT
import math
import argparse
ROOT.TH1.AddDirectory(False)
from ctypes import c_double
 
# note that the merging of bins requires an even numnber of phi_CP bins, and this number must be set to the specific value used in the dictionary below otherwise the method will give incorrect results!

cp_bins = {
        "tt_mva_fake": 1,
        "tt_mva_tau": 1,
        "tt_mva_higgs": 1,
        "tt_higgs_pia1" : 4,
        "tt_higgs_a11pra1" : 4,
        "tt_higgs_rhoa1" : 4,
        "tt_higgs_pia11pr" : 4,
        "tt_higgs_pirho" : 10,
        "tt_higgs_a1a1" : 4,
        "tt_higgs_pipi" : 4,
        "tt_higgs_rhoa11pr" : 4,
        "tt_higgs_rhorho" : 10,
}

def MergeXBins(hist, nxbins):
  histnew = hist.Clone()
  nbins = hist.GetNbinsX()
  chi2_total = 0.
  ndf_perbin = nxbins-1
  chi2_perbin=[]
  for i in range(1,nbins+1,nxbins):
    tot_err = c_double(0)
    tot = hist.IntegralAndError(i,i+nxbins-1,tot_err)
    tot_err = tot_err.value
    ave = tot / nxbins
    ave_err = tot_err / nxbins
    chi2 = 0.
    for j in range(i,i+nxbins):
      if hist.GetBinError(j) > 0: chi2 += (hist.GetBinContent(j)-ave)**2 / (hist.GetBinError(j)**2)
      histnew.SetBinContent(j,ave)
      histnew.SetBinError(j,ave_err)
    chi2_perbin.append(chi2)
  p_val_perbin = [ROOT.TMath.Prob(x, ndf_perbin) for x in chi2_perbin]
  chi2_total = sum(chi2_perbin)
  ndf_total = ndf_perbin * len(chi2_perbin)
  p_val_total = ROOT.TMath.Prob(chi2_total, ndf_total)

  return histnew, p_val_total, p_val_perbin

def Symmetrise(hist,nxbins):
  histnew=hist.Clone()
  nbins = hist.GetNbinsX()
  if nbins % 2:
    print('N X bins in 2D histogram is not even so cannot symmetrise!')
    return
  nybins = int(nbins/nxbins)
  chi2_perbin = [0 for i in range(nybins)]
  ndf_perbin = int(nxbins / 2)
  for i in range(1,int(nxbins/2)+1):
    lo_bin = i
    hi_bin = nxbins-i+1
    chi2 = 0.
    for j in range(1,nybins+1):
      lo_bin_ = lo_bin+(j-1)*nxbins
      hi_bin_ = hi_bin+(j-1)*nxbins
      c1 = hist.GetBinContent(lo_bin_)
      c2 = hist.GetBinContent(hi_bin_)
      e1 = hist.GetBinError(lo_bin_)
      e2 = hist.GetBinError(hi_bin_)
      cnew = (c1+c2)/2
      enew = math.sqrt(e1**2 + e2**2)/2
      histnew.SetBinContent(lo_bin_,cnew)
      histnew.SetBinContent(hi_bin_,cnew)
      histnew.SetBinError(lo_bin_,enew)
      histnew.SetBinError(hi_bin_,enew)

      chi2 = 0.
      if (e1**2 + e2**2) > 0: chi2 = (c1-c2)**2 / (e1**2 + e2**2)
      chi2_perbin[j-1] += chi2
    p_val_perbin = [ ROOT.TMath.Prob(x, ndf_perbin) for x in chi2_perbin]
    chi2_total = sum(chi2_perbin)
    ndf_total = ndf_perbin * nybins
    p_val_total = ROOT.TMath.Prob(chi2_total, ndf_total)

  return histnew, p_val_total, p_val_perbin

def ASymmetrise(hist,hsm,hps,nxbins):
  histnew=hist.Clone()
  hsub=hsm.Clone()
  hsub.Add(hps)
  hsub.Scale(0.5)
  for i in range(1,hsub.GetNbinsX()+1): 
    histnew.SetBinContent(i,histnew.GetBinContent(i)-hsub.GetBinContent(i))

  nbins = hist.GetNbinsX()
  if nbins % 2:
    print('N X bins in 2D histogram is not even so cannot symmetrise!')
    return
  nybins = int(nbins/nxbins)
  for i in range(1,int(nxbins/2)+1):
    lo_bin = i
    hi_bin = nxbins-i+1
    for j in range(1,nybins+1):
      lo_bin_ = lo_bin+(j-1)*nxbins
      hi_bin_ = hi_bin+(j-1)*nxbins

      mmi = hist.GetBinContent(lo_bin_)       
      mmj = hist.GetBinContent(hi_bin_)       
      smi = hsm.GetBinContent(lo_bin_) 
      smj = hsm.GetBinContent(hi_bin_)
      psi = hps.GetBinContent(lo_bin_)
      psj = hps.GetBinContent(hi_bin_)

      e_mmi = hist.GetBinError(lo_bin_)
      e_mmj = hist.GetBinError(hi_bin_)
      e_smi = hsm.GetBinError(lo_bin_)
      e_smj = hsm.GetBinError(hi_bin_)
      e_psi = hps.GetBinError(lo_bin_)
      e_psj = hps.GetBinError(hi_bin_) 

      c1_new = ( smj+psj-mmj + mmi)/2
      c2_new = ( smi+psi-mmi + mmj)/2

      e1_new = math.sqrt((e_smj+e_psj-e_mmj)**2 + e_mmi**2)/2 
      e2_new = math.sqrt((e_smi+e_psi-e_mmi)**2 + e_mmj**2)/2 

      histnew.SetBinContent(lo_bin_,c1_new)
      histnew.SetBinContent(hi_bin_,c2_new)
      histnew.SetBinError(lo_bin_,e1_new)
      histnew.SetBinError(hi_bin_,e2_new)

  return histnew

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname):
    directory = infile.Get(dirname)
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1D) or isinstance(histo,ROOT.TH1F):
            print('Processing:', dirname, histo.GetName()) 
            if dirname in cp_bins: nxbins = cp_bins[dirname]
            else: nxbins=1

            # we write all the old histograms to the output file
            outfile.cd()
            if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
            ROOT.gDirectory.cd(dirname)
            histo.Write()

            if nxbins== 1: continue

            # if data we skip
            if 'data_obs' in key.GetName(): continue

            # if mm signal then we only anti-symmetrise
            elif 'H_mm' in key.GetName() and 'htt125' in key.GetName() and nxbins>1:
                hsm = directory.Get(key.GetName().replace('H_mm_','H_sm_'))
                hps = directory.Get(key.GetName().replace('H_mm_','H_ps_'))
                histo_asym = ASymmetrise(histo,hsm,hps,nxbins)
                histo_asym.SetName(histo.GetName()+'_asym')
                histo_asym.Write()
                continue 

            # if not mm signal then we always symmetrise
            else:
                histo_sym, p_val_total, p_val_perbin = Symmetrise(histo,nxbins)
                histo_sym.SetName(histo.GetName()+'_sym')
                histo_sym.Write()

                # if not signal then we also flatten
                if 'htt125' not in key.GetName() or 'Higgs_flat' in key.GetName() or 'H_flat' in key.GetName():
                    histo_flat, p_val_total, p_val_perbin = MergeXBins(histo,nxbins)
                    histo_flat.SetName(histo.GetName()+'_flat')
                    histo_flat.Write()

        ROOT.gDirectory.cd('/')

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which we want to merge X bins')
parser.add_argument('--test', '-t', action='store_true', help= 'Run statistical tests on the histograms to check compatability of smoothed and unsmoothed histograms')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','-mergeXbins.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename,"RECREATE") 

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        getHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname)


##########

def rescale_chunks(h_toy, h_data, nxbins):
    nbins = h_toy.GetNbinsX()
    if nbins != h_data.GetNbinsX():
        raise ValueError("Histograms must have the same number of bins")
    
    for start_bin in range(1, nbins + 1, nxbins):
        end_bin = min(start_bin + nxbins - 1, nbins)

        # Compute total yield in this chunk
        sum_data = sum(h_data.GetBinContent(i) for i in range(start_bin, end_bin + 1))
        sum_toy  = sum(h_toy.GetBinContent(i) for i in range(start_bin, end_bin + 1))

        if sum_toy == 0:
            continue  # avoid division by zero

        scale_factor = sum_data / sum_toy

        # Rescale h_toy bins in this chunk
        for i in range(start_bin, end_bin + 1):
            content = h_toy.GetBinContent(i)
            error   = h_toy.GetBinError(i)

            h_toy.SetBinContent(i, content * scale_factor)
            h_toy.SetBinError(i, error * scale_factor)

def chi2_test(h_data, h_model, nxbins=1):
    h_data = h_data.Clone()
    h_model = h_model.Clone()

    if h_data.GetNbinsX() != h_model.GetNbinsX():
        raise ValueError("Histograms must have the same number of bins")
    
    chi2_obs = 0.0
    ndof = 0

    for i in range(1, h_data.GetNbinsX() + 1):  # bin 0 is underflow, skip it
        D = h_data.GetBinContent(i)
        M = h_model.GetBinContent(i)
        sigma = h_data.GetBinError(i)

        if sigma <= 0:
            continue  # skip bins with zero variance

        chi2_obs += (D - M)**2 / (sigma**2)
        ndof += 1

    p_value_old = ROOT.TMath.Prob(chi2_obs, ndof)

    # use toys to get p-value
    n_toys = 1000
    toys_chi2 = []

    N_err = c_double(0)
    N = h_data.IntegralAndError(-1,-1, N_err)
    N_eff = (N/N_err.value)**2

    for t in range(n_toys):

        toy_mode = 1

        h_toy = h_data.Clone()
        h_toy.Reset()
        if toy_mode == 1:
            h_toy.FillRandom(h_model, int(N_eff))
            h_toy.Scale(h_data.Integral()/h_toy.Integral())

        elif toy_mode == 2:
            for j in range(1, h_toy.GetNbinsX() + 1):
                central_value = h_model.GetBinContent(j)
                random_shift = ROOT.gRandom.Gaus(0, h_data.GetBinError(j)) # random shift based on the error
                #random_shift = ROOT.gRandom.Gaus(0, h_model.GetBinError(j)) # random shift based on the error
                h_toy.SetBinContent(j, central_value + random_shift)
                h_toy.SetBinError(j, h_data.GetBinError(j))

        chi2=0.0
        for j in range(1, h_toy.GetNbinsX() + 1):
            if h_toy.GetBinError(j)>0: chi2 += ((h_toy.GetBinContent(j) - h_model.GetBinContent(j))**2)/(h_toy.GetBinError(j)**2)
        toys_chi2.append(chi2)
    # p_value = number of times toy chi2 is greater than observed chi2 / total number of toys
    p_value = sum(1 for toy in toys_chi2 if toy >= chi2_obs) / n_toys
    print('!!!!', p_value, p_value_old)

    return p_value

if args.test:

    import os
    if not os.path.exists('test_results'):
        os.makedirs('test_results')

    def ZeroErrors(hist):
        for i in range(1,hist.GetNbinsX()+1):
            hist.SetBinError(i,0.)

    def LastNBinsHistogram(histo, nxbins):
        total_bins = histo.GetNbinsX()
        if nxbins > total_bins:
            raise ValueError(f"Requested {nxbins} bins, but histogram only has {total_bins}.")

        # Create a new histogram with x-axis from 0 to nxbins (arbitrary range)
        new_histo = ROOT.TH1F(histo.GetName() + "_lastBins", histo.GetTitle(), nxbins, 0, nxbins)

        # Fill the new histogram with the last nxbins from original
        for i in range(nxbins):
            orig_bin = total_bins - nxbins + i + 1
            new_histo.SetBinContent(i + 1, histo.GetBinContent(orig_bin))
            new_histo.SetBinError(i + 1, histo.GetBinError(orig_bin))

        return new_histo

    def PerformTests(histo, histo_smooth, nxbins=1):
        #ZeroErrors(histo_smooth)  # zero errors so they aren't counted twice in the test
        chi2 = histo.Chi2Test(histo_smooth, 'WWCHI2')
        if chi2 <= 0:
            chi2_pval = None
        else:
            ndf = histo.GetNbinsX()/2
            chi2_pval = histo.Chi2Test(histo_smooth, 'WW')

        chi2_last_pval = chi2_pval

        #chi2_pval = chi2_test(histo, histo_smooth, nxbins)
#
        #if nxbins>1: 
        #    # make a historam using only the last bins
        #    histo_last = LastNBinsHistogram(histo, nxbins)
        #    histo_smooth_last = LastNBinsHistogram(histo_smooth, nxbins)
        #    chi2_last_pval = chi2_test(histo_last, histo_smooth_last, nxbins)
        #    #ZeroErrors(histo_smooth_last)  # zero errors so they aren't counted twice in the test
        #    #chi2_last = histo_last.Chi2Test(histo_smooth_last, 'WWCHI2')
        #    #if chi2_last<=0: chi2_last_pval = None
        #    #else: chi2_last_pval = histo_last.Chi2Test(histo_smooth_last, 'WW')
        #else: chi2_last_pval = chi2_pval
#
        ## do chi2 tests for last BDT score bin (most sinal sensitive)
        ##return (chi2_pval, chi2_last_pval)
        return (chi2_pval, chi2_last_pval)

    test_results_sym = {}
    test_results_flat = {}
    test_results_asym = {}
    print('Running statistical tests on the histograms to check compatability of smoothed and unsmoothed histograms')
    for key in output_file.GetListOfKeys():
        if isinstance(output_file.Get(key.GetName()),ROOT.TDirectory):
            dirname=key.GetName()
            directory = output_file.Get(dirname)

            if dirname in cp_bins: nxbins = cp_bins[dirname]
            else: nxbins=1
 
            for subkey in directory.GetListOfKeys():
                if subkey.GetName().endswith('_sym') or subkey.GetName().endswith('_flat') or subkey.GetName().endswith('_asym'):
                    continue
                # skip systematics
                if subkey.GetName().endswith('Up') or subkey.GetName().endswith('Down'):
                    continue

                if 'unfiltered' in subkey.GetName(): continue

                if not ('ggH' in subkey.GetName() or 'qqH' in subkey.GetName() or 'ZTT' in subkey.GetName() or 'QCD' in subkey.GetName()): continue 
               
                if 'prod_ps' in subkey.GetName() or 'prod_mm' in subkey.GetName():
                    continue

                print('Performing tests on:', dirname, subkey.GetName())

                histo = directory.Get(subkey.GetName())
                histo_sym = directory.Get(subkey.GetName()+'_sym')
                histo_flat = directory.Get(subkey.GetName()+'_flat')
                histo_asym = directory.Get(subkey.GetName()+'_asym')

                if histo_sym and 'flat' not in subkey.GetName():
                    if 'H_ps' in subkey.GetName(): continue # since we use reweighting we don't do the test for both sm and ps
                    test_results = PerformTests(histo, histo_sym, nxbins)
                    if dirname not in test_results_sym:
                        test_results_sym[dirname] = {}
                    test_results_sym[dirname][subkey.GetName()] = test_results
    #            if histo_flat: 
    #                test_results = PerformTests(histo, histo_flat, nxbins)
    #                if dirname not in test_results_flat:
    #                    test_results_flat[dirname] = {}
    #                test_results_flat[dirname][subkey.GetName()] = test_results
#
    #            if histo_asym and 'flat' not in subkey.GetName():
    #                test_results = PerformTests(histo, histo_asym, nxbins)
    #                if dirname not in test_results_asym:
    #                    test_results_asym[dirname] = {}
    #                test_results_asym[dirname][subkey.GetName()] = test_results

    print('\n\nResults of the tests for symmetrised histograms:')
    for dirname, results in test_results_sym.items():
        print(f"Directory: {dirname}")
        for histo_name, (chi2_pval, chi2_last_pval) in results.items():
            print(f"  Histogram: {histo_name}, Chi2 p-value: {chi2_pval}, Last bin Chi2 p-value: {chi2_last_pval}")
    #print('\n\nResults of the tests for flattened histograms:')
    #for dirname, results in test_results_flat.items():
    #    print(f"Directory: {dirname}")
    #    for histo_name, (chi2_pval, chi2_last_pval) in results.items():
    #        print(f"  Histogram: {histo_name}, Chi2 p-value: {chi2_pval}, Last bin Chi2 p-value: {chi2_last_pval}")



