import ROOT
import math
import argparse
ROOT.TH1.AddDirectory(False)
from ctypes import c_double
import os
import numpy as np
ROOT.gROOT.SetBatch(True)
 
# note that the merging of bins requires an even numnber of phi_CP bins, and this number must be set to the specific value used in the dictionary below otherwise the method will give incorrect results!

cp_bins = {
        "tt_mva_fake": 1,
        "tt_mva_tau": 1,
        "tt_mva_higgs": 1,
        "tt_higgs_pia1" : 4,
        "tt_higgs_a11pra1" : 4,
        "tt_higgs_rhoa1" : 10,
        "tt_higgs_pia11pr" : 4,
        "tt_higgs_pirho" : 10,
        "tt_higgs_a1a1" : 4,
        "tt_higgs_pipi" : 4,
        "tt_higgs_rhoa11pr" : 10,
        "tt_higgs_rhorho" : 10,
        "mt_mva_tau": 1,
        "mt_mva_fake": 1,
        "mt_mva_higgs_mua1": 10,
        "mt_mva_higgs_mua11pr": 8,
        "mt_mva_higgs_mupi": 8,
        "mt_mva_higgs_murho": 10,
        "et_mva_tau": 1,
        "et_mva_fake": 1,
        "et_mva_higgs_ea1": 10,
        "et_mva_higgs_ea11pr": 8,
        "et_mva_higgs_epi": 8,
        "et_mva_higgs_erho": 10
}

test_results_sym = {}
test_results_flat = {}
test_results_asym = {}

ff_aiso_yields = {}
extra_hists = {}

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

def MergeYBins(hist, nxbins):
  histnew = hist.Clone()
  nbins = hist.GetNbinsX()

  nybins = int(nbins/nxbins)

  for i in range(1,nxbins+1):
    tot = 0.
    tot_err = 0.
    for j in range(1,nybins+1):
      bin_num = i+(j-1)*nxbins
      tot += hist.GetBinContent(bin_num)
      tot_err += hist.GetBinError(bin_num)**2
    tot_err = math.sqrt(tot_err)
    # now set all j bins to this value
    for j in range(1,nybins+1):
      bin_num = i+(j-1)*nxbins
      histnew.SetBinContent(bin_num, tot)
      histnew.SetBinError(bin_num, tot_err)

  return histnew


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

def GetDataAndEstFakes(dirname, infile):
    directory = infile.Get(dirname)
    # check if directory exists
    if not directory:
        print('Directory', dirname, 'not found in file', infile.GetName())
        return None, None
    data = directory.Get('data_obs')
    backgrounds = ['ZTT','ZL','TTT','VVT','JetFakesSublead']
    bkg = directory.Get(backgrounds[0])
    for b in backgrounds[1:]:
        bkg.Add(directory.Get(b))
    data.Add(bkg,-1)
    fake_est = directory.Get('JetFakes')
    if not data or not fake_est:
        print('Data or fake estimate not found in directory', dirname)
        return
    print(data, fake_est)
    return data, fake_est

def GetFFUncerts(dirname, infile):
    print('Getting FF uncertainties for:', dirname)
    dirname_aiso = dirname+'_aiso'

    extra_hists[dirname] = []

    data, fake_est = GetDataAndEstFakes(dirname_aiso, infile)

    ff_aiso_yields[dirname] = (data.Integral(), fake_est.Integral())

    if not data or not fake_est:
        return

    # now we estimate shape uncertainties for the BDT score and the aco-angle
    nxbins = cp_bins[dirname] if dirname in cp_bins else 1

    directory = infile.Get(dirname)
    fake_est_iso = directory.Get('JetFakes')
    rate = fake_est_iso.Integral()

    name_extra=""
    if 'higgs' in dirname:
        name_extra = '_higgs'
    elif 'tau' in dirname:
        name_extra = '_tau'
    elif 'fake' in dirname:
        name_extra = '_fake'

    for x in ['BDT', 'aco']:

        if x == 'BDT':
            # for the BDT shape uncertainties we rebin by nxbins so that the uncertainty does only depends on the BDT score and not the aco angle
            data_merged, _, _ = MergeXBins(data, nxbins)
            fake_est_merged, _, _ = MergeXBins(fake_est, nxbins)
            name_extra=""
            if 'higgs' in dirname:
                name_extra = '_higgs'
            elif 'tau' in dirname:
                name_extra = '_tau'
            elif 'fake' in dirname:
                name_extra = '_fake'
        elif x == 'aco':
            # for the aco angle shape uncertainties we rebin by nybins so that the uncertainty does only depends on the aco angle and not the BDT score
            data_merged = MergeYBins(data, nxbins)
            fake_est_merged = MergeYBins(fake_est, nxbins)
            name_extra = '_'+dirname
            
        uncert_up = fake_est_iso.Clone()
        if uncert_up.GetNbinsX() != data_merged.GetNbinsX():
            print('Number of bins in fake estimate and fake estimate iso do not match, cannot calculate uncertainties!')
            return

        uncert_up.Multiply(data_merged)
        uncert_up.Divide(fake_est_merged)
        uncert_up.Scale(rate/uncert_up.Integral())
        uncert_down = fake_est_iso.Clone()
        uncert_down.Multiply(fake_est_merged)
        uncert_down.Divide(data_merged)
        uncert_down.Scale(rate/uncert_down.Integral())

        uncert_up.SetName(f'JetFakes_ff_tt_syst_{x}shape{name_extra}Up')
        uncert_down.SetName(f'JetFakes_ff_tt_syst_{x}shape{name_extra}Down')
    
        extra_hists[dirname] += [uncert_up, uncert_down]

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname):
    directory = infile.Get(dirname)

    histos = []
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo, ROOT.TH1D) or isinstance(histo, ROOT.TH1F):
            histos.append(histo)
    # add extra_hists to histos
    if dirname in extra_hists:
        for histo in extra_hists[dirname]:
            if isinstance(histo, ROOT.TH1D) or isinstance(histo, ROOT.TH1F):
                histos.append(histo)


    for histo in histos:
        #print('Processing:', dirname, histo.GetName()) 
        if dirname in cp_bins: nxbins = cp_bins[dirname]
        else: nxbins = 1  
        # we write all the old histograms to the output file
        outfile.cd()
        if not ROOT.gDirectory.GetDirectory(dirname): ROOT.gDirectory.mkdir(dirname)
        ROOT.gDirectory.cd(dirname)
        histo.Write()    
        if nxbins== 1: continue    
        # if data we skip
        if 'data_obs' in histo.GetName(): continue    
        # if mm signal then we only anti-symmetrise
        elif 'H_mm' in histo.GetName() and 'htt125' in histo.GetName() and nxbins>1:
            hsm = directory.Get(histo.GetName().replace('H_mm_','H_sm_'))
            hps = directory.Get(histo.GetName().replace('H_mm_','H_ps_'))
            histo_asym = ASymmetrise(histo,hsm,hps,nxbins)
            histo_asym.SetName(histo.GetName()+'_sym')
            histo_asym.Write()
            continue    
        # if not mm signal then we always symmetrise
        else:
            histo_sym, p_val_total, p_val_perbin = Symmetrise(histo,nxbins)
            if dirname not in test_results_sym:
                test_results_sym[dirname] = {}
            test_results_sym[dirname][histo.GetName()] = (p_val_total, p_val_perbin[-1])
            histo_sym.SetName(histo.GetName()+'_sym')
            histo_sym.Write()    
            # if not signal then we also flatten
            if 'htt125' not in histo.GetName() or 'Higgs_flat' in histo.GetName() or 'H_flat' in histo.GetName():
                histo_flat, p_val_total, p_val_perbin = MergeXBins(histo,nxbins)
                if dirname not in test_results_flat:
                    test_results_flat[dirname] = {}
                test_results_flat[dirname][histo.GetName()] = (p_val_total, p_val_perbin[-1])
                histo_flat.SetName(histo.GetName()+'_flat')
                histo_flat.Write()

        ROOT.gDirectory.cd('/')

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which we want to merge bins')
parser.add_argument('--test', '-t', action='store_true', help= 'Run statistical tests on the histograms to check compatability of smoothed and unsmoothed histograms')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','-mergeXbins.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename,"RECREATE") 

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        if dirname in cp_bins and dirname.startswith('tt_') and not dirname.startswith("tt_mva_higgs"): GetFFUncerts(dirname, original_file)
        getHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname)

if 'tt_mva_higgs' in ff_aiso_yields and 'tt_mva_tau' in ff_aiso_yields and 'tt_mva_fake' in ff_aiso_yields:



    # now we estimate the normalization uncertainties for the FFs
    # get the correlated uncertainty from summing mva_higgs, mva_tau, and mva_fake
    data_total = 0.
    fake_total = 0.
    for x in ['tt_mva_higgs', 'tt_mva_tau', 'tt_mva_fake']:
        if x in ff_aiso_yields:
            data_total += ff_aiso_yields[x][0]
            fake_total += ff_aiso_yields[x][1]
    uncert_corr = data_total/fake_total

    # create a yml file to store the lnN uncertainties for the FFs
    with open('configs/ff_lnN_uncertainties.yml', 'w') as f:

        f.write('tt:')
        f.write(f'\n  correlated: {abs(1-uncert_corr)+1:.3f}\n')
    

        # now loop over 'tt' categories and get uncertainty per bin
        # in this case we scale the by the uncert_corr so that we don't double count the correlated uncertainty
        for dirname in ff_aiso_yields:
            if dirname.startswith('tt_'):
                data, fake_est = ff_aiso_yields[dirname]
                uncert = data / fake_est / (uncert_corr)
                uncert = abs(uncert - 1)+1
                f.write(f'  {dirname}: {uncert:.3f}\n')


if args.test:

    if not os.path.exists('test_results'):
        os.makedirs('test_results')

    # make ROOT plots of the p-values
    # 1D plot just showing the p-values to show they are uniform
    # 2D plot showing processes on the x-axis, categories on the y-axis and p-values as the z-axis
    gr_sym_id = ROOT.TGraph()
    gr_sym_last_id = ROOT.TGraph()
    gr_flat_id = ROOT.TGraph()
    gr_flat_last_id = ROOT.TGraph()

    hist_sym_1d = ROOT.TH1D('hist_sym_1d','Symmetrisation tests: all BDT bins',20,0,1)
    hist_flat_1d = ROOT.TH1D('hist_flat_1d','Flattening tests: all BDT bins',20,0,1)
    hist_sym_last_1d = ROOT.TH1D('hist_sym_last_1d','Symmetrisation tests: last BDT bin',20,0,1)
    hist_flat_last_1d = ROOT.TH1D('hist_flat_last_1d','Flattening tests: last BDT bin',20,0,1)

    N_categories = len(test_results_sym)
    hist_sym_2d = ROOT.TH2D('hist_sym_2d','Symmetrisation tests: all BDT bins',N_categories,0,N_categories,6,0,6)
    hist_flat_2d = ROOT.TH2D('hist_flat_2d','Flattening tests: all BDT bins',N_categories,0,N_categories,6,0,6)
    hist_sym_last_2d = ROOT.TH2D('hist_sym_last_2d','Symmetrisation tests: last BDT bin',N_categories,0,N_categories,6,0,6)
    hist_flat_last_2d = ROOT.TH2D('hist_flat_last_2d','Flattening tests: lst BDT bin',N_categories,0,N_categories,6,0,6)

    hist_sym_2d.GetZaxis().SetTitle('P-value')
    hist_flat_2d.GetZaxis().SetTitle('P-value')
    hist_sym_last_2d.GetZaxis().SetTitle('P-value (last BDT bin)')
    hist_flat_last_2d.GetZaxis().SetTitle('P-value (last BDT bin)')

    hist_sym_2d.GetYaxis().SetBinLabel(1,'ZTT')
    hist_sym_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    hist_sym_2d.GetYaxis().SetBinLabel(3,'qqH_sm_htt125')
    hist_sym_2d.GetYaxis().SetBinLabel(4,'ggH_sm_prod_sm_htt125')
    hist_sym_2d.GetYaxis().SetBinLabel(5,'ggH_ps_prod_ps_htt125')
    hist_sym_2d.GetYaxis().SetBinLabel(6,'ggH_flat_prod_mm_htt125')

    hist_flat_2d.GetYaxis().SetBinLabel(1,'ZTT')
    hist_flat_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    hist_flat_2d.GetYaxis().SetBinLabel(3,'qqH_flat_htt125')
    hist_flat_2d.GetYaxis().SetBinLabel(4,'ggH_flat_prod_sm_htt125')
    hist_flat_2d.GetYaxis().SetBinLabel(5,'ggH_flat_prod_ps_htt125')
    hist_flat_2d.GetYaxis().SetBinLabel(6,'ggH_flat_prod_mm_htt125')

    hist_sym_last_2d.GetYaxis().SetBinLabel(1,'ZTT')
    hist_sym_last_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    hist_sym_last_2d.GetYaxis().SetBinLabel(3,'qqH_sm_htt125')
    hist_sym_last_2d.GetYaxis().SetBinLabel(4,'ggH_sm_prod_sm_htt125')
    hist_sym_last_2d.GetYaxis().SetBinLabel(5,'ggH_ps_prod_ps_htt125')
    hist_sym_last_2d.GetYaxis().SetBinLabel(6,'ggH_flat_prod_mm_htt125')

    hist_flat_last_2d.GetYaxis().SetBinLabel(1,'ZTT')
    hist_flat_last_2d.GetYaxis().SetBinLabel(2,'JetFakes')
    hist_flat_last_2d.GetYaxis().SetBinLabel(3,'qqH_flat_htt125')
    hist_flat_last_2d.GetYaxis().SetBinLabel(4,'ggH_flat_prod_sm_htt125')
    hist_flat_last_2d.GetYaxis().SetBinLabel(5,'ggH_flat_prod_ps_htt125')
    hist_flat_last_2d.GetYaxis().SetBinLabel(6,'ggH_flat_prod_mm_htt125')

    hist_sym_2d.SetStats(0)
    hist_flat_2d.SetStats(0)
    hist_sym_last_2d.SetStats(0)
    hist_flat_last_2d.SetStats(0)

    hist_sym_1d.GetXaxis().SetTitle('P-value')
    hist_flat_1d.GetXaxis().SetTitle('P-value')
    hist_sym_last_1d.GetXaxis().SetTitle('P-value (last BDT bin)')
    hist_flat_last_1d.GetXaxis().SetTitle('P-value (last BDT bin)')

    hist_sym_1d.SetStats(0)
    hist_flat_1d.SetStats(0)
    hist_sym_last_1d.SetStats(0)
    hist_flat_last_1d.SetStats(0)
    # make sure 1D keeps track of errors
    hist_sym_1d.Sumw2()
    hist_flat_1d.Sumw2()
    hist_sym_last_1d.Sumw2()
    hist_flat_last_1d.Sumw2()

    #print('\n\nSymmetrisation test results:')
    i = 0
    for dirname, results in test_results_sym.items():
        i+=1
        #print(f'Directory: {dirname}')
        # set x-labels
        hist_sym_2d.GetXaxis().SetBinLabel(i, dirname)
        hist_sym_last_2d.GetXaxis().SetBinLabel(i, dirname)
        for hist_name, (p_val_total, p_val_last) in results.items():
            #if not(hist_name in ['JetFakes', 'ZTT','Higgs_flat_htt125'] or hist_name.startswith('ggH') or hist_name.startswith('qqH_sm')) or 'unfiltered' in hist_name: continue
            #print(f'  Histogram: {hist_name}, P-value total: {p_val_total:.4f}, P-value last bin: {p_val_last:.4f}')
            j = None
            if hist_name == 'ZTT': j = 1
            if hist_name == 'JetFakes': j = 2
            if hist_name == 'qqH_sm_htt125': j = 3
            if hist_name == 'ggH_sm_prod_sm_htt125': j = 4
            if hist_name == 'ggH_ps_prod_ps_htt125': j = 5
            if hist_name == 'ggH_flat_prod_mm_htt125': j = 6

            if j is None: continue

            hist_sym_2d.SetBinContent(i,j,p_val_total)
            hist_sym_last_2d.SetBinContent(i,j,p_val_last)

            if j: # only use statistically independent processes for 1D plot
              gr_sym_id.SetPoint(gr_sym_id.GetN(), p_val_total, 1)
              gr_sym_last_id.SetPoint(gr_sym_last_id.GetN(), p_val_last,1)
              hist_sym_1d.Fill(p_val_total)
              hist_sym_last_1d.Fill(p_val_last)

    #print('\n\nFlattening test results:')
    i=0
    for dirname, results in test_results_flat.items():
        i+=1
        #print(f'Directory: {dirname}')
        # set x-labels
        hist_flat_2d.GetXaxis().SetBinLabel(i, dirname)
        hist_flat_last_2d.GetXaxis().SetBinLabel(i, dirname)
        for hist_name, (p_val_total, p_val_last) in results.items():
            #if hist_name not in ['JetFakes', 'ZTT', 'qqH_flat_htt125', 'ggH_flat_prod_sm_htt125', 'ggH_flat_prod_ps_htt125', 'ggH_flat_prod_mm_htt125']: continue
            #print(f'  Histogram: {hist_name}, P-value total: {p_val_total:.4f}, P-value last bin: {p_val_last:.4f}')

            j = None
            if hist_name == 'ZTT': j = 1
            if hist_name == 'JetFakes': j = 2
            if hist_name == 'qqH_flat_htt125': j = 3
            if hist_name == 'ggH_flat_prod_sm_htt125': j = 4
            if hist_name == 'ggH_flat_prod_ps_htt125': j = 5
            if hist_name == 'ggH_flat_prod_mm_htt125': j = 6

            if j is None: continue

            print(hist_flat_2d)
            print(i,j,p_val_total)
            hist_flat_2d.SetBinContent(i,j,p_val_total)
            hist_flat_last_2d.SetBinContent(i,j,p_val_last)

            if j: # only use statistically independent processes for 1D plot
              hist_flat_1d.Fill(p_val_total)
              hist_flat_last_1d.Fill(p_val_last)
              gr_flat_id.SetPoint(gr_flat_id.GetN(), p_val_total, 1)
              gr_flat_last_id.SetPoint(gr_flat_last_id.GetN(), p_val_last, 1)

    fout = ROOT.TFile('test_results/test_results.root', 'RECREATE')
    fout.cd()
    hist_sym_2d.Write()
    hist_flat_2d.Write()
    hist_sym_last_2d.Write()
    hist_flat_last_2d.Write()
    gr_sym_id.Write('gr_sym_id')
    gr_sym_last_id.Write('gr_sym_last_id')
    gr_flat_id.Write('gr_flat_id')
    gr_flat_last_id.Write('gr_flat_last_id')
    hist_sym_1d.Write()
    hist_flat_1d.Write()
    hist_sym_last_1d.Write()
    hist_flat_last_1d.Write()

    c_2d = ROOT.TCanvas('', '', 1200, 600)
    # decrease number of sf shown in plot
    ROOT.gStyle.SetPaintTextFormat(".3g")
    ROOT.gStyle.SetPalette(ROOT.kBlackBody)

    #increase L and R margins
    c_2d.SetLeftMargin(0.15)
    c_2d.SetRightMargin(0.13)
    hist_sym_2d.GetZaxis().SetRangeUser(0,1)
    hist_sym_2d.Draw('COLZ TEXT')
    c_2d.Print('test_results/chi2tests_symmetrisation_2d.pdf')

    hist_sym_last_2d.GetZaxis().SetRangeUser(0,1)
    hist_sym_last_2d.Draw('COLZ TEXT')
    c_2d.Print('test_results/chi2tests_symmetrisation_last_2d.pdf')

    hist_flat_2d.GetZaxis().SetRangeUser(0,1)
    hist_flat_2d.Draw('COLZ TEXT')
    c_2d.Print('test_results/chi2tests_flattening_2d.pdf')

    hist_flat_last_2d.GetZaxis().SetRangeUser(0,1)
    hist_flat_last_2d.Draw('COLZ TEXT')
    c_2d.Print('test_results/chi2tests_flattening_last_2d.pdf')

    c_1d = ROOT.TCanvas('', '', 800, 600)
    hist_sym_1d.Draw()
    c_1d.Print('test_results/chi2tests_symmetrisation_1d.pdf')

    hist_flat_1d.Draw()
    c_1d.Print('test_results/chi2tests_flattening_1d.pdf')

    hist_sym_last_1d.Draw()
    c_1d.Print('test_results/chi2tests_symmetrisation_last_1d.pdf')

    hist_flat_last_1d.Draw()
    c_1d.Print('test_results/chi2tests_flattening_last_1d.pdf')
