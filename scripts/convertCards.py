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
  for i in range(1,nbins+1,nxbins):
    tot_err = c_double(0)
    tot = hist.IntegralAndError(i,i+nxbins-1,tot_err)
    tot_err = tot_err.value
    for j in range(i,i+nxbins):
      histnew.SetBinContent(j,tot/nxbins)
      histnew.SetBinError(j,tot_err/nxbins)
  return histnew

def Symmetrise(hist,nxbins):
  histnew=hist.Clone()
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
  return histnew

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
                histo_sym = Symmetrise(histo,nxbins)
                histo_sym.SetName(histo.GetName()+'_sym')
                histo_sym.Write()

                # if not signal then we also flatten
                if 'htt125' not in key.GetName() or 'Higgs_flat' in key.GetName() or 'H_flat' in key.GetName():
                    histo_flat = MergeXBins(histo,nxbins)
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

if args.test:
    print('Running statistical tests on the histograms to check compatability of smoothed and unsmoothed histograms')
    for key in output_file.GetListOfKeys():
        if isinstance(output_file.Get(key.GetName()),ROOT.TDirectory):
            dirname=key.GetName()
            directory = output_file.Get(dirname)
            for subkey in directory.GetListOfKeys():
                if subkey.endswith('_sym') or subkey.endswith('_flat') or subkey.endswith('_asym'):
                    continue
                # skip systematics
                if subkey.endswith('Up') or subkey.endswith('Down'):
                    continue
                
               
                histo = directory.Get(subkey.GetName())
                # check if asym, sym, or flat hitograms exist for this one
                