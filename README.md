## Setup Instructions:


### Combine Harvester

This package requires HiggsAnalysis/CombinedLimit to be in your local CMSSW area. We follow the release recommendations of the combine developers. The CombineHarvester framework is compatible with the CMSSW 14_1_X and 11_3_X series releases. The default branch, main, is for developments in the 14_1_X releases, and the current recommended tag is v3.0.0. The v2.1.0 tag should be used in CMSSW_11_3_X.

A new full release area can be set up and compiled in the following steps:

```
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
cd CombineHarvester
git checkout v3.0.0
scram b -j8
```

### CombineHarvester Repository for the Htautau CP in decay measurement

```
cd $CMSSW_BASE/src/CombineHarvester
git clone git@github.com:Ksavva1021/Combine_HtautauCP.git Combine_HtautauCP
scram b -j8
```

### The magic command

```
ulimit -c 0
```

### Produce txt datacards

Modify the options configs/harvestDatacards.yml as needed and then run
 
```
python3 scripts/harvestDatacards.py 
```


### Creating the workspaces

```
combineTool.py -m 125 -M T2W -P CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays -i outputs/cmb -o ws.root --parallel 8
```

### Run maximum likelihood fits

1D fit for alpha:

```
combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 --setParameterRanges alpha=-90,90 --points 21 --redefineSignalPOIs alpha  -d outputs/cmb/ws.root --algo grid -t -1 --there -n .alpha --alignEdges 1
```

### make plot of alpha scan

```
python3 scripts/plot1DScan.py --main=outputs/cmb/higgsCombine.alpha.MultiDimFit.mH125.root --POI=alpha --output=alpha_cmb --no-numbers --no-box --x-min=-90 --x-max=90 --y-max=8
```
