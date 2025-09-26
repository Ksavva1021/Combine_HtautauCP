# generate toys for alpha = 0 and 20

# generate toys
for mode in 0 1 2 3 4; do
  combineTool.py -m 125 -M GenerateOnly --setParameters muV=1,alpha=0,muggH=1,mutautau=1 -d outputs_mergemode${mode}/cmb/ws.root -t 50 --there -n .alpha0.mergemode${mode}Toys50 --saveToys -s 0:19:1
  combineTool.py -m 125 -M GenerateOnly --setParameters muV=1,alpha=20,muggH=1,mutautau=1 -d outputs_mergemode${mode}/cmb/ws.root -t 50 --there -n .alpha20.mergemode${mode}Toys50 --saveToys -s 0:19:1
done

# fit toys with different models
for mode1 in 0 1 2 3 4; do
  for mode2 in 0 1 2 3 4; do
      echo mode1 = ${mode1}, mode2 = ${mode2}
      for s in {0..19}; do
          combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=0,muggH=1,mutautau=1 -d outputs_mergemode${mode1}/cmb/ws.root -t 50 --there -n .alpha0.mergemode${mode2}Toys.mergemode${mode1}Fit --toysFile ../../outputs_mergemode${mode2}/cmb/higgsCombine.alpha0.mergemode${mode2}Toys50.GenerateOnly.mH125.${s}.root -s ${s} --algo singles --saveNLL --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10 --job-mode condor --task-name condor-alpha0_mergemode${mode2}Toys_mergemode${mode1}Fit_${s} --sub-opts='+MaxRuntime=10800'
          combineTool.py -m 125 -M MultiDimFit --setParameters muV=1,alpha=20,muggH=1,mutautau=1 -d outputs_mergemode${mode1}/cmb/ws.root -t 50 --there -n .alpha20.mergemode${mode2}Toys.mergemode${mode1}Fit --toysFile ../../outputs_mergemode${mode2}/cmb/higgsCombine.alpha20.mergemode${mode2}Toys50.GenerateOnly.mH125.${s}.root -s ${s} --algo singles --saveNLL --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1 --cminFallbackAlgo Minuit2,Migrad,0:1 --cminFallbackAlgo Minuit2,Migrad,0:2 --cminFallbackAlgo Minuit2,Migrad,0:4 --cminFallbackAlgo Minuit2,Migrad,0:10 --job-mode condor --task-name condor-alpha20_mergemode${mode2}Toys_mergemode${mode1}Fit_${s} --sub-opts='+MaxRuntime=10800'
      done
  done
done

## hadd results (after jobs have finished)
#for mode1 in 0 1 2 3 4; do
#    for mode2 in 0 1 2 3 4; do
#        hadd -f outputs_mergemode${mode1}/cmb/higgsCombine.alpha0.mergemode${mode2}Toys.mergemode${mode1}Fit.MultiDimFit.mH125.root outputs_mergemode${mode1}/cmb/higgsCombine.alpha0.mergemode${mode2}Toys.mergemode${mode1}Fit.MultiDimFit.mH125.*.root
#        hadd -f outputs_mergemode${mode1}/cmb/higgsCombine.alpha20.mergemode${mode2}Toys.mergemode${mode1}Fit.MultiDimFit.mH125.root outputs_mergemode${mode1}/cmb/higgsCombine.alpha20.mergemode${mode2}Toys.mergemode${mode1}Fit.MultiDimFit.mH125.*.root
#    done
#done
