import CombineHarvester.CombineTools.ch as ch
# Note CMS common systematics should be named following: https://gitlab.cern.ch/cms-analysis/general/systematics/-/blob/master/systematics_master.yml?ref_type=heads, analysis specific ones should eventually have the CADI number in the names

def AddSMRun3Systematics(cb):
    
    # define processes lists
    sig_procs = {}
    sig_procs['ggH'] = ['ggH_sm_htt','ggH_ps_htt','ggH_mm_htt']
    sig_procs['VBF'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt']
    sig_procs['ZH'] = ['ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']
    sig_procs['WH'] = ['WH_sm_htt','WH_ps_htt','WH_mm_htt']
    sig_procs['qqH'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt','WH_sm_htt','WH_ps_htt','WH_mm_htt','ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']
    
    dy_procs = ['ZTT', 'ZL']
    ttbar_procs = ['TTT']
    vv_procs = ['VVT']
    bkg_mc_procs = dy_procs + ttbar_procs + vv_procs
    
    mc_procs = bkg_mc_procs
    for p in sig_procs.values(): mc_procs+=p
   
    ###############################################
    # Luminosity
    ###############################################

    # lumi uncertainty from here: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3
    # 1.3% recomended for 2023 while 1.4% recomended for 2022 - used the largest number = 1.4% 
    cb.cp().process(mc_procs).AddSyst(cb, 'lumi_13p6TeV', 'lnN', ch.SystMap()(1.014))


    ###############################################
    # Pileup
    ###############################################

    # TODO: pileup? (not included for Run-2, but we could add it)

    ###############################################
    # Cross sections and BRs
    ###############################################

    # Cross-sections uncertainties - we keep naming consistent with Run-2 analysis for now

    # DY XS uncertainties from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV
    # Quadrature sum of scale, PDF, and difference between NLO additive vs multiplicative
    cb.cp().process(dy_procs).AddSyst(cb, 'CMS_htt_zjXsec_13p6TeV', 'lnN', ch.SystMap()((0.984,1.013)))
    
    # ttbar cross-section uncertainties from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    # Quadrature sum of scale, PDF, and mass uncerts
    cb.cp().process(ttbar_procs).AddSyst(cb, 'CMS_htt_tjXsec_13p6TeV', 'lnN', ch.SystMap()((0.949,1.044)))

    #For VV in principle can take NNLO numbers from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV
    #But since we mix together all VV + rare procs into this template, we take a conservative 5% uncertainty (same as for Run-2)
    cb.cp().process(ttbar_procs).AddSyst(cb, 'CMS_htt_vvXsec_13p6TeV', 'lnN', ch.SystMap()(1.05))

    # Higgs cross-section uncertainties from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap

    # QCD scale uncertainties
    cb.cp().process(sig_procs['ggH']).AddSyst(cb, 'QCDscale_ggH', 'lnN', ch.SystMap()(1.039))
    
    cb.cp().process(sig_procs['WH']).AddSyst(cb, 'QCDscale_VH', 'lnN', ch.SystMap()((0.993,1.004)))
    
    cb.cp().process(sig_procs['ZH']).AddSyst(cb, 'QCDscale_VH', 'lnN', ch.SystMap()((0.968,1.038)))
    
    # PDF uncertainties
    cb.cp().process(sig_procs['ggH']).AddSyst(cb, 'pdf_Higgs_gg', 'lnN', ch.SystMap()(1.032))
    
    cb.cp().process(sig_procs['VBF']).AddSyst(cb, 'pdf_Higgs_qqbar', 'lnN', ch.SystMap()(1.032))
    
    cb.cp().process(sig_procs['WH']).AddSyst(cb, 'pdf_Higgs_qqbar', 'lnN', ch.SystMap()(1.016))
    
    cb.cp().process(sig_procs['ZH']).AddSyst(cb, 'pdf_Higgs_qqbar', 'lnN', ch.SystMap()(1.013))
    
    # H->tautau BR uncertainties from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR#SM_Higgs_Branching_Ratios_and_To
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']+sig_procs['WH']+sig_procs['ZH']).AddSyst(cb, 'BR_htt_THU', 'lnN', ch.SystMap()((0.984,1.017)))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']+sig_procs['WH']+sig_procs['ZH']).AddSyst(cb, 'BR_htt_PU_mq', 'lnN', ch.SystMap()(1.010))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']+sig_procs['WH']+sig_procs['ZH']).AddSyst(cb, 'BR_htt_PU_alphas', 'lnN', ch.SystMap()(1.006))
    
    
    ###############################################
    # Shape/acceptance theory uncertainties
    ###############################################

    # TODO: signal theory uncertainties
    
    # TODO: DY shape uncertainty (e.g from pT/mass reweighting)
    
    # TODO: ttbar shape uncertainty (e.g from pT reweighting)
    
    ###############################################
    # Offline object identification
    ###############################################

    # TODO: muon ID
    
    # TODO: electron ID
    
    # TODO: tau ID (including for l fakes)
    
    # TODO: btag ID (including for fakes)
   
    ###############################################
    # Trigger
    ###############################################

    # TODO: muon trigger
    
    # TODO: electron trigger
    
    # TODO: tau trigger
    
    ###############################################
    # Lepton/Tau energy scales
    ###############################################

    # TODO: electron ES
    
    # TODO: muon ES ? (not included for Run-2 as it was small)
    
    # TODO: tau ES (including for l fakes)
    
    
    ###############################################
    # Jet/MET scale/resolutions 
    ###############################################

    # TODO: JES
    
    # TODO: JER
    
    # TODO: MET uncl
    
    # TODO: MET recoil
    
    ###############################################
    # jet->tau fake-factors
    ###############################################

    # TODO: FF uncertainties

    ###############################################
    # 4-vectors for CP angle reconstruction
    ###############################################
    
    # TODO: IP direction/scale (only significance included for Run-2)
    
    # TODO: pi0 direction/scale (not included for Run-2 but could add)
    
    # TODO: pi direction/scale (not included for Run-2 but could add)
    
    # TODO: SV vertex direction/scale/efficiency (only efficience added for Run-2)

    ###############################################
    # DM-migration uncertainties
    ###############################################

    #TODO: migration uncertainties for migrations between reco-decay mode bins (not included for Run-2 but we could add)

    return cb
