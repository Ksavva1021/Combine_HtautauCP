import CombineHarvester.CombineTools.ch as ch
import yaml
# Note CMS common systematics should be named following: https://gitlab.cern.ch/cms-analysis/general/systematics/-/blob/master/systematics_master.yml?ref_type=heads, analysis specific ones should eventually have the CADI number in the names

def AddSMRun3Systematics(cb):

    cats_tt = {
        1:'tt_mva_tau',
        2:'tt_mva_fake',
        3:'tt_higgs_rhorho',
        4:'tt_higgs_rhoa11pr',
        5:'tt_higgs_rhoa1',
        6:'tt_higgs_a1a1',
        7:'tt_higgs_pirho',
        8:'tt_higgs_pipi',
        9:'tt_higgs_pia1',
        10: 'tt_higgs_pia11pr',
        11: 'tt_higgs_a11pra1',
        }

    eras = ['2022', '2022EE', '2023', '2023BPix']
    
    # define processes lists
    sig_procs = {}
    sig_procs['ggH'] = ['ggH_sm_prod_sm_htt','ggH_ps_prod_sm_htt','ggH_mm_prod_sm_htt', 'ggH_sm_prod_ps_htt','ggH_ps_prod_ps_htt','ggH_mm_prod_ps_htt', 'ggH_sm_prod_mm_htt','ggH_ps_prod_mm_htt','ggH_mm_prod_mm_htt']
    sig_procs['VBF'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt']
    sig_procs['ZH'] = ['ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']
    sig_procs['WH'] = ['WH_sm_htt','WH_ps_htt','WH_mm_htt']
    sig_procs['qqH'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt','WH_sm_htt','WH_ps_htt','WH_mm_htt','ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']
    
    dy_procs = ['ZTT', 'ZL']
    ttbar_procs = ['TTT']
    vv_procs = ['VVT']
    bkg_mc_procs = dy_procs + ttbar_procs + vv_procs #+ ['JetFakesSublead']
    
    mc_procs = bkg_mc_procs
    for p in sig_procs.values(): mc_procs+=p

    recoil_procs = ['ZTT','ZL']
    for p in sig_procs.values(): recoil_procs+=p
   
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

    # signal theory uncertainties

    # QCD scale variations
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "QCDscale_ren_signal_ACCEPT", "shape", ch.SystMap()(1.0))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "QCDscale_fac_signal_ACCEPT", "shape", ch.SystMap()(1.0))

    # parton shower variations
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "ps_isr_signal", "shape", ch.SystMap()(1.0))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "ps_fsr_signal", "shape", ch.SystMap()(1.0))

    # DY shape uncertainty (e.g from pT/mass reweighting)
    cb.cp().process(dy_procs).AddSyst(cb, "dy_pt_reweighting", "shape", ch.SystMap()(1.0))

    # ttbar shape uncertainty (e.g from pT reweighting)
    cb.cp().process(ttbar_procs).AddSyst(cb, "top_pt_reweighting", "shape", ch.SystMap()(1.0))
    
    ###############################################
    # Offline object identification
    ###############################################

    # TODO: muon ID
    
    # TODO: electron ID
    
    # TODO: add uncerts for l->tau fakes as well
    # TODO: add specific bins for et and mt channels as well

    for era in eras:
        # statistical uncertainties from fitted function parameters
        for u in ['stat1','stat2']:
            for dm in ['0', '1', '2', '10']:
                cb.cp().process(mc_procs).process(['ZL'], False).bin_id([1,2]).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_{u}_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))
                cb.cp().process(['ZL']).channel(['tt']).bin_id([1,2]).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_{u}_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))

            cb.cp().process(mc_procs).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_{u}_DM0PNet_{era}', 'shape', ch.SystMap()(1.0))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([3,4,5,7]).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_{u}_DM1PNet_{era}', 'shape', ch.SystMap()(1.0))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([4,10,11]).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_{u}_DM2PNet_{era}', 'shape', ch.SystMap()(1.0))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_{u}_DM10PNet_{era}', 'shape', ch.SystMap()(1.0))

        # systematic that is correlated across decay modes but decorrelated across eras
        cb.cp().process(mc_procs).process(['ZL'], False).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_syst_{era}', 'shape', ch.SystMap()(1.0))
        cb.cp().process(['ZL']).channel(['tt']).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_syst_{era}', 'shape', ch.SystMap()(1.0))

    # systematic that is correlated across eras and decay modes
    cb.cp().process(mc_procs).process(['ZL'], False).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_syst_alleras', 'shape', ch.SystMap()(1.0))
    cb.cp().process(['ZL']).channel(['tt']).AddSyst(cb, f'CMS_eff_t_DeepTau2018v2p5_VSjet_dm_syst_alleras', 'shape', ch.SystMap()(1.0))

    # TODO: btag ID (including for fakes)
   
    ###############################################
    # Trigger
    ###############################################

    # TODO: muon trigger
    
    # TODO: electron trigger

    # statistical uncertainties
    for era in eras:
        # tau leg uncertainties
        for trig in ['ditau','ditaujet']:
            for dm in ['0', '1', '2', '10']:
                cb.cp().process(mc_procs).channel(['tt']).bin_id([1,2]).AddSyst(cb, f'CMS_trig_t_{trig}_VTight_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, f'CMS_trig_t_{trig}_VTight_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([3,4,5,7]).AddSyst(cb, f'CMS_trig_t_{trig}_VTight_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([4,10,11]).AddSyst(cb, f'CMS_trig_t_{trig}_VTight_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, f'CMS_trig_t_{trig}_VTight_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))

        # jet leg uncertainties
        cb.cp().process(mc_procs).channel(['tt']).AddSyst(cb, f'CMS_trig_j_ditaujet_{era}', 'shape', ch.SystMap()(1.0))

    # We also add a 3% systematic uncertainty due to the background modelling in the SF extraction based on the studies in https://indico.cern.ch/event/1263107/contributions/5306043/attachments/2606862/4503028/tautriggerSF_checks.pdf
    cb.cp().process(mc_procs).channel(['tt']).AddSyst(cb, "CMS_trig_t_ditau_syst", "lnN", ch.SystMap()(1.03))
    
    ###############################################
    # Lepton/Tau energy scales
    ###############################################

    # TODO: electron ES
    
    # TODO: muon ES ? (not included for Run-2 as it was small)

    for era in eras:
        for dm in ['0', '1', '2', '10']:
            cb.cp().process(mc_procs).process(['ZL'], False).bin_id([1,2]).AddSyst(cb, f'CMS_scale_t_DeepTau2018v2p5_DM{dm}PNet_{era}_genTau', 'shape', ch.SystMap()(1.0))

        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, f'CMS_scale_t_DeepTau2018v2p5_DM0PNet_{era}_genTau', 'shape', ch.SystMap()(1.0))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([3,4,5,7]).AddSyst(cb, f'CMS_scale_t_DeepTau2018v2p5_DM1PNet_{era}_genTau', 'shape', ch.SystMap()(1.0))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([4,10,11]).AddSyst(cb, f'CMS_scale_t_DeepTau2018v2p5_DM2PNet_{era}_genTau', 'shape', ch.SystMap()(1.0))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, f'CMS_scale_t_DeepTau2018v2p5_DM10PNet_{era}_genTau', 'shape', ch.SystMap()(1.0))


    #TODO: add for mt and et specific bins as well         
    #TODO: add the ones for the l->tau fakes for mt and et channels (ZL only)

    
    
    ###############################################
    # Jet/MET scale/resolutions 
    ###############################################

    # TODO: Might need to split by eras or even into all the sources

    cb.cp().process(mc_procs).AddSyst(cb, 'CMS_scale_j_13p6TeV', 'shape', ch.SystMap()(1.0))

    cb.cp().process(mc_procs).AddSyst(cb, 'CMS_res_j_13p6TeV', 'shape', ch.SystMap()(1.0))
    
    # TODO: MET uncl
    
    # MET recoil uncertainties

    cb.cp().process(recoil_procs).AddSyst(cb,'CMS_scale_met', 'shape', ch.SystMap()(1.0))
    cb.cp().process(recoil_procs).AddSyst(cb,'CMS_res_met', 'shape', ch.SystMap()(1.0))
    
    ###############################################
    # jet->tau fake-factors
    ###############################################

    # TODO: FF uncertainties - need to be added for et and mt channels as well

    # tt channel statistical uncertainties
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,7,8,9,10]).AddSyst(cb, "ff_tt_stat_dm0", "shape", ch.SystMap()(1.0))
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,3,4,5,7]).AddSyst(cb, "ff_tt_stat_dm1", "shape", ch.SystMap()(1.0))
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,4,10,11]).AddSyst(cb, "ff_tt_stat_dm2", "shape", ch.SystMap()(1.0))
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,5,6,9,11]).AddSyst(cb, "ff_tt_stat_dm10", "shape", ch.SystMap()(1.0))

    # add lnN FF uncertainty from yml file located in configs/ff_lnN_uncertainties.yml
    # open yml file and read uncertainties
    with open('configs/ff_lnN_uncertainties.yml', 'r') as f:
        ff_lnN_uncertainties = yaml.safe_load(f)
        cb.cp().process(['JetFakes']).channel(['tt']).AddSyst(cb, "ff_tt_syst_lnN", "lnN", ch.SystMap()(ff_lnN_uncertainties['tt']['correlated']))
        for i in range(1, 12):
            # add lnN uncertainty for each decay mode
            cb.cp().process(['JetFakes']).channel(['tt']).bin_id([i]).AddSyst(cb, "ff_tt_syst_lnN_$BIN", "lnN", ch.SystMap()(ff_lnN_uncertainties['tt'][cats_tt[i]]))

    # add shape uncertainties for BDT score, this is decorrelated between tau, fake, and higgs categories
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1]).AddSyst(cb, "ff_tt_syst_BDTshape_tau", "shape", ch.SystMap()(1.0))
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([2]).AddSyst(cb, "ff_tt_syst_BDTshape_fake", "shape", ch.SystMap()(1.0))
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2], False).AddSyst(cb, "ff_tt_syst_BDTshape_higgs", "shape", ch.SystMap()(1.0))

    # add shape uncertainties for aco angle, this is decorrelated between categories
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2], False).AddSyst(cb, "ff_tt_syst_acoshape_$BIN", "shape", ch.SystMap()(1.0))

    # uncertainty due to subtracted real taus in the tt channel
    cb.cp().process(['JetFakes']).channel(['tt']).AddSyst(cb, "ff_tt_sub_syst", "shape", ch.SystMap()(1.0))
    # lnN uncertainty for the JetFakesSublead as it is estimated from MC
    cb.cp().process(['JetFakesSublead']).channel(['tt']).AddSyst(cb, "ff_tt_mc", "lnN", ch.SystMap()(1.3))

    ###############################################
    # 4-vectors for CP angle reconstruction
    ###############################################
    
    # IP direction/scale
    cb.cp().process(mc_procs).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, "CMS_res_IP", "shape", ch.SystMap()(1.0))
    
    # TODO: pi0 direction/scale (not included for Run-2 but could add)
    
    # TODO: pi direction/scale (not included for Run-2 but could add)
    
    # SV vertex resolution uncertainty
    cb.cp().process(mc_procs).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, "CMS_res_SV", "shape", ch.SystMap()(1.0))


    ###############################################
    # DM-migration uncertainties
    ###############################################

    #TODO: migration uncertainties for migrations between reco-decay mode bins (not included for Run-2 but we could add)

    return cb
