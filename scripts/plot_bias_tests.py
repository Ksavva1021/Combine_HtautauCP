import ROOT
from array import array
ROOT.gROOT.SetBatch(True)


def plot_toys(tree1,tree2,title='',plot_name='bias_test.pdf'):
    #c1 = ROOT.TCanvas("c1","c1",800,600)
    c1 = ROOT.TCanvas()
    h1 = ROOT.TH1F("h1","h1",100,-90,90)
    h2 = ROOT.TH1F("h2","h2",100,-90,90)

    h1.Sumw2()
    h2.Sumw2()

    tree1.Draw("alpha>>h1(100,-90,90)","quantileExpected<-0.5")
    h1 = tree1.GetHistogram()
    tree2.Draw("alpha>>h2(100,-90,90)","quantileExpected<-0.5")
    h2 = tree2.GetHistogram()

    h1.SetStats(0)
    h2.SetStats(0)
    h1.SetTitle(title)
    h2.SetTitle(title)
    h1.GetXaxis().SetTitle('#phi_{CP} (#circ)')

    h1.SetLineColor(ROOT.kRed)
    h2.SetLineColor(ROOT.kBlue)
    h1.SetMaximum(1.2*max(h1.GetMaximum(),h2.GetMaximum()))
    h1.Draw("HIST")
    h2.Draw("HIST SAME")
    mean1 = h1.GetMean()
    mean2 = h2.GetMean()
    rms1 = h1.GetRMS()
    rms2 = h2.GetRMS()
    leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
    leg.AddEntry(h1,f"fit=toys mean={mean1:.1f}, rms={rms1:.1f}","l")
    leg.AddEntry(h2,f"fit#neq toys mean={mean2:.1f}, rms={rms2:.1f}","l")
    leg.Draw()
    c1.SaveAs(plot_name)

def plot_pull(tree1,tree2,title='',plot_name='bias_test_pull.pdf'):
    #c1 = ROOT.TCanvas("c1","c1",800,600)
    c1 = ROOT.TCanvas()
    h1_pull = ROOT.TH1F("h1_pull","h1_pull",100,-90,90)
    h2_pull = ROOT.TH1F("h2_pull","h2_pull",100,-90,90)

    h1_pull.Sumw2()
    h2_pull.Sumw2()

    print(tree1.GetEntries(), tree2.GetEntries())

    tree1.Draw("alpha/alpha_e>>h1_pull(100,-4,4)")
    h1_pull = tree1.GetHistogram()
    tree2.Draw("alpha/alpha_e>>h2_pull(100,-4,4)")
    h2_pull = tree2.GetHistogram()

    h1_pull.SetStats(0)
    h2_pull.SetStats(0)
    h1_pull.SetTitle(title)
    h2_pull.SetTitle(title)
    h1_pull.GetXaxis().SetTitle('#phi_{CP}/#sigma_{#phi_{CP}}')

    h1_pull.SetLineColor(ROOT.kRed)
    h2_pull.SetLineColor(ROOT.kBlue)
    h1_pull.SetMaximum(1.2*max(h1_pull.GetMaximum(),h2_pull.GetMaximum()))
    h1_pull.Draw("HIST")
    h2_pull.Draw("HIST SAME")
    mean1 = h1_pull.GetMean()
    mean2 = h2_pull.GetMean()
    rms1 = h1_pull.GetRMS()
    rms2 = h2_pull.GetRMS()
    leg = ROOT.TLegend(0.6,0.7,0.9,0.9)
    leg.AddEntry(h1_pull,f"fit=toys mean={mean1:.1f}, rms={rms1:.1f}","l")
    leg.AddEntry(h2_pull,f"fit#neq toys mean={mean2:.1f}, rms={rms2:.1f}","l")
    leg.Draw()
    c1.SaveAs(plot_name)


def ConvertTree(intree):

    outtree = ROOT.TTree("new"+intree.GetName(), "tree with alpha uncertainties")

    outtree.alpha_up   = array('f', [0.])
    outtree.alpha_down = array('f', [0.])
    outtree.alpha_e    = array('f', [0.])
    outtree.alpha_nom  = array('f', [0.])

    outtree.Branch("alpha_up",   outtree.alpha_up,   "alpha_up/F")
    outtree.Branch("alpha_down", outtree.alpha_down, "alpha_down/F")
    outtree.Branch("alpha_e",    outtree.alpha_e,    "alpha_e/F")
    outtree.Branch("alpha",      outtree.alpha_nom,  "alpha/F")


    # Helper dict to group events by (iToy, iSeed)
    events = {}

    for entry in intree:
        key = (entry.iToy, entry.iSeed)
        if key not in events:
            events[key] = {}
        if abs(entry.quantileExpected + 1) < 1e-6:  # nominal
            events[key]["nominal"] = entry.alpha
        elif abs(entry.quantileExpected - 0.32) < 1e-6:
            events[key]["up"] = entry.alpha
        elif abs(entry.quantileExpected + 0.32) < 1e-6:
            events[key]["down"] = entry.alpha

    # Loop again to fill tree
    for key, vals in events.items():
        if "nominal" in vals and "up" in vals and "down" in vals:

            outtree.alpha_up[0] = vals["up"] - vals["nominal"]
            outtree.alpha_down[0] = vals["down"] - vals["nominal"]
            outtree.alpha_e[0] = abs(vals["up"] - vals["down"])/2
            outtree.alpha_nom[0] = vals["nominal"]
            outtree.Fill()

    return outtree

# test fit mode with different toys

for mode1 in [0,1,2,3,4]: # loop over fit models
    file1 = ROOT.TFile(f"outputs_mergemode{mode1}/cmb/higgsCombine.alpha0.mergemode{mode1}Toys.mergemode{mode1}Fit.MultiDimFit.mH125.root")
    tree1 = file1.Get("limit")
    for mode2 in [0,1,2,3,4]: # loop over toy models
        if mode1 == mode2: continue
        print(f"fit mode={mode1}, toys mode={mode2}")
        file2 = ROOT.TFile(f"outputs_mergemode{mode1}/cmb/higgsCombine.alpha0.mergemode{mode2}Toys.mergemode{mode1}Fit.MultiDimFit.mH125.root")
        tree2 = file2.Get("limit")
        title = f"fit mode={mode1}, toys mode={mode2}"
        plot_name = f"bias_test_plots/bias_test_fitmode{mode1}_toymode{mode2}.pdf"
        plot_toys(tree1,tree2,title,plot_name)

        newtree1 = ConvertTree(tree1)
        newtree2 = ConvertTree(tree2)
        newtree1.SetName(f"tree1_fitmode{mode1}_toymode{mode2}")
        newtree2.SetName(f"tree2_fitmode{mode1}_toymode{mode2}")
        #fout = ROOT.TFile(f"bias_test_plots/converted_trees_fitmode{mode1}_toymode{mode2}.root","RECREATE")
        #newtree1.Write('tree1')
        #newtree2.Write('tree2')
        plot_name = f"bias_test_plots/bias_test_pull_fitmode{mode1}_toymode{mode2}.pdf"
        plot_pull(newtree1,newtree2,title,plot_name)
