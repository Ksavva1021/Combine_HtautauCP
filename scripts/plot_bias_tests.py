import ROOT
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
