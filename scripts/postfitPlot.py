import CombineHarvester.CombineTools.plotting as plot
import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)

# these can eventually be set by command line arguments:
lumi = '61.9 fb^{-1} (13.6 TeV)'
bin_name = 'htt_tt_1_13p6TeV'
channel = 'tt'
if channel == '':  channel=binname[4:6]
input_file = ROOT.TFile('shapes_output.root')
input_file_ps = 'shapes_output_ps.root'
if input_file_ps:
    input_file_ps = ROOT.TFile(input_file_ps)

extra_pad = 0.15
norm_bins = False
autoblind = 0.02
ratio_range = '0.5,1.5'
signal_scale = 100.0

if norm_bins:
    y_title = 'Events / bin width'
else: 
    y_title = 'Events'
x_title = 'BDT output'

def getHistogram(fname, histname, dirname='', allowEmpty=False):
    outname = fname.GetName()
    for key in fname.GetListOfKeys():
        histo = fname.Get(key.GetName())
        dircheck = False
        if dirname == '' : dircheck=True
        elif dirname in key.GetName(): dircheck=True
        if isinstance(histo,ROOT.TH1F) and key.GetName()==histname:
            return [histo,outname]
        elif isinstance(histo,ROOT.TDirectory) and dircheck:
            return getHistogram(histo,histname, allowEmpty=allowEmpty)
    print('Failed to find histogram with name %(histname)s in file %(fname)s '%vars())
    if allowEmpty:
        return [ROOT.TH1F('empty', '', 1, 0, 1), outname]
    else:
        return None

def backgroundComp(leg,plots,colour):
    return dict([('leg_text',leg),('plot_list',plots),('colour',colour)])

def createAxisHists(n,src,xmin=0,xmax=499):
    result = []
    for i in range(0,n):
        res = src.Clone()
        res.Reset()
        res.SetTitle("")
        res.SetName("axis%(i)d"%vars())
        res.SetAxisRange(xmin,xmax)
        res.SetStats(0)
        result.append(res)
    return result

def PositionedLegendUnrolled(width, height, pos, offset):
    o = offset
    w = width
    h = height
    l = ROOT.gPad.GetLeftMargin()
    t = ROOT.gPad.GetTopMargin()
    b = ROOT.gPad.GetBottomMargin()
    r = ROOT.gPad.GetRightMargin()
    if pos == 1:
        return ROOT.TLegend(l + o, 1 - t - o - h, l + o + w, 1 - t - o, '', 'NBNDC')
    if pos == 2:
        c = l + 0.5 * (1 - l - r)
        return ROOT.TLegend(c - 0.5 * w, 1 - t - o - h, c + 0.5 * w, 1 - t - o, '', 'NBNDC')
    if pos == 3:
        return ROOT.TLegend(1 - r - o - w, 1 - t - o - h, 1 - r - o, 1 - t - o, '', 'NBNDC')
    if pos == 4:
        return ROOT.TLegend(l + o, b + o, l + o + w, b + o + h, '', 'NBNDC')
    if pos == 5:
        c = l + 0.5 * (1 - l - r)
        return ROOT.TLegend(c - 0.5 * w, b + o, c + 0.5 * w, b + o + h, '', 'NBNDC')
    if pos == 6:
        return ROOT.TLegend(1 - r - o - w, b + o, 1 - r - o, b + o + h, '', 'NBNDC')
    if pos == 7:
        return ROOT.TLegend(1 - o - w, 1 - t - o - h, 1 - o, 1 - t - o, '', 'NBNDC')

def DrawTitleUnrolled(pad, text, align, scale=1):
    pad_backup = ROOT.gPad
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()

    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if pad_ratio < 1.:
        pad_ratio = 1.

    textSize = 0.6
    textOffset = 0.2

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(textSize * t * pad_ratio * scale)

    y_off = 1 - t + textOffset * t + 0.01
    if align == 1:
        latex.SetTextAlign(11)
        latex.DrawLatex(l, y_off, text)
    if align == 2:
        latex.SetTextAlign(21)
        latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text)
    if align == 3:
        latex.SetTextAlign(31)
        latex.DrawLatex(1 - r, y_off, text)
    pad_backup.cd()

is_1d_bin = bin_name in ['htt_mt_1_13p6TeV', 'htt_mt_2_13p6TeV', 'htt_tt_1_13p6TeV', 'htt_tt_2_13p6TeV']
if is_1d_bin:
    plot.ModTDRStyle(width=900, height=800, r=0.3, l=0.15)
#else:
    #plot.ModTDRStyle(width=1800, height=700, r=0.4, l=0.16, t=0.12,b=0.15)

canv = ROOT.TCanvas()    
pads=plot.TwoPadSplit(0.35,0.01,0.01)

bkghist = getHistogram(input_file,'TotalProcs',bin_name)[0]
bkgonlyhist = getHistogram(input_file,'TotalBkg',bin_name)[0]

sighist = getHistogram(input_file,'TotalSig',bin_name)[0]

# if pseudoscalar file is provided then prdoduce a signal histogram from it and average it with the scalar one for the autoblinding
# if it isn't provided then just use the scalar one
if input_file_ps:
    sighist_ps = getHistogram(input_file_ps,'TotalSig',bin_name)[0]
    ave_sig = sighist.Clone()
    ave_sig.Add(sighist_ps)
    ave_sig.Scale(0.5)
else :
    sighist_ps = None
    ave_sig = sighist.Clone()

datahist = getHistogram(input_file,"data_obs",bin_name)[0]
datahist.SetMarkerStyle(20)

if autoblind:
    # get rid of any data bins when the ave_sig/bkghist fraction is > autoblind by settinfg the yield to -0.1 and the error to 0
    for i in range(1, datahist.GetNbinsX()+1):
        if bkghist.GetBinContent(i) == 0: blind = True
        elif ave_sig.GetBinContent(i) / bkghist.GetBinContent(i) > autoblind: blind = True
        else: blind = False
        if blind:
            datahist.SetBinContent(i, -0.1)
            datahist.SetBinError(i, 0)

bkghist.SetMarkerSize(0)
bkghist.SetFillColor(2001)
bkghist.SetLineColor(0)
bkghist.SetLineWidth(1)
bkghist.SetFillColor(plot.CreateTransparentColor(12,0.4))
bkgonlyhist.SetMarkerSize(0)
bkgonlyhist.SetFillColor(2001)
bkgonlyhist.SetLineColor(0)
bkgonlyhist.SetLineWidth(1)
bkgonlyhist.SetFillColor(plot.CreateTransparentColor(12,0.4))
if norm_bins: 
    bkghist.Scale(1.0,"width")
    bkgonlyhist.Scale(1.0,"width")
    sighist.Scale(1.0,"width")
    datahist.Scale(1.0,"width")

background_schemes = {
    'mt':[
            backgroundComp("H#rightarrow#tau#tau",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
            backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)),
            backgroundComp("Jet#rightarrow#tau_{h} fakes",["JetFakes"],ROOT.TColor.GetColor(192,232,100)),
            backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(248,206,104)),
            ],
    'tt':[
            backgroundComp("H#rightarrow#tau#tau",["TotalSig"],ROOT.TColor.GetColor(51,51,230)),
            backgroundComp("Others",["ZL","VVT","TTT"],ROOT.TColor.GetColor(100,192,232)),
            backgroundComp("Jet#rightarrow#tau_{h} fakes",['JetFakes','JetFakesSublead'],ROOT.TColor.GetColor(192,232,100)),
            backgroundComp("Z#rightarrow#tau#tau",["ZTT"],ROOT.TColor.GetColor(248,206,104)),
            ],
}
background_schemes['et'] = background_schemes['mt']

stack_histos = []

for t in background_schemes[channel]:
    plots = t['plot_list']
    isHist = False
    h = ROOT.TH1F()
    for k in plots:
        if h.GetEntries()==0 and getHistogram(input_file,k, bin_name,False) is not None:
            isHist = True
            h = getHistogram(input_file,k, bin_name)[0]
            h.SetName(k)
        else:
            if getHistogram(input_file,k, bin_name,False) is not None:
                isHist = True
                h.Add(getHistogram(input_file,k, bin_name)[0])

    h.SetFillColor(t['colour'])
    h.SetLineColor(ROOT.kBlack)
    h.SetMarkerSize(0)
    
    if norm_bins: h.Scale(1.0,"width")
    if isHist:
        stack_histos.append(h)

stack = ROOT.THStack("hs","")
for hists in stack_histos:
    stack.Add(hists)

signal_stack = ROOT.THStack("hs_sig","")
if sighist.GetEntries() > 0:
    sighist.SetFillStyle(0)
    sighist.SetLineColor(ROOT.kRed)
    sighist.SetLineWidth(2)
    if norm_bins: sighist.Scale(1.0,"width")
    sighist.Scale(signal_scale)
    signal_stack.Add(sighist)
if sighist_ps and sighist_ps.GetEntries() > 0:
    sighist_ps.SetFillStyle(0)
    sighist_ps.SetLineColor(ROOT.kBlue)
    sighist_ps.SetLineWidth(2)
    if norm_bins: sighist_ps.Scale(1.0,"width")
    sighist_ps.Scale(signal_scale)
    signal_stack.Add(sighist_ps)

axish = createAxisHists(2,bkghist,bkghist.GetXaxis().GetXmin(),bkghist.GetXaxis().GetXmax()-0.01)
axish[0].SetMaximum((extra_pad+1)*bkghist.GetMaximum())
axish[0].GetYaxis().SetTitle(y_title)
axish[0].GetXaxis().SetLabelSize(0)
axish[0].GetXaxis().SetTitleSize(0)
axish[1].GetXaxis().SetTitle(x_title)


pads[0].cd()
pads[0].SetTicks(1)
axish[0].Draw("hist")

stack.Draw("histsame")
bkghist.Draw("e2same")
signal_stack.Draw("histsamenostack")
datahist.Draw("P Z 0 same")
pads[0].RedrawAxis()

# now we produce the ratio plot
# For now the ratio just shows the ratio of the observed data to the signal+background expectation
# this can be changed in teh future to show the ratio to background only and signal+background etc seperatly 

pads[1].cd()
pads[1].SetTicks(1)
pads[1].SetGrid(0,1)
axish[1].Draw("axis")
axish[1].SetMinimum(float(ratio_range.split(',')[0]))
axish[1].SetMaximum(float(ratio_range.split(',')[1]))
axish[1].GetYaxis().SetTitle("Obs./Exp.")
# draw dashed line at 1
line = ROOT.TLine(axish[1].GetXaxis().GetXmin(), 1, axish[1].GetXaxis().GetXmax(), 1)
line.SetLineStyle(2)
line.SetLineColor(ROOT.kBlack)
line.Draw("same")

ratio_bkghist = plot.MakeRatioHist(bkghist,bkghist,True,False)
ratio_bkghist.Draw("e2same")
ratio_datahist = plot.MakeRatioHist(datahist,bkghist,True,False)
ratio_datahist.Draw("P Z 0 same")
#sbhist = bkghist.Clone()
#sbhist.SetLineColor(ROOT.kRed)
#sbhist.SetLineWidth(2)
#sbhist.SetFillStyle(0)
#ratio_sighist = plot.MakeRatioHist(sbhist,bkghist,True,False)
#ratio_sighist.Draw("histsame")
pads[1].RedrawAxis()

# draw the title, labels, and legend

#CMS and lumi labels
extra='Preliminary'
if is_1d_bin:
    plot.DrawCMSLogo(pads[0], 'CMS', extra, 0, 0.07, -0.0, 2.0, '', 0.85)
    plot.DrawTitle(pads[0], lumi, 3, textSize=0.6)
else:
    plot.DrawCMSLogo(pads[0], 'CMS', extra, 0, 0.07, -0.0, 2.0, '', 0.6)
    DrawTitleUnrolled(pads[0], lumi, 3, scale=0.7)

legend = PositionedLegendUnrolled(0.25,0.5,7,0.)
legend.SetTextSize(0.03) 
legend.SetTextFont(42)
legend.SetFillStyle(0)

legend.AddEntry(datahist,"Observed","PEl")
# loop backwards over the stack_histos by index to add them in the right order
for legi in range(len(stack_histos)-1,-1,-1):
    hists = stack_histos[legi]
    legend.AddEntry(hists,background_schemes[channel][legi]['leg_text'],"f")

legend.AddEntry(bkghist,"Uncertainty","f")
legend.AddEntry(sighist,"H#rightarrow#tau#tau (#alpha^{H#tau#tau}=0^{#circ})","l")
if sighist_ps:
    legend.AddEntry(sighist_ps,"H#rightarrow#tau#tau (#alpha^{H#tau#tau}=90^{#circ})","l")
legend.Draw("same")

canv.Print('test.pdf')
