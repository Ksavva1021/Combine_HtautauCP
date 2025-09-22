import ROOT
import math

ROOT.gErrorIgnoreLevel = ROOT.kFatal

def compare_histograms(h1: ROOT.TH1, h2: ROOT.TH1, metric: str = "chi2") -> float:
    """
    Compare two histograms using a chosen metric.

    Parameters:
        h1 (ROOT.TH1): First histogram
        h2 (ROOT.TH1): Second histogram
        metric (str): Similarity metric. Options:
            - "chi2"   : Chi-squared distance (ROOT built-in)
            - "kolmogorov" : Kolmogorov-Smirnov test (ROOT built-in)
            - "integral"   : Relative difference in total integrals
            - "bhattacharyya" : Bhattacharyya coefficient (manual implementation)

    Returns:
        float: similarity/distance measure
    """

    if metric == "chi2":
        # ROOT Chi2Test (option "CHI2" returns reduced chi2)
        return h1.Chi2Test(h2, "WW")

    elif metric == "kolmogorov":
        # ROOT Kolmogorov-Smirnov test
        return h1.KolmogorovTest(h2,"WW")

    elif metric == "integral":
        # Compare total number of entries (normalized)
        i1, i2 = h1.Integral(), h2.Integral()
        if i1 == 0 and i2 == 0:
            return 0.0
        return abs(i1 - i2) / max(i1, i2)

    elif metric == "bhattacharyya":
        # Normalize histograms
        i1, i2 = h1.Integral(), h2.Integral()
        if i1 == 0 or i2 == 0:
            return float("inf")

        bc = 0.0
        for bin in range(1, h1.GetNbinsX() + 1):
            p = max(h1.GetBinContent(bin) / i1,0.)
            q = max(h2.GetBinContent(bin) / i2,0.)
            bc += math.sqrt(p * q)

        return bc # Returns the Bhattacharyya coefficient
        # Bhattacharyya distance
        #return -math.log(bc) if bc > 0 else float("inf")

    else:
        raise ValueError(f"Unknown metric: {metric}")


if __name__ == "__main__":

    for i in range(1,10):

        bin_name = f'htt_tt_{i}_13p6TeV'
        print(f"Comparing histograms in bin: {bin_name}")

        f1 = ROOT.TFile(f'merge_datacards/outputs_mergemode_0/cmb/common/{bin_name}_input.root')
        f2 = ROOT.TFile(f'merge_datacards/outputs_mergemode_2/cmb/common/{bin_name}_input.root')

        # get directory with bin_name
        d1 = f1.Get(bin_name)
        d2 = f2.Get(bin_name)

        # loop over all histograms in d1 and find matching ones in d2
        for key in d1.GetListOfKeys():
            obj1 = key.ReadObj()
            if not isinstance(obj1, ROOT.TH1):
                continue

            name = obj1.GetName()
            if '_bbb_' in name: # skip bbb histograms
                continue
            obj2 = d2.Get(name)
            if not obj2:
                #print(f"Histogram {name} not found in second file")
                continue

            chi2 = compare_histograms(obj1, obj2, metric="chi2")
            ks = compare_histograms(obj1, obj2, metric="kolmogorov")
            integral_diff = compare_histograms(obj1, obj2, metric="integral")
            bh_coeff = compare_histograms(obj1, obj2, metric="bhattacharyya")

            # we only print if the histograms differ significantly
            if integral_diff > 0.001 or bh_coeff < 0.95:  
                print(f"Histogram: {name}")
                print(f"  Chi2: {chi2:.4f}")
                print(f"  Kolmogorov: {ks:.4f}")
                print(f"  Integral difference: {integral_diff:.7f}")
                print(f"  Bhattacharyya coefficient: {bh_coeff:.7f}")
                print()

        print("--------------------------------------------------\n")