def green(string,**kwargs):
    '''Displays text in green text inside a black background'''
    return kwargs.get('pre',"")+"\x1b[0;32;40m%s\033[0m"%string

def NegativeBins(p):
    '''Replaces negative bins in hists with 0'''
    hist = p.shape()
    has_negative = False
    for i in range(1,hist.GetNbinsX()+1):
      if hist.GetBinContent(i) < 0:
         has_negative = True
         print("(Process, Channel, Bin) = (%s, %s, %s) has negative bins." % (p.process(), p.channel(), p.bin()))
    if (has_negative):
      for i in range(1,hist.GetNbinsX()+1):
         if hist.GetBinContent(i) < 0:
            hist.SetBinContent(i,0)
    p.set_shape(hist,False)


def NegativeYields(p):
    '''If process has negative yield then set to 0'''
    if p.rate()<0:
        print("(Process, Channel, Bin) = (%s, %s, %s) has a negative yield" % (p.process(), p.channel(), p.bin()))
        p.set_rate(0.)
