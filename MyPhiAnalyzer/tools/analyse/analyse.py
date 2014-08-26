from __future__ import division
import glob
import ROOT
from ROOT import TFile, TH1, TGraph, TGraphErrors, TGraphAsymmErrors, TDirectory
import scipy
import numpy
from numpy import array
from numpy import where
from r2numpy import *
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.legend as pyleg
import matplotlib.gridspec as gridspec
#from matplotlib import rc
#rc('text', usetex=True)

filelist = []
for i,f in enumerate(["Histos_minbias_old.root","Histos_minbias_new.root","Histos_minbias_new_overlapfix.root"]):
    filelist.append( [ROOT.TFile(f), ["k","b","r","g"][i], ["old", "new", "overlap-fixed"][i]] )

print "File list: ", filelist

histolist = []
histolist.append(["phi_det","linear", 'sector', '$1/N <E_{\\rm det}>$ / GeV', "Detector $\phi$-profile"])
histolist.append(["phi_gen","linear", '$\phi$', '$1/N <E_{\\rm gen}>$ / GeV', "Generated $\phi$-profile"])
histolist.append(["z","log", 'module', '$1/N <E>$ / GeV', "$z$-profile / GeV"])
histolist.append(["z_near","log", 'module', '$1/N <E>$ / GeV', "$z_{\\rm near}$-profile"])
histolist.append(["z_far","log", 'module', '$1/N <E>$ / GeV', "$z_{\\rm far}$-profile"])
histolist.append(["Efac_over_eta_gen", "linear", "$\eta_{\\rm GEN}$", "$E_{\\rm tot,gen}/E_{\\rm tot,det}$", "correction factor"])
histolist.append(["E_over_eta_gen", "linear", "$\eta_{\\rm GEN}$", "$E_{\\rm tot,det}$", "$\eta$ boundaries of CASTOR"])
#histolist.append(["hEnergy_res","linear", '?', '$\sigma E/E$', "Energy resolution"])

histolist.append(["etaBinned/z_eta_n5.3","log", 'module', '$E$ / arb.u.', "$z$-profile @ $\\eta=5.3$"])
histolist.append(["etaBinned/z_eta_n5.7","log", 'module', '$E$ / arb.u.', "$z$-profile @ $\\eta=5.7$"])
histolist.append(["etaBinned/z_eta_n6.1","log", 'module', '$E$ / arb.u.', "$z$-profile @ $\\eta=6.1$"])
histolist.append(["etaBinned/z_eta_n6.5","log", 'module', '$E$ / arb.u.', "$z$-profile @ $\\eta=6.5$"])

histolist.append(["etaBinned/phi_det_eta_n5.3","linear", 'sector', '$E$ / arb.u. ', "$\\phi$-profile @ $\\eta=5.3$"])
histolist.append(["etaBinned/phi_det_eta_n5.7","linear", 'sector', '$E$ / arb.u. ', "$\\phi$-profile @ $\\eta=5.7$"])
histolist.append(["etaBinned/phi_det_eta_n6.1","linear", 'sector', '$E$ / arb.u. ', "$\\phi$-profile @ $\\eta=6.1$"])
histolist.append(["etaBinned/phi_det_eta_n6.5","linear", 'sector', '$E$ / arb.u. ', "$\\phi$-profile @ $\\eta=6.5$"])


histolist = [["demo/" + h[0],h[1],h[2],h[3],h[4]] for h in histolist]

for h in histolist:
    print "Creating histogram: {histo} ...".format(histo=h[4] + ' (' + h[0] + ')')
    fig = plt.figure(figsize=[8,8])
    ax = fig.gca()
    ax.axis('off')
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) #used to change size of subplots
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1],sharex=ax1)
    ax1.label_outer()
    ax2.autoscale(tight=True)
    plt.tight_layout(1.1)
    ax.xaxis.set_visible(False) #don't show ticks or axes for main plot. only subplots
    ax.yaxis.set_visible(False)
    ax.set_title(h[4],fontsize=20)
    ax2.set_xlabel(h[2],fontsize=20)
    ax1.set_ylabel(h[3],fontsize=20)
    
    f=filelist[0] #first file
    hist = f[0].Get(h[0])
    if type(hist) != ROOT.TH1D:
        print 'skipping', h[0], "of type", type(hist)
        continue
    xold=hist.GetXaxis().getBinEdges()
    x=[(a + b) / 2 for a, b in zip(xold[::1], xold[1::1])]
    yref=hist.getBuffer()
    yeref=hist.getErrBuffer()
    #print type(yref), yref
    #print type(yeref), yeref
    if type(yeref) != type(None):
        yeref = numpy.sqrt(yeref) 
    else:
        print "No errors found for histogram. Setting to 1."
        yeref = numpy.ones_like(yref)
    
    #plotting
    ax1.plot(x,yref,label=str(f[2]),color=f[1])
    positive = yref - yeref > 0 
    ax1.fill_between(x, yref-yeref, yref+yeref,where=positive,facecolor=f[1],alpha=0.2)
    
    for f in filelist[1:]:
        #obtaining values
        hist = f[0].Get(h[0])
        y=hist.getBuffer()
        ye=hist.getErrBuffer()
        if type(ye) != type(None):
            ye = numpy.sqrt(ye)
        else:
            ye = numpy.ones_like(yref)
                    
        #print zip(x,y)
       
        
        #plotting
        ax1.plot(x,y,label=str(f[2]),color=f[1])
        positive = y - ye > 0
        ax1.fill_between(x, y-ye, y+ye,where=positive,facecolor=f[1],alpha=0.2)
        print "len", (y>0).all()
        if (y>0).all():
            ax1.set_yscale(h[1])
        
        #ratio plot
        ydiff  = array([y[i]/yref[i] if yref[i] != 0 else 0 for i in range(len(y))])
        ydiffe = array([sqrt(ye[i]**2/yref[i]**2 + yeref[i]**2*y[i]**2/yref[i]**4) if yref[i] != 0 else 0 for i in range(len(y))])
        ax2.plot(x,ydiff,label=str(f),color=f[1])
        positive = ydiff - ydiffe > 0
        ax2.fill_between(x, ydiff-ydiffe, ydiff+ydiffe,where=positive,facecolor=f[1],alpha=0.2)
        ax2.axhline(y=1,color=filelist[0][1])
        ax2.set_ylim(max(1/3,numpy.nanmin(ydiff - ydiffe)),min(3,numpy.nanmax(ydiff+ydiffe)))
        #ax2.set_yscale(h[1])
        
    leg = ax1.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    
    print "saving ", str(h[0]), '\n'
    fig.savefig(str(h[0])+".png", layout='tight')
    fig.savefig(str(h[0])+".pdf", layout='tight')

