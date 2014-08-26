#!/usr/bin/env python
# encoding: utf-8

from mpl_tools import *
from r2numpy import *

def drawGraph( g, *args, **kwargs ):
    """Plot TGraph"""
    from pylab import plot
    x, y = g.getBuffers()
    return plot( x, y, *args, **kwargs )
##end def 
from ROOT import TGraph
setattr( TGraph, 'draw', drawGraph )

def drawGraphErr( g, *args, **kwargs ):
    """Plot TGraphErrors"""
    from pylab import errorbar
    x, y = g.getBuffers()
    ex, ey = g.getErrBuffers()
    from numpy import all
    ex = ex if all( ex!=0.0 ) else None
    ey = ey if all( ey!=0.0 ) else None
    return errorbar( x, y, ey, ex, *args, **kwargs )
##end def 
from ROOT import TGraphErrors
setattr( TGraphErrors, 'draw', drawGraphErr )

def drawGraphAsymmErr( g, *args, **kwargs ):
    """Plot TGraphErrors"""
    from pylab import errorbar
    x, y = g.getBuffers()
    exl, exh, eyl, eyh = g.getErrBuffers()
    from numpy import array
    ex = array( ( exl, exh ) )
    ey = array( ( eyl, eyh ) )
    from numpy import all
    ex = ex if all( ex!=0.0 ) else None
    ey = ey if all( ey!=0.0 ) else None
    return errorbar( x, y, ey, ex, *args, **kwargs )
##end def 
from ROOT import TGraphAsymmErrors
setattr( TGraphAsymmErrors, 'draw', drawGraphAsymmErr )

def drawHist1Dbar( h, *args, **kwargs ):
    """Plot 1-dimensinal histogram using pyplot.bar"""
    height = h.getBuffer()
    ax = h.GetXaxis()
    lims, fixed = ax.getBinEdges( type=True )
    width=None
    left  = lims[:-1]
    if fixed: width = ax.GetBinWidth( 1 )
    else: width = lims[1:] - left        

    from pylab import bar
    return bar( left, height, width, *args, **kwargs )
##end drawHist1D
from ROOT import TH1
setattr( TH1, 'drawBar', drawHist1Dbar )

def drawHist1Dline( h, *args, **kwargs ):
    """Plot 1-dimensinal histogram using pyplot.plot"""
    lims=h.GetXaxis().getBinEdges()
    height=h.getBuffer()
    return plot_hist( lims, height, *args, **kwargs )
##end drawHist1D
from ROOT import TH1
setattr( TH1, 'drawLine', drawHist1Dline )

def drawHist1D( h, *args, **kwargs ):
    """Plot 1-dimensinal histogram using pyplot.plot and pyplot.bar 
       baroptions are passed as baropts=dict(opt=value)
    """
    ax = h.GetXaxis()
    height = h.getBuffer()
    lims, fixed = ax.getBinEdges( type=True )

    from mpl_tools import plot_histbar
    return plot_histbar( lims, height, *args, fixedwidth=fixed, **kwargs )
##end drawHist1D
from ROOT import TH1
setattr( TH1, 'draw', drawHist1D )

def drawHist2Dmesh( h, *args, **kwargs ):
    """Plot TH2 using matplotlib.pcolormesh"""
    mask, colz = [ kwargs.pop(x) if x in kwargs else None for x in ['mask', 'colz'] ]
    baropts = baropts or {}

    # get bin edges first
    x1 = h.GetXaxis().getBinEdges()
    y1 = h.GetYaxis().getBinEdges()

    # make a 2D mesh
    from numpy import meshgrid
    x, y = meshgrid( x1, y1 )
    # print 'mesh x', x
    # print 'mesh y', y

    # get data bufer w/o underflow/overflow bins
    buf = h.getBuffer( mask=mask )

    # plot
    from pylab import pcolormesh
    res = pcolormesh( x, y, buf, *args, **kwargs )
    if colz:
        from pylab import colorbar
        cbar = colorbar()
        # return res, cbar 
    ##end if colz
    return res
##end def drawHist2D
from ROOT import TH2
setattr( TH2, 'drawMesh', drawHist2Dmesh )

def drawHist2D( h, *args, **kwargs ):
    """Plot TH2 using matplotlib.pcolorfast"""
    mask, colz = [ kwargs.pop(x) if x in kwargs else None for x in ['mask', 'colz'] ]

    # get bin edges first
    xax = h.GetXaxis()
    yax = h.GetYaxis()
    if xax.GetXbins().GetSize()>0 or yax.GetXbins().GetSize()>0:
        print 'Can not draw 2D a histogram with variable bin widths'
        print 'Use drawMesh method or draweHist2Dmesh function instead'
        return
    ##end if
    x = [ xax.GetXmin(), xax.GetXmax() ]
    y = [ yax.GetXmin(), yax.GetXmax() ]

    # get data bufer w/o underflow/overflow bins
    buf = h.getBuffer( mask=mask )

    # plot
    from pylab import axes
    ax = axes()
    res = ax.pcolorfast( x, y, buf, *args, **kwargs )
    if colz:
        from pylab import gcf
        cbar = gcf().colorbar( res )
        # return res, cbar
    ##end if colz

    return res
##end def drawHist2D
from ROOT import TH2
setattr( TH2, 'draw', drawHist2D )

def drawTMatrix( m, *args, **kwargs ):
    """Plot TMatrixD using matplotlib.pcolorfast"""
    mask, colz, limits = [ kwargs.pop(x) if x in kwargs else None for x in ['mask', 'colz', 'limits'] ]
    x, y = limits!=None and limits or ([ 0.0, m.GetNcols() ], [ 0.0, m.GetNrows() ])

    buf = m.getBuffer()
    if mask:
        from numpy import ma
        buf = ma.array( buf, mask = buf==mask )
    ##end if

    # plot
    from pylab import axes
    ax = axes()
    res = ax.pcolorfast( x, y, buf, *args, **kwargs )
    if colz:
        from pylab import gcf
        cbar = gcf().colorbar( res )
        # return res, cbar
    ##end if colz

    return res
##end def drawHist2D
from ROOT import TMatrixD, TMatrixF
setattr( TMatrixD, 'draw', drawTMatrix )
setattr( TMatrixF, 'draw', drawTMatrix )

def drawTMatrixDiag( m, *args, **kwargs ):
    """Plot TH2 diagoanl ad TH1"""
    assert m.GetNcols()==m.GetNrows(), 'Matrix is not square'
    limits, = [ kwargs.pop(x) if x in kwargs else None for x in ['limits'] ]

    buf = m.getBuffer()
    from numpy import diagonal
    bins = diagonal( buf )
    lims = None
    if limits:
        from numpy import linspace
        lims = linspace( limits[0], limits[1], m.GetNcols()+1 )
    else: 
        from numpy import arange
        lims = arange( 0.0, m.GetNcols()+1, 1.0 )
    ##end if

    from mpl_tools import plot_histbar
    return plot_histbar( lims, bins, *args, **kwargs )
##end def drawHist2D
setattr( TMatrixD, 'drawDiag', drawTMatrixDiag )
setattr( TMatrixF, 'drawDiag', drawTMatrixDiag )

def drawHist2Ddiag( h, *args, **kwargs ):
    """Plot TH2 diagoanl ad TH1"""
    assert h.GetNbinsX()==h.GetNbinsY(), 'Histogram is not square'
    buf = h.getBuffer()
    from numpy import diagonal
    bins = diagonal( buf )
    lims, fixedwidth = h.GetXaxis().getBinEdges( type=True )

    from mpl_tools import plot_histbar
    return plot_histbar( lims, bins, *args, fixedwidth=fixedwidth, **kwargs )
##end def drawHist2D
from ROOT import TH2
setattr( TH2, 'drawDiag', drawHist2Ddiag )

def drawTF1( f, x=None, *args, **kwargs ):
    """Plot TF1
       if x is an array-like it's used as array of x
       if x is integer it is used as N of points in function range
       if x is float it is used as step
       if x is None, the TF1->Npx is used as number of points
    """
    import numpy
    tp = type(x)
    if x==None:
        x = numpy.linspace( f.GetXmin(), f.GetXmax(), f.GetNpx() )
    elif tp==int:
        x = numpy.linspace( f.GetXmin(), f.GetXmax(), x )
    elif tp==float:
        x = numpy.arange( f.GetXmin(), f.GetXmax(), x )
    ##end

    return drawFun( f.Eval, x, *args, **kwargs )
##end def drawTF1
from ROOT import TF1
setattr( TF1, 'draw', drawTF1 )

def importTitle( o, what='' ):
    """Import object titles to canvas/axis labels"""
    if what=='': return
    from pylab import axes
    ax = axes()
    if 'x' in what: ax.set_xlabel( o.GetXaxis().GetTitle() )
    if 'y' in what: ax.set_ylabel( o.GetYaxis().GetTitle() )
    if 't' in what: ax.set_title( o.GetTitle() )
    if 'n' in what: ax.set_title( o.GetName() )
    if 'f' in what: ax.set_title( o.GetExpFormula() )
##end def importTitle
from ROOT import TF1, TH1, TGraph
setattr( TF1,    'showTitle', importTitle )
setattr( TH1,    'showTitle', importTitle )
setattr( TGraph, 'showTitle', importTitle )

