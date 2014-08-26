#!/usr/bin/env python
# encoding: utf-8

def savefig( name, *args, **kwargs ):
    """Save fig and print output filename"""
    if not name: return
    from pylab import savefig
    savefig( name, *args, **kwargs )
    print 'Save figure', name
##end def savefig

def set_title( t ):
    """Set window title"""
    from pylab import canvas
    canvas.set_window_title( t )
##end if

def plot_hist( lims, height, *args, **kwargs ):
    """Plot histogram with lines. Like bar(), but without lines between bars."""
    from numpy import ravel, empty, vstack
    y = empty( len(height)*2+2 )
    y[1:-1] = vstack( ( height, height ) ).ravel( order='F' ) 
    y[0]=0.0
    y[-1]=0.0
    x = vstack( ( lims, lims ) ).ravel( order='F' )

    from pylab import plot
    return plot( x, y, *args, **kwargs )
##end def plot_hist

def plot_histbar( lims, height, *args, **kwargs ):
    fixedwidth, baropts = [ kwargs.pop(x) if x in kwargs else None\
                            for x in ['fixedwidth', 'baropts'] ]
    left  = lims[:-1]
    width = None
    if fixedwidth: width = lims[1]-lims[0]
    else: width = lims[1:] - left        

    baropts = baropts or {}
    if not 'linewidth' in baropts: baropts['linewidth']=0

    from pylab import bar
    return plot_hist( lims, height, *args, **kwargs ), bar( left, height, width, **baropts )
##end def plot_histbar

def drawFun( f, x, *args, **kwargs ):
    """Draw a function values for x"""
    from numpy import frompyfunc
    fun = frompyfunc( f, 1, 1 )

    from pylab import plot
    return plot( x, fun(x), *args, **kwargs )
##end def drawFun

