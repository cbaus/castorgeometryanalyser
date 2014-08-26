#!/usr/bin/env python
# encoding: utf-8

def tGraphGetBuffers( g ):
    """Get TGraph x and y buffers"""
    npoints = g.GetN()
    if npoints==0: return None, None
    from numpy import frombuffer, double
    return frombuffer( g.GetX(), double, npoints)\
         , frombuffer( g.GetY(), double, npoints)
##end def function
from ROOT import TGraph
setattr( TGraph, 'getBuffers', tGraphGetBuffers )

def tGraphGetErrBuffers( g ):
    """Get TGraph x and y buffers"""
    npoints = g.GetN()
    if npoints==0: return None, None
    from numpy import frombuffer, double
    return frombuffer( g.GetEX(), double, npoints)\
         , frombuffer( g.GetEY(), double, npoints)
##end def function
from ROOT import TGraphErrors
setattr( TGraphErrors, 'getErrBuffers', tGraphGetErrBuffers )

def tGraphGetAsymmErrBuffers( g ):
    """Get TGraph x and y buffers"""
    npoints = g.GetN()
    if npoints==0: return None, None, None, None
    from numpy import frombuffer, double
    return frombuffer( g.GetEXlow(), double, npoints)\
         , frombuffer( g.GetEXhigh(), double, npoints)\
         , frombuffer( g.GetEYlow(), double, npoints)\
         , frombuffer( g.GetEYhigh(), double, npoints)
##end def function
from ROOT import TGraphAsymmErrors
setattr( TGraphAsymmErrors, 'getErrBuffers', tGraphGetAsymmErrBuffers )

def histConvert1( h, cls, newname=None, addbins=None, multbins=None, noerr=False ):
    """
    Convert TH1X to TH1Y
    addbins=(n, x2a) will add n bins to the end of the histogram increasing x2 to x2a
                     new bins will be empty
    multbins=n will rebin the resulting histogram increasing the number of bins n times
             each bin will be split in n bins with n times smaller height 
             and n times smaller err^2
    addbins and multbins can not be applied in the same time
    """
    ax = h.GetXaxis()
    xbins = ax.GetXbins()
    newh = None
    if newname==None: newname=h.GetName()
    n = h.GetNbinsX() 
    if xbins.GetSize()==0:
        x1, x2 = ax.GetXmin(), ax.GetXmax()
        if addbins:
            assert multbins==None, 'Can not add/mult bins in the same time'
            na, x2a = addbins
            newh = cls( newname, h.GetTitle(), na, x1, x2a )
        elif multbins:
            na = n*multbins
            newh = cls( newname, h.GetTitle(), na, x1, x2 )
        else:
            newh = cls( newname, h.GetTitle(), n, x1, x2 )
        ##end if
    else:
        assert addbins==None, 'Can not add bins to variable bins histogram'
        newh = cls( newname, h.GetTitle(), h.GetNbinsX(), xbins.GetArray() )
    ##end if
    err = h.getErrBuffer( flows=True ) if not noerr else None
    b = h.getBuffer( flows=True )
    newb = newh.getBuffer( flows=True )

    if multbins:
        if err!=None:
            newh.Sumw2()
            newerr = newh.getErrBuffer( flows=True )
            newerr[0], newerr[-1] = err[0], err[-1]
            newerr1=newerr[1:-1].reshape( ( multbins, n ), order='F' )
            newerr1[:]=err[1:-1]/float(multbins)
        ##end if err

        newb[0], newb[-1] = b[0], b[-1]
        newb1=newb[1:-1].reshape( ( multbins, n ), order='F' )
        newb1[:]=b[1:-1]/float(multbins)
    else:
        if err!=None:
            newh.Sumw2()
            newerr = newh.getErrBuffer( flows=True )
            newerr[:n+2] = err[:]
        ##end if err

        newb[:n+2] = b[:]
    ##end if multbins

    newh.SetEntries( h.GetEntries() )

    return newh
##end def histConvert
from ROOT import TH1
setattr( TH1, 'convert', histConvert1 )

def vectortGetBuffer( self ):
    from numpy import frombuffer, dtype
    buf = self.GetMatrixArray()
    return frombuffer( buf, dtype( buf.typecode ), self.GetNoElements() )
##end def vectortGetBuffer
from ROOT import TVectorD, TVectorF
setattr( TVectorD, 'getBuffer', vectortGetBuffer )
setattr( TVectorF, 'getBuffer', vectortGetBuffer )

def matrixGetBuffer( m ):
    from numpy import frombuffer, dtype
    buf = frombuffer( cbuf, dtype( cbuf.typecode ), m.GetNoElements() )
    return buf.reshape( m.GetNrows(), m.GetNcols() )
##end def
from ROOT import TMatrixD, TMatrixF
setattr( TMatrixD, 'getBuffer', matrixGetBuffer )
setattr( TMatrixF, 'getBuffer', matrixGetBuffer )

def histGetBuffer1( h, flows=False ):
    """Return histogram data buffer
    if flows=False, exclude underflow and overflow
    """
    from numpy import frombuffer, dtype
    buf = h.GetArray()
    buf = frombuffer(buf , dtype( buf.typecode ), h.GetNbinsX()+2 )
    if not flows: buf = buf[1:-1] 
    return buf
##end def histGetBuffer
from ROOT import TH1
setattr( TH1, 'getBuffer', histGetBuffer1 )

def histGetErrBuffer1( h, flows=False ):
    """Return histogram error buffer
    if flows=False, exclude underflow and overflow
    """
    sw2 = h.GetSumw2()
    if sw2.GetSize()==0: return None 

    from numpy import frombuffer, double
    buf = frombuffer( sw2.GetArray(), double, h.GetNbinsX()+2 )
    if not flows: buf = buf[1:-1] 
    return buf
##end def histGetErrBuffer1
from ROOT import TH1
setattr( TH1, 'getErrBuffer', histGetErrBuffer1 )

def histGetBuffer2( h, flows=False, mask=None ):
    """Return histogram data buffer
    if flows=False, exclude underflow and overflow
    if mask=0.0 than bins with 0.0 content will be white, but not colored
    NOTE: buf[biny][binx] is the right access signature
    """
    from numpy import frombuffer, dtype
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    buf = h.GetArray()
    buf = frombuffer( buf, dtype( buf.typecode ), (nx+2)*(ny+2) ).reshape( ( ny+2, nx+2 ) )
    if mask!=None:
        from numpy import ma
        buf = ma.array( buf, mask = buf==mask )
    ##end if mask!=None
    if flows: return buf

    buf = buf[1:ny+1,1:nx+1]
    return buf
##end def histGetBuffer
from ROOT import TH2
setattr( TH2, 'getBuffer', histGetBuffer2 )

def axisGetBinEdges( ax, type=False ):
    """Get the array with bin edges"""
    xbins = ax.GetXbins()
    n = xbins.GetSize()
    lims=None
    fixed = False
    if n>0:
        from numpy import frombuffer, double
        lims = frombuffer( xbins.GetArray(), double, n )
        fixed = False
    else:
        from numpy import linspace
        lims = linspace( ax.GetXmin(), ax.GetXmax(), ax.GetNbins()+1 )
        fixed = True
    ##end if 

    if type: return lims, fixed
    return lims
##end def
from ROOT import TAxis
setattr( TAxis, 'getBinEdges', axisGetBinEdges )

