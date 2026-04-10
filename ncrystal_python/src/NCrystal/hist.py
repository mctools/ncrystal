
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

"""A few simple Histogram classes to support other NCrystal features.

Hist1D: A histogram class which is primarily intended as an objected-oriented
        representation of the histogram data returned by the MiniMC framework,
        and which allows for convenient plotting/dumping of contents,
        serialisation, and statistical comparisons. It can be instantiated
        either from the JSON data produced by the NCrystal C++ layer as a result
        of MiniMC simulations, or via the HistFiller1D class.

HistFiller1D: A class which can be used to efficiently histogram data which is
              in the form of numpy arrays, making it suitable for usage in for
              instance MiniMC callback functions. Once all data has been filled
              into the histogram, it can be converted to a Hist1D object by
              calling its .to_hist1d() method.
"""

__all__ = ['Hist1D','HistFiller1D']

from ._numpy import _np, _ensure_numpy, _np_linspace
from .exceptions import NCBadInput, NCCalcError
import math

class Hist1D:

    """Histogram class which is primarily intended as an objected-oriented
    representation of the histogram data returned by the MiniMC framework.

    The histograms contents are created on the C++ side (using a class declared
    the NCHists.hh C++ header file), and are usually exported as JSON data. This
    data can be loaded via the present Python API in order to facilitate
    plotting and data analysis. Users should not normally instantiate these
    objects directly themselves.

    """

    def __init__( self, data ):
        """Initialise histogram. This is intended to be using data provided by
        the C++ histogram's toJSON() method, or from the .to_dict() or
        .to_json() of this Python Hist1D class itself. Most users should not
        have a need to call this directly.
        """

        if data=='_no_init_':
            return
        _ensure_numpy()
        if isinstance(data,str):
            import json
            data = json.loads(data)

        self.__stat_integral = None
        self.__stat_rms = None
        self.__stat_mean = None
        self.__stat_minfilled = None
        self.__stat_maxfilled = None

        def float_or_None( x ):
            return float(x) if x is not None else None
        if 'stats' in data:
            _ = data['stats']
            self.__stat_integral = float_or_None(_.get('integral'))
            self.__stat_rms = float_or_None(_.get('rms'))
            self.__stat_mean = float_or_None(_.get('mean'))
            self.__stat_minfilled = float_or_None(_.get('minfilled'))
            self.__stat_maxfilled = float_or_None(_.get('maxfilled'))

        ns = sum(int(e is not None) for e in [self.__stat_rms,
                                              self.__stat_mean,
                                              self.__stat_minfilled,
                                              self.__stat_maxfilled])
        if ns not in (0,4):
            raise NCBadInput('Inconsistent input: all of rms/mean/minfilled'
                             '/maxfilled stats must be available if any is.')

        self.__title = data.get('title')
        hb = data['bindata']
        self.__flow_under = float_or_None(hb.get('underflow'))
        self.__flow_over = float_or_None(hb.get('overflow'))
        self.__flow_under_errorsq = float_or_None(hb.get('underflow_errorsq'))
        self.__flow_over_errorsq = float_or_None(hb.get('overflow_errorsq'))

        if (self.__flow_under is None) != (self.__flow_over is None):
            raise NCBadInput('Inconsistent input: both under and overflow'
                             ' must be available if one of them is.')
        if (self.__flow_under_errorsq is None) != (self.__flow_over_errorsq is None):
            raise NCBadInput('Inconsistent input: both under and overflow'
                             ' error^2 info must be available if one of them is.')
        self.__xmin = hb['xmin']
        self.__xmax = hb['xmax']
        self.__nbins = hb['nbins']
        self.__y = _np.asarray(hb['content'],dtype=float)
        esq = hb.get('errorsq')
        self.__yerrsq = ( _np.asarray(esq,dtype=float)
                          if esq is not None else None )
        self.__yerr = None

        assert self.__nbins == len(self.__y)
        assert self.__yerrsq is None or self.__nbins == len(self.__yerrsq)

        if self.__stat_integral is None:
            self.__stat_integral = self.integrate_bins()[0]

    @classmethod
    def objectify_data( cls, data ):
        """Traverse data and replace any dictionaries inside it that are really
        the result of a Hist1D.to_dict() (or histogram JSON data) with a
        corresponding Hist1D object.  Returns a new object with the result of
        the transformation.
        """
        import copy
        def o( d ):
            if isinstance(d,list):
                return [o(e) for e in d]
            if isinstance(d,tuple):
                return tuple([o(e) for e in d])
            if isinstance(d,dict):
                if d.get('datatype')=='NCrystalHist1D_v1':
                    return cls(d)
                else:
                    return dict( (copy.deepcopy(k),o(v))
                                 for k,v in d.items() )
            #something else, just pass through:
            return copy.deepcopy( d )
        return o( data )

    def _to_json_compat_object( self ):
        return self.to_dict( json_compat = True )

    def to_dict( self, json_compat = False ):
        """Serialise as dictionary. The returned string can be used to
        initialise a Hist1D object again. If json_compat is True, the returned
        dictionary will contain lists rather than numpy arrays.
        """
        return dict( datatype = 'NCrystalHist1D_v1',
                     stats = self.stats,
                     title = self.__title,
                     bindata = self.bindata( json_compat = json_compat ) )

    def to_json( self ):
        """Serialise as JSON string. The returned string can be used to
        initialise a Hist1D object again."""
        import json
        return json.dumps( self.to_dict( json_compat = True ) )

    def clone( self, rebin_factor = 1 ):
        """Return clone of object. Optionally, the resulting histogram can be
        rebinning by a given factor, which must be a divisor in the current
        number of bins.
        """
        c = Hist1D('_no_init_')
        c.__stat_integral = self.__stat_integral
        c.__stat_rms = self.__stat_rms
        c.__stat_mean = self.__stat_mean
        c.__stat_minfilled = self.__stat_minfilled
        c.__stat_maxfilled = self.__stat_maxfilled
        c.__flow_under = self.__flow_under
        c.__flow_over = self.__flow_over
        c.__flow_under_errorsq = self.__flow_under_errorsq
        c.__flow_over_errorsq = self.__flow_over_errorsq
        c.__title = self.__title
        c.__xmin = self.__xmin
        c.__xmax = self.__xmax
        c.__nbins = self.__nbins
        c.__y = self.__y.copy()
        c.__yerrsq = self.__yerrsq.copy() if self.__yerrsq is not None else None
        c.__yerr = self.__yerr.copy() if self.__yerr is not None else None
        if rebin_factor > 1:
            c.__rebin( rebin_factor )
        return c

    def integrate( self, xlow, xhigh, tolerance = 1e-5 ):
        """Returns integrated contents of the histogram over the area
        [xlow,xhigh] along with the error of that value in a tuple
        (content,error).

        This is done translating xlow and xhigh to exact bin edges and then
        calling integrate_bins. If that is not possible within the
        tolerance, an exception is raised.

        If xlow is None it will taken to be xmin and underflow content will be
        included. Likewise, if xhigh is None it will be taken to be xmax and
        overflow content will be included.

        """
        include_underflow = False
        include_overflow = False
        if xlow is None:
            xlow = self.xmin
            include_underflow = True
        if xhigh is None:
            xhigh = self.xmax
            include_overflow = True

        if not ( xhigh >= xlow ):
            raise NCBadInput('Invalid integration range requested.')

        bw = self.binwidth
        def _findedge(x):
            if x <= self.__xmin:
                return 0
            if x >= self.__xmax:
                return self.nbins
            r = ( x - self.__xmin ) / bw
            ir = int(r+0.5)
            if abs(r-ir) > tolerance:
                raise NCBadInput(f'Value {x} does not correspond exactly'
                                 ' to a bin edge within the tolerance.')
            return ir
        e_low = _findedge(xlow)
        e_high = _findedge(xhigh)
        if e_low == e_high:
            return ( 0.0, 0.0 )
        assert e_low >= 0
        assert e_low < e_high <= self.nbins
        assert not include_underflow or e_low==0
        assert not include_overflow or e_high==self.nbins
        return self.integrate_bins( None if include_underflow else e_low,
                                    None if include_overflow else e_high )

    @property
    def integral( self ):
        """Returns integral, including values in under and overflow bins."""
        return self.__stat_integral

    @property
    def has_stats( self ):
        """Whether stats (mean/rms/minfilled/maxfilled) are available."""
        return self.__stat_mean is not None

    @property
    def has_flow( self ):
        """Whether under/overflow contents are available."""
        return self.__flow_under is not None

    @property
    def mean( self ):
        """Returns mean value or None if not available."""
        return self.__stat_mean

    @property
    def rms( self ):
        """Returns root-mean-square value or None if not available."""
        return self.__stat_rms

    @property
    def minfilled( self ):
        """Returns smallest x-value filled, or None if not available."""
        return self.__stat_minfilled

    @property
    def maxfilled( self ):
        """Returns largest x-value filled, or None if not available."""
        return self.__stat_maxfilled

    @property
    def underflow( self ):
        """Returns underflow content or None if not available."""
        return self.__flow_under

    @property
    def overflow( self ):
        """Returns overflow content or None if not available."""
        return self.__flow_over

    @property
    def underflow_error( self ):
        """Returns error on underflow content or None if not available."""
        if self.__flow_under is not None:
            return math.sqrt( self.__flow_under_errorsq
                              if self.__flow_under_errorsq is not None
                              else self.__flow_under)

    @property
    def overflow_error( self ):
        """Returns error on overflow content or None if not available."""
        if self.__flow_over is not None:
            return math.sqrt( self.__flow_over_errorsq
                              if self.__flow_over_errorsq is not None
                              else self.__flow_over)

    @property
    def empty( self ):
        """Returns integral, including values in under and overflow bins."""
        return not (self.integral>0.0)

    def integrate_bins( self, bin_low = None, bin_up = None ):
        """Returns integrated contents of the bins [bin_low,bin_up[ along with
        the error of that value in a tuple (content,error).

        If bin_low is None the integration will start at the first bin and
        include the underflow bin.

        If bin_up is None the integration will end at the last bin and include
        the overflow bin. To include the last regular bin but NOT the overflow
        content, bin_up should equal nbins.

        """

        add_overflow, add_underflow = False, False
        if bin_low is None:
            add_underflow = True
            bin_low = 0
            underflow_c = self.__flow_under
            underflow_e2 = self.__flow_under_errorsq
            assert ( (underflow_c is None) ==
                     (underflow_e2 is None) ), 'Inconsistent underflow info.'
            if underflow_c is None:
                add_underflow = False

        if bin_up is None:
            add_overflow = True
            bin_up = self.__nbins
            overflow_c = self.__flow_over
            overflow_e2 = self.__flow_over_errorsq
            assert ( (overflow_c is None) ==
                     (overflow_e2 is None) ), 'Inconsistent overflow info.'
            if overflow_c is None:
                add_overflow = False

        bin_low, bin_up = int(bin_low), int(bin_up)
        if bin_up < bin_low or bin_low<0 or bin_up > self.__nbins:
            raise NCBadInput('Invalid bin range requested.')
        content_integral = self.__y[bin_low:bin_up].sum()

        if self.__yerrsq is None:
            #unweighted, just base errors on contents:
            errorsq_integral = content_integral
        else:
            errorsq_integral = self.__yerrsq[bin_low:bin_up].sum()

        if add_underflow:
            content_integral += underflow_c
            errorsq_integral += underflow_e2
        if add_overflow:
            content_integral += overflow_c
            errorsq_integral += overflow_e2

        return ( float(content_integral), math.sqrt(errorsq_integral) )

    def adopt_contents( self, other_hist, keep_title = False ):
        """Adopt the contents (including binning and statistics) of other
        histogram in place of the current contents. If keep_title is True, the
        title will not be updated. Returns self.
        """
        o = other_hist
        self.__xmin == o.__xmin
        self.__xmax == o.__xmax
        self.__nbins == o.__nbins
        self.__stat_integral = o.__stat_integral
        self.__stat_rms = o.__stat_rms
        self.__stat_mean = o.__stat_mean
        self.__stat_minfilled = o.__stat_minfilled
        self.__stat_maxfilled = o.__stat_maxfilled
        self.__flow_under = o.__flow_under
        self.__flow_over = o.__flow_over
        self.__flow_under_errorsq = o.__flow_under_errorsq
        self.__flow_over_errorsq = o.__flow_over_errorsq
        self.__yerr = None
        self.__y = o.__y.copy()
        self.__yerrsq = ( o.__yerrsq.copy()
                          if o.__yerrsq is not None else None )
        if not keep_title:
            self.__title = o.__title
        return self

    def add_contents( self, other_hist ):
        """Add the contents of other histogram to this histogram. The two
        histograms should have compatible binnings and over/underflow
        settings. Statistics will be updated if available on both histograms,
        otherwise dropped. The title will be kept. Returns self.

        """
        assert isinstance(other_hist,Hist1D)
        o = other_hist
        if self.binning != o.binning:
            raise NCBadInput("incompatible binning (%i,%g,%g) vs. (%i,%g,%g)."
                             %( *self.binning, *o.binning))

        def flow_sig( obj ):
            return ( obj.__flow_under is None,
                     obj.__flow_under_errorsq is None,
                     obj.__flow_over is None,
                     obj.__flow_over_errorsq is None )

        if flow_sig(self) != flow_sig(o):
            raise NCBadInput("incompatible under/overflow settings.")

        assert self.__xmin == o.__xmin
        assert self.__xmax == o.__xmax
        assert self.__nbins == o.__nbins
        if not ( o.integral > 0.0 ):
            #No contents in other histogram
            return self # noop
        if not self.integral:
            #simply adopt contents of other histogram
            return self.adopt_contents( o, keep_title = True )

        if self.mean is None or o.mean is None:
            self.__stat_integral += o.__stat_integral
            self.__stat_rms = None
            self.__stat_mean = None
            self.__stat_minfilled = None
            self.__stat_maxfilled = None
        else:
            sumw = self.integral
            assert sumw > 0.0
            sumwx = self.mean * sumw
            o_sumw = o.integral
            o_sumwx = o.mean * o_sumw
            sumwx2 = (self.rms**2 + self.mean**2)*sumw
            o_sumwx2 = (o.rms**2 + o.mean**2)*o_sumw
            sumw += o_sumw
            sumwx += o_sumwx
            sumwx2 += o_sumwx2
            self.__stat_integral = sumw
            self.__stat_mean = sumwx / sumw
            rms2 = sumwx2 / sumw - (self.__stat_mean)**2
            self.__stat_rms = math.sqrt(abs(rms2))#abs() vs num. issues.
            if o.__stat_minfilled < self.__stat_minfilled:
                self.__stat_minfilled = o.__stat_minfilled
            if o.__stat_maxfilled > self.__stat_maxfilled:
                self.__stat_maxfilled = o.__stat_maxfilled

        if self.__flow_under is not None:
            self.__flow_under += o.__flow_under
        if self.__flow_under_errorsq is not None:
            self.__flow_under_errorsq += o.__flow_under_errorsq
        if self.__flow_over is not None:
            self.__flow_over += o.__flow_over
        if self.__flow_over_errorsq is not None:
            self.__flow_over_errorsq += o.__flow_over_errorsq

        self.__yerr = None
        if self.__yerrsq is None:
            if o.__yerrsq is None:
                pass#done
            else:
                self.__yerrsq = self.__y.copy()
                self.__yerrsq += o.__yerrsq
        else:
            if o.__yerrsq is None:
                self.__yerrsq += o.__y
            else:
                self.__yerrsq += o.__yerrsq
        self.__y += o.__y

    def __rebin( self, rebin_factor ):
        if rebin_factor == 1:
            return
        if not rebin_factor>=1 or not self.__nbins % rebin_factor == 0:
            raise NCBadInput( f'Can not rebin nbins={self.__nbins} with '
                              f'rebin_factor={rebin_factor} (rebin_factor'
                              ' must be a divisor of nbins).' )
        def _dorebin(x):
            return _np.sum( x.reshape( len(x)//rebin_factor, rebin_factor ),
                            axis=1)
        self.__y = _dorebin(self.__y)
        self.__yerr = None
        if self.__yerrsq is not None:
            self.__yerrsq = _dorebin(self.__yerrsq)
        self.__nbins = self.__nbins // rebin_factor
        assert self.__yerrsq is None or len(self.__yerrsq)==len(self.__y)
        assert self.__nbins==len(self.__y)

    @property
    def stats( self ):
        """Return statistics of unbinned filled values as array. This include
        integral, rms, mean, as well as the range of filled values."""
        return dict ( integral = self.__stat_integral,
                      rms = self.__stat_rms,
                      mean = self.__stat_mean,
                      minfilled = self.__stat_minfilled,
                      maxfilled = self.__stat_maxfilled )

    def bindata( self, json_compat = False ):
        """Return bindata as dictionary. This include binnings, contents,
        errors, as well as under/overflow information. If json_compat is True,
        the returned dictionary will contain lists rather than numpy arrays."""
        def export_array( a ):
            return ( a.tolist()
                     if (json_compat and hasattr(a,'tolist')) else a.copy() )

        d = dict( xmin = self.xmin,
                  xmax = self.xmax,
                  nbins = self.nbins,
                  content = export_array(self.__y) )
        if self.__yerrsq is not None:
            d['errorsq'] = export_array( self.__yerrsq )
        def add_if_not_None( key, val ):
            if val is not None:
                d[key] = val
        add_if_not_None('underflow_errorsq',self.__flow_under_errorsq)
        add_if_not_None('overflow_errorsq',self.__flow_over_errorsq)
        add_if_not_None('underflow',self.__flow_under)
        add_if_not_None('overflow',self.__flow_over)
        return d


    @property
    def contents( self ):
        """Return numpy array of bin contents."""
        return self.__y

    @property
    def errors( self ):
        """Return numpy array of bin errors."""
        if self.__yerr is None:
            self.__yerr = _np.sqrt( self.errors_squared )
        return self.__yerr

    @property
    def errors_squared( self ):
        """Return numpy array of bin errors squared."""
        return ( self.__yerrsq
                 if self.__yerrsq is not None
                 else self.__y )


    def set_title( self, title ):
        """Sets the title and returns self."""
        assert isinstance(title,str)
        self.__title = title
        return self

    @property
    def title( self ):
        """The title assigned to the histogram, or None if no title was
        assigned. The title can be modified via the set_title method
        """
        return self.__title

    @property
    def xmin( self ):
        """The lower edge of the bin range."""
        return self.__xmin

    @property
    def xmax( self ):
        """The upper edge of the bin range."""
        return self.__xmax

    @property
    def binning( self ):
        """Tuple of (nbins,xmin,xmax)"""
        return ( self.__nbins, self.__xmin, self.__xmax )

    @property
    def binwidth( self ):
        """Binwidth."""
        return (self.__xmax-self.__xmin)/self.__nbins

    @property
    def nbins( self ):
        """Number of bins."""
        return self.__nbins

    @property
    def bincenters( self ):
        """Create and return a numpy array of bin centers."""
        halfbw = 0.5*self.binwidth
        return _np_linspace(self.__xmin+halfbw,
                            self.__xmax-halfbw,
                            self.__nbins)

    @property
    def binedges( self ):
        """Create and return a numpy array of bin edges."""
        return _np_linspace(self.__xmin,
                            self.__xmax,
                            self.__nbins+1)


    def dump( self,
              prefix = '',
              contents = 'auto',
              highres = False,
              do_print = True, ):
        """Creates an overview of the histogram in printable string form,
        optionally prefixing each line with provided prefix. If do_print is
        True, this will be printed directly, and otherwise it will be
        returned. By default, the contents of the histogram will only be
        included in case of 50 or fewer bins. Setting contents to True or False
        can be used to explicitly override this.  If highres is True, all
        floating point numbers will be printed with full
        precision. Alternatively it can be a fmt string like '%.10g'.

        """
        if contents == 'auto':
            contents = self.nbins <= 50

        result = []
        def pr(s):
            result.append(s)

        p = prefix
        def fmt_highres( x ):
            s = '%.14g'%x
            return s if float(s)==x else '%.19g'%x
        def fmt( x ):
            return '%.6g'%x
        if highres:
            if isinstance(highres,str):
                def fmt_highres( x ):
                    return highres % x
            fmt = fmt_highres

        pr("%sHist1D(nbins=%i,xmin=%s,xmax=%s):"%( p,
                                                   self.nbins,
                                                   fmt_highres(self.xmin),
                                                   fmt_highres(self.xmax)))
        pr("%s  title     : %s"%(p, self.title or "<none>"))
        pr("%s  integral  : %s"%(p,fmt(self.integral)))
        if self.has_stats:
            pr("%s  mean      : %s"%(p,fmt(self.mean)))
            pr("%s  rms       : %s"%(p,fmt(self.rms)))
            pr("%s  minfilled : %s"%(p,fmt_highres(self.minfilled)))
            pr("%s  maxfilled : %s"%(p,fmt_highres(self.maxfilled)))
        else:
            pr("%s  mean      : <n/a>"%p)
            pr("%s  rms       : <n/a>"%p)
            pr("%s  minfilled : <n/a>"%p)
            pr("%s  maxfilled : <n/a>"%p)
        if self.has_flow:
            pr("%s  underflow : %s +- %s"%(p,fmt(self.underflow),
                                           fmt(self.underflow_error)))
            pr("%s  overflow  : %s +- %s"%(p,fmt(self.overflow),
                                           fmt(self.overflow_error)))

        if contents:
            c,e = self.contents, self.errors
            for ibin in range(self.nbins):
                pr("%s  content[ibin=%i] : %s +- %s"%(p,ibin,
                                                      fmt(c[ibin]),
                                                      fmt(e[ibin])))

        result.append('')
        result = '\n'.join(result)
        if do_print:
            from ._common import print as ncprint
            ncprint(result,end='')
        else:
            return result

    def _hist_curve( self ):
        be = self.binedges
        y = self.contents
        #if error_offset:
        #    y += self.errors * error_offset
        cx = _np.empty(self.__nbins*2)
        cy = _np.empty(self.__nbins*2)
        i = 0
        for ibin in range(self.__nbins):
            cx[i] = be[ibin]
            cy[i] = y[ibin]
            i+=1
            cx[i] = be[ibin+1]
            cy[i] = y[ibin]
            i+=1
        return cx,cy

    def errorbar_args( self, style = True, **kwargs ):
        """Return a kwargs dictionary suitable for using with the matplotlib
        errorbar function. If style is False, no parameters related to plotting
        style or colouring are included. Any excess kwargs are simply passed
        along to the returned dictionary.
        """
        d = {'x':self.bincenters,
             'y':self.contents,'xerr':0.5*self.binwidth,
             'yerr':self.errors }
        if style:
            d.update({'fmt':'.',#dont connect with line
                      'mec':'black','mfc':'black',
                      #'ms':4,'mew':1,
                      'ecolor':'black','elinewidth':1.0 })
        d.update(kwargs)
        return d

    def bar_args( self, **kwargs ):
        """Return a kwargs dictionary suitable for using with the matplotlib bar
        function. Any excess kwargs are simply passed along to the returned
        dictionary.
        """
        d = {'x' : self.binedges[:-1],
             'height': self.contents,
             'width': self.binwidth,
             'align':'edge'
             }
        d.update(kwargs)
        return d

    def plot( self, label=None, show_errors=True, set_xlim = True,
              logy = False, error_bands = None, alpha = None, color = None,
              title = True, do_grid = True, do_legend=False, **kw_plot ):
        """Produce a matplotlib plot of the histogram. If plt is None,
        matplotlib.pyplot is used. If axis is None, plt.gca() is used. Unless
        do_show is False, plt.show() wil be called ultimately. The show_errors,
        label, set_xlim, alpha, and color are all options whose affects are
        hopefully self explanatory to users with experience in matplotlib.

        If error bands is set to a positive value, the curve and errors are
        instead shown as a coloured band spanning y_i +- error_bands*yerr_i.

        The title parameter can be True, False, or a string. If not False, a
        title will be added to the plot.
        """
        from .plot import PlotContext
        pctx = PlotContext(**kw_plot).check_unused()
        if error_bands:
            #need to repeat last entry for proper error band visualisation of
            #the last bin:
            def repeat_last( x ):
                return _np.concatenate( (x, _np.asarray( [x[-1]],
                                                         dtype=x.dtype ) ) )
            cc = repeat_last(self.contents)
            ee = repeat_last(self.errors)*error_bands
            fill_between_args = dict( x = self.binedges, step = 'post',
                                      y1 = cc - ee, y2 = cc + ee )
            pctx.axis.fill_between(**fill_between_args,alpha=alpha,color=color,
                                   label=label)
        else:
            if color != 'none':
                pctx.axis.bar( **self.bar_args( label = label ),
                               alpha=alpha, color=color )
            if show_errors:
                pctx.axis.errorbar(**self.errorbar_args(),alpha=alpha,
                                   label = label if color=='none' else None)
        xmin,xmax,binwidth = self.xmin, self.xmax, self.binwidth
        if set_xlim:
            pctx.axis.set_xlim(xmin-1e-6*binwidth,xmax+1e-6*binwidth)
        if logy:
            pctx.axis.semilogy()
        if title:
            t = title if isinstance(title,str) else self.title
            if t:
                pctx.axis.set_title(t)
        return pctx.finalise(do_grid=do_grid,do_legend=do_legend)

    def scale( self, factor ):
        """Scale contents by a positive factor and return self. This also
        affects errors and under/overflow statistics."""
        if factor == 1.0:
            return self
        if not ( 0.0 < factor < 1e150 ):
            raise NCBadInput('histogram scale factor out of range.')
        self.__yerr = None
        #Force error array from this point:
        if self.__yerrsq is None:
            self.__yerrsq = self.__y.copy()
        self.__yerrsq *= factor**2
        self.__y *= factor
        self.__stat_integral *= factor

        #NB: rms/mean/minfilled/maxfilled: no update needed!

        if self.__flow_under is not None:
            assert self.__flow_over is not None
            self.__flow_under *= factor
            self.__flow_over *= factor

        if self.__flow_under_errorsq is not None:
            assert self.__flow_over_errorsq is not None
            self.__flow_under_errorsq *= factor**2
            self.__flow_over_errorsq *= factor**2

        return self

    def norm( self ):
        """Scale histogram with 1/integral and return self. Does nothing if
        integral is not positive."""
        f = self.integral
        if f > 0.0:
            self.scale( 1.0 / f )
        return self

    def __eq__(self, o ):
        if id(self) == id(o):
            return True
        if ( self.binning != o.binning
             or self.stats != o.stats
             or self.title != o.title
             or not _np.array_equal(self.__y,o.__y) ):
            return False
        if ( self.__flow_under != o.__flow_under
             or self.__flow_under_errorsq != o.__flow_under_errorsq
             or self.__flow_over != o.__flow_over
             or self.__flow_over_errorsq != o.__flow_over_errorsq ):
            return False
        return _np.array_equal( self.errors_squared, o.errors_squared )

    def check_compat( self, other_hist,
                      threshold = 0.05,
                      force_norm = True,
                      return_pval = False,
                      check=False ):
        """Finds the p-value based on a chi-square comparison between the
        contents of this and another histogram, which must have identical
        binning. If force_norm is False, an error will be produced if histograms
        do not have identical statistics (integral). A high p-value indicates
        that the histogram contents are likely to come from the same underlying
        distribution. By default, the calculated p-value is compared against the
        threshold and True is returned if it is higher, False if it is
        lower. If return_pval is True, the p-value is instead returned. If check
        is True, a CalcError exception is raised in case the p-value fails to
        exceed threshold.

        Note, if you are running N histogram comparisons (e.g. in a test suite),
        and you want to have a p=5% chance of false positive in case you
        e.g. change the random seed of the test suite, you should use a
        threshold of roughly p/N=0.05/N (this is called Bonferroni
        correction and was confirmed to work well with simple Monte Carlo test).
        """
        chisq, k = self.chisquare_dist( other_hist,
                                        force_norm = force_norm )
        p_value = 1.0 - _chisq_cdf(chisq, k)
        assert -1e-6 <= p_value <= 1.0+1e-6
        p_value = min( 1.0, max( 0.0, p_value ) )
        ok = ( p_value >= threshold )
        if check and not ok:
            #fmt for unit test reproducibility
            fmt='p-value=%g'%p_value if p_value>=1e-6 else 'p-value'
            raise NCCalcError(f'check_compat failed: {fmt}'
                             f' is not greater than {threshold}.')
        if return_pval:
            return p_value
        return ok

    def chisquare_dist( self, other, force_norm = True ):
        """Estimate chi-square difference to contents of other histogram, which
        must have identical binning. If force_norm is False, an error will be
        produced if histograms do not have identical statistics
        (integral). Returns both the chi-square value and the estimated degrees
        of freedom k in a tuple. The degrees of freedom is simply estimated as
        one less than the number of bins (not counting bins empty in both
        histograms), but at least 1.

        """

        if self.binning != other.binning:
            raise NCBadInput('chisquare_dist: incompatible binnings '
                             '(%s vs %s).'%(self.binning,other.binning))
        if self.has_flow != other.has_flow:
            raise NCBadInput('chisquare_dist: '
                             'incompatible overflow settings.')

        s, o = self, other
        if force_norm:
            s, o = s.clone().norm(), o.clone().norm()
        else:
            if s.integral != o.integral:
                raise NCBadInput('chisquare_dist: '
                                 'incompatible histogram integrals.')

        c1,e1 = s.contents, s.errors
        c2,e2 = o.contents, o.errors

        #use mask to avoid zero division warnings in bins where both datasets
        #have no content.
        mask = (e1>0.0) | (e2>0.0)
        c1, e1, c2, e2 = c1[mask], e1[mask], c2[mask], e2[mask]

        n = len(c1)
        if n > 0:
            c1sum = c1.sum()+(s.underflow or 0.0)+(s.overflow or 0.0)
            c2sum = c2.sum()+(o.underflow or 0.0)+(o.overflow or 0.0)
            assert ( c1sum==0.0
                     or c2sum==0.0
                     or (c1sum-c2sum)<1e-5*(c1sum+c2sum) ),"must be normalised"
            chi_squared = ( (c1-c2)**2 / (e1**2 + e2**2) ).sum()
        else:
            chi_squared = 0.0
        if s.has_flow:
            n += 2
            a, b = 0.0, 0.0
            if s.underflow>0.0 or o.underflow>0.0:
                a = ( (s.underflow - o.underflow)**2
                      / ( s.underflow_error**2 + o.underflow_error**2 ) )
            if s.overflow>0.0 or o.overflow>0.0:
                b = ( (s.overflow - o.overflow)**2
                      / ( s.overflow_error**2 + o.overflow_error**2 ) )
            chi_squared += ( a + b )

        return ( chi_squared, max( 1, n - 1) )

def _chisq_cdf( x, k ):
    """Chi-squared CDF function. Similar to scipy scipy.stats.chi2.cdf, but with
    a slightly lower precision."""
    if x<=1e-10*k:
        return 0.0
    if k==1:
        return math.erf(math.sqrt(x/2))
    if k==2:
        return 1.0 - math.exp(-0.5*x)
    #Approx. of Luisa Canal (2005), https://math.stackexchange.com/a/4915053
    #We use it all the way down to k=3, but it seems to be OK for our purposes.
    u = x / k
    t = math.fsum([u**(1/6),-(1/2)*u**(1/3),+(1/3)*u**(1/2),
                   -5/6,(1/(9*k)),(7/(648*k*k)),-(25/(2187*k*k*k))])
    s = math.fsum([(2/(18*k)),(2/(162*k*k)),-(74/(11664*k*k*k))])
    return 0.5*(1.0+math.erf(t/math.sqrt(s)))

class HistFiller1D:
    """Mutable utility class which is intended for filling numpy arrays into 1D
    histograms, and ultimately either simply accessing the bin data with the
    .bindata() method, or converting the results to a Hist1D object with the
    to_hist1d() method, which contains more functionality for things like
    persistification, plotting, and statistical comparisons.

    To prepare a histogram of 100 bins from xmin=0.0 to xmax=1.0, with a title
    'My histogram' you can use either of the calls:

      h = HistFiller1D( 100, 0.0, 1.0, 'My histogram' )
      h = HistFiller1D( (100, 0.0, 1.0), title = 'My histogram' )

    The tuple (100, 0.0, 1.0) is also what h.binning returns afterwards. To
    create a new empty histogram with the same binning (optionally overriding
    the title) you can use h.clone_empty('new title').

    To fill array data into the histogram, use the fill method which accepts
    numpy arrays of data (if weights are not supplied, they will be assumed to
    be unity):

      h.fill( values, weights )

    """

    def __init__(self, nbins_or_binning, xmin=None, xmax=None, title = None):
        """Initialise a HistFiller1D object. If xmin=xmax=None, the
        nbins_or_binning argument must be a tuple of (nbins,xmin,xmax) values,
        otherwise it must be nbins.
        """
        _ensure_numpy()
        if xmin is None and xmax is None:
            nbins, xmin, xmax = nbins_or_binning
        else:
            nbins = nbins_or_binning
        del nbins_or_binning
        xmin, xmax = float(xmin), float(xmax)
        assert int(nbins)==nbins and nbins >= 1
        assert xmin < xmax
        self.__nbins = nbins
        self.__xmin, self.__xmax = float(xmin), float(xmax)
        self.__title = None if title is None else str(title)
        self.__invbw = nbins/(xmax-xmin)
        self.__bin_edges = _np_linspace(xmin, xmax, nbins+1)
        self.__sumw = _np.zeros(nbins, dtype=float)
        self.__sumw2 = _np.zeros(nbins, dtype=float)
        self.__count = 0

    @property
    def binning( self ):
        """Tuple of (nbins,xmin,xmax)"""
        return ( self.__nbins, self.__xmin, self.__xmax )

    @property
    def integral( self ):
        """Calculates total sum of weights filled within [xmin,xmax]."""
        return self.__sumw.sum()

    @property
    def count( self ):
        """Total number of fills (an integer) done within the bin range."""
        return self.__count

    def clone_empty( self, title = None ):
        """Create a new HistFiller1D object with the same binning. Unless a new
        title is supplied, the title of the resulting object will also be the
        same.
        """
        return HistFiller1D( self.binning,
                             title = self.__title if title is None else title )

    def clone( self, title = None ):
        """Create a new HistFiller1D object with the same binning and
        contents. Unless a new title is supplied, the title of the resulting
        object will also be the same.
        """
        h = self.clone_empty(title=title)
        h.__sumw = self.__sumw.copy()
        h.__sumw2 = self.__sumw2.copy()
        h.__count = self.__count
        return h

    def reset( self ):
        """Reset histogram contents, as if no data was never filled into it."""
        self.__sumw.fill(0.0)
        self.__sumw2.fill(0.0)
        self.__count = 0
        return self

    def fill(self, values, weights = None):
        """Fill histogram with values (scalar or arrays). Unit weights are
        assumed unless explicitly provided."""
        def as_1darray( v ):
            v = _np.asarray(v,dtype=float)
            ns = len(v.shape)
            assert ns<=1, "fill works only with 1-dimensional arrays"
            return v if ns==1 else _np.array(v,ndmin=1)
        x = as_1darray(values)
        mask = (x >= self.__xmin) & (x <= self.__xmax)
        x = x[mask]
        n = len(x)
        if n==0:
            return
        self.__count += n
        #For fixed bins, _np.bincount is a lot faster than _np.histogram:
        idxs = ((x - self.__xmin)*self.__invbw).astype(int)
        idxs = _np.clip(idxs, 0, self.__nbins - 1)
        w = None
        if weights is not None:
            w = as_1darray(weights)
            assert len(mask)==len(w), ("length of values and weight arrays"
                                       " are unequal")
            w = w[mask]
        sw = _np.bincount(idxs, weights=w, minlength=self.__nbins)
        self.__sumw += sw
        self.__sumw2 += ( sw if w is None
                          else _np.bincount(idxs, weights=w**2,
                                            minlength=self.__nbins) )

    def add_contents( self, other_hist ):
        """Add the contents of other histogram to this histogram. The two
        histograms should have compatible binnings. The title will be
        kept. Returns self.
        """
        assert isinstance(other_hist,HistFiller1D)
        if self.binning != other_hist.binning:
            raise NCBadInput('Incompatible binnings')
        self.__sumw += other_hist.__sumw
        self.__sumw2 += other_hist.__sumw2
        self.__count += other_hist.__count
        return self

    def bindata( self, json_compat = False ):
        """Return bindata as dictionary. This include binnings, contents, and
        errors-squared. If json_compat is True, the returned dictionary will
        contain lists rather than Numpy arrays.
        """
        def export_array( a ):
            return ( a.tolist()
                     if (json_compat and hasattr(a,'tolist'))
                     else a.copy() )
        return dict( xmin = self.__xmin,
                     xmax = self.__xmax,
                     nbins = self.__nbins,
                     content = export_array(self.__sumw),
                     errorsq = export_array(self.__sumw2) )

    def to_dict(self, json_compat = False ):
        """Serialise as dictionary. If json_compat is True, the returned
        dictionary will contain lists rather than numpy arrays.
        """
        return dict ( title = self.__title,
                      bindata = self.bindata(json_compat=json_compat),
                      count = self.__count )

    def to_hist1d(self):
        """Initialise and return a Hist1D object representing this one."""
        return Hist1D( self.to_dict() )
