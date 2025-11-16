
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

# In this internal file we provide a histogram class which can be used to load
# the JSON data representing the (likewise internal) histograms from
# "NCrystal/internal/utils/NCHists.hh".

__all__ = ['Hist1D']

from ._numpy import _np, _ensure_numpy, _np_linspace, _np_trapezoid
from .exceptions import NCBadInput, NCCalcError
import math

class Hist1D:

    def __init__( self, data ):
        """Initialise histogram. This is intended to be using data provided by
        the C++ histogram's toJSON() method, or from the .to_dict() or
        .to_json() of this Python Hist1D class itself
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

    def to_dict( self, json_compat = False ):
        """Serialise as dictionary. The returned string can be used to
        initialise a Hist1D object again. If json_compat is True, the returned
        dictionary will contain lists rather than numpy arrays.
        """
        return dict( stats = self.stats,
                     title = self.__title,
                     bindata = self.bindata( json_compat = json_compat ) )

    def to_json( self ):
        """Serialise as JSON string. The returned string can be used to
        initialise a Hist1D object again."""
        import json
        return json.dumps( self.to_dict( json_compat = True ) )

    def clone( self, rebin_factor = 1 ):
        """Clone object. This is useful to keep the original histogram intact in
        case of subsequent operations on the histogram is desired, such as
        rebinning."""
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
            c.rebin( rebin_factor )
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

    def _flow_signature( self ):
        return ( self.__flow_under is None,
                 self.__flow_under_errorsq is None,
                 self.__flow_over is None,
                 self.__flow_over_errorsq is None )

    def add_contents( self, other_hist ):
        """Add the contents of other histogram to this histogram. The two
        histograms should have compatible binnings and over/underflow
        settings. Statistics will be updated if available on both histograms,
        otherwise dropped. The title will be kept. Returns self.

        """
        o = other_hist
        if self.binning != o.binning:
            raise NCBadInput("incompatible binning (%i,%g,%g) vs. (%i,%g,%g)."
                             %( *self.binning, *o.binning))

        if self._flow_signature() != o._flow_signature():
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

    def rebin( self, rebin_factor ):
        """Reduce granularity of binnings in histogram by provided factor, which
        must be a divisor in the current number of bins. Returns self."""
        assert self.__nbins % rebin_factor == 0
        if rebin_factor == 1:
            return self
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
        return self

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
            return a.tolist() if (json_compat and hasattr(a,'tolist')) else a

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
    def content( self ):
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
        """Returns an overview of the histogram in printable string form,
        optionally prefixing each line with provided prefix. If do_print is
        True, this will also be printed directly. By default, the contents of
        the histogram will only be included in case of 50 or fewer bins. Setting
        contents to True or False can be used to explicitly override this.  If
        highres is True, all floating point numbers will be printed with full
        precision.
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
            c,e = self.content, self.errors
            for ibin in range(self.nbins):
                pr("%s  content[ibin=%i] : %s +- %s"%(p,ibin,
                                                      fmt(c[ibin]),
                                                      fmt(e[ibin])))

        result.append('')
        result = '\n'.join(result)
        if do_print:
            from ._common import print as ncprint
            ncprint(result,end='')
        return result

    def _hist_curve( self ):
        be = self.binedges
        y = self.content
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
             'y':self.content,'xerr':0.5*self.binwidth,
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
             'height': self.content,
             'width': self.binwidth,
             'align':'edge'
             }
        d.update(kwargs)
        return d

    def plot_hist( self, *args, **kwargs ):
        """Alias for the .plot() method."""
        return self.plot_hist(*args,**kwargs)

    def plot( self, plt=None, axis=None, label=None,
                   show_errors=True, do_show = True, set_xlim = True ):
        """Produce a matplotlib plot of the histogram. If plt is None,
        matplotlib.pyplot is used. If axis is None, plt.gca() is used. Unless
        do_show is False, plt.show() wil be called ultimately. The show_errors,
        label, and set_xlim are all boolean flags whose affects are hopefully
        self explanatory.
        """
        if not plt and not axis:
            from .plot import _import_matplotlib_plt
            plt = _import_matplotlib_plt()
        if not axis:
            axis = plt.gca()
        axis.bar(**self.bar_args(label=label))
        if show_errors:
            axis.errorbar(**self.errorbar_args())
        xmin,xmax,binwidth = self.xmin, self.xmax, self.binwidth
        if set_xlim:
            axis.set_xlim(xmin-1e-6*binwidth,xmax+1e-6*binwidth)
        if do_show and plt:
            plt.show()

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

    def check_compat( self, other_hist,
                      p_value_threshold = 0.05,
                      force_norm = True,
                      return_pval = False,
                      check=False ):
        """Finds the p-value based on a chi-square comparison between the
        contents of this and another histogram, which must have identical
        binning. If force_norm is False, an error will be produced if histograms
        do not have identical statistics (integral). A high p-value indicates
        that the histogram contents are likely to come from the same underlying
        distribution. By default, the calculated p-value is compared against the
        p_value_threshold and True is returned if it is higher, False if it is
        lower. If return_pval is True, the p-value is instead returned. If check
        is True, a CalcError exception is raised in case the p-value fails to
        exceed p_value_threshold.

        Note, if you are running N histogram comparisons (e.g. in a test suite),
        and you want to have a p=5% chance of false positive in case you
        e.g. change the random seed of the test suite, you should use a
        p_value_threshold of roughly p/N=0.05/N (this is called Bonferroni
        correction and was confirmed to work well with simple Monte Carlo test).
        """
        def chisq_cdf( x, k ):
            #Could have used scipy.stats.chi2.cdf if we used scipy.
            if x<=0:
                return 0.0
            xx = _np_linspace(0.0,x,10000)
            c = (1.0 / (2.0 ** (0.5*k) * math.gamma(0.5*k)))
            yy = c * ( (xx ** (0.5*k - 1.0)) * _np.exp(- 0.5 * xx) )
            return _np_trapezoid( yy, xx )

        chisq, k = self.chisquare_dist( other_hist,
                                        force_norm = force_norm )
        p_value = 1.0 - chisq_cdf(chisq, k)
        assert -1e-6 <= p_value <= 1.0+1e-6
        p_value = min( 1.0, max( 0.0, p_value ) )
        ok = ( p_value >= p_value_threshold )
        if check and not ok:
            #fmt for unit test reproducibility
            fmt='p-value=%g'%p_value if p_value>=1e-6 else 'p-value'
            raise NCCalcError(f'check_compat failed: {fmt}'
                             f' is not greater than {p_value_threshold}.')
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
            raise NCBadInput('chisquare_dist: incompatible binnings.')
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

        c1,e1 = s.content, s.errors
        c2,e2 = o.content, o.errors

        #use mask to avoid zero division warnings in bins with no content in
        #either dataset:
        mask = (e1>0.0) | (e2>0.0)
        c1, e1, c2, e2 = c1[mask], e1[mask], c2[mask], e2[mask]

        if len(c1) > 0:
            c1sum = c1.sum()+(s.underflow or 0.0)+(s.overflow or 0.0)
            c2sum = c2.sum()+(o.underflow or 0.0)+(o.overflow or 0.0)
            assert (c1sum-c2sum)<1e-5*(c1sum+c2sum),"must be normalised"
            chi_squared = ( (c1-c2)**2 / (e1**2 + e2**2) ).sum()
        else:
            chi_squared = 0.0
        n = len(c1)
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
