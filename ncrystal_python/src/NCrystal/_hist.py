
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
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

from ._numpy import _np, _ensure_numpy, _np_linspace

class Hist1D:

    def __init__( self, data ):
        if data=='_no_init_':
            return
        _ensure_numpy()
        self.__stats = data.get('stats')
        self.__title = data['title']
        hb = data['bindata']
        self.__xmin = hb['xmin']
        self.__xmax = hb['xmax']
        self.__nbins = hb['nbins']
        self.__y = _np.asarray(hb['content'],dtype=float)
        esq = hb.get('errorsq')
        self.__yerrsq = ( _np.asarray(esq,dtype=float)
                          if esq is not None
                          else None )
        self.__yerr = None
        assert self.__nbins == len(self.__y)
        assert self.__yerrsq is None or self.__nbins == len(self.__yerrsq)

    def clone( self, rebin_factor = 1 ):
        c = Hist1D('_no_init_')
        c.__stats = self.__stats
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
        """
        Returns integrated contents of the histogram over the area
        [xlow,xhigh] along with the error of that value in a tuple
        (content,error).

        This is done translating xlow and xhigh to exact bin edges and then
        calling integrate_bins. If that is not possible within the
        tolerance, an exception is raised.

        """
        if not ( xhigh >= xlow ):
            from .exceptions import NCBadInput
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
                from .exceptions import NCBadInput
                raise NCBadInput(f'Value {x} does not correspond exactly'
                                 ' to a bin edge within the tolerance.')
            return ir
        e_low = _findedge(xlow)
        e_high = _findedge(xhigh)
        if e_low == e_high:
            return ( 0.0, 0.0 )
        assert e_low >= 0
        assert e_low < e_high <= self.nbins
        return self.integrate_bins( e_low, e_high - 1 )

    def integrate_bins( self, bin_low = None, bin_up = None ):
        """
        Returns integrated contents of the bins [bin_low,bin_up[ along with
        the error of that value in a tuple (content,error).

        If bin_low is None the integration will start at the first bin and
        include the underflow bin.

        If bin_up is None the integration will end at the last bin and
        include the overflow bin.
        """

        add_overflow, add_underflow = False, False
        if bin_low is None:
            add_underflow = True
            bin_low = 0
            underflow_c = self.stats.get('underflow')
            underflow_e2 = self.stats.get('underflow_errorsq')
            if bool(underflow_c is None) != bool(underflow_e2 is None):
                from .exceptions import NCBadInput
                raise NCBadInput('Inconsistent underflow info')
            if underflow_c is None:
                add_underflow = False

        if bin_up is None:
            add_overflow = True
            bin_up = self.__nbins
            overflow_c = self.stats.get('overflow')
            overflow_e2 = self.stats.get('overflow_errorsq')
            if bool(overflow_c is None) != bool(overflow_e2 is None):
                from .exceptions import NCBadInput
                raise NCBadInput('Inconsistent overflow info')
            if overflow_c is None:
                add_overflow = False

        bin_low, bin_up = int(bin_low), int(bin_up)
        if bin_up < bin_low or bin_low<0 or bin_up > self.__nbins:
            from .exceptions import NCBadInput
            raise NCBadInput('Invalid bin range requested')
        content_integral = self.__y[bin_low:bin_up].sum()
        if add_underflow:
            content_integral += underflow_c
        if add_overflow:
            content_integral += overflow_c
        if self.__yerrsq is None:
            #unweighted, just base erros on contents:
            return ( content_integral, _np.sqrt(content_integral) )
        errorsq_integral = self.__yerrsq[bin_low:bin_up].sum()
        if add_underflow:
            errorsq_integral += underflow_e2
        if add_overflow:
            errorsq_integral += overflow_e2
        return ( content_integral, _np.sqrt(errorsq_integral) )

    def add_contents( self, other_hist ):
        o = other_hist
        assert self.__xmin == o.__xmin
        assert self.__xmax == o.__xmax
        assert self.__nbins == o.__nbins
        self.__stats = {}
        self.__title = '<edited>'
        self.__y += o.__y
        self.__yerr = None
        if self.__yerrsq is None:
            if o.__yerrsq is None:
                pass#done
            else:
                self.__yerrsq = self.__y + o.__yerrsq
        else:
            if o.__yerrsq is None:
                self.__yerrsq += o.__y
            else:
                self.__yerrsq += o.__yerrsq

    def rebin( self, rebin_factor ):
        assert self.__nbins % rebin_factor == 0
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
        return self.__stats or {}

    @property
    def errors( self ):
        if self.__yerr is None:
            self.__yerr = _np.sqrt( self.errors_squared )
        return self.__yerr

    @property
    def errors_squared( self ):
        return ( self.__yerrsq
                 if self.__yerrsq is not None
                 else self.__y )

    @property
    def content( self ):
        return self.__y

    @property
    def title( self ):
        return self.__title

    @property
    def xmin( self ):
        return self.__xmin

    @property
    def xmax( self ):
        return self.__xmax

    @property
    def binwidth( self ):
        return (self.__xmax-self.__xmin)/self.__nbins

    @property
    def nbins( self ):
        return self.__nbins

    @property
    def bincenters( self ):
        halfbw = 0.5*self.binwidth
        return _np_linspace(self.__xmin+halfbw,
                            self.__xmax-halfbw,
                            self.__nbins)

    @property
    def binedges( self ):
        return _np_linspace(self.__xmin,
                            self.__xmax,
                            self.__nbins+1)


    def _hist_curve( self, error_offset = 0.0 ):
        be = self.binedges
        y = self.content
        if error_offset:
            y += self.errors * error_offset
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

    def bar_args( self, style = True, **kwargs ):
        d = {'x' : self.binedges[:-1],
             'height': self.content,
             'width': self.binwidth,
             'align':'edge'
             }
        d.update(kwargs)
        return d

    def plot_hist( self, plt=None, axis=None, style=True, label=None,
                   show_errors=True, do_show = True, set_xlim = True ):
        if not plt and not axis:
            from .plot import _import_matplotlib_plt
            plt = _import_matplotlib_plt()
        if not axis:
            axis = plt.gca()
        axis.bar(**self.bar_args(label=label))
        if show_errors:
            ( error_markers,
              ecaplines,
              ebarlinecols ) = axis.errorbar(**self.errorbar_args())

        xmin,xmax,binwidth = self.xmin, self.xmax, self.binwidth
        if set_xlim:
            axis.set_xlim(xmin-1e-6*binwidth,xmax+1e-6*binwidth)
        if do_show and plt:
            plt.show()
