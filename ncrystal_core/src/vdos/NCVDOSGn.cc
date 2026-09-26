////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <deque>

namespace NC=NCrystal;
namespace NCV=NCrystal::VDOS;

namespace NCRYSTAL_NAMESPACE {

  namespace VDOS {

    namespace {
      static std::atomic<bool> s_verbose_vdosgn( ncgetenv_bool("DEBUG_PHONON") );

      struct CfgDecoded {
        int minThinOrder = 4;//Below this order, no thinning takes place
        //(0=always thin, -1=never thin)
        unsigned thinNBins = 1000;//double binwidth whenever number of bins
        //exceeds this value (0 disables)
        int minThinAgressiveOrder = 20;//same but for more agressive thinning
        unsigned thinAgressiveNBins = 300;//same but for more agressive thinning
        int minTruncOrder = 0;//Below this order, no truncation takes place
        //(0=always, -1=never)
        double truncationThreshold = 1e-13;//trim ranges to remove negligible
        //noise at edges (0 disables)
        bool legacyConvolve = false;
        bool directConvolve = false;//use FastConvolve::convolveDirect
        //instead of convolve, for every order (see VDOSGn::Cfg::MaxLux).
      };
      CfgDecoded decodeCfg( VDOSGn::Cfg choice ) {
        CfgDecoded res;
        if (choice == VDOSGn::Cfg::Legacy) {
          //NCrystal 2.x-4.x behaviour:
          res.minThinOrder = res.minTruncOrder = 5;
          res.thinNBins = res.thinAgressiveNBins = 1000;
          res.minThinAgressiveOrder = 50000;
          res.truncationThreshold = 1e-14;
          res.legacyConvolve = true;
        } else if (choice == VDOSGn::Cfg::MaxLux) {
          res.directConvolve = true;
        } else {
          nc_assert( choice == VDOSGn::Cfg::Default );
        }
        return res;
      }

      class VDOSGnData : private MoveOnly {
      public:
        //The egrid_lower value is the energy of the first point in the
        //spectrum, and startIdx is the integer number such that egrid_lower =
        //startIdx*egrid_binwidth (only used in the non-legacy mode where this
        //is how we ensure that the energy grids of all orders are on lattices
        //which are anchored at energy=0, independently of rounding errors and
        //where truncation happens to cut).
        //Fixme: Safer with std::int64_t for startIdx?
        VDOSGnData( const VectD &spec,
                    double egrid_lower,
                    double egrid_binwidth,
                    unsigned long thinFactor,
                    long startIdx = 0 );
        long getStartIdx() const { return m_startIdx; }
        double interpolateDensity(double energy) const;
        void interpolateDensityMany(Span<const double>, VectD&, VectD&) const;

        const VectD& getSpectrum() const { return m_spec; }
        double getEGridLower() const {return m_egrid_lower;}
        double getEGridUpper() const {return m_egrid_upper;}
        double getEGridBinwidth() const {return m_egrid_binwidth;}
        double maxDensity() const { return m_specMaxVal; }
        unsigned long getThinFactor() const { return m_thinFactor; }
        VDOSGnData( VDOSGnData&& ) = default;
        VDOSGnData& operator=( VDOSGnData&& ) = default;
      private:
        VectD m_spec;
        std::size_t m_spec_size_minus_2;
        double m_egrid_lower, m_egrid_upper;
        double m_egrid_binwidth, m_egrid_invbinwidth, m_specMaxVal;
        unsigned long m_thinFactor;//NB: binwidth is G1's binwidth multiplied
                                   //    by m_thinFactor
        long m_startIdx;
      };
    }
  }
}

NCV::VDOSGnData::VDOSGnData( const VectD &spec,
                             double egrid_lower,
                             double egrid_binwidth,
                             unsigned long thinFactor,
                             long startIdx )
  : m_spec(spec.begin(),spec.end()),
    m_thinFactor(thinFactor),
    m_startIdx(startIdx)
{
  m_egrid_lower = egrid_lower;
  m_egrid_binwidth = egrid_binwidth;
  nc_assert(m_egrid_binwidth>0.0);
  m_egrid_invbinwidth = 1.0/egrid_binwidth;
  nc_assert(m_spec.size()>3);
  m_spec_size_minus_2 = m_spec.size() - 2;
  m_egrid_upper = equidistantGridPoint( m_egrid_lower, m_egrid_binwidth,
                                        m_spec.size()-1 );
  nc_assert_always(!m_spec.empty());
  //The assymetric Gn functions are constructed/required to have unit area, so
  //must normalise.  NB: In principle higher-order spectra are guaranteed to be
  //normalised when lower order spectra are. However, this re-normalisation acts
  //as a safeguard in the face of thinning, truncation and numerical errors in
  //general.
  double asym_area = 0.;
  for(unsigned i=0;i<m_spec.size();i++)
    asym_area +=  m_spec[i];
  asym_area *= egrid_binwidth;
  double inv_asym_area(1.0/asym_area);
  for(unsigned i=0;i<m_spec.size();i++)
    m_spec[i] *= inv_asym_area;
  m_specMaxVal = *std::max_element(m_spec.begin(),m_spec.end());
}

double NCV::VDOSGnData::interpolateDensity(double energy) const
{
  if (!valueInInterval(m_egrid_lower,m_egrid_upper,energy))
    return 0.0;
  double a = (energy-m_egrid_lower)*m_egrid_invbinwidth;
  double floor_a = std::floor(a);
  //clamp to safe-guard against numerical errors.:
  std::size_t index = std::min<std::size_t>(m_spec_size_minus_2,
                                            static_cast<std::size_t>(floor_a));
  double f = a - floor_a;//a-index instead would mix int and double => slower.
  nc_assert( index+1 < m_spec.size() );
  const double * valptr = &m_spec[index];
  return nclerp( *valptr, *(valptr+1), f );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Batched form of interpolateDensity's nclerp(spec[ix],spec[ix+1],f)
    //call above; must match it bit for bit (see the debug assertion in
    //interpolateDensityMany below). See doc/devel_fma_attribute.md for
    //the audit this requires:
    NCRYSTAL_FMADISPATCH_ATTR
    void vdosGnInterpolateDensityRun( double* out, const double* f,
                                      const double* ix_as_dbl,//indices, as double
                                      const double* spec, std::size_t count )
    {
      for ( std::size_t k = 0; k < count; ++k ) {
        const std::size_t ix = static_cast<std::size_t>(ix_as_dbl[k]);
        out[k] = nclerp( spec[ix], spec[ix+1], f[k] );
      }
    }
  }
}

void NCV::VDOSGnData::interpolateDensityMany( Span<const double> energy,
                                              VectD& out,
                                              VectD& workbuf) const
{
#ifndef NDEBUG
  nc_assert(energy.size() <
            static_cast<std::size_t>(1000000000));

  const std::size_t npts = m_spec.size();

  nc_assert(npts >= 2);
  nc_assert(m_spec_size_minus_2 == npts - 2);
  nc_assert(std::isfinite(m_egrid_lower));
  nc_assert(std::isfinite(m_egrid_upper));
  nc_assert(std::isfinite(m_egrid_binwidth));
  nc_assert(m_egrid_binwidth > 0.0);
  nc_assert(std::isfinite(m_egrid_invbinwidth));
  nc_assert(m_egrid_upper >= m_egrid_lower);

  nc_assert( nc_is_grid( energy ) );
  for (std::size_t i = 0; i < energy.size(); ++i) {
    nc_assert(std::isfinite(energy[i]));
    if (i != 0)
      nc_assert(energy[i - 1] < energy[i]);
  }

  for (std::size_t i = 0; i < npts; ++i)
    nc_assert(std::isfinite(vectAt(m_spec, i)));

  const auto no_overlap = [](const double* a, std::size_t an,
                             const double* b, std::size_t bn) {
    const std::uintptr_t ab =
      reinterpret_cast<std::uintptr_t>(a);
    const std::uintptr_t ae = ab + an * sizeof(double);
    const std::uintptr_t bb =
      reinterpret_cast<std::uintptr_t>(b);
    const std::uintptr_t be = bb + bn * sizeof(double);

    return ae <= bb || be <= ab;
  };

  nc_assert(no_overlap(
                       energy.data(), energy.size(), out.data(), out.size()));
  nc_assert(no_overlap(
                       energy.data(), energy.size(), workbuf.data(), workbuf.size()));
  nc_assert(no_overlap(
                       out.data(), out.size(), workbuf.data(), workbuf.size()));
#endif

  const std::size_t n = energy.size();

  out.assign(n, 0.0);
  workbuf.resize(2 * n);

  const double* ncrestrict ep = energy.data();
  double* ncrestrict op = out.data();
  double* ncrestrict buf_f = workbuf.data();
  double* ncrestrict buf_ix = buf_f + n;

  const std::size_t beg = fastLowerBoundIdx( energy.data(), n, m_egrid_lower );
  const std::size_t end = beg + fastUpperBoundIdx( energy.data()+beg, n-beg,
                                                    m_egrid_upper );

  for (std::size_t i = beg; i < end; ++i) {
    const double a =
      (ep[i] - m_egrid_lower) * m_egrid_invbinwidth;
    const double fa = std::floor(a);
    const std::size_t ix = ncmin(
                                 m_spec_size_minus_2, static_cast<std::size_t>(fa));

    buf_f[i] = a - fa;
    buf_ix[i] = static_cast<double>(ix);
  }

#ifndef NDEBUG
  //vdosGnInterpolateDensityRun below indexes spec[ix]/spec[ix+1] for every
  //ix in buf_ix[beg,end) without any bounds check of its own (it only gets
  //a count, not the spectrum's size) -- verify here, from the caller's
  //side, that every value it will read is exactly what it is assumed to
  //be: an index in [0,m_spec_size_minus_2], leaving spec[ix+1] safely
  //within m_spec (size m_spec_size_minus_2+2):
  for (std::size_t i = beg; i < end; ++i) {
    nc_assert( buf_ix[i] >= 0.0 );
    nc_assert( buf_ix[i] <= static_cast<double>(m_spec_size_minus_2) );
  }
#endif

  vdosGnInterpolateDensityRun( op + beg, buf_f + beg, buf_ix + beg,
                               m_spec.data(), end - beg );

#ifndef NDEBUG
  for (std::size_t i = 0; i < energy.size(); ++i)
    nc_assert(op[i] == interpolateDensity(energy[i]));
#endif
}

struct NCV::VDOSGn::Impl {
  Impl(const VDOSEval& vde, const VDOSGn::Cfg& );
  std::deque<VDOSGnData> m_gndata;//deque, to not change address of VDOSGnData
                                  //objects when growMaxOrder() is called.
  std::vector<VDOSGnData> m_mt_pending_gndata;//if in MT mode, we might
                                              //precalculate next orders
                                              //concurrently and same them here.
  Optional<FactoryJobs> m_mt_jobs;
  static_assert( std::is_nothrow_default_constructible<Optional<VDOSGnData>>::value, "");
  SmallVector<Optional<VDOSGnData>,10> m_mt_buffer;
  CfgDecoded m_cfg;
  SmallVector<FastConvolve,4> m_fastConvolve;
  int m_nmaxconcurrent = 0;

  void produceNewOrderByConvolution(Order);
  VDOSGnData produceNewOrderByConvolutionImpl( Order, FastConvolve& ) const;
  VDOSGnData& accessAtOrder(Order n)
  {
    nc_assert(n.value()<=m_gndata.size());
    return m_gndata[n.value()-1];
  }
  const VDOSGnData& accessAtOrder(Order n) const
  {
    nc_assert(n.value()<=m_gndata.size());
    return m_gndata[n.value()-1];
  }

};

NCV::VDOSGn::Impl::Impl(const VDOSEval& vde,
                        const VDOSGn::Cfg& cfg )
  : m_cfg(decodeCfg(cfg)),
    m_nmaxconcurrent(ncgetenv_int("VDOSGN_CONCURRENT",4))
{
  auto gridinfo = vde.getGridInfo();
  nc_assert(gridinfo.npts>1);

  //egrid starting from 0.0:
  unsigned long nbins = gridinfo.npts_extended - 1;

  //Thicken if too few bins for robust numerical integration (not really tuned,
  //but seems sensible to increase very low number of pts a bit). We apply the
  //factor to nbins, not npts, since we want e.g. thicken_factor=2 to correspond
  //to the placement of 1 extra point in the middle of all existing bins.:
  constexpr unsigned long min_nbins = 400;
  const unsigned long thicken_factor = static_cast<unsigned long>(std::ceil(double(min_nbins)/nbins));

  if ( s_verbose_vdosgn && thicken_factor != 1 )
    NCRYSTAL_MSG("VDOSGn Thickening provided VDOS egrid for G1 by a"
                 " factor of "<<thicken_factor<<" resulting in number of grid"
                 " points for [-emax,emax] increasing "<<nbins*2+1
                 <<" -> "<<nbins*thicken_factor*2+1);

  nbins *= thicken_factor;
  nc_assert_always( nbins < 10000000);

  auto egrid = linspace(0.0,gridinfo.emax,nbins+1);
  const double binwidth = egrid.back() / nbins;

  //Initialise G1 array on the egrid, from -emax to +emax:
  VectD G1spectrum(egrid.size()*2-1,0.0);

  const double gamma0 = vde.calcGamma0();

  if ( s_verbose_vdosgn )
    NCRYSTAL_MSG("VDOSGn ctor: entering G1-fill loop (nbins="<<nbins
                 <<", G1spectrum.size()="<<G1spectrum.size()<<")");

  for (auto e: enumerate(egrid) ) {
    nc_assert(e.val>=0.0);
    auto g1_vals = vde.evalG1AsymmetricAtEPair( e.val, gamma0 );
    //Fill at +e.val:
    vectAt(G1spectrum,nbins+e.idx) = g1_vals.second;
    //Fill at -e.val:
    vectAt(G1spectrum,nbins-e.idx) = g1_vals.first;
  }

  if ( s_verbose_vdosgn )
    NCRYSTAL_MSG("VDOSGn ctor: G1-fill loop done");

  nc_assert_always( valueInInterval(0.0,0.1,m_cfg.truncationThreshold) );
  nc_assert_always( m_cfg.minThinOrder >= -1 );
  nc_assert_always( m_cfg.minThinAgressiveOrder >= -1 );
  nc_assert_always( m_cfg.minTruncOrder >= -1 );

  //Discard excess zeroes at edges for G1, keeping at most a single entry with 0
  //at each edge. This might in particular happen at very low energies where the
  //detailed balance factor can become numerically identical to zero for even
  //relatively small positive arguments:
  double actual_edgelower = -gridinfo.emax;
  auto itB = G1spectrum.begin();
  auto itE = G1spectrum.end();
  auto itFirst = itB;
  auto itLast = std::prev(itE);
  while ( itFirst != itLast && !(*itFirst>0.0) && !(*std::next(itFirst)>0.0) )
    ++itFirst;
  while ( itLast != itB && !(*itLast>0.0) && !(*std::prev(itLast)>0.0) )
    --itLast;
  if ( itFirst >= itLast || std::distance(itFirst,itLast) < 3 )
    NCRYSTAL_THROW(CalcError,"Too few non-zero pts in G1 spectrum.");
  //Index of first point, in units of binwidth (point 'nbins' is at energy 0):
  nc_assert( std::distance( itB, itFirst ) >= 0 );
  const auto distBF = std::distance( itB, itFirst );
#ifndef NDEBUG
  nc_assert( static_cast<std::size_t>(nbins)+1
             < static_cast<std::size_t>(std::numeric_limits<long>::max())/4 );
  nc_assert( static_cast<std::size_t>(distBF)+1
             < static_cast<std::size_t>(std::numeric_limits<long>::max())/4 );
#endif
  const long g1StartIdx
    = ( static_cast<long>( distBF ) - static_cast<long>( nbins ) );
  if ( itFirst != itB || itLast != std::prev(itE) ) {
    actual_edgelower += std::distance( itB, itFirst ) * binwidth;
    VectD( itFirst, std::next(itLast) ).swap( G1spectrum );
  }

  if ( s_verbose_vdosgn )
    NCRYSTAL_MSG("VDOSGn ctor: edge-trim done (G1spectrum.size()="
                 <<G1spectrum.size()<<", g1StartIdx="<<g1StartIdx<<")");

  //Place G1:
  if ( m_cfg.legacyConvolve ) {
    m_gndata.emplace_back( G1spectrum, actual_edgelower, binwidth, 1 );
  } else {
    m_gndata.emplace_back( G1spectrum,
                           static_cast<double>(g1StartIdx) * binwidth,
                           binwidth, 1, g1StartIdx );
  }

  if (s_verbose_vdosgn)
    NCRYSTAL_MSG("VDOSGn constructed (input spectrum size: "<<G1spectrum.size()
                 <<", thinning with minOrder="<<m_cfg.minThinOrder
                 <<" and thinNBins="<<m_cfg.thinNBins
                 <<" and more agressive thinning with minOrder="
                 <<m_cfg.minThinAgressiveOrder
                 <<" and thinNBins="<<m_cfg.thinAgressiveNBins
                 <<", truncation with minOrder="<<m_cfg.minTruncOrder
                 <<" and truncationThreshold="<<m_cfg.truncationThreshold
                 <<(m_cfg.legacyConvolve?", mode=legacy":"")
                 <<(m_cfg.directConvolve?", mode=directConvolve":"")
                 <<")");
}

NCV::VDOSGn::~VDOSGn() {
  //A moved-from VDOSGn (e.g. the husk left behind when GnExpansion's
  //implicit move constructor is actually invoked, rather than elided via
  //NRVO -- confirmed to happen on MSVC for expandVDOSToGnFcts's "return
  //res;" where GCC/Clang instead reliably apply NRVO here) has a null
  //m_impl and must not be dereferenced below:
  if ( !m_impl )
    return;
  if ( m_impl->m_mt_jobs.has_value() ) {
    //End running jobs, so they don't write to suddenly non-existent buffers:
    m_impl->m_mt_jobs.value().waitAll();
  }
  if (s_verbose_vdosgn)
    NCRYSTAL_MSG("VDOSGn destructed (final max order: "
                 <<maxOrder().value()<<")")

      }

NCV::VDOSGn::VDOSGn( const VDOSEval& vde, Cfg cfg )
  : m_impl(vde,cfg),
    m_kT(vde.kT())
{
}

NCV::VDOSGn::Order NCV::VDOSGn::maxOrder() const
{
  return static_cast<unsigned>( m_impl->m_gndata.size() );
}

void NCV::VDOSGn::growMaxOrder( Order target_n )
{
  Order n = maxOrder();
  ++n;
  for ( ; n <= target_n; ++n )
    m_impl->produceNewOrderByConvolution(n);

  nc_assert( maxOrder().value() == target_n.value() );
}

double NCV::VDOSGn::eval( Order n, double energy ) const
{
  return m_impl->accessAtOrder(n).interpolateDensity(energy);
}

void NCV::VDOSGn::evalMany( Order n, Span<const double> egrid,
                            VectD& out, VectD& workbuf ) const
{
  m_impl->accessAtOrder(n).interpolateDensityMany(egrid,out,workbuf);
}

const NC::VectD& NCV::VDOSGn::getRawSpectrum( Order n ) const
{
  return m_impl->accessAtOrder(n).getSpectrum();
}

double NCV::VDOSGn::binWidth( Order n) const
{
  return m_impl->accessAtOrder(n).getEGridBinwidth();
}

NC::PairDD NCV::VDOSGn::eRange( Order n, double relthreshold ) const
{
  nc_assert(relthreshold>0.0&&relthreshold<1.0);
  const auto& p = m_impl->accessAtOrder(n);
  return estimateGnErange( p.getEGridLower(), p.getEGridBinwidth(),
                           p.getSpectrum(), relthreshold );
}

NC::PairDD NCV::VDOSGn::eRange( Order n ) const
{
  const auto& p = m_impl->accessAtOrder(n);
  return { p.getEGridLower(), p.getEGridUpper() };
}

void NCV::VDOSGn::enableVerboseOutput(bool status)
{
  s_verbose_vdosgn = status;
}

bool NCV::VDOSGn::verboseOutputEnabled()
{
  return s_verbose_vdosgn;
}

void NCV::VDOSGn::Impl::produceNewOrderByConvolution( Order order )
{
  const unsigned current_maxorder = static_cast<unsigned>( m_gndata.size() );
  nc_assert_always( order.value() == current_maxorder + 1 );

  if ( m_mt_jobs.has_value() ) {
    m_mt_jobs.value().waitAll();
    //Transfer concurrently generated results:
    for ( auto i : ncrange( m_mt_buffer.size() ) )
      m_mt_pending_gndata
        .emplace_back( std::move( m_mt_buffer.at( m_mt_buffer.size()-1-i ).value() ) );
    m_mt_buffer.clear();
    m_mt_jobs.reset();
  }

  if ( !m_mt_pending_gndata.empty() ) {
    //Easy, we already calculated that order previously (hopefully taking
    //advantage of multithreading):
    m_gndata.emplace_back( std::move( m_mt_pending_gndata.back() ) );
    m_mt_pending_gndata.pop_back();
    return;
  }

  //We have to do actual work. In case we want to take advantage of
  //multi-threading, we need to go ahead and produce more than just the
  //requested order. For that, we note that if we already have order N, then at
  //most we can produce up to order 2N ("G_2n = G_n(x)G_n"), so we can never
  //produce more new orders from existing data than we already have:
  const int nconcurrent = std::min<int>((int)m_gndata.size(),m_nmaxconcurrent);

  //We always need at least one convolver:
  if ( m_fastConvolve.empty() )
    m_fastConvolve.emplace_back();

  if ( nconcurrent<=1 ) {
    //Produce one without concurrency:
    m_gndata.emplace_back( produceNewOrderByConvolutionImpl( order, m_fastConvolve.back() ) );
    return;
  }

  //We want concurrent production, if there is a thread-pool available:
  m_mt_jobs.emplace();
  if ( !m_mt_jobs.value().isMT() ) {
    //Abort, we don't actually have a thread pool anyway:
    m_mt_jobs.reset();
    m_gndata.emplace_back( produceNewOrderByConvolutionImpl( order, m_fastConvolve.back() ) );
    return;
  }

  while ( m_fastConvolve.size() < (std::size_t)(nconcurrent) )
    m_fastConvolve.emplace_back();

  //Do nconcurrent-1 in the thread pool, and 1 here in the current thread (nb:
  //we do not invoke m_mt_jobs.waitAll() yet).
  m_mt_buffer.resize(nconcurrent-1);
  unsigned max_new_order = current_maxorder + (unsigned)nconcurrent;
  auto itFC = m_fastConvolve.begin();
  auto itBuffer = m_mt_buffer.begin();
  //NCRYSTAL_MSG("nconcurrent = "<<nconcurrent);
  for ( unsigned target_order = current_maxorder + 2;
        target_order <= max_new_order;
        ++target_order ) {
    FastConvolve *fcptr = &(*itFC++);
    Optional<VDOSGnData>* resbufptr = &(*itBuffer++);
    m_mt_jobs.value().queue( [fcptr,resbufptr,target_order,this]()
    {
      resbufptr->emplace( this->produceNewOrderByConvolutionImpl( Order{target_order}, *fcptr ) );
    });
  }
  m_gndata.emplace_back( this->produceNewOrderByConvolutionImpl( Order{current_maxorder+1}, *itFC ) );
}

NCV::VDOSGnData
NCV::VDOSGn::Impl::produceNewOrderByConvolutionImpl( Order order,
                                                     FastConvolve& fastConvolve ) const
{
  Order order2 = order.value()/2;
  Order order1 = order.value()-order2.value();

  const auto& p1 = accessAtOrder(order1);
  const auto& p2 = accessAtOrder(order2);

  //Function which can thin a vector (i.e. increase binwidth by merging bins),
  //used two places below:
  auto thinVector = [](unsigned thinFactor, const VectD& v)
  {
    nc_assert(thinFactor>1);
    VectD vt;
    auto newsize = ( v.size() + thinFactor - 1 ) / thinFactor;
    vt.reserve( newsize );
    for ( std::size_t i = 0; i<v.size(); i+= thinFactor )
      vt.push_back(vectAt(v,i));
    nc_assert_always( vt.size() == newsize  );
    return vt;
  };

  //In non-legacy mode, all grids are anchored at energy 0. I.e. the energy of
  //point j in a spectrum of G1 is (startIdx+j)*binwidth, with integer startIdx
  //(in units of the binwidth of that spectrum). Thinning then keeps the
  //points where startIdx+j is a multiple of the thinning factor, and
  //truncation just adds to startIdx. This makes the energy grid of every order
  //independent of rounding errors, and of where truncation happens to cut
  //(which could otherwise shift the phase of a thinned grid by a fraction of
  //a bin, whenever tiny noise in the spectrum changed the cut by one bin).
  const bool anchored = !m_cfg.legacyConvolve;
  auto thinVectorAnchored = [](unsigned long thinFactor, const VectD& v,
                               long& startIdx)
  {
    nc_assert(thinFactor>1);
    const long f = static_cast<long>(thinFactor);
    long r = startIdx % f;
    if ( r < 0 )
      r += f;
    const std::size_t skip = static_cast<std::size_t>( r ? f - r : 0 );
    nc_assert_always( v.size() > skip );
    VectD vt;
    vt.reserve( ( v.size() - skip + thinFactor - 1 ) / thinFactor );
    for ( std::size_t i = skip; i < v.size(); i += thinFactor )
      vt.push_back(vectAt(v,i));
    nc_assert( ( startIdx + static_cast<long>(skip) ) % f == 0 );
    startIdx = ( startIdx + static_cast<long>(skip) ) / f;
    return vt;
  };
  long startIdx1 = p1.getStartIdx();
  long startIdx2 = p2.getStartIdx();

  double dt;
  double dt1 = p1.getEGridBinwidth();
  double dt2 = p2.getEGridBinwidth();
  unsigned long thinFactor1 = p1.getThinFactor();
  unsigned long thinFactor2 = p2.getThinFactor();

  VectD vtmp;
  const VectD* input1_spec = &p1.getSpectrum();
  const VectD* input2_spec = &p2.getSpectrum();
  bool dt_mismatch = false;
  if ( thinFactor1 == thinFactor2 ) {
    nc_assert(floateq(dt1,dt2));
    //This happens most of the time:
    dt = dt1;
  } else {
    dt_mismatch = true;
    //But at certain transition points, we might need this on-demand thinning of
    //one spectrum:
    dt = std::max<double>(dt1,dt2);
    if ( thinFactor1 > thinFactor2 ) {
      //thin dt2
      nc_assert(thinFactor1%thinFactor2==0);
      unsigned long thinFactor = thinFactor1 / thinFactor2;
      nc_assert_always( floateq(dt,dt2*thinFactor) );
      vtmp = ( anchored
               ? thinVectorAnchored(thinFactor,*input2_spec,startIdx2)
               : thinVector(thinFactor,*input2_spec) );
      input2_spec = &vtmp;
      thinFactor2 *= thinFactor;
    } else {
      //thin dt1
      nc_assert(thinFactor2%thinFactor1==0);
      unsigned long thinFactor = thinFactor2 / thinFactor1;
      nc_assert_always( floateq(dt,dt1*thinFactor) );
      vtmp = ( anchored
               ? thinVectorAnchored(thinFactor,*input1_spec,startIdx1)
               : thinVector(thinFactor,*input1_spec) );
      input1_spec = &vtmp;
      thinFactor1 *= thinFactor;
    }
  }
  nc_assert_always(thinFactor1==thinFactor2);

  if ( s_verbose_vdosgn )
    NCRYSTAL_MSG("VDOSGn convolve: order="<<order.value()
                 <<" pre-thinning done (n1="<<input1_spec->size()
                 <<", n2="<<input2_spec->size()<<", dt_mismatch="
                 <<(dt_mismatch?"yes":"no")<<")");

  VectD phonon_spe;
  double start_energy = p1.getEGridLower() + p2.getEGridLower();
  long startIdx = startIdx1 + startIdx2;//(only used when anchored)
  if ( m_cfg.legacyConvolve )
    fastConvolve.convolveLegacy( *input1_spec, *input2_spec, phonon_spe, dt );
  else if ( m_cfg.directConvolve )
    fastConvolve.convolveDirect( *input1_spec, *input2_spec, phonon_spe, dt );
  else
    fastConvolve.convolve( *input1_spec, *input2_spec, phonon_spe, dt );
  auto orig_npts_result = phonon_spe.size();

  if ( s_verbose_vdosgn )
    NCRYSTAL_MSG("VDOSGn convolve: order="<<order.value()
                 <<" fastConvolve done (phonon_spe.size()="
                 <<orig_npts_result<<")");

  if ( m_cfg.minTruncOrder >= 0
       && m_cfg.truncationThreshold > 0.0
       && order.value() >= static_cast<unsigned>(m_cfg.minTruncOrder) ) {
    // => do truncation
    const double spec_max = *std::max_element(phonon_spe.begin(),phonon_spe.end());
    const double spec_cutoff = m_cfg.truncationThreshold * spec_max;
    //FastConvolve's own (size-dependent) FFT round-off noise floor can, for
    //large/high-order spectra, exceed the fixed spec_cutoff above (already
    //used below to gate whether the edge-crossing itself needs the more
    //careful taper treatment). Values kept above spec_cutoff but still below
    //this floor are not reliable signal -- they are FFT round-off, which
    //genuinely differs (not just last-ULP) across platforms/compilers, and
    //were observed doing so on real production data (Li2O, vdoslux=2004)
    //feeding into VDOS::determineAlphaBetaGridFromGn's grid-point selection.
    //Computed once here and reused both for that existing gate and for the
    //final value-cleanup below, so the same criterion decides both whether
    //the discrete edge needs refining and which retained values are trusted
    //enough to keep as-is rather than snapped to exactly 0.0. See
    //docs/claude_session_vdos_fma_reprod.md.
    constexpr double noiseFloorSafetyFactor = 8.0;
    const double noiseFloorGate = VDOS::estimateFFTConvolutionNoiseFloor(
      spec_max, phonon_spe.size(), noiseFloorSafetyFactor );
    const double cleanupThreshold = ncmax( spec_cutoff, noiseFloorGate );
    std::size_t ifront(0), iback(phonon_spe.size()-1);
    for (;ifront<iback;++ifront) {
      if (phonon_spe.at(ifront)>spec_cutoff)
        break;
    }
    for (;iback>ifront;--iback) {
      if (phonon_spe.at(iback)>spec_cutoff)
        break;
    }
    //Refine the two edges with a windowed quadratic-log-fit crossing
    //estimate (same tool as VDOSGn::eRange's estimateGnErange), then smoothly
    //taper the spectrum around that estimate (applyCrossingTaper) before
    //truncating a safe margin beyond it: this way, neither a single point
    //which only clears spec_cutoff by chance round-off, nor a genuine
    //near-tie in the continuous crossing estimate itself (confirmed on real
    //production data: the estimate can shift by several tenths of an index
    //unit between an -mfma and a plain build, near an array edge where the
    //fit window is asymmetric), can flip a discrete array-length/content
    //decision by more than a numerically negligible, already-tapered-near-
    //zero amount. Only attempted where FastConvolve's own (size-dependent)
    //noise floor is not safely below spec_cutoff -- for most orders (thinned
    //to a few hundred/thousand points) the plain crossing is already far
    //more reliable than that floor, and refining those too was found to
    //relocate near-ties into orders that never needed it, without reducing
    //the total count. Skipped for legacyConvolve (must reproduce NCrystal
    //2.x-4.x bit-for-bit) and directConvolve (VDOSGn::Cfg::MaxLux;
    //convolveDirect has no FFT noise floor for this estimate to describe).
    //See docs/claude_session_vdos_fma_reprod.md.
    if ( iback > ifront && !m_cfg.legacyConvolve && !m_cfg.directConvolve ) {
      constexpr std::size_t crossingNExtra = 6;
      constexpr double noiseFloorGateFactor = 2.0;
      if ( noiseFloorGate > noiseFloorGateFactor*spec_cutoff ) {
        //Taper half-width matches the crossing estimate's own fit window,
        //comfortably covering the largest xcross shift observed between an
        //-mfma and a plain build in production data (~0.75 index units, at
        //an array edge where the fit window is asymmetric/less averaged --
        //see applyCrossingTaper's doc comment):
        constexpr double taperHalfWidth = double(crossingNExtra);
        if ( ifront > 0 ) {
          const double x = VDOS::estimateSpectrumCrossing( 0.0, 1.0, phonon_spe,
                                                           ifront-1, ifront,
                                                           spec_cutoff, crossingNExtra );
          VDOS::applyCrossingTaper( phonon_spe, x, true, taperHalfWidth );
          //Truncate well past the taper's own zero-plateau (floor, not
          //round), so a residual last-ULP-level shift in x cannot flip
          //which integer gets chosen for something that still matters --
          //both plausible choices there already hold an equally negligible,
          //tapered value:
          const double cut = x - taperHalfWidth - 1.0;
          ifront = ( cut > 0.0 ? static_cast<std::size_t>(std::floor(cut)) : 0 );
          ifront = ncmin( ifront, iback );
        }
        if ( iback+1 < phonon_spe.size() ) {
          const double x = VDOS::estimateSpectrumCrossing( 0.0, 1.0, phonon_spe,
                                                           iback, iback+1,
                                                           spec_cutoff, crossingNExtra );
          VDOS::applyCrossingTaper( phonon_spe, x, false, taperHalfWidth );
          const double cut = x + taperHalfWidth + 1.0;
          iback = ( cut < double(phonon_spe.size()-1)
                    ? static_cast<std::size_t>(std::ceil(cut))
                    : phonon_spe.size()-1 );
          iback = ncmax( iback, ifront );
        }
      }
    }
    if (iback>ifront) {
      VectD truncated_spec(phonon_spe.begin()+ifront,phonon_spe.begin()+iback+1);
      truncated_spec.swap(phonon_spe);
    }
    //Remove non-cross-platform-reproducible noise from the FFT alg by snapping
    //tiny noise to 0.0 (also internally, not just at the edges). Uses
    //cleanupThreshold (not the bare spec_cutoff) so that values which only
    //cleared the fixed spec_cutoff by virtue of FastConvolve's own
    //size-dependent round-off floor (see above) are cleaned up too:
    //fixme: with non-legacy convolve we MUST do this, or we can get negative values in the spectra
    if ( !m_cfg.legacyConvolve ) {
      for ( auto&e : phonon_spe) {
        if ( e < cleanupThreshold ) {
          nc_assert( e > -1e-12*spec_max );
          e = 0.0;
        }
      }
    }
    start_energy = equidistantGridPoint( start_energy, dt, ifront );
    startIdx += static_cast<long>( ifront );
  }
  if ( s_verbose_vdosgn )
    NCRYSTAL_MSG("VDOSGn convolve: order="<<order.value()
                 <<" truncation/taper block done (phonon_spe.size()="
                 <<phonon_spe.size()<<")");

  int minThinOrder = m_cfg.minThinOrder;
  unsigned thinNBins = m_cfg.thinNBins;
  nc_assert( order.value() < 65000u );
  const int order_int = static_cast<int>(order.value());
  if ( m_cfg.minThinAgressiveOrder >= 0
       && order_int >= m_cfg.minThinAgressiveOrder ) {
    minThinOrder = m_cfg.minThinAgressiveOrder;
    thinNBins = m_cfg.thinAgressiveNBins;
  }

  unsigned long extraThinFactor = 1;
  if ( minThinOrder >= 0 && thinNBins > 0
       && order_int >= minThinOrder
       && phonon_spe.size() > static_cast<std::size_t>(thinNBins) ) {
    // => do thinning
    while ( phonon_spe.size() > thinNBins*extraThinFactor)
      extraThinFactor *= 2;//always orders of 2, allows for on-demand thinning
                           //later (above) without incompatible fractions of
                           //thinFactors.
    if ( extraThinFactor >= 8
         && order.value() <= static_cast<unsigned>(minThinOrder*2) ) {
      //Make brutal thinning slightly less brutal for orders between
      //minThinOrder and (minThinOrder-1)*2:
      extraThinFactor /= 2;
    }

    if ( anchored )
      phonon_spe = thinVectorAnchored( extraThinFactor, phonon_spe, startIdx );
    else
      phonon_spe = thinVector( extraThinFactor, phonon_spe );
    dt *= extraThinFactor;
  }
  if ( s_verbose_vdosgn )
    NCRYSTAL_MSG("VDOSGn convolve: order="<<order.value()
                 <<" extra-thinning block done (phonon_spe.size()="
                 <<phonon_spe.size()<<", extraThinFactor="
                 <<extraThinFactor<<")");

  if ( anchored )
    start_energy = static_cast<double>(startIdx) * dt;

  if (s_verbose_vdosgn) {
    std::ostringstream msg;
    msg<<"VDOSGn Convolved G"<<order1.value()<<"(x)G"<<order2.value()
       <<" -> G"<<order.value()
       <<" ("<<(dt_mismatch?" one input spectrum had to be thinned,":"")
       <<" resulting npts="<<orig_npts_result;
    if (orig_npts_result!=phonon_spe.size())
      msg<<" -> "<<phonon_spe.size()<<" after thinning/truncation";
    msg<<" )";
    NCRYSTAL_MSG(msg.str());
  }

  return VDOSGnData{ phonon_spe, start_energy, dt,
                     thinFactor1*extraThinFactor, startIdx };
}

NCV::VDOSGn::VDOSGn( VDOSGn&& ) = default;
NCV::VDOSGn& NCV::VDOSGn::operator=( VDOSGn&& ) = default;
