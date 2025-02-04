////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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
#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <deque>

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace {
    static std::atomic<bool> s_verbose_vdosgn( ncgetenv_bool("DEBUG_PHONON") );

    class VDOSGnData : private MoveOnly {
    public:
      VDOSGnData( const VectD &spec,
                  double egrid_lower,
                  double egrid_binwidth,
                  unsigned long thinFactor );
      double interpolateDensity(double energy) const;
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
      double m_egrid_lower, m_egrid_upper, m_egrid_binwidth, m_egrid_invbinwidth, m_specMaxVal;
      unsigned long m_thinFactor;//binwidth is G1's binwidth multiplied by m_thinFactor
    };
  }
}

NC::VDOSGnData::VDOSGnData( const NC::VectD &spec,
                            double egrid_lower,
                            double egrid_binwidth,
                            unsigned long thinFactor )
  : m_spec(spec.begin(),spec.end()),
    m_thinFactor(thinFactor)
{
  m_egrid_lower = egrid_lower;
  m_egrid_binwidth = egrid_binwidth;
  nc_assert(m_egrid_binwidth>0.0);
  m_egrid_invbinwidth = 1.0/egrid_binwidth;
  nc_assert(m_spec.size()>3);
  m_spec_size_minus_2 = m_spec.size() - 2;
  m_egrid_upper = m_egrid_lower+(m_spec.size()-1)*m_egrid_binwidth;
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

double NC::VDOSGnData::interpolateDensity(double energy) const
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
  return (*valptr) * (1.0-f) +  f * (*(valptr+1));
}

struct NC::VDOSGn::Impl {
  Impl(const VDOSEval& vde, TruncAndThinningParams);
  std::deque<VDOSGnData> m_gndata;//deque, to not change address of VDOSGnData
                                  //objects when growMaxOrder() is called.
  std::vector<VDOSGnData> m_mt_pending_gndata;//if in MT mode, we might
                                              //precalculate next orders
                                              //concurrently and same them here.
  Optional<FactoryJobs> m_mt_jobs;
  static_assert( std::is_nothrow_default_constructible<Optional<VDOSGnData>>::value, "");
  SmallVector<Optional<VDOSGnData>,10> m_mt_buffer;
  TruncAndThinningParams m_ttpars;
  SmallVector<FastConvolve,4> m_fastConvolve;
  int m_nmaxconcurrent = 0;

  void produceNewOrderByConvolution(Order);
  VDOSGnData produceNewOrderByConvolutionImpl( Order, FastConvolve& ) const;
  VDOSGnData& accessAtOrder(Order n) { nc_assert(n.value()<=m_gndata.size()); return m_gndata[n.value()-1]; }
  const VDOSGnData& accessAtOrder(Order n) const { nc_assert(n.value()<=m_gndata.size()); return m_gndata[n.value()-1]; }

};

NC::VDOSGn::TruncAndThinningParams::TruncAndThinningParams(TruncAndThinningChoices choice)
  : TruncAndThinningParams()
{
  if (choice == TruncAndThinningChoices::Disabled)
    minOrder = -1;
}

NC::VDOSGn::Impl::Impl(const VDOSEval& vde, const TruncAndThinningParams ttpars)
  : m_ttpars(ttpars),
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

  for (auto e: enumerate(egrid) ) {
    nc_assert(e.val>=0.0);
    auto g1_vals = vde.evalG1AsymmetricAtEPair( e.val, gamma0 );
    //Fill at +e.val:
    vectAt(G1spectrum,nbins+e.idx) = g1_vals.second;
    //Fill at -e.val:
    vectAt(G1spectrum,nbins-e.idx) = g1_vals.first;
  }

  nc_assert_always( valueInInterval(0.0,0.1,m_ttpars.truncationThreshold) );
  nc_assert_always( m_ttpars.minOrder >= -1 );

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
  if ( itFirst != itB || itLast != std::prev(itE) ) {
    actual_edgelower += std::distance( itB, itFirst ) * binwidth;
    VectD( itFirst, std::next(itLast) ).swap( G1spectrum );
  }

  //Place G1:
  m_gndata.emplace_back( G1spectrum, actual_edgelower, binwidth, 1 );

  if (s_verbose_vdosgn)
    NCRYSTAL_MSG("VDOSGn constructed (input spectrum size: "<<G1spectrum.size()
                 <<", truncation/thinning with minOrder="<<ttpars.minOrder
                 <<" thinNBins="<<ttpars.thinNBins
                 <<" truncationThreshold="<<ttpars.truncationThreshold
                 <<")");
}

NC::VDOSGn::~VDOSGn() {
  if ( m_impl->m_mt_jobs.has_value() ) {
    //End running jobs, so they don't write to suddenly non-existent buffers:
    m_impl->m_mt_jobs.value().waitAll();
  }
  if (s_verbose_vdosgn)
    NCRYSTAL_MSG("VDOSGn destructed (final max order: "
                 <<maxOrder().value()<<")")

}

NC::VDOSGn::VDOSGn( const NC::VDOSEval& vde, NC::VDOSGn::TruncAndThinningParams ttpars )
  : m_impl(vde,ttpars),
    m_kT(vde.kT())
{
}

NC::VDOSGn::Order NC::VDOSGn::maxOrder() const
{
  return static_cast<unsigned>( m_impl->m_gndata.size() );
}

void NC::VDOSGn::growMaxOrder( Order target_n )
{
  Order n = maxOrder();
  ++n;
  for ( ; n <= target_n; ++n )
    m_impl->produceNewOrderByConvolution(n);

  nc_assert( maxOrder().value() == target_n.value() );
}

double NC::VDOSGn::eval( Order n, double energy ) const
{
  return m_impl->accessAtOrder(n).interpolateDensity(energy);
}

const NC::VectD& NC::VDOSGn::getRawSpectrum( NC::VDOSGn::Order n ) const
{
  return m_impl->accessAtOrder(n).getSpectrum();
}

double NC::VDOSGn::binWidth( NC::VDOSGn::Order n) const
{
  return m_impl->accessAtOrder(n).getEGridBinwidth();
}

NC::PairDD NC::VDOSGn::eRange( NC::VDOSGn::Order n, double relthreshold ) const
{
  nc_assert(relthreshold>0.0&&relthreshold<1.0);
  const auto& p = m_impl->accessAtOrder(n);
  const auto& spec = p.getSpectrum();
  const double spec_max = p.maxDensity();
  const double threshold = relthreshold * spec_max;
  PairDD erange(p.getEGridLower(), p.getEGridUpper());

  for ( auto e :  enumerate(spec) ) {
    if ( e.val >= threshold ) {
      erange.first = p.getEGridLower() + e.idx * p.getEGridBinwidth();
      break;
    }
  }

  for (std::size_t i = spec.size(); i>0; --i) {
    if ( vectAt(spec,i-1) >= threshold ) {
      erange.second = ncmin(erange.second,p.getEGridLower() + (i-1) * p.getEGridBinwidth());
      break;
    }
  }
  nc_assert( erange.second >= erange.first );
  return erange;
}

NC::PairDD NC::VDOSGn::eRange( NC::VDOSGn::Order n ) const
{
  const auto& p = m_impl->accessAtOrder(n);
  return { p.getEGridLower(), p.getEGridUpper() };
}

void NC::VDOSGn::enableVerboseOutput(bool status)
{
  s_verbose_vdosgn = status;
}

bool NC::VDOSGn::verboseOutputEnabled()
{
  return s_verbose_vdosgn;
}

void NC::VDOSGn::Impl::produceNewOrderByConvolution( Order order )
{
  const unsigned current_maxorder = static_cast<unsigned>( m_gndata.size() );
  nc_assert_always( order.value() == current_maxorder + 1 );

  if ( m_mt_jobs.has_value() ) {
    m_mt_jobs.value().waitAll();
    //Transfer concurrently generated results:
    for ( auto i : ncrange( m_mt_buffer.size() ) )
      m_mt_pending_gndata.emplace_back( std::move( m_mt_buffer.at( m_mt_buffer.size()-1-i ).value() ) );
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

NC::VDOSGnData NC::VDOSGn::Impl::produceNewOrderByConvolutionImpl( Order order, FastConvolve& fastConvolve ) const
{
  Order order2 = order.value()/2;
  Order order1 = order.value()-order2.value();

  //NCRYSTAL_MSG("Convolve G"<<order.value()<<" = conv( G"<<order2.value()<<", G"<<order.value()<<" )");

  const auto& p1 = accessAtOrder(order1);
  const auto& p2 = accessAtOrder(order2);

  //Function which can thin a vector (i.e. increase binwidth by merging bins),
  //used two places below:
  auto thinVector = [](unsigned thinFactor, const NC::VectD& v)
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
      vtmp = thinVector(thinFactor,*input2_spec);
      input2_spec = &vtmp;
      thinFactor2 *= thinFactor;
    } else {
      //thin dt1
      nc_assert(thinFactor2%thinFactor1==0);
      unsigned long thinFactor = thinFactor2 / thinFactor1;
      nc_assert_always( floateq(dt,dt1*thinFactor) );
      vtmp = thinVector(thinFactor,*input1_spec);
      input1_spec = &vtmp;
      thinFactor1 *= thinFactor;
    }
  }
  nc_assert_always(thinFactor1==thinFactor2);

  VectD phonon_spe;
  double start_energy = p1.getEGridLower() + p2.getEGridLower();
  fastConvolve.convolve( *input1_spec, *input2_spec, phonon_spe, dt );
  auto orig_npts_result = phonon_spe.size();

  unsigned long extraThinFactor = 1;
  if ( m_ttpars.minOrder >= 0 && order.value() >= static_cast<unsigned>(m_ttpars.minOrder) ) {
    //We should do truncation and/or thinning at this order.
    if (m_ttpars.truncationThreshold > 0 ) {
      // => do truncation
      const double spec_max = *std::max_element(phonon_spe.begin(),phonon_spe.end());
      const double spec_cutoff = m_ttpars.truncationThreshold * spec_max;
      std::size_t ifront(0), iback(phonon_spe.size()-1);
      for (;ifront<iback;++ifront) {
        if (phonon_spe.at(ifront)>spec_cutoff)
          break;
      }
      for (;iback>ifront;--iback) {
        if (phonon_spe.at(iback)>spec_cutoff)
          break;
      }
      if (iback>ifront) {
        VectD truncated_spec(phonon_spe.begin()+ifront,phonon_spe.begin()+iback+1);
        truncated_spec.swap(phonon_spe);
      }
      start_energy += ifront*dt;
    }
    if ( m_ttpars.thinNBins > 0 && phonon_spe.size() > m_ttpars.thinNBins ) {
      // => do thinning
      while ( phonon_spe.size() > m_ttpars.thinNBins*extraThinFactor)
        extraThinFactor *= 2;//always orders of 2, allows for on-demand thinning
                             //later (above) without incompatible fractions of
                             //thinFactors.
      if ( extraThinFactor >= 8 && order.value() <= static_cast<unsigned>(m_ttpars.minOrder*2) ) {
        //Make brutal thinning slightly less brutal for orders between minOrder
        //and (minOrder-1)*2:
        extraThinFactor /= 2;
      }

      phonon_spe = thinVector( extraThinFactor, phonon_spe );
      dt *= extraThinFactor;
    }
  }

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

  return VDOSGnData{ phonon_spe, start_energy, dt, thinFactor1*extraThinFactor };
}
