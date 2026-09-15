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

#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCStableDbl.hh"
#include "NCrystal/internal/utils/NCTinyVector.hh"
#include <complex>
namespace NC = NCrystal;

//Temporarily uncomment the following define to test with safer but slower code:
//#define NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH

namespace NCRYSTAL_NAMESPACE {

  namespace {

    //Independent of data, we need the same tables of W factors and swap
    //patterns. We make sure we can reuse calculations for these if needed, by
    //keeping a common global cache of these in a dedicated manager class:

    //For maximal precision, we use StableDbl's for our phase calculations and a
    //lookup table for values of exp(i*2pi/2^n):

    using ComplexSD = std::pair<StableDbl, StableDbl>;

    inline ComplexSD complexMult( const ComplexSD& a, const ComplexSD& b)
    {
      return { a.first * b.first - a.second * b.second,
               a.first * b.second + a.second * b.first };
    }

    class FastConvolveCacheMgr final : private NoCopyMove {
    public:
      struct SwapPatternCache final {
        std::vector<std::pair<unsigned long,unsigned long>> pattern;
        int output_log_size = 0;
      };
      using WTable = std::vector< std::complex<double> >;

      static shared_obj<WTable> defaultWTable()
      {
        static auto wt = makeSO<WTable>();
        return wt;
      }

      static shared_obj<SwapPatternCache> defaultSwapPattern()
      {
        static auto spc = makeSO<SwapPatternCache>();
        return spc;
      }

      shared_obj<WTable> getWTable( unsigned long k ) const
      {
        {
          NCRYSTAL_LOCK_GUARD(m_w_cache_mutex);
          auto it = m_w_cache.find(k);
          if ( it != m_w_cache.end() )
            return it->second;
        }
        //Create from scratch:
        auto res = makeSO<WTable>();
        initWTable( k, res );
        {
          //insert (unless someone beat us to it):
          NCRYSTAL_LOCK_GUARD(m_w_cache_mutex);
          auto it = m_w_cache.find(k);
          if ( it != m_w_cache.end() )
            return it->second;
          m_w_cache.insert( {k,res} );
        }
        return res;
      }

      shared_obj<SwapPatternCache> getSwapPattern( int k ) const
      {
        {
          NCRYSTAL_LOCK_GUARD(m_swap_cache_mutex);
          auto it = m_swap_cache.find(k);
          if ( it != m_swap_cache.end() )
            return it->second;
        }
        //Create from scratch:
        auto res = makeSO<SwapPatternCache>();
        initSwapPattern( k, res );
        {
          //insert (unless someone beat us to it):
          NCRYSTAL_LOCK_GUARD(m_swap_cache_mutex);
          auto it = m_swap_cache.find(k);
          if ( it != m_swap_cache.end() )
            return it->second;
          m_swap_cache.insert( {k, res} );
        }
        return res;
      }

      void clearCaches()
      {
        NCRYSTAL_LOCK_GUARD(m_w_cache_mutex);
        NCRYSTAL_LOCK_GUARD(m_swap_cache_mutex);
        m_w_cache.clear();
        m_swap_cache.clear();
      }
      static ComplexSD calcPhaseSD(unsigned long k, unsigned long n);
    private:
      void initWTable( unsigned long, WTable& ) const;
      void initSwapPattern( int, SwapPatternCache& ) const;

      mutable std::map<int,shared_obj<WTable>> m_w_cache;
      mutable std::map<int,shared_obj<SwapPatternCache>> m_swap_cache;
      mutable std::mutex m_w_cache_mutex;
      mutable std::mutex m_swap_cache_mutex;
    };

    FastConvolveCacheMgr& getFastConvolveCacheMgr()
    {
      static FastConvolveCacheMgr mgr;
      static const bool dummy = []()
      {
        registerCacheCleanupFunction([]()
        {
          getFastConvolveCacheMgr().clearCaches();
        });
        return false;
      }();
      (void) dummy;
      return mgr;
    }
  }

  struct FastConvolve::Impl {
    using WTable = FastConvolveCacheMgr::WTable;
    using SwapPatternCache = FastConvolveCacheMgr::SwapPatternCache;

    Impl()
      : m_w(FastConvolveCacheMgr::defaultWTable()),
        m_swap(FastConvolveCacheMgr::defaultSwapPattern())
    {
    }

    shared_obj<WTable> m_w;
    shared_obj<SwapPatternCache> m_swap;

    //The actual fast-fourier transform algorithm:
    template<bool is_forward>
    void fft( std::vector<std::complex<double> > &inout,
              unsigned long minimum_output_size );
    void applySwaps( const SwapPatternCache&,
                     std::vector<std::complex<double>>& ) const;
    void convolve( const VectD& a1, const VectD& a2,
                   VectD& y, double dt,
                   bool useLegacyBehaviour );

  };
}

NC::FastConvolve::FastConvolve( FastConvolve&& o ) noexcept
  : m_impl( std::move(o.m_impl) )
{
  static_assert( std::is_nothrow_move_constructible<decltype(m_impl)>::value, "" );
}

NC::FastConvolve& NC::FastConvolve::operator=( FastConvolve&& o ) noexcept
{
  static_assert( std::is_nothrow_move_assignable<decltype(m_impl)>::value, "" );
  m_impl = std::move(o.m_impl);
  return *this;
}

void NC::FastConvolveCacheMgr::initWTable( unsigned long n_size_raw,
                                           std::vector< std::complex<double> >& wtable ) const
{
  static_assert( sizeof(unsigned long) >= sizeof(std::uint32_t), "");
  nc_assert_always(n_size_raw<=1000000000u);
  //round n_size up to next power of 2, and also:
  unsigned long nsize = 1;
  unsigned long log2_nsize = 0;
  while ( nsize < n_size_raw ) {
    nsize <<= 1;
    ++log2_nsize;
  }

  wtable.clear();
  wtable.reserve(nsize);

  ComplexSD phaseval{1.0,0.0};
  const ComplexSD phase1n = calcPhaseSD(1, log2_nsize);
  for( unsigned long i = 0; i < nsize; ++i ) {
    if ( i%2==1 ) {
      //for odd i, we anyway would end up calculating it like this inside
      //calcPhase - but this way we can reuse the previous phaseval:
      phaseval = complexMult( phase1n, phaseval);
    } else {
      phaseval = calcPhaseSD(i, log2_nsize);
    }
    wtable.emplace_back(std::complex<double>(phaseval.first.value(),
                                             phaseval.second.value()));
  }
}

void NC::FastConvolve::convolve( const VectD& a1, const VectD& a2,
                                 VectD& y, double dt )
{
  m_impl->convolve(a1,a2,y,dt,false);
}

void NC::FastConvolve::convolveLegacy( const VectD& a1, const VectD& a2,
                                       VectD& y, double dt )
{
  m_impl->convolve(a1,a2,y,dt,true);
}


void NC::FastConvolve::Impl::convolve( const VectD& a1, const VectD& a2,
                                       VectD& y, double dt,
                                       bool useLegacyBehaviour )
{
  const int minimum_out_size = a1.size() + a2.size() - 1;

  //Note: We could calculate the next two fft calls concurrently, but it was
  //attempted and did not work as well as instead employing concurrency in
  //NCVDOSGn.cc, so we keep that for now.

  std::vector<std::complex<double> > b1(a1.begin(),a1.end());
  fft<true>(b1,minimum_out_size);
  std::vector<std::complex<double> > b2(a2.begin(),a2.end());
  fft<true>(b2,minimum_out_size);

  nc_assert(b1.size()==b2.size());
  std::vector<std::complex<double> >::iterator
    itb1(b1.begin()), itb1E(b1.end()), itb2(b2.begin());
  while (itb1!=itb1E)
    *itb1++ *= *itb2++;

  fft<false>(b1,minimum_out_size);

  y.resize(minimum_out_size);
  const double k = dt/b1.size();
  nc_assert(b1.size()==b2.size());
  nc_assert(y.size()<=b1.size());
  VectD::iterator ity(y.begin()), ityE(y.end());
  itb1 = b1.begin();

  //LegacyBehaviour: Wrongly used std::abs. However, result is supposed to be
  //real already, with only tiny numerical fluctuations in the imaginary
  //component. Correct behaviour is instead take the real part of result.
  //
  //For a use case like VDOS convolutions where we expect output to be purely
  //positive, the calling code should itself clamp output to positive values if
  //needed.

  if ( useLegacyBehaviour ) {
    //Wrongly used std::abs. However, result is supposed to be real already,
    //with only tiny numerical fluctuations in the imaginary component.
#   ifdef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
    for(;ity!=ityE;++ity,++itb1) {
      //NB: Use std::abs which calls std::hypot behind the
      //    scenes (expensive but can avoid overflows)
      *ity = k * std::abs(*itb1);
    }
#   else
    //Naive and simple, avoids std::hypot, more easily vectorisable:
    for(;ity!=ityE;++ity,++itb1) {
      double a(itb1->real());
      double b(itb1->imag());
      *ity = a*a+b*b;
    }
    for(ity = y.begin();ity!=ityE;++ity)
      *ity = std::sqrt(*ity);
    for(ity = y.begin();ity!=ityE;++ity)
      *ity *= k;
#   endif
  } else {
    //Result is the real part of the calculations:
    for(;ity!=ityE;++ity,++itb1)
      *ity = k * itb1->real();
  }
}

void NC::FastConvolveCacheMgr::initSwapPattern( int output_log_size, SwapPatternCache& swap_cache ) const
{
  nc_assert( swap_cache.output_log_size != output_log_size );

  decltype(swap_cache.pattern) pattern;
  pattern.swap(swap_cache.pattern);
  swap_cache.output_log_size = 0;
  pattern.clear();
  pattern.reserve(65536);

  const int output_size = 1 << output_log_size;
  const int output_size_m1 = output_size-1;
  for(int j=1;j<output_size_m1;++j)
    {
#if 1
      int tmp = j;
      int i = tmp & 1;
      tmp >>= 1;
      for( int k = 1; k < output_log_size; ++k )
        {
          i = (i<<1) | (tmp&1);
          tmp >>= 1;
        }
#else
      int i=0;
      for(int k=1,tmp=j;
          k<output_size;
          i=(i<<1)|(tmp&1),k<<=1,tmp>>=1)
        {
        }
#endif
      if(j<i) {
        nc_assert( j < output_size );
        nc_assert( i < output_size );
        pattern.emplace_back(i<<1,j<<1);//<<1 is *2
      }
    }

  swap_cache.output_log_size = output_log_size;
  pattern.swap(swap_cache.pattern);
}

void NC::FastConvolve::Impl::applySwaps( const SwapPatternCache& swapcache,
                                         std::vector<std::complex<double>> &data ) const
{
  //Avoid direct std::complex usage. This next cast is actually OK by the c++11
  //standard (https://stackoverflow.com/questions/69591371):
  double * rawdata = reinterpret_cast<double*>( data.data() );
  for ( auto& e : swapcache.pattern ) {
    double* it1 = std::next(rawdata,(e.first));//We already had a *2 applied to the index
    double* it2 = std::next(rawdata,(e.second));//We already had a *2 applied to the index
    std::swap(*it1++, *it2++);
    std::swap(*it1, *it2);
  }
}

template<bool is_forward>
void NC::FastConvolve::Impl::fft( std::vector<std::complex<double>> &data,
                                  unsigned long minimum_output_size )
{
  const double output_log_size_fp = std::ceil(std::log2(minimum_output_size));
  static_assert( sizeof(int) == sizeof(std::int32_t), "" );//otherwise we need
                                                           //to update the code
                                                           //in this class!!
  nc_assert_always(output_log_size_fp<32);
  const int output_log_size = output_log_size_fp;
  const int output_size = ( 1 << output_log_size );//this is now
                                                   //minimum_output_size rounded
                                                   //up to next power of 2
  nc_assert_always( data.size() <= (std::size_t)output_size );

  if ( m_w->size() < (std::size_t)output_size )
    m_w = getFastConvolveCacheMgr().getWTable( output_size );

  const auto& wtable = *m_w;

  if( data.size() != (size_t)output_size )
    data.resize(output_size,std::complex<double>());
#if 1
  if ( output_log_size != m_swap->output_log_size )
    m_swap = getFastConvolveCacheMgr().getSwapPattern( output_log_size );
  nc_assert( data.size() == (std::size_t)( 1 << output_log_size ) );
  nc_assert( output_log_size == m_swap->output_log_size );
  applySwaps( m_swap, data );
#else
  //Old, without cached swaps:
  const int output_size_m1 = output_size-1;
  for(int j=1;j<output_size_m1;++j) {
    int i=0;
    //NB: output_size is power of two (2**output_log_size to be exact)
    for(int k=1,tmp=j;
        k<output_size;
        i=(i<<1)|(tmp&1),k<<=1,tmp>>=1)
      {
      }
    if(j<i) {
      std::swap(data[i],data[j]);
    }
  }
#endif

  nc_assert_always(wtable.size()%output_size==0);
  const int jump = wtable.size()/output_size;

#ifndef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
  double * rawdata = reinterpret_cast<double*>( data.data() );
  const double * raww = reinterpret_cast<const double*>( wtable.data() );
#endif
  for(int i=0;i<output_log_size;++i){
    int z=0;
    const int i1 = (1<<i);
#ifndef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
    auto twoi1 = i1*2;
#endif

    const int i1m1 = i1-1;
    const int i2 = 1<<(output_log_size-i-1);

    for(int j=0;j<output_size;++j){

#if 0
      //original (integer division very expensive!):
      if((j/i1)%2) {
#else
      //NB: j/i1=j>>i since i1=(1<<i):
      if( (j>>i)%2 ) {
        //Todo: Perhaps we can figure out a bit mask for the check (j>>i)%2
        //instead of doing it each time?
        nc_assert( j/i1 == j>>i );
#endif
#ifdef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
        std::complex<double>& data_j = data[j];
        std::complex<double>& data_sympos = data[j-i1];
        //std::complex<> multiplication is slow since it takes care of proper inf/nan/overflow
        data_j *= ((!is_forward)?std::conj(wtable[z*jump]):wtable[z*jump]);
        //and the -=,+= operators seems to carry significant overhead for some reason:
        std::complex<double> temp = data_sympos;
        data_sympos += data_j;
        temp -= data_j;
        data_j = temp;
#else
        //naive and simple is faster:
#  if 1
        //v2:
        auto twoj = j*2;
        double* rawdata_j = std::next(rawdata,twoj);
        double* rawdata_sympos = std::next(rawdata,twoj-twoi1);
        const double* rawdata_wzjump = std::next(raww,z*jump*2);
        //const auto& w_zjump = wtable[z*jump];
        double& a = *rawdata_j;// accessReal( data_j );
        double& b = *std::next(rawdata_j); //accessImag( data_j );
        const double& c = *rawdata_wzjump;//accessReal( w_zjump );
        const double d = ( is_forward
                           ? (*std::next(rawdata_wzjump))//accessImag( w_zjump );
                           : -(*std::next(rawdata_wzjump)) );//accessImag( w_zjump );
        //const double d = convfact * (*std::next(rawdata_wzjump));//accessImag( w_zjump );
        const double jr(a*c-b*d);
        const double ji(a*d+b*c);
        double& sr = *rawdata_sympos;//accessReal(data_sympos);
        double& si = *std::next(rawdata_sympos);//accessImag(data_sympos);
        a = sr - jr;
        b = si - ji;
        sr += jr;
        si += ji;
#  else
        //v1:
        const std::complex<double>& w_zjump = wtable[z*jump];
        const double a(data_j.real()), b(data_j.imag()), c(w_zjump.real());
        const double d = ( is_forward ? w_zjump.imag() : - w_zjump.imag() );
        const double jr(a*c-b*d);
        const double ji(a*d+b*c);
        const double sr(data_sympos.real());//no by-ref access to real/imag parts
        const double si(data_sympos.imag());
        data_j.real( sr - jr );
        data_j.imag( si - ji );
        data_sympos.real( sr + jr );
        data_sympos.imag( si + ji );
#  endif
#endif
        z += i2;
      } else {
        z = 0;
        //will be the same result for the next i1-1 loops, so skip ahead:
        j += i1m1;
      }
    }
  }
}

NC::PairDD NC::FastConvolve::calcPhase(unsigned long k, unsigned long n)
{
  auto p = FastConvolveCacheMgr::calcPhaseSD(k, n);
  return { p.first.value(), p.second.value() };
}

NC::ComplexSD NC::FastConvolveCacheMgr::calcPhaseSD(unsigned long k,
                                                    unsigned long n)
{
  //Calculate exp(i*2*pi*k/2^n). Aim is to do so fast, but more importantly as
  //precise as possible with results reproducible on all platforms.
  //
  //Idea, use exp(a*b)=exp(a)*exp(b) and a small cache of
  //exp(i2pi/2^N), N=0,1,2,. when k!=1, one can combine as in these examples:
  //
  // 13/32 = 1/32 + 12/32 = 1/32 + 3/8 = 1/32 + 1/8 + 1/4
  //
  // 17/32 = 1/32 + 1/2
  //
  // 11/32 = 1/32 + 5/16 = 1/32 + 1/16 + 1/4
  //
  // Additional idea: when k/2^n > 1/2 we can initially lower the effective
  // value by evaluating with q=2^n-k instead of k, which gives us the
  // conjugated value of the result we are after, since:
  //
  // exp(i*2*pi*k/2^n) = exp(-i*2*pi*(q/2^n) = conjugate(exp(i*2*pi*(q/2^n)))
  //
  // That gives us the advantage of less numbers to multiply, with less
  // accumulated error.

  constexpr unsigned long ncache = 33;
  static_assert( sizeof(unsigned long) >= sizeof(std::uint32_t), "" );
  nc_assert( n <= 30u );//We could go up to ncache-1, but at ~1e9 2^30 already
                        //consumes too much memory.

  //Trivial case:
  if ( n == 0 || k == 0 ) {
    return { 1.0, 0.0 };
  }

  //Eliminate common factors of 2 in fraction k/2^n:
  while ( k%2==0 ) {
    nc_assert( n>=1 );
    n -= 1;
    k /= 2;
  }
  const unsigned long two_to_n = 1u<<n;
  nc_assert(two_to_n<=1073741824);
  nc_assert( k < two_to_n );
  if ( n == 1 && k == 1 )
    return { -1.0, 0.0 };

  nc_assert(n>=2);//since (n,k) = (0,0),(1,0), (1,1) already dealt with
  const unsigned long two_to_nminus1 = 1u<<(n-1u);

  if ( k >= two_to_nminus1 ) {
    if ( k==two_to_nminus1 )
      return { -1.0, 0.0 };
    //i2pi*k/2^n is in 3rd of 4th quadrant => map to 1st or 2nd quadrant.
    auto p = calcPhaseSD(two_to_n-k, n);
    return { p.first, -p.second };
  }
  const unsigned long two_to_nminus2 = 1u<<(n-2u);
  if ( k >= two_to_nminus2 ) {
    if ( k == two_to_nminus2 )
      return { 0.0, 1.0 };
    //i2pi*k/2^n is in 2nd quadrant => map to 1st quadrant.
    auto p = calcPhaseSD( two_to_nminus1-k, n );
    return { -p.first, p.second };
  }

  //Ok, we now have k < 2^(n-2), i.e. we are in the 1st quadrant.

  if ( k == 1 ) {
    //Fundamental form, trigonometric argument is 2pi/2^n. Since there are only
    //~32 possible n values, we use high-res lookup table for these to ensure
    //maximal precision and cross-platform consistency.
    nc_assert(n>=1);//k>0, so it follows from k<=n-1 that n>1
    //Cache high-precision numbers (from Python's mpmath module) for efficiency,
    //accuracy and reproducability. The lookup tables contains 37 significant
    //digits, so would even be OK for float128):
    static std::array<double, ncache> cosvals = {
      1.0,  // cos(2pi/2^0)
      -1.0, // cos(2pi/2^1)
      0.0,  // cos(2pi/2^2)
      7.071067811865475244008443621048490393e-1, // cos(2pi/2^3)
      9.238795325112867561281831893967882868e-1, // cos(2pi/2^4)
      9.80785280403230449126182236134239037e-1, // cos(2pi/2^5)
      9.951847266721968862448369531094799216e-1, // cos(2pi/2^6)
      9.987954562051723927147716047591006944e-1, // cos(2pi/2^7)
      9.996988186962042201157656496661721969e-1, // cos(2pi/2^8)
      9.999247018391445409216464911963832244e-1, // cos(2pi/2^9)
      9.999811752826011426569904377285677162e-1, // cos(2pi/2^10)
      9.999952938095761715115801257001198996e-1, // cos(2pi/2^11)
      9.99998823451701909929025710171526019e-1, // cos(2pi/2^12)
      9.999997058628822191602282177387656771e-1, // cos(2pi/2^13)
      9.999999264657178511447314807073878569e-1, // cos(2pi/2^14)
      9.999999816164292938083469154029097145e-1, // cos(2pi/2^15)
      9.99999995404107312890971933139606149e-1, // cos(2pi/2^16)
      9.999999988510268275626733077945541084e-1, // cos(2pi/2^17)
      9.999999997127567068494139722186417761e-1, // cos(2pi/2^18)
      9.999999999281891767097750958838504903e-1, // cos(2pi/2^19)
      9.999999999820472941772826241477841074e-1, // cos(2pi/2^20)
      9.999999999955118235443105841729973244e-1, // cos(2pi/2^21)
      9.999999999988779558860770165517525365e-1, // cos(2pi/2^22)
      9.999999999997194889715192147947195845e-1, // cos(2pi/2^23)
      9.999999999999298722428798012397287368e-1, // cos(2pi/2^24)
      9.999999999999824680607199501562477367e-1, // cos(2pi/2^25)
      9.999999999999956170151799875294566562e-1, // cos(2pi/2^26)
      9.999999999999989042537949968817638342e-1, // cos(2pi/2^27)
      9.999999999999997260634487492204034379e-1, // cos(2pi/2^28)
      9.999999999999999315158621873050985144e-1, // cos(2pi/2^29)
      9.99999999999999982878965546826274482e-1, // cos(2pi/2^30)
      9.999999999999999957197413867065686114e-1, // cos(2pi/2^31)
      9.999999999999999989299353466766421523e-1 // cos(2pi/2^32)
    };
    static std::array<double, ncache> sinvals = {
      0.0, // sin(2pi/2^0)
      0.0, // sin(2pi/2^1)
      1.0, // sin(2pi/2^2)
      7.071067811865475244008443621048490393e-1, // sin(2pi/2^3)
      3.826834323650897717284599840303988668e-1, // sin(2pi/2^4)
      1.950903220161282678482848684770222409e-1, // sin(2pi/2^5)
      9.801714032956060199419556388864184586e-2, // sin(2pi/2^6)
      4.906767432741801425495497694268265831e-2, // sin(2pi/2^7)
      2.454122852291228803173452945928292507e-2, // sin(2pi/2^8)
      1.227153828571992607940826195100321214e-2, // sin(2pi/2^9)
      6.135884649154475359640234590372580917e-3, // sin(2pi/2^10)
      3.067956762965976270145365490919842519e-3, // sin(2pi/2^11)
      1.53398018628476561230369715026407908e-3, // sin(2pi/2^12)
      7.669903187427045269385683579485766431e-4, // sin(2pi/2^13)
      3.834951875713955890724616811813812634e-4, // sin(2pi/2^14)
      1.917475973107033074399095619890009335e-4, // sin(2pi/2^15)
      9.587379909597734587051721097647635119e-5, // sin(2pi/2^16)
      4.793689960306688454900399049465887275e-5, // sin(2pi/2^17)
      2.396844980841821872918657716502182009e-5, // sin(2pi/2^18)
      1.19842249050697064215215615969889848e-5, // sin(2pi/2^19)
      5.99211245264242784287971180889086173e-6, // sin(2pi/2^20)
      2.996056226334660750454812808357059812e-6, // sin(2pi/2^21)
      1.498028113169011228854278846155361121e-6, // sin(2pi/2^22)
      7.490140565847157211304985667306556372e-7, // sin(2pi/2^23)
      3.745070282923841239031691790846331774e-7, // sin(2pi/2^24)
      1.872535141461953448688245765935636171e-7, // sin(2pi/2^25)
      9.362675707309808279906728668088562019e-8, // sin(2pi/2^26)
      4.681337853654909269511551813854009696e-8, // sin(2pi/2^27)
      2.340668926827455275950549341903484404e-8, // sin(2pi/2^28)
      1.17033446341372771812462135032381038e-8, // sin(2pi/2^29)
      5.851672317068638690809790100834139694e-9, // sin(2pi/2^30)
      2.925836158534319357928230469068955902e-9, // sin(2pi/2^31)
      1.46291807926715968052953216186596371e-9 // sin(2pi/2^32)
    };

    nc_assert( n < ncache );
    return { cosvals[n], sinvals[n] };
  }

  //Non-fundamental form, must combine results from several fundamental forms
  //using multiplication of complex numbers.
  nc_assert(k%2==1);//must be odd at this point

  //NB: This turns k-1 into an even value, meaning the second factor will
  //actually be evaluated at ((k-1)/2,n-1):
  return complexMult( calcPhaseSD( 1, n ),
                      calcPhaseSD( k-1, n ) );
}

//Here due to pimpl:
NC::FastConvolve::FastConvolve() = default;
NC::FastConvolve::~FastConvolve() = default;
