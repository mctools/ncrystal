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

#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <complex>
namespace NC = NCrystal;

//Temporarily uncomment the following define to test with safer but slower code:
//#define NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH

namespace NCRYSTAL_NAMESPACE {

  namespace {

    //Independent of data, we need the same tables of W factors and swap
    //patterns. We make sure we can reuse calculations for these if needed, by
    //keeping a common global cache of these in a dedicated manager class:

    class FastConvolveCacheMgr : NoCopyMove {
    public:
      struct SwapPatternCache {
        std::vector<std::pair<unsigned,unsigned>> pattern;
        int output_log_size = 0;
      };
      using WTable = std::vector< std::complex<double> >;

      static shared_obj<WTable> defaultWTable() { static auto wt = makeSO<WTable>(); return wt; }
      static shared_obj<SwapPatternCache> defaultSwapPattern() { static auto spc = makeSO<SwapPatternCache>(); return spc; }

      shared_obj<WTable> getWTable( unsigned k ) const
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

      static PairDD calcPhase(unsigned k, unsigned n);
    private:
      void initWTable( unsigned, WTable& ) const;
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
        registerCacheCleanupFunction([](){ getFastConvolveCacheMgr().clearCaches(); });
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
    void fft( std::vector<std::complex<double> > &inout,unsigned minimum_output_size );
    void applySwaps( const SwapPatternCache&, std::vector<std::complex<double>>& ) const;
    void convolve( const NC::VectD& a1, const NC::VectD& a2, NC::VectD& y, double dt );

  };

  namespace {

    inline PairDD detail_mult_imag( PairDD a, PairDD b ) noexcept
    {
      // (a+ib)*(c+id) = (a*c-b*d) + i* (a*d+b*c)
      return { (a.first*b.first-a.second*b.second),
               (a.first*b.second+a.second*b.first) };
    }

    // //According to cppreference.com, the layout of std::complex<double> MUST be
    // //the same as a double[2] array. Hence, faster access is provided with:
    // inline double* complexAsArray( std::complex<double>& v ) noexcept { return reinterpret_cast<double*>(&v); }
    // inline const double* complexAsArray( const std::complex<double>& v ) noexcept { return reinterpret_cast<const double*>(&v); }
    // inline double& accessReal( std::complex<double>& v ) noexcept { return reinterpret_cast<double(&)[2]>(v)[0]; }
    // inline const double& accessReal( const std::complex<double>& v ) noexcept { return reinterpret_cast<const double(&)[2]>(v)[0]; }
    // inline double& accessImag( std::complex<double>& v ) noexcept { return reinterpret_cast<double(&)[2]>(v)[1]; }
    // inline const double& accessImag( const std::complex<double>& v ) noexcept { return reinterpret_cast<const double(&)[2]>(v)[1]; }

  }
}

NC::FastConvolve::FastConvolve() = default;
NC::FastConvolve::~FastConvolve() = default;

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

void NC::FastConvolveCacheMgr::initWTable( unsigned n_size_raw, std::vector< std::complex<double> >& wtable ) const
{
  //round n_size up to next power of 2, and also:
  unsigned nsize = 1;
  unsigned log2_nsize = 0;
  while ( nsize < n_size_raw ) {
    nsize <<= 1;
    ++log2_nsize;
  }

  wtable.clear();
  wtable.reserve(nsize);

#if 1
  PairDD phaseval{1.0,0.0};
  const PairDD phase1n = ::NC::FastConvolve::calcPhase(1, log2_nsize);
  for( unsigned i = 0; i < nsize; ++i ) {
    if ( i%2==1 ) {
      //for odd i, we anyway would end up calculating it like this inside
      //calcPhase - but this way we can reuse the previous phaseval:
      phaseval = detail_mult_imag(phase1n,phaseval);
    } else {
      phaseval = calcPhase(i, log2_nsize);
    }
    wtable.emplace_back(std::complex<double>(phaseval.first,phaseval.second));
  }
#else
  for( unsigned i = 0; i < nsize; ++i ) {
    PairDD phaseval = calcPhase(i, log2_nsize);
    wtable.emplace_back(std::complex<double>(phaseval.first,phaseval.second));
  }
#endif
}

void NC::FastConvolve::convolve( const NC::VectD& a1, const NC::VectD& a2,
                                 NC::VectD& y, double dt )
{
  m_impl->convolve(a1,a2,y,dt);
}


void NC::FastConvolve::Impl::convolve( const NC::VectD& a1, const NC::VectD& a2,
                                       NC::VectD& y, double dt )
{
  const int minimum_out_size = a1.size() + a2.size() - 1;

  //Note: We could calculate the next two fft calls concurrently, but it was
  //attempted and did not work as well, as instead employing concurrency in
  //NCVDOSGn.cc, so we keep that for now.

  std::vector<std::complex<double> > b1(a1.begin(),a1.end());
  fft<true>(b1,minimum_out_size);
  std::vector<std::complex<double> > b2(a2.begin(),a2.end());
  fft<true>(b2,minimum_out_size);

  nc_assert(b1.size()==b2.size());
  std::vector<std::complex<double> >::iterator itb1(b1.begin()), itb1E(b1.end()), itb2(b2.begin());
  while (itb1!=itb1E)
    *itb1++ *= *itb2++;

  fft<false>(b1,minimum_out_size);

  y.resize(minimum_out_size);
  const double k = dt/b1.size();
  nc_assert(b1.size()==b2.size());
  nc_assert(y.size()<=b1.size());
  VectD::iterator ity(y.begin()), ityE(y.end());
  itb1 = b1.begin();
#ifdef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
  for(;ity!=ityE;++ity,++itb1) {
    //use std::abs which calls std::hypot behind the scenes (expensive but can avoid overflows)
    *ity = k * std::abs(*itb1);
  }
#else
  //naive and simple, avoids std::hypot, more easily vectorisable:
  for(;ity!=ityE;++ity,++itb1) {
    double a(itb1->real());
    double b(itb1->imag());
    *ity = a*a+b*b;
  }
  for(ity = y.begin();ity!=ityE;++ity)
    *ity = std::sqrt(*ity);
  for(ity = y.begin();ity!=ityE;++ity)
    *ity *= k;
#endif
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

void NC::FastConvolve::Impl::applySwaps( const SwapPatternCache& swapcache, std::vector<std::complex<double>> &data ) const
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
                                  unsigned minimum_output_size )
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
  //const double convfact = (is_forward?1.0:-1.0);
#endif

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

NC::PairDD NC::FastConvolve::calcPhase(unsigned k, unsigned n)
{
  return FastConvolveCacheMgr::calcPhase(k,n);
}

NC::PairDD NC::FastConvolveCacheMgr::calcPhase(unsigned k, unsigned n)
{
  //Calculate exp(i*2*pi*k/2^n) where n must a nonzero number n=1,2,3,4,...
  //and k must be a number in 0...n-1.
  //
  //NB: Using PairDD for (real,imag) parts, rather than std::complex<double> => for
  //efficiency (we don't need the safety checks for obscure cases...).
  //
  //Idea, use exp(a*b)=exp(a)*exp(b) and a small cache of
  //exp(i2pi/2^N), N=0,1,2,. when k!=1, one can combine as in these examples:
  //
  // 13/32 = 1/32 + 12/32 = 1/32 + 3/8 = 1/32 + 1/8 + 1/4
  //
  // 17/32 = 1/32 + 1/2
  //
  // 11/32 = 1/32 + 5/16 = 1/32 + 1/16 + 1/4

  //Trivial case:
  if ( k == 0 )
    return { 1.0, 0.0 };

  //Eliminate common factors of 2 in fraction k/2^n:
  while ( k%2==0 ) {
    nc_assert( n>=1 );
    n -= 1;
    k /= 2;
  }

  if ( k == 1 ) {
    //Fundamental form.
    nc_assert(n>=1);//k>0, so it follows from k<=n-1 that n>1

    //Cache high-precision numbers (from Python's mpmath module) (for
    //efficiency, accuracy and reproducability):
    constexpr unsigned ncache = 21;
    static std::array<double, ncache> cosvals = { 1.0, // cos(2pi/2^0)
                                                  -1.0, // cos(2pi/2^1)
                                                  0.0, // cos(2pi/2^2)
                                                  0.707106781186547524401, // cos(2pi/2^3)
                                                  0.923879532511286756128, // cos(2pi/2^4)
                                                  0.980785280403230449126, // cos(2pi/2^5)
                                                  0.995184726672196886245, // cos(2pi/2^6)
                                                  0.998795456205172392715, // cos(2pi/2^7)
                                                  0.999698818696204220116, // cos(2pi/2^8)
                                                  0.999924701839144540922, // cos(2pi/2^9)
                                                  0.999981175282601142657, // cos(2pi/2^10)
                                                  0.999995293809576171512, // cos(2pi/2^11)
                                                  0.999998823451701909929, // cos(2pi/2^12)
                                                  0.99999970586288221916, // cos(2pi/2^13)
                                                  0.999999926465717851145, // cos(2pi/2^14)
                                                  0.999999981616429293808, // cos(2pi/2^15)
                                                  0.999999995404107312891, // cos(2pi/2^16)
                                                  0.999999998851026827563, // cos(2pi/2^17)
                                                  0.999999999712756706849, // cos(2pi/2^18)
                                                  0.99999999992818917671, // cos(2pi/2^19)
                                                  0.999999999982047294177 }; // cos(2pi/2^20)

    static std::array<double, ncache> sinvals = { 0.0, // sin(2pi/2^0)
                                                  0.0, // sin(2pi/2^1)
                                                  1.0, // sin(2pi/2^2)
                                                  0.707106781186547524401, // sin(2pi/2^3)
                                                  0.382683432365089771728, // sin(2pi/2^4)
                                                  0.195090322016128267848, // sin(2pi/2^5)
                                                  0.0980171403295606019942, // sin(2pi/2^6)
                                                  0.049067674327418014255, // sin(2pi/2^7)
                                                  0.0245412285229122880317, // sin(2pi/2^8)
                                                  0.0122715382857199260794, // sin(2pi/2^9)
                                                  0.00613588464915447535964, // sin(2pi/2^10)
                                                  0.00306795676296597627015, // sin(2pi/2^11)
                                                  0.0015339801862847656123, // sin(2pi/2^12)
                                                  0.000766990318742704526939, // sin(2pi/2^13)
                                                  0.000383495187571395589072, // sin(2pi/2^14)
                                                  0.00019174759731070330744, // sin(2pi/2^15)
                                                  0.0000958737990959773458705, // sin(2pi/2^16)
                                                  0.000047936899603066884549, // sin(2pi/2^17)
                                                  0.0000239684498084182187292, // sin(2pi/2^18)
                                                  0.0000119842249050697064215, // sin(2pi/2^19)
                                                  0.00000599211245264242784288 }; // sin(2pi/2^20)

    if ( n < ncache ) {
      nc_assert( n < cosvals.size() );
      nc_assert( n < sinvals.size() );
      return { cosvals[n], sinvals[n] };
    }

    //Since n>20 (ncache=21), 2pi/2^n must be less than 3e-6, so we can safely
    //use the fast sincos_mpi256pi256 which requires arguments to be less than
    //pi/256~=1e-2. However, note that we usually do not get here at all, since
    //n=20 covers W tables of a size up to ~1e6.
    double cosval, sinval;
    sincos_mpi256pi256( k2Pi/std::exp2(n) ,cosval,sinval);
    return { cosval, sinval };
  }

  //Non-fundamental form, must combine results from several fundamental forms
  //using multiplication of complex numbers.
  nc_assert(k%2==1);//must be odd at this point
  return detail_mult_imag( calcPhase(1, n), calcPhase(k-1, n) );
}
