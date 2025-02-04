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

#include "NCrystal/internal/utils/NCHists.hh"

namespace NC = NCrystal;

void NC::Hists::RunningStats1D::merge( const RunningStats1D& o )
{
  if (!o.m_sumw) {
    //other has no data
    return;
  }
  if (!m_sumw) {
    //we don't have data yet, so simply copy over stats:
    m_rms_state = o.m_rms_state;
    m_sumw = o.m_sumw;
    m_sumwx = o.m_sumwx;
    m_minfilled = o.m_minfilled;
    m_maxfilled = o.m_maxfilled;
    return;
  }

  //Both had something, do a proper merge:
#if NCRYSTAL_HIST_ROOT_STYLE_RMS
  m_rms_state += o.m_rms_state;
#else
  const double w1(m_sumw);//s1
  const double w2(o.m_sumw);//s2
  const double wx1(m_sumwx);//b1
  const double wx2(o.m_sumwx);//b2
  const double k(w2*wx1-w1*wx2);
  assert((w1+w2)>0);
  m_rms_state += o.m_rms_state + k*k/(w1*w2*(w1+w2));
#endif
  m_sumw += o.m_sumw;
  m_sumwx += o.m_sumwx;
  m_minfilled = std::min<double>(m_minfilled,o.m_minfilled);
  m_maxfilled = std::max<double>(m_maxfilled,o.m_maxfilled);
}


void NC::Hists::RunningStats1D::registerNValues( double val, std::size_t N )
{
  update_maxmin(val);
#if NCRYSTAL_HIST_ROOT_STYLE_RMS
  m_rms_state += N * val*val;
#else
  const double dN = static_cast<double>(N);
  const double d1 = m_sumw * val - m_sumwx;
  const double d2 = m_sumw * ( dN + m_sumw );
  if ( d2 )
    m_rms_state += dN * ( d1 * d1  / d2 );
#endif
  m_sumw += dN;
  m_sumwx += dN * val;
}

void NC::Hists::RunningStats1D::toJSON( std::ostringstream& os ) const
{
  streamJSONDictEntry(os, "integral", m_sumw,JSONDictPos::FIRST);
  if ( hasData() ) {
    streamJSONDictEntry(os, "rms", calcRMS() );
    streamJSONDictEntry(os, "mean", calcMean() );
    streamJSONDictEntry(os, "minfilled", m_minfilled );
    streamJSONDictEntry(os, "maxfilled", m_maxfilled , JSONDictPos::LAST);
  } else {
    auto none = json_null_t{};
    streamJSONDictEntry(os, "rms", none );
    streamJSONDictEntry(os, "mean", none );
    streamJSONDictEntry(os, "minfilled", none );
    streamJSONDictEntry(os, "maxfilled", none, JSONDictPos::LAST);
  }
}

double NC::Hists::RunningStats1D::calcRMSSq() const
{
  if (!hasData())
    NCRYSTAL_THROW(CalcError,"RMS not well defined in empty histograms");
#ifdef NCRYSTAL_HIST_ROOT_STYLE_RMS
  //Unstable calculation:
  return m_rms_state / m_sumw - (m_sumwx*m_sumwx)/(m_sumw*m_sumw);
#else
  double rms2 = m_rms_state / m_sumw;
  nc_assert( std::isfinite(rms2) && rms2 >= 0.0 );
  return rms2;
#endif
}
