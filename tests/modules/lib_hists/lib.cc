
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

#include "NCTestUtils/NCTestModUtils.hh"
#include "NCrystal/internal/utils/NCHists.hh"

NCTEST_CTYPE_DICTIONARY
{
  return
    "uint nctest_hist1d_book( uint, double, double, int, int );"
    "uint nctest_hist1d_fill( uint, double );"
    "uint nctest_hist1d_fillw( uint, double, double );"
    "const char * nctest_hist1d_tojson( uint );"
    ;
}

namespace {
      // AllowWeights ALLOW_WEIGHTS = AllowWeights::YES,
      // OverflowHandling OF_HANDLING = OverflowHandling::Record,

  using AllowWeights = NC::Hists::AllowWeights;
  using OverflowHandling = NC::Hists::OverflowHandling;
  using Hist1D_W = NC::Hists::Hist1D<AllowWeights::YES,OverflowHandling::Clamp>;
  using Hist1D_WO = NC::Hists::Hist1D<AllowWeights::YES,OverflowHandling::Record>;
  using Hist1D = NC::Hists::Hist1D<AllowWeights::NO,OverflowHandling::Clamp>;
  using Hist1D_O = NC::Hists::Hist1D<AllowWeights::NO,OverflowHandling::Record>;

  class GenericHist : NC::NoCopyMove {
  public:
    GenericHist( unsigned nbins, double xmin, double xmax,
                 AllowWeights W, OverflowHandling O )
    {
      if ( W == AllowWeights::YES ) {
        if ( O == OverflowHandling::Record )
          m_h_WO = std::make_unique<Hist1D_WO>(nbins,xmin,xmax);
        else
          m_h_W = std::make_unique<Hist1D_W>(nbins,xmin,xmax);
      } else {
        if ( O == OverflowHandling::Record )
          m_h_O = std::make_unique<Hist1D_O>(nbins,xmin,xmax);
        else
          m_h = std::make_unique<Hist1D>(nbins,xmin,xmax);
      }
    }
    void fill( double x )
    {
      if ( m_h_W != nullptr )
        m_h_W->fill(x);
      else if ( m_h_WO != nullptr )
        m_h_WO->fill(x);
      else if ( m_h_O != nullptr )
        m_h_O->fill(x);
      else if ( m_h != nullptr )
        m_h->fill(x);
      else
        nc_assert_always(false);
    }

    void fill( double x, double w )
    {
      if ( m_h_W != nullptr )
        m_h_W->fill(x,w);
      else if ( m_h_WO != nullptr )
        m_h_WO->fill(x,w);
      else
        nc_assert_always(false);
    }

    std::string toJSON() const
    {
      std::ostringstream os;
      if ( m_h_W != nullptr )
        m_h_W->toJSON(os);
      else if ( m_h_WO != nullptr )
        m_h_WO->toJSON(os);
      else if ( m_h_O != nullptr )
        m_h_O->toJSON(os);
      else if ( m_h != nullptr )
        m_h->toJSON(os);
      else
        nc_assert_always(false);
      return os.str();
    }

  private:
    std::unique_ptr<Hist1D_W> m_h_W;
    std::unique_ptr<Hist1D_WO> m_h_WO;
    std::unique_ptr<Hist1D_O> m_h_O;
    std::unique_ptr<Hist1D> m_h;
  };

  struct HistDB {
    std::mutex mtx;
    std::map<unsigned, std::shared_ptr<GenericHist>> id2hist;
  };
  HistDB& getHistDB()
  {
    static HistDB db;
    return db;
  }
  unsigned addHistToDB( std::shared_ptr<GenericHist> hist )
  {
    auto& db = getHistDB();
    NCRYSTAL_LOCK_GUARD(db.mtx);
    unsigned newid  = db.id2hist.size() + 1;
    db.id2hist[newid] = std::move(hist);
    return newid;
  }
  GenericHist& findHistInDB( unsigned id )
  {
    auto& db = getHistDB();
    NCRYSTAL_LOCK_GUARD(db.mtx);
    auto it = db.id2hist.find(id);
    nc_assert_always( it!=db.id2hist.end() );
    std::shared_ptr<GenericHist>& sp_h = it->second;
    GenericHist * h = sp_h.get();
    nc_assert_always( h != nullptr );
    return *h;
  }
}

NCTEST_CTYPES unsigned nctest_hist1d_book( unsigned nbins,
                                           double xmin,
                                           double xmax,
                                           int allow_weights,
                                           int overflow_record )
{
  unsigned id = 0;
  try {
    auto O = ( overflow_record
               ? OverflowHandling::Record
               : OverflowHandling::Clamp );
    auto W = ( allow_weights
               ? AllowWeights::YES
               : AllowWeights::NO );
    auto h = std::make_shared<GenericHist>( static_cast<std::size_t>( nbins ),
                                            xmin, xmax, W, O );
    id = addHistToDB(std::move(h));
  } NCCATCH;
  return id;
}

NCTEST_CTYPES unsigned nctest_hist1d_fill( unsigned id, double x )
{
  try {
    auto& h = findHistInDB( id );
    h.fill(x);
  } NCCATCH;
  return id;
}

NCTEST_CTYPES unsigned nctest_hist1d_fillw( unsigned id, double x, double w )
{
  try {
    auto& h = findHistInDB( id );
    h.fill(x,w);
  } NCCATCH;
  return id;
}

NCTEST_CTYPES const char * nctest_hist1d_tojson( unsigned id )
{
  static char buf[1048576];//1MB, plenty for tests
  buf[0] = '\0';
  try {
    auto& h = findHistInDB( id );
    std::string res_str = h.toJSON();
    if ( !res_str.empty() )
      std::strncat(buf,res_str.data(),res_str.size());
  } NCCATCH;
  return &buf[0];
}
