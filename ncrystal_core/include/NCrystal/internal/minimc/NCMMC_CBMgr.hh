#ifndef NCrystal_MMC_CBMgr_hh
#define NCrystal_MMC_CBMgr_hh

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

#include "NCrystal/internal/minimc/NCMMC_Basket.hh"

#ifndef NCRYSTAL_DISABLE_THREADS
#  include <thread>
#  include <queue>
#  include <condition_variable>
#endif

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    namespace CB {

      //Data structure and manager related to getting results in the form of
      //lists of event data back via a call-back function. All configuration is
      //handled via an CBMgrInput object.

      class DataArea;
      using CallBackFct = std::function<void(const DataArea&)>;

      struct CBMgrInput final {
        //The actual function that should receive callbacks. Should not be a
        //nullptr.
        CallBackFct callbackfct;

        //Maximum number of neutrons provided in each callback. Note that memory
        //will always be allocated for this number in advance, so it is
        //important that it is not too large (values smaller than basket_N are
        //rounded up to basket_N):
        std::size_t max_neutrons_per_callback = 262144;

        //Maximum number of caches to keep around. Too small, and threads might
        //stall while waiting for an available cache. Too big, and the memory
        //usage might get too large. The actual number of caches used, though,
        //will not be above the number of threads used for the simulation:
        std::size_t max_cache_areas = 8;//fixme: test with a value of 1

        //About memory usage: Total number of values per neutron is mostly
        //likely around 100bytes, and more specifically from 69 to 160 bytes. So
        //the default values above will lead to a maximum memory overhead of the
        //callback caches between 0.14-0.42 gigabyte.
      };

      class DataArea final : NoCopyMove {
        //Data area where neutrons are collated into long arrays of double
        //precision values.

        //Fixme: just have one global fieldstype/baskettype in NCDefs.hh
      public:
        enum class FieldsType : unsigned long {
          NOT_INITIALISED = 0,
          BASIC = 1,// [x,y,z,ux,uy,uz,ekin,w,nscat,sawinelas]//fixme: sawinelas -> nscat_inelas?
          WITH_INITIAL = 2,// BASIC + [x0,y0,z0,ux0,uy0,uz0,ekin0,w0]
        };
        constexpr static unsigned nfieldsmax = 18;//BASIC 10, WITH_INITIAL 18
      private:
        std::unique_ptr<double[]> m_memholder;
        std::size_t m_size = 0;
        std::size_t m_capacity;
        FieldsType m_fieldsType = FieldsType::NOT_INITIALISED;
        double * m_datacache[nfieldsmax];
      public:
        DataArea( std::size_t capacity );
        std::size_t capacity() const { return m_capacity; }
        std::size_t size() const { return m_size; }
        const double * const * view_data() const { return m_datacache; }
        unsigned nFields() const;
        unsigned long fieldsTypeForCallBack() const
        {
          nc_assert_always(m_fieldsType != FieldsType::NOT_INITIALISED);
          return static_cast<unsigned long>(m_fieldsType);
        }

        //For internal usage:
        class Mutable;
        friend class Mutable;
      };

      class CBMgr : NoCopyMove {
        //Callback function manager. Class which collects data from worker
        //threads, collates them, and provides them to the provided callback
        //function.
      public:
        CBMgr( CBMgrInput );
        ~CBMgr();

        //Function called by the worker threads:
        void registerData( const BasketView& );

        //Function called after all worker threads have finished, but before
        //CBMgr is destructed:
        void flush();

      private:
        using DataAreaPtr = std::unique_ptr<DataArea>;
#ifndef NCRYSTAL_DISABLE_THREADS
        DataAreaPtr threadAcquireCache();
        void threadReturnCache( DataAreaPtr );
        void fireCallback( const DataArea& );
        std::mutex m_cachemtx;
        std::condition_variable m_cachecondvar;
        std::queue<DataAreaPtr> m_caches;
        std::size_t m_ncaches_created = 0;
        std::size_t m_nmax;
        std::size_t m_nmax_caches;
        std::mutex m_cbmtx;
        CallBackFct m_callback;
#else
        DataAreaPtr m_cache;
        std::size_t m_nmax;
        CallBackFct m_callback;
#endif
      };
    }
  }
}
#endif
