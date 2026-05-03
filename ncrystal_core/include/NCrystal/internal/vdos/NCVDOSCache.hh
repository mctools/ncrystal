#ifndef NCrystal_VDOSCache_hh
#define NCrystal_VDOSCache_hh

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

#include "NCrystal/interfaces/NCInfoTypes.hh"
#include "NCrystal/interfaces/NCSABData.hh"

namespace NCRYSTAL_NAMESPACE {

  //Infrastructure which is meant to facilitate caching outcomes of processing
  //VDOSData and energy grid, by shared storage and consistent UIDs.

  class VDOSDataHashPtr final {
  public:
    //Class which keeps a VDOSData object in a shared pointer, as well as a hash
    //of it for fast comparisons. The UID is contained in the VDOSData object
    //itself.
    VDOSDataHashPtr( shared_obj<const VDOSData> );
    const shared_obj<const VDOSData>& dataSO() const { return m_data; }
    const VDOSData& data() const { return m_data; }
    bool operator<(const VDOSDataHashPtr&o) const
    {
      return ( m_hash != o.m_hash ? m_hash < o.m_hash : *m_data<*o.m_data );
    }
    bool operator==(const VDOSDataHashPtr&o) const
    {
      return m_hash == o.m_hash && *m_data == *o.m_data;
    }
  private:
    shared_obj<const VDOSData> m_data;
    std::size_t m_hash;
  };

  class EnergyGridPtr final : public MoveOnly {
    //Class which keeps an energy grid vector in a shared pointer, as well as a
    //a UID. Both of these are created based on the egrid caches in NCSABFactory.hh.
    //
    //The energy grid can be empty or unset (always kept as a nullptr), or have
    //a length >=3 (3 with the special meaning of (emin,emax,npoints).
  public:
    EnergyGridPtr( std::shared_ptr<const VectD> );
    EnergyGridPtr( const VectD& );
    EnergyGridPtr( NullOptType ) {}

    const std::shared_ptr<const VectD>& dataShPtr() const { return m_data; }
    ncconstexpr17 UniqueIDValue getUniqueID() const noexcept { return m_uid; }

    EnergyGridPtr clone() const {
      EnergyGridPtr v( NullOpt );
      v.m_data = m_data;
      v.m_uid = m_uid;
      return v;
    }

  private:
    std::shared_ptr<const VectD> m_data = nullptr;
    UniqueIDValue m_uid = {0};
  };

  //Returns a VDOSDataHashPtr based on the provided VDOSData. Crucially, the
  //returned object might preferably return a previously returned object to
  //which the argument is identical (except possibly its UniqueID and
  //address). This means that passing two VDOSData objects with the same
  //contents through this function, will *normally* result in the same object
  //being returned both times => a good foundation on which to build caching of
  //derived objects. This might fail if the internal cache size is exceeded.
  VDOSDataHashPtr getCachedVDOSDataHashPtr( VDOSData&& );

  //Specialisation of DI_VDOS which supports keeping and accessing the VDOSData
  //as a VDOSDataHashPtr, and the energy grid as an EnergyGridPtr:
  class DI_VDOSShPtr : public DI_VDOS {
  public:
    virtual const VDOSDataHashPtr& vdosDataHashPtr() const = 0;
    virtual const EnergyGridPtr& energyGridPtr() const = 0;
    using DI_VDOS::DI_VDOS;
    virtual ~DI_VDOSShPtr();
  };

}
#endif
