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

  class EnergyGridHashPtr final {
    //Class which keeps an energy grid vector in a shared pointer, as well as a
    //hash of it for fast comparisons, and a UID. The energy grid can be unset
    //(nullptr), or have a length >=3 (3 with the special meaning of
    //(emin,emax,npoints).
    EnergyGridHashPtr( std::shared_ptr<const VectD> v,
                       std::size_t hash,
                       UniqueIDValue uid )
      : m_data(std::move(v)), m_hash(hash), m_uid(uid) {}
    friend class EnergyGridHashPtrFact;
  public:
    //NB: Empty initial vector -> null shared ptr:
    const std::shared_ptr<const VectD>& dataShPtr() const { return m_data; }
    std::size_t hash() const { return m_hash; }

    EnergyGridHashPtr clone() const { return { m_data, m_hash, m_uid }; }

    ncconstexpr17 UniqueIDValue getUniqueID() const noexcept { return m_uid; }


  private:
    std::shared_ptr<const VectD> m_data;
    std::size_t m_hash = 0;
    UniqueIDValue m_uid;
  };

  //Returns a VDOSDataHashPtr based on the provided VDOSData. Crucially, the
  //returned object might preferably return a previously returned object to
  //which the argument is identical (except possibly its UniqueID and
  //address). This means that passing two VDOSData objects with the same
  //contents through this function, will *normally* result in the same object
  //being returned both times => a good foundation on which to build caching of
  //derived objects. This might fail if the internal cache size is exceeded.
  VDOSDataHashPtr getCachedVDOSDataHashPtr( VDOSData&& );

  //Same for energy grid:
  EnergyGridHashPtr getCachedEnergyGridHashPtr( VectD&& );

  //Specialisation of DI_VDOS which supports keeping and accessing the VDOSData
  //as a VDOSDataHashPtr, and the energy grid as an EnergyGridHashPtr:
  class DI_VDOSShPtr : public DI_VDOS {
  public:
    virtual const VDOSDataHashPtr& vdosDataHashPtr() const = 0;
    virtual const EnergyGridHashPtr& energyGridHashPtr() const = 0;
    using DI_VDOS::DI_VDOS;
    virtual ~DI_VDOSShPtr();
  };

}
#endif
