////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/core/NCImmutBuf.hh"
#include "NCrystal/core/NCSmallVector.hh"
#include <iostream>
namespace NC = NCrystal;

template<class TImmutableBuffer>
void inspectBufferType() {
  using B = TImmutableBuffer;
  static_assert( std::is_same<
                 B,
                 NC::ImmutableBuffer<B::buffer_local_requested_size,B::buffer_alignment,typename B::metadata_type>
                 >::value, "" );
  std::cout<<"  ImmutableBuffer<"<<B::buffer_local_requested_size<<","<<B::buffer_alignment;
  if ( B::has_metadata )
    std::cout<< ",TMetaData(size="<<sizeof(typename B::metadata_type)<<",align="<<alignof(typename B::metadata_type)<<")";
  std::cout<<">\n";
  static_assert(sizeof(B)==B::object_size,"");
  static_assert(alignof(B)==B::object_alignment,"");
  std::cout<<"     |-object size: "<<sizeof(B)<<"\n";
  std::cout<<"     |-object_alignment: "<<B::object_alignment<<"\n";
  std::cout<<"     |-unused_trailing_bytes: "<<B::unused_trailing_bytes<<"\n";
  std::cout<<"     |\n";
  //  std::cout<<"     |-_detail_min_total: "<<B::_detail_min_total<<"\n";
  std::cout<<"     |-metadata_size: "<<B::metadata_size<<"\n";

  std::cout<<"     |\n";
  std::cout<<"     |-buffer_local_size: "<<B::buffer_local_size<<"\n";
  std::cout<<"     |-buffer_alignment: "<<B::buffer_alignment<<"\n";
#if 0
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
  std::cout<<"     |-: "<<B::<<"\n";
#endif
  std::cout<<"     \\----------------------------------------------\n";
}


template<unsigned LOCALBUF_MINSIZE, unsigned BUF_ALIGN, class TMetaData = NC::NullOptType>
void test_specific_immutbuf()
{
  //Test data. Note the null char inside:
  constexpr const char data1[] = "lksdjvn09\0jv;lskdmf[0j20[mi;kmv;lsdfp2ij[nw;lksdfsd;lfjsldkflskjdfljsdl;fj;sldjf;jsd;fj[p0439f0004jfjweljflksdjf-94-934j9gsig=s=0ifs=-if=2k-j4wpfndslkmv;smc=k2=k=2k=dk3=kd=k=kcds[m[k=0i=k=i=skv[om2[pk[okq";

  std::cout<<"Testing:\n";
  using B = NC::ImmutableBuffer<LOCALBUF_MINSIZE,BUF_ALIGN,TMetaData>;
  inspectBufferType<B>();

  auto testdata = [](const char* thedata, std::size_t len) {
    std::cout<<"  ==> Testing with data of length "<<len<<std::endl;
    B b(thedata,len,TMetaData());
    static_assert(B::buffer_local_size >= LOCALBUF_MINSIZE,"");
    static_assert(alignof(B)%BUF_ALIGN == 0,"");
    static_assert(sizeof(B)%BUF_ALIGN == 0,"");
    static_assert(sizeof(B)%alignof(std::shared_ptr<void>) == 0,"");
    nc_assert_always(0==std::memcmp(thedata,b.data(),len));
    B b2(b);//copy construct
    nc_assert_always(0==std::memcmp(thedata,b.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b2.data(),len));
    nc_assert_always( ( len > B::buffer_local_size ) == ( b.data() == b2.data() ) );//same address if and only if remote
    B b3(std::move(b2));//move construct [ b2 object is now moved-from! ]
    nc_assert_always(0==std::memcmp(thedata,b.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b3.data(),len));
    nc_assert_always( ( len > B::buffer_local_size ) == ( b.data() == b3.data() ) );//same address if and only if remote
    const char* dummy = "a";
    B b4(dummy,1,TMetaData());
    b4 = b;//copy assign over previous local mode
    nc_assert_always(0==std::memcmp(thedata,b.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b3.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b4.data(),len));
    nc_assert_always( ( len > B::buffer_local_size ) == ( b.data() == b4.data() ) );//same address if and only if remote
    B b5(dummy,1,TMetaData());
    b5 = std::move(b4);//move assign over previous local mode [ b4 object is now moved-from! ]
    nc_assert_always(0==std::memcmp(thedata,b.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b3.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b5.data(),len));
    nc_assert_always( ( len > B::buffer_local_size ) == ( b.data() == b5.data() ) );//same address if and only if remote
    std::vector<char> dummy_large;
    dummy_large.resize(1234567,'!');
    B b4alt(&dummy_large[0],dummy_large.size(),TMetaData());
    b4alt = b;//copy assign over previous remote mode
    nc_assert_always(0==std::memcmp(thedata,b.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b3.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b5.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b4alt.data(),len));
    nc_assert_always( ( len > B::buffer_local_size ) == ( b.data() == b4alt.data() ) );//same address if and only if remote
    B b5alt(&dummy_large[0],dummy_large.size(),TMetaData());
    b5alt = std::move(b4alt);//move assign over previous remote mode [ b4alt object is now moved-from! ]
    nc_assert_always(0==std::memcmp(thedata,b.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b3.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b5.data(),len));
    nc_assert_always(0==std::memcmp(thedata,b5alt.data(),len));
    nc_assert_always( ( len > B::buffer_local_size ) == ( b.data() == b5alt.data() ) );//same address if and only if remote
  };

  testdata( data1, 1 );
  testdata( data1, ( LOCALBUF_MINSIZE > 2 ? LOCALBUF_MINSIZE-2 : 1 ) );
  testdata( data1, ( LOCALBUF_MINSIZE > 1 ? LOCALBUF_MINSIZE-1 : 1 ) );
  testdata( data1, LOCALBUF_MINSIZE );
  testdata( data1, LOCALBUF_MINSIZE+1 );
  testdata( data1, sizeof(data1) );
  testdata( data1, strlen(data1) );
  testdata( data1, B::buffer_local_size-1 );
  testdata( data1, B::buffer_local_size );
  testdata( data1, B::buffer_local_size+1 );
  std::vector<char> v;
  v.resize(2000000,char(17));
  testdata( &v[0], v.size() );
}

void test_immutbuf()
{
  test_specific_immutbuf<1,1>();
  test_specific_immutbuf<1,32>();
  test_specific_immutbuf<1,1>();
  test_specific_immutbuf<1,8>();
  test_specific_immutbuf<16,32>();
  test_specific_immutbuf<16,1>();
  //test_specific_immutbuf<16,3>();
  test_specific_immutbuf<16,8>();
  test_specific_immutbuf<15,32>();
  test_specific_immutbuf<15,1>();
  test_specific_immutbuf<15,8>();

  test_specific_immutbuf<18,8>();
  test_specific_immutbuf<18,8,unsigned>();
  test_specific_immutbuf<27,8,uint32_t>();
  test_specific_immutbuf<19,8,uint32_t>();
  test_specific_immutbuf<1,1,uint32_t>();
  test_specific_immutbuf<1,1,uint32_t>();

  test_specific_immutbuf<19,8,uint64_t>();
  // test_specific_immutbuf<27,8,uint_fast8_t>();
  // test_specific_immutbuf<27,8,uint_fast16_t>();
  //  test_specific_immutbuf<27,8,uint_fast64_t>();



#if 0
  using MDType = uint32_t;
  static constexpr unsigned LBSIZE = 17;
  static constexpr unsigned BA = 8;
  using BM = NC::ImmutableBufferWithMetaData< uint32_t, LBSIZE, BA>;
  static_assert( alignof(BM::buffer_type)%BA == 0,"");
  static_assert( BM::buffer_type::buffer_local_size>=LBSIZE ,"");
  std::cout<<"Requested min local buffer size: "<<LBSIZE<<std::endl;
  std::cout<<"Requested min alignment: "<<BA<<std::endl;
  std::cout<<"alignof(BM): "<<alignof(BM)<<std::endl;
  std::cout<<"alignof(BM::buffer_type): "<<alignof(BM::buffer_type)<<std::endl;
  std::cout<<"alignof(MDType): "<<alignof(MDType)<<std::endl;
  std::cout<<"sizeof(BM): "<<sizeof(BM)<<std::endl;
  std::cout<<"BM::buffer_type::object_size: "<<BM::buffer_type::object_size<<std::endl;
  std::cout<<"local buffer useful size: "<<BM::buffer_type::buffer_local_size<<std::endl;
  std::cout<<"sizeof(MDType): "<<sizeof(MDType)<<std::endl;
  std::cout<<"internal overhead in BM::buffer_type (ideally 1): "<<sizeof(BM::buffer_type[100])/100-BM::buffer_type::buffer_local_size<<std::endl;
  std::cout<<"Padding in BM: "<<sizeof(BM)-sizeof(BM::buffer_type)-sizeof(MDType)<<std::endl;
  std::cout<<"Total overhead (ideally 1): "<<sizeof(BM[100])/100-BM::buffer_type::buffer_local_size-sizeof(MDType)<<std::endl;
  struct Dummy {
    alignas(8) char dummy[24];
  };
  static_assert(alignof(char[24])==1,"");
  static_assert(alignof(Dummy)==8,"");
  struct Dummy2 {
    Dummy dummy;
    uint32_t dummy2;
  };
  //  static_assert(alignof(char[24])==1,"");
  static_assert(alignof(Dummy2)==8,"");
  static_assert(sizeof(Dummy2)==32,"");
  static_assert(sizeof(Dummy2)==32||sizeof(Dummy2)==28,"");
#endif
}

int main() {
  using B = NC::ImmutableBuffer< 23, 16, uint32_t>;
  inspectBufferType<B>();
  std::cout<<sizeof(B)<<std::endl;
  std::cout<<sizeof(NC::SmallVector<B,8,NC::SVMode::FASTACCESS>)<<std::endl;
  test_immutbuf();
  return 0;
}
