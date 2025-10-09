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

#include "NCrystal/internal/extn_utils/NCExtnBC2025.hh"
namespace NCBC = NCrystal::Extn::BC2025;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///                       Test code for BC2025 Recipes                       ///
///           Test function to call: bc2015_test_implementation()            ///
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <stdexcept>

void bc2015_test_implementation()
{
  auto test_primary = []( double x, double sinth, double refval)
  {
    const double val = NCBC::y_primary(x,sinth);
    const double val_lux = NCBC::y_primary_lux(x,sinth);
    if (!(std::fabs(val-refval) <= 0.001*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_primary("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
    if (!(std::fabs(val_lux-refval) <= 1e-06*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_primary_lux("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
  };

  auto test_scndgauss = []( double x, double sinth, double refval)
  {
    const double val = NCBC::y_scndgauss(x,sinth);
    const double val_lux = NCBC::y_scndgauss_lux(x,sinth);
    if (!(std::fabs(val-refval) <= 0.001*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_scndgauss("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
    if (!(std::fabs(val_lux-refval) <= 1e-06*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_scndgauss_lux("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
  };

  auto test_scndlorentz = []( double x, double sinth, double refval)
  {
    const double val = NCBC::y_scndlorentz(x,sinth);
    const double val_lux = NCBC::y_scndlorentz_lux(x,sinth);
    if (!(std::fabs(val-refval) <= 0.001*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_scndlorentz("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
    if (!(std::fabs(val_lux-refval) <= 1e-06*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_scndlorentz_lux("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
  };

  auto test_scndfresnel = []( double x, double sinth, double refval)
  {
    const double val = NCBC::y_scndfresnel(x,sinth);
    const double val_lux = NCBC::y_scndfresnel_lux(x,sinth);
    if (!(std::fabs(val-refval) <= 0.001*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_scndfresnel("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
    if (!(std::fabs(val_lux-refval) <= 1e-06*std::fmin(refval,1.0-refval))) {
      std::ostringstream ss;
      ss << "bc2015_test_implementation: "
         << "Failure in bc2025_y_scndfresnel_lux("<<x<<","<<sinth<<")\n";
      throw std::runtime_error(ss.str());
    }
  };

  std::cout << "Testing bc2025_y_primary + bc2025_y_primary_lux"
            << std::endl;
  test_primary(1e-20,0,1);
  test_primary(1e-20,0.37,1);
  test_primary(1e-20,0.7071,1);
  test_primary(1e-20,0.93,1);
  test_primary(1e-20,1,1);
  test_primary(1e-13,0,0.999999999999905742);
  test_primary(1e-13,0.37,0.999999999999905742);
  test_primary(1e-13,0.7071,0.999999999999905742);
  test_primary(1e-13,0.93,0.999999999999905742);
  test_primary(1e-13,1,0.999999999999905742);
  test_primary(1e-12,0,0.999999999999057088);
  test_primary(1e-12,0.37,0.999999999999057088);
  test_primary(1e-12,0.7071,0.999999999999057088);
  test_primary(1e-12,0.93,0.999999999999057088);
  test_primary(1e-12,1,0.999999999999057088);
  test_primary(0.0001,0,0.999905722489449977);
  test_primary(0.0001,0.37,0.999905722830989552);
  test_primary(0.0001,0.7071,0.999905724213736913);
  test_primary(0.0001,0.93,0.999905725910026466);
  test_primary(0.0001,1,0.99990572659042809);
  test_primary(0.001,0,0.999057962697086133);
  test_primary(0.001,0.37,0.999057996797014347);
  test_primary(0.001,0.7071,0.999058134748896043);
  test_primary(0.001,0.93,0.99905830389750816);
  test_primary(0.001,1,0.999058371729258488);
  test_primary(0.009,0,0.991580310627700912);
  test_primary(0.009,0.37,0.991583034240936856);
  test_primary(0.009,0.7071,0.991593979472062736);
  test_primary(0.009,0.93,0.991607341081938376);
  test_primary(0.009,1,0.991612688397355946);
  test_primary(0.01,0,0.990652882165304338);
  test_primary(0.01,0.37,0.990656238776665754);
  test_primary(0.01,0.7071,0.990669716639168807);
  test_primary(0.01,0.93,0.990686161056473114);
  test_primary(0.01,1,0.99069274044635991);
  test_primary(0.011,0,0.989727059403622267);
  test_primary(0.011,0.37,0.989731113812797214);
  test_primary(0.011,0.7071,0.989747380102478536);
  test_primary(0.011,0.93,0.989767215916285248);
  test_primary(0.011,1,0.989775150196034903);
  test_primary(0.09,0,0.921378643745266457);
  test_primary(0.09,0.37,0.921615887415132673);
  test_primary(0.09,0.7071,0.922510563580389631);
  test_primary(0.09,0.93,0.923558752275423278);
  test_primary(0.09,1,0.923970372765745229);
  test_primary(0.1,0,0.913359883597570255);
  test_primary(0.1,0.37,0.913647965770169601);
  test_primary(0.1,0.7071,0.914726352615906402);
  test_primary(0.1,0.93,0.915983847452636613);
  test_primary(0.1,1,0.916476606013764128);
  test_primary(0.11,0,0.905473627109827595);
  test_primary(0.11,0.37,0.905816514756760061);
  test_primary(0.11,0.7071,0.907090710502909459);
  test_primary(0.11,0.93,0.908569675300849533);
  test_primary(0.11,1,0.909148001400050698);
  test_primary(0.5,0,0.677397998871731022);
  test_primary(0.5,0.37,0.681375531144782887);
  test_primary(0.5,0.7071,0.693094594967365474);
  test_primary(0.5,0.93,0.70493150458645526);
  test_primary(0.5,1,0.709285545104949011);
  test_primary(1.5,0,0.427059012343876054);
  test_primary(1.5,0.37,0.438829477369566623);
  test_primary(1.5,0.7071,0.463370141438253225);
  test_primary(1.5,0.93,0.484182069478211519);
  test_primary(1.5,1,0.491318844971881885);
  test_primary(4,0,0.256478305440050736);
  test_primary(4,0.37,0.271450489435057285);
  test_primary(4,0.7071,0.294058279041176651);
  test_primary(4,0.93,0.311164652224688087);
  test_primary(4,1,0.316803783227092128);
  test_primary(10,0,0.160692615730259891);
  test_primary(10,0.37,0.172334879934710278);
  test_primary(10,0.7071,0.187275015600959244);
  test_primary(10,0.93,0.19822718194603256);
  test_primary(10,1,0.201806603479064706);
  test_primary(30,0,0.0927496567337003086);
  test_primary(30,0.37,0.0995863832091555823);
  test_primary(30,0.7071,0.108079350860056767);
  test_primary(30,0.93,0.114326427593817989);
  test_primary(30,1,0.116371315352427424);
  test_primary(999,0,0.0160653610398071058);
  test_primary(999,0.37,0.0172401058436648311);
  test_primary(999,0.7071,0.0187047403164292345);
  test_primary(999,0.93,0.0197827741536802069);
  test_primary(999,1,0.0201357305432466258);
  test_primary(1000,0,0.0160573261407892312);
  test_primary(1000,0.37,0.0172314830286544526);
  test_primary(1000,0.7071,0.0186953847519196886);
  test_primary(1000,0.93,0.0197728792897164975);
  test_primary(1000,1,0.0201256591124670571);
  test_primary(1000.1,0,0.0160565233346921486);
  test_primary(1000.1,0.37,0.0172306215191156964);
  test_primary(1000.1,0.7071,0.0186944500527839452);
  test_primary(1000.1,0.93,0.0197718907198941307);
  test_primary(1000.1,1,0.020124652904976368);
  test_primary(10000,0,0.00507777237370555211);
  test_primary(10000,0.37,0.00544907338330845401);
  test_primary(10000,0.7071,0.00591199975492481909);
  test_primary(10000,0.93,0.00625273344550764977);
  test_primary(10000,1,0.00636429222075187481);
  test_primary(1000000000000,0,5.07777237370555256e-07);
  test_primary(1000000000000,0.37,5.44907338330845463e-07);
  test_primary(1000000000000,0.7071,5.91199975492481964e-07);
  test_primary(1000000000000,0.93,6.2527334455076501e-07);
  test_primary(1000000000000,1,6.36429222075187487e-07);
  test_primary(1e+99,0,1.60573261407892309e-50);
  test_primary(1e+99,0.37,1.72314830286544515e-50);
  test_primary(1e+99,0.7071,1.86953847519196873e-50);
  test_primary(1e+99,0.93,1.97728792897164961e-50);
  test_primary(1e+99,1,2.01256591124670573e-50);
  test_primary(1e+200,0,5.0777723737055523e-101);
  test_primary(1e+200,0.37,5.44907338330845441e-101);
  test_primary(1e+200,0.7071,5.91199975492481873e-101);
  test_primary(1e+200,0.93,6.25273344550764931e-101);
  test_primary(1e+200,1,6.36429222075187476e-101);

  std::cout << "Testing bc2025_y_scndgauss + bc2025_y_scndgauss_lux"
            << std::endl;
  test_scndgauss(1e-20,0,1);
  test_scndgauss(1e-20,0.37,1);
  test_scndgauss(1e-20,0.7071,1);
  test_scndgauss(1e-20,0.93,1);
  test_scndgauss(1e-20,1,1);
  test_scndgauss(1e-13,0,0.999999999999893974);
  test_scndgauss(1e-13,0.37,0.999999999999893974);
  test_scndgauss(1e-13,0.7071,0.999999999999893974);
  test_scndgauss(1e-13,0.93,0.999999999999893974);
  test_scndgauss(1e-13,1,0.999999999999893974);
  test_scndgauss(1e-12,0,0.999999999998939293);
  test_scndgauss(1e-12,0.37,0.999999999998939293);
  test_scndgauss(1e-12,0.7071,0.999999999998939293);
  test_scndgauss(1e-12,0.93,0.999999999998939293);
  test_scndgauss(1e-12,1,0.999999999998939293);
  test_scndgauss(0.0001,0,0.999893943219759662);
  test_scndgauss(0.0001,0.37,0.999893943604313717);
  test_scndgauss(0.0001,0.7071,0.999893945161208886);
  test_scndgauss(0.0001,0.93,0.999893947071134992);
  test_scndgauss(0.0001,1,0.999893947837228736);
  test_scndgauss(0.001,0,0.998940262922392908);
  test_scndgauss(0.001,0.37,0.998940301317097323);
  test_scndgauss(0.001,0.7071,0.998940456643843677);
  test_scndgauss(0.001,0.93,0.998940647096726808);
  test_scndgauss(0.001,1,0.998940723471927972);
  test_scndgauss(0.009,0,0.990528399718734187);
  test_scndgauss(0.009,0.37,0.990531466451105014);
  test_scndgauss(0.009,0.7071,0.99054379074675869);
  test_scndgauss(0.009,0.93,0.990558836030893675);
  test_scndgauss(0.009,1,0.990564857181265301);
  test_scndgauss(0.01,0,0.98948511172574416);
  test_scndgauss(0.01,0.37,0.9894888912145966);
  test_scndgauss(0.01,0.7071,0.989504067325894687);
  test_scndgauss(0.01,0.93,0.989522583996211957);
  test_scndgauss(0.01,1,0.989529992530914226);
  test_scndgauss(0.011,0,0.98844363174074501);
  test_scndgauss(0.011,0.37,0.988448196954850977);
  test_scndgauss(0.011,0.7071,0.988466512940368802);
  test_scndgauss(0.011,0.93,0.988488848519800345);
  test_scndgauss(0.011,1,0.988497782751834531);
  test_scndgauss(0.09,0,0.911562637444368273);
  test_scndgauss(0.09,0.37,0.91182984491342689);
  test_scndgauss(0.09,0.7071,0.912837666638713241);
  test_scndgauss(0.09,0.93,0.914018527518092538);
  test_scndgauss(0.09,1,0.914482268329535386);
  test_scndgauss(0.1,0,0.902543729477368939);
  test_scndgauss(0.1,0.37,0.902868207385230104);
  test_scndgauss(0.1,0.7071,0.904083032834019806);
  test_scndgauss(0.1,0.93,0.90549977608786314);
  test_scndgauss(0.1,1,0.906054964136642105);
  test_scndgauss(0.11,0,0.893674079978205804);
  test_scndgauss(0.11,0.37,0.894060300442503064);
  test_scndgauss(0.11,0.7071,0.895495779046981633);
  test_scndgauss(0.11,0.93,0.897162136668675281);
  test_scndgauss(0.11,1,0.897813774121218655);
  test_scndgauss(0.5,0,0.637345276720855503);
  test_scndgauss(0.5,0.37,0.641830710434428853);
  test_scndgauss(0.5,0.7071,0.655054692256161419);
  test_scndgauss(0.5,0.93,0.668416681710815808);
  test_scndgauss(0.5,1,0.673332496976859263);
  test_scndgauss(1.5,0,0.357597444689184474);
  test_scndgauss(1.5,0.37,0.370896435653487344);
  test_scndgauss(1.5,0.7071,0.39864416210881054);
  test_scndgauss(1.5,0.93,0.422176335657189539);
  test_scndgauss(1.5,1,0.430245061090626379);
  test_scndgauss(4,0,0.17297395239130367);
  test_scndgauss(4,0.37,0.189879012656458984);
  test_scndgauss(4,0.7071,0.215258501298703303);
  test_scndgauss(4,0.93,0.234348559295702613);
  test_scndgauss(4,1,0.240621677545467538);
  test_scndgauss(10,0,0.0807180267555285547);
  test_scndgauss(10,0.37,0.0935745735759661579);
  test_scndgauss(10,0.7071,0.109430399854542268);
  test_scndgauss(10,0.93,0.120737870379589413);
  test_scndgauss(10,1,0.12438630443554409);
  test_scndgauss(30,0,0.0308030026998234036);
  test_scndgauss(30,0.37,0.0374331506022029178);
  test_scndgauss(30,0.7071,0.0446630365746027921);
  test_scndgauss(30,0.93,0.0496682757865696753);
  test_scndgauss(30,1,0.0512666946970163212);
  test_scndgauss(999,0,0.00122143170251923856);
  test_scndgauss(999,0.37,0.00157928392333777571);
  test_scndgauss(999,0.7071,0.00192676609247600461);
  test_scndgauss(999,0.93,0.00216076771498466182);
  test_scndgauss(999,1,0.00223476218498307828);
  test_scndgauss(1000,0,0.00122028414965030048);
  test_scndgauss(1000,0.37,0.00157781375102613859);
  test_scndgauss(1000,0.7071,0.00192497830984177045);
  test_scndgauss(1000,0.93,0.00215876529765583443);
  test_scndgauss(1000,1,0.00223269181037356315);
  test_scndgauss(1000.1,0,0.00122017030814190763);
  test_scndgauss(1000.1,0.37,0.00157766655522962587);
  test_scndgauss(1000.1,0.7071,0.00192479872672215134);
  test_scndgauss(1000.1,0.93,0.00215856390431820762);
  test_scndgauss(1000.1,1,0.00223248352035886272);
  test_scndgauss(10000,0,0.000142383928135902927);
  test_scndgauss(10000,0.37,0.00018410082586283293);
  test_scndgauss(10000,0.7071,0.00022460832045571348);
  test_scndgauss(10000,0.93,0.000251886811028229865);
  test_scndgauss(10000,1,0.000260512627627712404);
  test_scndgauss(1000000000000,0,4.8917130164293988e-12);
  test_scndgauss(1000000000000,0.37,6.32493019401069932e-12);
  test_scndgauss(1000000000000,0.7071,7.71659736570022946e-12);
  test_scndgauss(1000000000000,0.93,8.65377158999020815e-12);
  test_scndgauss(1000000000000,1,8.95011837497809759e-12);
  test_scndgauss(1e+99,0,3.29959753303634066e-93);
  test_scndgauss(1e+99,0.37,4.26634268909302845e-93);
  test_scndgauss(1e+99,0.7071,5.2050612016247022e-93);
  test_scndgauss(1e+99,0.93,5.83721148274434557e-93);
  test_scndgauss(1e+99,1,6.03710569513684424e-93);
  test_scndgauss(1e+200,0,1.92957191898721976e-187);
  test_scndgauss(1e+200,0.37,2.49491490014387209e-187);
  test_scndgauss(1e+200,0.7071,3.04386817807524513e-187);
  test_scndgauss(1e+200,0.93,3.4135433941631538e-187);
  test_scndgauss(1e+200,1,3.53043954744815707e-187);

  std::cout << "Testing bc2025_y_scndlorentz + bc2025_y_scndlorentz_lux"
            << std::endl;
  test_scndlorentz(1e-20,0,1);
  test_scndlorentz(1e-20,0.37,1);
  test_scndlorentz(1e-20,0.7071,1);
  test_scndlorentz(1e-20,0.93,1);
  test_scndlorentz(1e-20,1,1);
  test_scndlorentz(1e-13,0,0.999999999999899969);
  test_scndlorentz(1e-13,0.37,0.999999999999899969);
  test_scndlorentz(1e-13,0.7071,0.999999999999899969);
  test_scndlorentz(1e-13,0.93,0.999999999999899969);
  test_scndlorentz(1e-13,1,0.999999999999899969);
  test_scndlorentz(1e-12,0,0.999999999999000022);
  test_scndlorentz(1e-12,0.37,0.999999999999000022);
  test_scndlorentz(1e-12,0.7071,0.999999999999000022);
  test_scndlorentz(1e-12,0.93,0.999999999999000022);
  test_scndlorentz(1e-12,1,0.999999999999000022);
  test_scndlorentz(0.0001,0,0.999900010665679129);
  test_scndlorentz(0.0001,0.37,0.999900011109701836);
  test_scndlorentz(0.0001,0.7071,0.999900012907317559);
  test_scndlorentz(0.0001,0.93,0.999900015112513807);
  test_scndlorentz(0.0001,1,0.999900015997037594);
  test_scndlorentz(0.001,0,0.999001065679801958);
  test_scndlorentz(0.001,0.37,0.999001109992181391);
  test_scndlorentz(0.001,0.7071,0.999001289216636779);
  test_scndlorentz(0.001,0.93,0.999001508937245819);
  test_scndlorentz(0.001,1,0.999001597042950551);
  test_scndlorentz(0.009,0,0.9910856851515224);
  test_scndlorentz(0.009,0.37,0.991089210630832884);
  test_scndlorentz(0.009,0.7071,0.991103348821331553);
  test_scndlorentz(0.009,0.93,0.99112058472810427);
  test_scndlorentz(0.009,1,0.991127478159255615);
  test_scndlorentz(0.01,0,0.99010568685861311);
  test_scndlorentz(0.01,0.37,0.990110029591393381);
  test_scndlorentz(0.01,0.7071,0.990127426844987224);
  test_scndlorentz(0.01,0.93,0.990148621204925639);
  test_scndlorentz(0.01,1,0.990157095078302407);
  test_scndlorentz(0.011,0,0.989127763578492414);
  test_scndlorentz(0.011,0.37,0.989133006571792373);
  test_scndlorentz(0.011,0.7071,0.989153988224505398);
  test_scndlorentz(0.011,0.93,0.989179531663385569);
  test_scndlorentz(0.011,1,0.989189741100344522);
  test_scndlorentz(0.09,0,0.917968747806152341);
  test_scndlorentz(0.09,0.37,0.918264757046562075);
  test_scndlorentz(0.09,0.7071,0.91936110921278158);
  test_scndlorentz(0.09,0.93,0.920631113580382632);
  test_scndlorentz(0.09,1,0.92112729294202611);
  test_scndlorentz(0.1,0,0.9097528189826527);
  test_scndlorentz(0.1,0.37,0.910110749489652537);
  test_scndlorentz(0.1,0.7071,0.91142447519078873);
  test_scndlorentz(0.1,0.93,0.912937647317605849);
  test_scndlorentz(0.1,1,0.913527312148028403);
  test_scndlorentz(0.11,0,0.901699446207167576);
  test_scndlorentz(0.11,0.37,0.90212370842963574);
  test_scndlorentz(0.11,0.7071,0.903667060105980635);
  test_scndlorentz(0.11,0.93,0.905434839585329443);
  test_scndlorentz(0.11,1,0.906121994198515912);
  test_scndlorentz(0.5,0,0.679484857014460153);
  test_scndlorentz(0.5,0.37,0.683835026258809542);
  test_scndlorentz(0.5,0.7071,0.695952597070977896);
  test_scndlorentz(0.5,0.93,0.707838551395702509);
  test_scndlorentz(0.5,1,0.712159062730127967);
  test_scndlorentz(1.5,0,0.450909553995101353);
  test_scndlorentz(1.5,0.37,0.46204677668140548);
  test_scndlorentz(1.5,0.7071,0.484151009069781701);
  test_scndlorentz(1.5,0.93,0.502678312036372255);
  test_scndlorentz(1.5,1,0.509016994374947407);
  test_scndlorentz(4,0,0.287152109145697998);
  test_scndlorentz(4,0.37,0.300558728640546835);
  test_scndlorentz(4,0.7071,0.321257651180918202);
  test_scndlorentz(4,0.93,0.337268045522034154);
  test_scndlorentz(4,1,0.342597462329693514);
  test_scndlorentz(10,0,0.183933260109467284);
  test_scndlorentz(10,0.37,0.195084924510890434);
  test_scndlorentz(10,0.7071,0.210345782278526683);
  test_scndlorentz(10,0.93,0.221803740614884309);
  test_scndlorentz(10,1,0.225580067220738728);
  test_scndlorentz(30,0,0.106766732862608862);
  test_scndlorentz(30,0.37,0.114087479074691439);
  test_scndlorentz(30,0.7071,0.123518815342400348);
  test_scndlorentz(30,0.93,0.130506961130937377);
  test_scndlorentz(30,1,0.132800000000000001);
  test_scndlorentz(999,0,0.0185489882266850793);
  test_scndlorentz(999,0.37,0.0199022500866932808);
  test_scndlorentz(999,0.7071,0.0215914499946390122);
  test_scndlorentz(999,0.93,0.0228350617321228556);
  test_scndlorentz(999,1,0.0232422609737816636);
  test_scndlorentz(1000,0,0.0185397128633782056);
  test_scndlorentz(1000,0.37,0.0198923006742199558);
  test_scndlorentz(1000,0.7071,0.0215806574956853257);
  test_scndlorentz(1000,0.93,0.0228236482946756088);
  test_scndlorentz(1000,1,0.0232306441938155799);
  test_scndlorentz(1000.1,0,0.0185387859472531669);
  test_scndlorentz(1000.1,0.37,0.0198913061337761578);
  test_scndlorentz(1000.1,0.7071,0.0215795785437312634);
  test_scndlorentz(1000.1,0.93,0.0228225071978424253);
  test_scndlorentz(1000.1,1,0.0232294827487135452);
  test_scndlorentz(10000,0,0.00586277198137972414);
  test_scndlorentz(10000,0.37,0.00629049780314381614);
  test_scndlorentz(10000,0.7071,0.00682440310903509924);
  test_scndlorentz(10000,0.93,0.00721747131257928064);
  test_scndlorentz(10000,1,0.00734617471654232766);
  test_scndlorentz(1000000000000,0,5.86277198137972435e-07);
  test_scndlorentz(1000000000000,0.37,6.29049780314381553e-07);
  test_scndlorentz(1000000000000,0.7071,6.82440310903509955e-07);
  test_scndlorentz(1000000000000,0.93,7.2174713125792808e-07);
  test_scndlorentz(1000000000000,1,7.34617471654232877e-07);
  test_scndlorentz(1e+99,0,1.85397128633782046e-50);
  test_scndlorentz(1e+99,0.37,1.98923006742199557e-50);
  test_scndlorentz(1e+99,0.7071,2.15806574956853264e-50);
  test_scndlorentz(1e+99,0.93,2.28236482946756082e-50);
  test_scndlorentz(1e+99,1,2.32306441938155811e-50);
  test_scndlorentz(1e+200,0,5.86277198137972341e-101);
  test_scndlorentz(1e+200,0.37,6.29049780314381608e-101);
  test_scndlorentz(1e+200,0.7071,6.8244031090350985e-101);
  test_scndlorentz(1e+200,0.93,7.21747131257928036e-101);
  test_scndlorentz(1e+200,1,7.34617471654232763e-101);

  std::cout << "Testing bc2025_y_scndfresnel + bc2025_y_scndfresnel_lux"
            << std::endl;
  test_scndfresnel(1e-20,0,1);
  test_scndfresnel(1e-20,0.37,1);
  test_scndfresnel(1e-20,0.7071,1);
  test_scndfresnel(1e-20,0.93,1);
  test_scndfresnel(1e-20,1,1);
  test_scndfresnel(1e-13,0,0.999999999999899969);
  test_scndfresnel(1e-13,0.37,0.999999999999899969);
  test_scndfresnel(1e-13,0.7071,0.999999999999899969);
  test_scndfresnel(1e-13,0.93,0.999999999999899969);
  test_scndfresnel(1e-13,1,0.999999999999899969);
  test_scndfresnel(1e-12,0,0.999999999999000022);
  test_scndfresnel(1e-12,0.37,0.999999999999000022);
  test_scndfresnel(1e-12,0.7071,0.999999999999000022);
  test_scndfresnel(1e-12,0.93,0.999999999999000022);
  test_scndfresnel(1e-12,1,0.999999999999000022);
  test_scndfresnel(0.0001,0,0.999900008799360918);
  test_scndfresnel(0.0001,0.37,0.999900009165697434);
  test_scndfresnel(0.0001,0.7071,0.999900010648836823);
  test_scndfresnel(0.0001,0.93,0.999900012468281973);
  test_scndfresnel(0.0001,1,0.99990001319808286);
  test_scndfresnel(0.001,0,0.999000879361239824);
  test_scndfresnel(0.001,0.37,0.999000915936696932);
  test_scndfresnel(0.001,0.7071,0.999001063902838493);
  test_scndfresnel(0.001,0.93,0.999001245329921517);
  test_scndfresnel(0.001,1,0.999001318085486423);
  test_scndfresnel(0.009,0,0.991070816626705819);
  test_scndfresnel(0.009,0.37,0.991073737786429154);
  test_scndfresnel(0.009,0.7071,0.991085476503148044);
  test_scndfresnel(0.009,0.93,0.99109980646752649);
  test_scndfresnel(0.009,1,0.991105541260750567);
  test_scndfresnel(0.01,0,0.990087364760898203);
  test_scndfresnel(0.01,0.37,0.990090964803397422);
  test_scndfresnel(0.01,0.7071,0.990105419594910052);
  test_scndfresnel(0.01,0.93,0.990123055546527797);
  test_scndfresnel(0.01,1,0.990130111588998552);
  test_scndfresnel(0.011,0,0.989105635014708051);
  test_scndfresnel(0.011,0.37,0.989109983428291484);
  test_scndfresnel(0.011,0.7071,0.989127428567477374);
  test_scndfresnel(0.011,0.93,0.989148701360755611);
  test_scndfresnel(0.011,1,0.989157210324545733);
  test_scndfresnel(0.09,0,0.916686686396509076);
  test_scndfresnel(0.09,0.37,0.916940982241070679);
  test_scndfresnel(0.09,0.7071,0.917899652712336533);
  test_scndfresnel(0.09,0.93,0.919022575437036138);
  test_scndfresnel(0.09,1,0.9194634997882567);
  test_scndfresnel(0.1,0,0.908198193066665715);
  test_scndfresnel(0.1,0.37,0.908506958401812503);
  test_scndfresnel(0.1,0.7071,0.909662351988239859);
  test_scndfresnel(0.1,0.93,0.911009322804470889);
  test_scndfresnel(0.1,1,0.911537084737879644);
  test_scndfresnel(0.11,0,0.899851686707168641);
  test_scndfresnel(0.11,0.37,0.900219164896514989);
  test_scndfresnel(0.11,0.7071,0.901584198877591336);
  test_scndfresnel(0.11,0.93,0.90316818524855913);
  test_scndfresnel(0.11,1,0.903787502974892076);
  test_scndfresnel(0.5,0,0.659409131196108778);
  test_scndfresnel(0.5,0.37,0.663659802995088977);
  test_scndfresnel(0.5,0.7071,0.676160714686383946);
  test_scndfresnel(0.5,0.93,0.688772024068108335);
  test_scndfresnel(0.5,1,0.69340829954302019);
  test_scndfresnel(1.5,0,0.399391537851936262);
  test_scndfresnel(1.5,0.37,0.411881935845156089);
  test_scndfresnel(1.5,0.7071,0.437771759058127041);
  test_scndfresnel(1.5,0.93,0.459644997052124471);
  test_scndfresnel(1.5,1,0.467132755649424358);
  test_scndfresnel(4,0,0.228676356271608677);
  test_scndfresnel(4,0.37,0.244302243726306212);
  test_scndfresnel(4,0.7071,0.26752883192663951);
  test_scndfresnel(4,0.93,0.284943660965013057);
  test_scndfresnel(4,1,0.290661423460214396);
  test_scndfresnel(10,0,0.13889515089169982);
  test_scndfresnel(10,0.37,0.150636185020350977);
  test_scndfresnel(10,0.7071,0.165214272401134205);
  test_scndfresnel(10,0.93,0.175733804134711902);
  test_scndfresnel(10,1,0.179149979699797152);
  test_scndfresnel(30,0,0.0791555029478179517);
  test_scndfresnel(30,0.37,0.085571907221417845);
  test_scndfresnel(30,0.7071,0.093235749158265338);
  test_scndfresnel(30,0.93,0.0988114821782814789);
  test_scndfresnel(30,1,0.100629854949566708);
  test_scndfresnel(999,0,0.0136388653535980068);
  test_scndfresnel(999,0.37,0.0146400998180710069);
  test_scndfresnel(999,0.7071,0.0158858871922350139);
  test_scndfresnel(999,0.93,0.016802473622760522);
  test_scndfresnel(999,1,0.017102530942043951);
  test_scndfresnel(1000,0,0.0136320418730609214);
  test_scndfresnel(1000,0.37,0.0146327712142333838);
  test_scndfresnel(1000,0.7071,0.0158779327680368988);
  test_scndfresnel(1000,0.93,0.0167940591461938557);
  test_scndfresnel(1000,1,0.0170939659039612904);
  test_scndfresnel(1000.1,0,0.013631360322083166);
  test_scndfresnel(1000.1,0.37,0.0146320396305409921);
  test_scndfresnel(1000.1,0.7071,0.0158771389309357858);
  test_scndfresnel(1000.1,0.93,0.0167932195062090207);
  test_scndfresnel(1000.1,1,0.0170931112697631227);
  test_scndfresnel(10000,0,0.00431083014776604666);
  test_scndfresnel(10000,0.37,0.00462728855171251618);
  test_scndfresnel(10000,0.7071,0.00502104320820185614);
  test_scndfresnel(10000,0.93,0.00531074780615552796);
  test_scndfresnel(10000,1,0.00540558665017767687);
  test_scndfresnel(1000000000000,0,4.3108301477660464e-07);
  test_scndfresnel(1000000000000,0.37,4.6272885517125163e-07);
  test_scndfresnel(1000000000000,0.7071,5.02104320820185633e-07);
  test_scndfresnel(1000000000000,0.93,5.31074780615552798e-07);
  test_scndfresnel(1000000000000,1,5.40558665017767726e-07);
  test_scndfresnel(1e+99,0,1.36320418730609202e-50);
  test_scndfresnel(1e+99,0.37,1.46327712142333832e-50);
  test_scndfresnel(1e+99,0.7071,1.58779327680368993e-50);
  test_scndfresnel(1e+99,0.93,1.67940591461938551e-50);
  test_scndfresnel(1e+99,1,1.70939659039612907e-50);
  test_scndfresnel(1e+200,0,4.31083014776604605e-101);
  test_scndfresnel(1e+200,0.37,4.62728855171251623e-101);
  test_scndfresnel(1e+200,0.7071,5.02104320820185612e-101);
  test_scndfresnel(1e+200,0.93,5.31074780615552798e-101);
  test_scndfresnel(1e+200,1,5.40558665017767687e-101);
  std::cout<<"All tests completed"<<std::endl;
}



int main() {
  bc2015_test_implementation();
  return 0;
}
