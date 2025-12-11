#include "NCrystal/ncrystal.h"
#include <iostream>
#include <fstream>

int main( int argc, char** argv )
{
  if ( argc!=6 ) {
    std::cout
      <<"Please provide args: [cfgstr] [geomcfg] [srccfg] [enginecfg] [outfile]"
      <<std::endl;
    std::cout
      <<"(outfile should be a json file, \"none\" or \"stdout\")."
      <<std::endl;
    return 1;
  }
  const std::string dest(argv[5]);

  std::string res;
  {
    //Invoke MiniMC:
    char* res_raw = ncrystal_minimc( argv[1], argv[2], argv[3], argv[4] );
    res = std::string(res_raw);
    ncrystal_dealloc_string(res_raw);
  }

  if ( dest == "none " ) {
    //nothing
  } else if ( dest == "stdout" ) {
    std::cout<<res<<std::endl;
  } else {
    {
      std::ofstream fh(dest);
      fh << res;
      fh.close();
    }
    std::cout<<"Wrote: "<<dest<<std::endl;
  }
  return 0;
}
