#include <cmath>
#include <limits>
#include <iostream>

void test( double * __restrict a,
           const double * __restrict b,
           std::size_t n )
{
  for ( std::size_t i = 0; i < n*128; ++i ) {
    a[i]=(std::max<double>
          (0.0,1.0-std::max<double>
           (0.0,b[i]-std::numeric_limits<double>::max())));
  }

}


int main()
{
  double a[128];
  double b[128];
  for ( int i =0;i<128;++i)
    b[i]=1e300;
  b[0] = std::numeric_limits<double>::infinity();
  test(a,b,1);
  std::cout<<a[0]<<std::endl;
  std::cout<<a[1]<<std::endl;
  std::cout<<a[2]<<std::endl;
  std::cout<<a[3]<<std::endl;

  return 0;

}
