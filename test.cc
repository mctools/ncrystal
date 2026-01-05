#include <limits>
#include <iostream>

int main()

{
  constexpr double a = std::numeric_limits<double>::denorm_min();
  std::cout<< (1e-15 / std::max<double>(0.0*a,0.0) ) <<std::endl;
  return 0;
}
