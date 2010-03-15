#include <iostream>

#include "xgstats.H"

using namespace std;

int main(int argc, const char *argv[])
{
  cout << log(0.0) << endl;
  test_mixture();
  test_finite_distn();
  test_mixture_b();
  test_uniform_mixture();
  return 0;
}
