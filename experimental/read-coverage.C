#include <iostream>

#include "xgstats.H"

using namespace std;


void demo_poisson_log_pmf()
{
  int observed_n = 60;
  int expected_n = 20;
  cout << poisson_log_pmf(observed_n, expected_n) << endl;
}

void demo_geometric_log_pmf()
{
  int obs = 5;
  double pr = 0.1;
  cout << geometric_log_pmf(obs, pr) << endl;
}


int main(int argc, const char *argv[])
{
  cout << log(0.0) << endl;
  test_mixture();
  test_finite_distn();
  test_mixture_b();
  test_uniform_mixture();
  demo_poisson_log_pmf();
  demo_geometric_log_pmf();
  return 0;
}

