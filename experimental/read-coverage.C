#include <iostream>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>

#include "xgstats.H"
#include "dgrpstats.H"

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

vector<vector<int> > get_test_observations()
{
  vector<vector<int> > observations;
  vector<int> obs;
  /* a heterozygous observation */
  obs.clear();
  obs.push_back(20);
  obs.push_back(20);
  obs.push_back(1);
  obs.push_back(0);
  observations.push_back(obs);
  /* a homozyous non-reference observation */
  obs.clear();
  obs.push_back(2);
  obs.push_back(40);
  obs.push_back(1);
  obs.push_back(0);
  observations.push_back(obs);
  /* a bizarre observation */
  obs.clear();
  obs.push_back(5);
  obs.push_back(25);
  obs.push_back(18);
  obs.push_back(2);
  observations.push_back(obs);
  /* return the observations */
  return observations;
}

void demo_zygosity_models()
{
  // define some model parameters
  int low = 2;
  int med = 20;
  int high = 1000;
  double x = 0.1;
  double y = 0.01;
  double z = .000001;
  double seqerr = .01;
  int nomcoverage = 20;
  int kmulticoverages = 4;
  // create the models
  HMMGarbage garbage(low, med, high);
  HMMRecent recent(x, y, z, seqerr, nomcoverage, kmulticoverages);
  HMMAncient ancient(x, y, z, seqerr, nomcoverage, kmulticoverages);
  // evaluate several observations
  vector<vector<int> > test_obs = get_test_observations();
  vector<vector<int> >::iterator test_it;
  for (test_it = test_obs.begin(); test_it != test_obs.end(); ++test_it)
  {
    vector<int> obs = *test_it;
    // show the observation
    vector<int>::iterator it;
    cout << "observation: ";
    for (it = obs.begin(); it != obs.end(); ++it)
    {
      cout << *it << " ";
    }
    cout << endl;
    // compute a likelihood for each model
    cout << "recent likelihood: " << recent.get_lik(obs) << endl;
    cout << "ancient likelihood: " << ancient.get_lik(obs) << endl;
    cout << "other likelihood: " << garbage.get_lik(obs) << endl;
    cout << endl;
  }
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
  demo_zygosity_models();
  return 0;
}

