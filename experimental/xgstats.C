#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>

#include "xgstats.H"

using namespace std;


double binomial_log_pmf(int observed_n, int max_n, double p_success)
{
  // TODO special cases
  double accum = 0;
  accum += lgamma(max_n + 1);
  accum -= lgamma(observed_n + 1);
  accum -= lgamma((max_n - observed_n) + 1);
  accum += observed_n * log(p_success);
  accum += (max_n - observed_n) * log1p(-p_success);
  return accum;
}

/*
 * @param observed_n: the number of completed events
 * @param pr: the probability of quitting
 */
double geometric_log_pmf(int observed_n, double pr)
{
  if (pr == 0)
    return -numeric_limits<double>::infinity();
  if (pr == 1.0)
  {
      if (observed_n)
        return -numeric_limits<double>::infinity();
      else
        return log(pr);
  }
  return observed_n * log1p(-pr) + log(pr);
}

double poisson_log_pmf(int observed_n, int expected_n)
{
  if (expected_n == 0 && observed_n == 0)
    return 0;
  if (expected_n == 0 && observed_n != 0)
    return -numeric_limits<double>::infinity();
  double accum = 0;
  accum += observed_n * log(expected_n);
  accum -= expected_n;
  accum -= lgamma(observed_n + 1);
  return accum;
}

void FiniteDistn::set_distn(const vector<double> &distn)
{
  distn_ = distn;
  log_distn_.clear();
  transform(distn_.begin(), distn_.end(),
      back_inserter(log_distn_), (double (*)(double)) log);
}

void test_uniform_mixture()
{
  // define the mixture components
  vector<Model<int> *> models;
  models.push_back(new FairD6());
  models.push_back(new LoadedD6());
  // create the model
  UniformMixture<int> m(models);
  // show some likelihoods
  int i;
  cout << "uniform mixture likelihoods:" << endl;
  for (i=0; i<6; i++)
  {
    cout << i << endl;
    cout << m.get_lik(i) << endl;
    cout << m.get_log_lik(i) << endl;
    cout << endl;
  }
  vector<double> post = m.get_posterior_distn(5);
  cout << "uniform mixture posterior distribution:" << endl;
  vector<double>::iterator it;
  for (it=post.begin(); it != post.end(); ++it)
    cout << *it << endl;
  cout << endl;
  // cleanup
  for_each(models.begin(), models.end(), delete_object());
}

void test_mixture_b()
{
  // define the mixture parameters
  vector<double> distn;
  distn.push_back(.4);
  distn.push_back(.6);
  // define the mixture components
  vector<Model<int> *> models;
  models.push_back(new FairD6());
  models.push_back(new LoadedD6());
  // create the model
  Mixture<int> m(distn, models);
  // show some likelihoods
  int i;
  cout << "non-uniform mixture likelhoods:" << endl;
  for (i=0; i<6; i++)
  {
    cout << i << endl;
    cout << m.get_lik(i) << endl;
    cout << m.get_log_lik(i) << endl;
    cout << endl;
  }
  vector<double> post = m.get_posterior_distn(5);
  cout << "non-uniform mixture posterior distribution:" << endl;
  vector<double>::iterator it;
  for (it=post.begin(); it != post.end(); ++it)
    cout << *it << endl;
  cout << endl;
  // cleanup
  for_each(models.begin(), models.end(), delete_object());
}

void test_mixture()
{
  /* define a mixture model and an observation */
  cout << "testing mixture..." << endl;
  Mixture<vector<int> > m;
  int arr[3] = {3, 5, 8};
  vector<int> v(arr, arr+3);
  /* compute a posterior distribution */
  cout << "computing posterior distribution..." << endl;
  vector<double> post = m.get_posterior_distn(v);
  cout << "posterior distribution:" << endl;
  vector<double>::iterator it;
  for (it=post.begin(); it != post.end(); ++it)
    cout << *it << endl;
  cout << endl;
}

void test_finite_distn()
{
  cout << "finite distribution:" << endl;
  vector<double> v(6, 1.0/6.0);
  FiniteDistn d(v);
  int i;
  for (i=0; i<6; i++)
  {
    cout << i << endl;
    cout << d.get_lik(i) << endl;
    cout << d.get_log_lik(i) << endl;
    cout << endl;
  }
}
