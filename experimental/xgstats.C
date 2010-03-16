#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>

#include "xgstats.H"

using namespace std;

template <typename T>
double logsumexp(T begin, T end)
{
  // logarithm of empty sum is -inf
  if (begin == end)
    return -numeric_limits<double>::infinity();
  // if sum is not empty then do the logsumexp trick
  T max_it = max_element(begin, end);
  double max_value = *max_it;
  double accum = 0;
  T it;
  for (it=begin; it!=end; it++)
  {
    if (it != max_it)
    {
      accum += exp(*it - max_value);
    }
  }
  return log1p(accum) + max_value;
}


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

/*
 * This should be in scipy.stats but it isn't.
 * @param distribution: the distribution over classes
 * @param counts: the observed counts over classes
 */
template <typename T_distn, typename T_counts>
double multinomial_log_pmf(T_distn distn_begin, T_distn distn_end,
    T_counts counts_begin, T_counts counts_end)
{
  // handle special cases
  T_distn d_it = distn_begin;
  T_counts c_it = counts_begin;
  for (; d_it != distn_end; ++d_it, ++c_it)
    if (*d_it == 0 && *c_it != 0)
      return -numeric_limits<double>::infinity();
  // handle less unusual cases
  int n = accumulate(counts_begin, counts_end, 0);
  // initialize the log probability mass
  int accum = 0;
  // add the contribution of n to the multinomial coefficient
  if (n > 1)
      accum += lgamma(n + 1);
  // add the contribution of the counts to the multinomial coefficient
  for (c_it = counts_begin; c_it != counts_end; ++c_it)
  {
    if (*c_it > 1)
      accum -= lgamma(*c_it + 1);
  }
  // add the contribution of probabilities
  d_it = distn_begin;
  c_it = counts_begin;
  for (; d_it != distn_end; ++d_it, ++c_it)
  {
    if (*c_it > 0)
      accum += *c_it * log(*d_it);
  }
  return accum;
}


template<typename T>
double UniformMixture<T>::get_log_lik(const T& obs) const
{
  if (!nmodels_)
    return numeric_limits<double>::signaling_NaN();
  // get the log likelihoods
  int i;
  vector<double> log_liks;
  for (i=0; i<nmodels_; i++)
  {
    log_liks.push_back(models_[i]->get_log_lik(obs));
  }
  // deal with an unexplainable observation
  vector<double>::iterator it = max_element(log_liks.begin(), log_liks.end());
  if (*it == -numeric_limits<double>::infinity())
    return -numeric_limits<double>::infinity();
  // return the log likelihood
  return logsumexp(log_liks.begin(), log_liks.end()) - log_nmodels_;
}

template<typename T>
vector<double> UniformMixture<T>::get_posterior_distn(const T& obs)
{
  if (!nmodels_)
    return vector<double>();
  int i;
  // get the log likelihoods
  vector<double> post;
  for (i=0; i<nmodels_; i++)
  {
    post.push_back(models_[i]->get_log_lik(obs));
  }
  // an unexplainable observation is a bad problem
  vector<double>::iterator it = max_element(post.begin(), post.end());
  if (*it == -numeric_limits<double>::infinity())
    return vector<double>(nmodels_, numeric_limits<double>::signaling_NaN());
  // compute the posterior distribution
  double log_lik_sum = logsumexp(post.begin(), post.end());
  for (i=0; i<nmodels_; i++)
  {
    post[i] = exp(post[i] - log_lik_sum);
  }
  return post;
}

void FiniteDistn::set_distn(const vector<double> &distn)
{
  distn_ = distn;
  log_distn_.clear();
  transform(distn_.begin(), distn_.end(),
      back_inserter(log_distn_), (double (*)(double)) log);
}

template<typename T>
Mixture<T>::Mixture(const vector<double> &distn,
    const vector<Model<T> *> &models) : Model<T>(),
    distn_(distn), models_(models)
{
  update_info();
}

/*
 * Compute the posterior distribution.
 * This is proportional to the vector of elementwise products
 * of priors with likelihoods.
 * When likelihoods are tiny, this normalization factor
 * cannot be accurately represented with double precision.
 * Therefore we do tricksy things in log space.
 */
template<typename T>
vector<double> Mixture<T>::get_posterior_distn(const T& obs)
{
  if (!nmodels_)
    return vector<double>();
  int i;
  // get the log likelihoods
  vector<double> post;
  for (i=0; i<nmodels_; i++)
  {
    post.push_back(models_[i]->get_log_lik(obs));
  }
  // an unexplainable observation is a bad problem
  vector<double>::iterator it = max_element(post.begin(), post.end());
  if (*it == -numeric_limits<double>::infinity())
    return vector<double>(nmodels_, numeric_limits<double>::signaling_NaN());
  // get the weighted log likelihoods
  for (i=0; i<nmodels_; i++)
  {
    post[i] += log_distn_[i];
  }
  // compute the posterior distribution
  double obs_log_lik = logsumexp(post.begin(), post.end());
  for (i=0; i<nmodels_; i++)
  {
    post[i] = exp(post[i] - obs_log_lik);
  }
  return post;
}

template<typename T>
double Mixture<T>::get_log_lik(const T& obs) const
{
  if (!nmodels_)
    return numeric_limits<double>::signaling_NaN();
  // get the log likelihoods
  int i;
  vector<double> post;
  for (i=0; i<nmodels_; i++)
  {
    post.push_back(models_[i]->get_log_lik(obs));
  }
  // deal with an unexplainable observation
  vector<double>::iterator it = max_element(post.begin(), post.end());
  if (*it == -numeric_limits<double>::infinity())
    return -numeric_limits<double>::infinity();
  // get the weighted log likelihoods
  for (i=0; i<nmodels_; i++)
  {
    post[i] += log_distn_[i];
  }
  return logsumexp(post.begin(), post.end());
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
  Mixture<vector<int> > m;
  int arr[3] = {3, 5, 8};
  vector<int> v(arr, arr+3);
  /* compute a posterior distribution */
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
