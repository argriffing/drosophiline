#ifndef XGSTATS_H
#define XGSTATS_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

typedef vector<double> vdouble_t;

/*
 * In this file, 'lik' is an abbreviation of 'likelihood'.
 * Also, 'distn' is an abbreviation of 'distribution'.
 */

// a deletion functor
struct delete_object
{
  template <typename T>
  void operator()(T *ptr){ delete ptr;}
};

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
  for (it=begin; it!=end; ++it)
  {
    if (it != max_it)
    {
      accum += exp(*it - max_value);
    }
  }
  return log1p(accum) + max_value;
}

/*
 * The following few functions are directly ported from python.
 * They should be parameterized the same way as their scipy counterparts,
 * but they have additional checking for special cases,
 * and they can somewhat deal with hugely negative log likelihoods.
 * Nevertheless, they should probably all be rewritten.
 */

double binomial_log_pmf(int observed_n, int max_n, double p_success);

double geometric_log_pmf(int observed_n, double pr);

double poisson_log_pmf(int observed_n, int expected_n);

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

// statistical distribution templated base class
template <typename T>
class Model
{
  public:
    virtual ~Model() {}
    virtual double get_lik(const T& obs) const = 0;
    virtual double get_log_lik(const T& obs) const = 0;
};

template <typename T>
class Mixture: public Model<T>
{
  public:
    Mixture() {}
    Mixture(const vector<double> &distn, const vector<Model<T> *> &models) :
        distn_(distn), models_(models) {update_info();}
    double get_lik(const T& obs) const {return exp(get_log_lik(obs));}
    double get_log_lik(const T& obs) const;
    vector<double> get_posterior_distn(const T& obs);
  protected:
    void update_info()
    {
      nmodels_ = models_.size();
      log_distn_.clear();
      transform(distn_.begin(), distn_.end(),
          back_inserter(log_distn_), (double (*)(double)) log);
    }
    vector<double> distn_;
    vector<Model<T> *> models_;
    int nmodels_;
    vector<double> log_distn_;
};

/*
 * The name suggests that this is a subclass of Mixture.
 * But this is not the case.
 */
template <typename T>
class UniformMixture: public Model<T>
{
  public:
    UniformMixture() {}
    UniformMixture(const vector<Model<T> *> &models)
    {
      models_ = models;
      update_info();
    }
    double get_lik(const T&obs) const {return exp(get_log_lik(obs));}
    double get_log_lik(const T& obs) const
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
      vector<double>::iterator it;
      it = max_element(log_liks.begin(), log_liks.end());
      if (*it == -numeric_limits<double>::infinity())
        return -numeric_limits<double>::infinity();
      // return the log likelihood
      return logsumexp(log_liks.begin(), log_liks.end()) - log_nmodels_;
    }
    vector<double> get_posterior_distn(const T& obs)
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
        return vdouble_t(nmodels_, numeric_limits<double>::signaling_NaN());
      // compute the posterior distribution
      double log_lik_sum = logsumexp(post.begin(), post.end());
      for (i=0; i<nmodels_; i++)
      {
        post[i] = exp(post[i] - log_lik_sum);
      }
      return post;
    }
  protected:
    void update_info()
    {
      nmodels_ = models_.size();
      log_nmodels_ = log(nmodels_);
    }
    vector<Model<T> *> models_;
    int nmodels_;
    double log_nmodels_;
};

class FiniteDistn: public Model<int>
{
  public:
    // override base class member functions
    FiniteDistn() {}
    FiniteDistn(const vector<double> &distn) {set_distn(distn);}
    double get_lik(const int &obs) const {return distn_[obs];}
    double get_log_lik(const int& obs) const {return log_distn_[obs];}
    // add some new functions
    void set_distn(const vector<double>&);
  private:
    vector<double> distn_;
    vector<double> log_distn_;
};

class FairD6: public FiniteDistn
{
  public:
    FairD6() : FiniteDistn(vector<double>(6, 1.0/6.0)) {}
};

class LoadedD6: public FiniteDistn
{
  public:
    LoadedD6()
    {
      vector<double> v(6, 0.1);
      v[5] = .5;
      set_distn(v);
    }
};

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
  if (models_.empty())
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

void test_uniform_mixture();

void test_mixture_b();

void test_mixture();

void test_finite_distn();

#endif
