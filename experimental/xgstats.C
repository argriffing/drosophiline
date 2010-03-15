#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>

using namespace std;

/*
 * In this file, 'lik' is an abbreviation of 'likelihood'.
 * Also, 'distn' is an abbreviation of 'distribution'.
 */

/*
 * This is a functor.
 */
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
  for (it=begin; it!=end; it++)
  {
    if (it != max_it)
    {
      accum += exp(*it - max_value);
    }
  }
  return log1p(accum) + max_value;
}

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
    // override base class member functions
    Mixture() : Model<T>(), nmodels_(0) {}
    Mixture(const vector<double>&, const vector<Model<T> *>&);
    double get_lik(const T& obs) const {return exp(get_log_lik(obs));}
    double get_log_lik(const T& obs) const;
    // add some new functions
    void set_distn_and_models(const vector<double>&,
        const vector<Model<T> *>&);
    vector<double> get_posterior_distn(const T& obs);
  private:
    int nmodels_;
    vector<Model<T> *> models_;
    vector<double> distn_;
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
    UniformMixture() :
      Model<T>(), nmodels_(0), log_nmodels_(log(nmodels_)) {}
    UniformMixture(const vector<Model<T> *> &models) :
      Model<T>(), nmodels_(models.size()), models_(models),
      log_nmodels_(log(nmodels_)) {}
    double get_lik(const T&obs) const {return exp(get_log_lik(obs));}
    double get_log_lik(const T&obs) const;
    vector<double> get_posterior_distn(const T& obs);
  private:
    int nmodels_;
    vector<Model<T> *> models_;
    double log_nmodels_;
};

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

class FiniteDistn: public Model<int>
{
  public:
    // override base class member functions
    FiniteDistn() : Model<int>() {}
    FiniteDistn(const vector<double> &distn) : Model<int>() {set_distn(distn);}
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
    FairD6() : FiniteDistn() {set_distn(vector<double>(6, 1.0/6.0));}
};

class LoadedD6: public FiniteDistn
{
  public:
    LoadedD6() : FiniteDistn()
    {
      vector<double> v(6, 0.1);
      v[5] = .5;
      set_distn(v);
    }
};

void FiniteDistn::set_distn(const vector<double> &distn)
{
  distn_ = distn;
  log_distn_.clear();
  transform(distn_.begin(), distn_.end(),
      back_inserter(log_distn_), (double (*)(double)) log);
}

template<typename T>
Mixture<T>::Mixture(const vector<double> &distn,
    const vector<Model<T> *> &models) : Model<T>()
{
  set_distn_and_models(distn, models);
}

template<typename T>
void Mixture<T>::set_distn_and_models(const vector<double> &distn,
    const vector<Model<T> *> &models)
{
  distn_ = distn;
  models_ = models;
  /* precompute some stuff */
  nmodels_ = models.size();
  log_distn_.clear();
  transform(distn_.begin(), distn_.end(),
      back_inserter(log_distn_), (double (*)(double)) log);
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
