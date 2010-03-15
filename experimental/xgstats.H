#ifndef XGSTATS_H
#define XGSTATS_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

using namespace std;

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

// a template function for stable numerics
template <typename T>
double logsumexp(T begin, T end);

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

void test_uniform_mixture();

void test_mixture_b();

void test_mixture();

void test_finite_distn();

#endif
