#ifndef DGRPSTATS_H
#define DGRPSTATS_H

#include <iostream>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>

#include "xgstats.H"

using namespace std;

/*
 * Get lists of RABC distributions.
 * R is the reference count,
 * and ABC are the three non-reference counts.
 */

vector<vector<double> > get_RR_distns(double r);
vector<vector<double> > get_RA_distns(double r);
vector<vector<double> > get_AA_distns(double r);
vector<vector<double> > get_AB_distns(double r);

class WrappedPoisson: public Model<int>
{
  public:
    WrappedPoisson(int expectation) : expectation_(expectation) {}
    double get_lik(const int &obs) const {return exp(get_log_lik(obs));}
    double get_log_lik(const int &obs) const
    {
      return poisson_log_pmf(obs, expectation_);
    }
  private:
    double expectation_;
};

class GoodMultiCoverage: public UniformMixture<int>
{
  public:
    GoodMultiCoverage(int nomcoverage, int kmulticoverages)
    {
      // create the models
      int i;
      for (i=0; i<kmulticoverages; i++)
      {
        int e = nomcoverage * (i + 1);
        models_.push_back(new WrappedPoisson(e));
      }
      update_info();
    }
    ~GoodMultiCoverage()
    {
      // destroy the models
      for_each(models_.begin(), models_.end(), delete_object());
    }
};

/*
 * Defines a nucleotide distribution and a coverage distribution.
 * Together, these give a distribution over count vectors.
 */
class SinglePatternState: public Model<vector<int> >
{
  public:
    SinglePatternState(const vector<double> &distn,
      const Model<int> *pcoverage_model) : 
      distn_(distn), pcoverage_model_(pcoverage_model) {}
    double get_lik(const vector<int> &obs) const
    {
      return exp(get_log_lik(obs));
    }
    double get_log_lik(const vector<int> &obs) const;
  private:
    vector<double> distn_;
    const Model<int> *pcoverage_model_;
};

class ZygosityState: public UniformMixture<vector<int> >
{
  public:
    ZygosityState(const vector<vector<double> >&, const Model<int>*);
    ~ZygosityState()
    {
      for_each(models_.begin(), models_.end(), delete_object());
    }
  private:
    const Model<int> *pcoverage_model_;
};

/*
 * This is supposed to be a somewhat flat distribution.
 * Each of the counts is sampled independently
 * according to a geometric distribution.
 * The distribution of the sum of counts is negative binomial.
 */
class FlatState: public Model<vector<int> >
{
  public:
    FlatState(int nstates, int expected_coverage);
    double get_log_lik(const vector<int>&) const;
    double get_lik(const vector<int> &obs) const
    {
      return exp(get_log_lik(obs));
    }
  private:
    int nstates_;
    double mu_;
    double pr_;
    double log_pr_;
    double log_not_pr_;
};

class GoodState: public Mixture<vector<int> >
{
  public:
    GoodState(double dref, double dchild, double seqerr,
        int nomcoverage, int kmulticoverages);
    ~GoodState()
    {
      // clean up
      for_each(models_.begin(), models_.end(), delete_object());
      delete pcoverage_model_;
    }
  private:
    GoodMultiCoverage *pcoverage_model_;
};

/*
 * This region has mostly homozygous states.
 */
class HMMRecent : public GoodState
{
  public:
    HMMRecent (double x, double y, double z,
        double seqerr, int nomcoverage, int kmulticoverages) :
        GoodState(x+y, y, seqerr, nomcoverage, kmulticoverages) {}
};

/*
 * This region has many heterozygous states.
 */
class HMMAncient: public GoodState
{
  public:
    HMMAncient(double x, double y, double z,
        double seqerr, int nomcoverage, int kmulticoverages) :
        GoodState(x, y+z, seqerr, nomcoverage, kmulticoverages) {}
};

/*
 * This region has states with ill-defined zygosity.
 */
class HMMGarbage: public UniformMixture<vector<int> >
{
  public:
    HMMGarbage(int low, int med, int high)
    {
      models_.push_back(new FlatState(4, low));
      models_.push_back(new FlatState(4, med));
      models_.push_back(new FlatState(4, high));
      update_info();
    }
    ~HMMGarbage()
    {
      for_each(models_.begin(), models_.end(), delete_object());
    }
};

#endif
