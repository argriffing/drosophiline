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

vector<vector<double> > get_RR_distns(double r)
{
  vector<vector<double> > distns;
  vector<double> d(4, r/4);
  d[0] = 1 - 3*r/4;
  distns.push_back(d);
  return distns;
}

vector<vector<double> > get_RA_distns(double r)
{
  vector<vector<double> > distns;
  int i;
  for (i = 1; i < 4; ++i)
  {
    vector<double> d(4, r/4);
    d[0] = .5 - r/4;
    d[i] = .5 - r/4;
    distns.push_back(d);
  }
  return distns;
}

vector<vector<double> > get_AA_distns(double r)
{
  vector<vector<double> > distns;
  int i;
  for (i = 1; i < 4; ++i)
  {
    vector<double> d(4, r/4);
    d[i] = 1 - 3*r/4;
    distns.push_back(d);
  }
  return distns;
}

vector<vector<double> > get_AB_distns(double r)
{
  vector<vector<double> > distns;
  int i, j;
  for (i = 1; i < 4; ++i)
  {
    for (j = i+1; j < 4; ++j)
    {
      vector<double> d(4, r/4);
      d[i] = .5 - r/4;
      d[j] = .5 - r/4;
      distns.push_back(d);
    }
  }
  return distns;
}

/*
 * @param d: the expected number of substitutions along a branch
 * @return: the probability that the endpoints are different
 */
double jc69_distance_to_probability(double d)
{
  return (3.0 / 4.0) * (1.0 - exp(-(4.0 / 3.0) * d));
}

/*
 * This is based on the Jukes-Cantor model on a three taxon tree.
 * @param ref_length: length of the reference taxon branch
 * @param child_length: length of each child taxon branch
 * @return: the distribution (RR, RA, AA, AB)
 */
vector<double> get_zygosity_distn(double ref_length, double child_length)
{
  double p_ref_change = jc69_distance_to_probability(ref_length);
  double p_child_change = jc69_distance_to_probability(child_length);
  // For now sum over all possibilities of non-reference nodes.
  // This could be done more efficiently using Felsenstein pruning,
  // but I am ignoring this for now.
  double p_RR = 0.0;
  double p_RA = 0.0;
  double p_AA = 0.0;
  double p_AB = 0.0;
  const int ref = 0;
  int c1, c2, c12;
  double p1, p2, p12;
  for (c12 = 0; c12 < 4; ++c12)
  {
    if (c12 == ref)
      p12 = 1.0 - p_ref_change;
    else
      p12 = p_ref_change / 3.0;
    for (c1 = 0; c1 < 4; ++c1)
    {
      if (c1 == c12)
        p1 = p12 * (1.0 - p_child_change);
      else
        p1 = p12 * (p_child_change / 3.0);
      for (c2 = 0; c2 < 4; ++c2)
      {
        if (c2 == c12)
          p2 = p1 * (1.0 - p_child_change);
        else
          p2 = p1 * (p_child_change / 3.0);
        // Classify the joint distribution
        // and add weight to the appropriate state.
        if (c1 == ref && c2 == ref)
            p_RR += p2;
        else if (c1 == ref || c2 == ref)
            p_RA += p2;
        else if (c1 == c2)
            p_AA += p2;
        else
            p_AB += p2;
      }
    }
  }
  vector<double> distn;
  distn.push_back(p_RR);
  distn.push_back(p_RA);
  distn.push_back(p_AA);
  distn.push_back(p_AB);
  return distn;
}

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

double SinglePatternState::get_log_lik(const vector<int> &obs) const
{
  if (distn_.empty())
    return numeric_limits<double>::signaling_NaN();
  int n = accumulate(obs.begin(), obs.end(), 0);
  double accum = 0;
  accum += pcoverage_model_->get_log_lik(n);
  accum += multinomial_log_pmf(distn_.begin(), distn_.end(),
      obs.begin(), obs.end());
  return accum;
}

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

ZygosityState::ZygosityState(const vector<vector<double> > &distns,
    const Model<int> *pcoverage_model)
{
  vector<vector<double> >::const_iterator it;
  for (it = distns.begin(); it != distns.end(); ++it)
  {
    models_.push_back(new SinglePatternState(*it, pcoverage_model));
  }
}

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
    ~FlatState();
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

FlatState::FlatState(int nstates, int expected_coverage)
{
  nstates_ = nstates;
  mu_ = expected_coverage / (double) nstates_;
  pr_ = 1 / (mu_ + 1);
  log_pr_ = log(pr_);
  log_not_pr_ = log(1.0 - pr_);
}

double FlatState::get_log_lik(const vector<int> &obs) const
{
  if (obs.empty())
    return numeric_limits<double>::signaling_NaN();
  if (pr_ == 0.0)
    return -numeric_limits<double>::infinity();
  int obs_sum = accumulate(obs.begin(), obs.end(), 0);
  if (pr_ == 1.0)
  {
    if (obs_sum == 0)
      return 0;
    else
      return -numeric_limits<double>::infinity();
  }
  // handle non-special cases
  double accum = 0;
  accum += obs_sum * log_not_pr_;
  accum += nstates_ * log_pr_;
  return accum;
}


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
 * @param dref: a branch length parameter
 * @param dchild: a branch length parameter
 * @param r: probability of sequence randomization
 * @param nomcoverage: nominal coverage
 * @param kmulticoverages: allowed multiples of nominal coverage
 */
GoodState::GoodState(double dref, double dchild, double r,
    int nomcoverage, int kmulticoverages)
{
  // create the coverage model
  pcoverage_model_ = new GoodMultiCoverage(nomcoverage, kmulticoverages);
  // create the mixture components
  models_.push_back(new ZygosityState(get_RR_distns(r), pcoverage_model_));
  models_.push_back(new ZygosityState(get_RA_distns(r), pcoverage_model_));
  models_.push_back(new ZygosityState(get_AA_distns(r), pcoverage_model_));
  models_.push_back(new ZygosityState(get_AB_distns(r), pcoverage_model_));
  // define the mixture parameters
  distn_ = get_zygosity_distn(dref, dchild);
}
/*
 * A predominantly homozygous region.
 */
class HMMRecent: public GoodState
{
  public:
    HMMRecent(double x, double y, double z,
        double seqerr, int nomcoverage, int kmulticoverages) :
        GoodState(x+y, z, seqerr, nomcoverage, kmulticoverages) {}
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

void demo_zygosity_models()
{
  // define an observation of RABC counts
  vector<int> obs;
  obs.push_back(20);
  obs.push_back(2);
  obs.push_back(1);
  obs.push_back(0);
  // define some model parameters
  int low = 2;
  int med = 20;
  int high = 1000;
  double x = 0.1;
  double y = 0.01;
  double z = .0001;
  double seqerr = .1;
  int nomcoverage = 20;
  int kmulticoverages = 4;
  // create the models
  HMMGarbage garbage(low, med, high);
  HMMRecent recent(x, y, z, seqerr, nomcoverage, kmulticoverages);
  HMMAncient ancient(x, y, z, seqerr, nomcoverage, kmulticoverages);
  Model<vector<int> > *pgarbage = &garbage;
  Model<vector<int> > *precent = &recent;
  Model<vector<int> > *pancient = &ancient;
  // show the observation
  vector<int>::iterator it;
  cout << "observation: ";
  for (it = obs.begin(); it != obs.end(); ++it)
  {
    cout << *it << " ";
  }
  cout << endl;
  // compute a likelihood for each model
  cout << "recent likelihood: " << precent->get_lik(obs) << endl;
  cout << "ancient likelihood: " << pancient->get_lik(obs) << endl;
  cout << "other likelihood: " << pgarbage->get_lik(obs) << endl;
}

int main(int argc, const char *argv[])
{
  //cout << log(0.0) << endl;
  //test_mixture();
  //test_finite_distn();
  //test_mixture_b();
  //test_uniform_mixture();
  demo_poisson_log_pmf();
  demo_geometric_log_pmf();
  demo_zygosity_models();
  return 0;
}

