#include <iostream>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>

#include "xgstats.H"
#include "phylo.H"
#include "dgrpstats.H"

using namespace std;

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

ZygosityState::ZygosityState(const vector<vector<double> > &distns,
    const Model<int> *pcoverage_model)
{
  vector<vector<double> >::const_iterator it;
  for (it = distns.begin(); it != distns.end(); ++it)
  {
    models_.push_back(new SinglePatternState(*it, pcoverage_model));
  }
  update_info();
}

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
  // update cached info
  update_info();
}
