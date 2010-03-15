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
    double get_lik(const T& obs) const;
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
double Mixture<T>::get_lik(const T& obs) const
{
  return exp(get_log_lik(obs));
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

int main(int argc, const char *argv[])
{
  cout << log(0.0) << endl;
  test_mixture();
  test_finite_distn();
  test_mixture_b();
  test_uniform_mixture();
  return 0;
}



/*

class Flat: public Model
{
  protected:
    double m_mu; // geometric distribution parameter
    double m_pr; // geometric distribution parameter
    double m_log_pr; // precalculated logarithm
    double m_log_not_pr; // precalculated logarithm
  public:
    int m_nstates;
    int expected_coverage;
    Flat(int nstates, int expected_coverage);
    get_likelihood(const Obs *pobs);
    get_log_likelihood(const Obs *pobs);
}

// Each state has a geometrically distributed count.
Flat::Flat(int nstates, int expected_coverage) :
  m_nstates(nstates), m_expected_coverage(expected_coverage)
{
  // precalculate part of the log likelihood
  m_mu = m_expected_coverage / (double) m_nstates;
  m_pr = 1.0 / (m_mu + 1.0);
  if (m_pr != 0.0)
  {
    m_log_pr = log(m_pr);
  }
  if (m_pr != 1.0)
  {
    m_log_not_pr = log(1.0 - m_pr);
  }
}

double Flat::get_likelihood(const Obs *pobs)
{
  return exp(get_log_likelihood(pobs));
}

Flat::get_log_likelihood(const Obs *pobs)
{
  if len(observation) != self.nstates:
      raise ValueError('expected a vector of %d integers' % self.nstates)
  if self.pr == 0.0:
      return float('-inf')
  if self.pr == 1.0:
      if any(observation):
          return float('-inf')
      else:
          return 0
  return sum(observation) * self.log_not_pr + self.nstates * self.log_pr
}




"""
This module defines some distributions related to resequencing.

The defined distributions are for hidden Markov models.
For each distribution observations can be sampled,
and observations can be assigned a likelihood.
Observations are ordered sequences of four integers
corresponding to reads of A, C, G, and T.
"""

def get_homozygous_distributions(randomization_rate):
    """
    Each distribution is over four states.
    @param randomization_rate: the probability that a read is randomized
    @return: a list of four distributions
    """
    distributions = []
    for favorite in range(4):
        d = [randomization_rate/4.0]*4
        d[favorite] = 1 - 3*randomization_rate/4.0
        distributions.append(d)
    return distributions

def get_heterozygous_distributions(randomization_rate):
    """
    Each distribution is over four states.
    @param randomization_rate: the probability that a read is randomized
    @return: a list of six distributions
    """
    distributions = []
    for first_index in range(4):
        for second_index in range(first_index):
            d = [randomization_rate/4.0]*4
            d[first_index] = .5 - randomization_rate/4.0
            d[second_index] = .5 - randomization_rate/4.0
            distributions.append(d)
    return distributions


class Mixture:
    """
    This allows sampling and likelihood calculations for a mixture model.
    This class can act as a HMM hidden state.
    """

    def __init__(self, states, distribution):
        """
        @param states: a sequence of hidden states
        @param distribution: the distribution of the hidden states
        """
        # do some validation
        if not len(states):
            raise ValueError('no states were specified')
        if len(states) != len(distribution):
            msg = 'the number of states should match the distribution length'
            raise ValueError(msg)
        if not np.allclose(sum(distribution), 1):
            raise ValueError('expected the distribution to sum to 1.0')
        if min(distribution) < 0:
            msg = 'expected the distribution to be a stochastic vector'
            raise ValueError(msg)
        # store the arguments, leaving out states with zero probability
        self.states = [state for state, d in zip(states, distribution) if d]
        self.distribution = [d for d in distribution if d]
        # precompute part of the likelihood
        self.log_distribution = [math.log(p) if p else float('-inf')
                for p in self.distribution]

    def get_posterior_distribution(self, observation):
        log_likelihoods = [state.get_log_likelihood(observation)
                for state in self.states]
        weighted_lls = [ll + log_p
                for ll, log_p in zip(log_likelihoods, self.log_distribution)]
        obs_ll = scipy.maxentropy.logsumexp(weighted_lls)
        return [math.exp(ll - obs_ll) for ll in weighted_lls]

    def sample_observation(self):
        """
        @return: an observation sampled from the distribution
        """
        state = self.states[xgstats.random_weighted_int(self.distribution)]
        return state.sample_observation()

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))

    def get_log_likelihood(self, observation):
        log_likelihoods = [state.get_log_likelihood(observation)
                for state in self.states]
        if all(ll==float('-inf') for ll in log_likelihoods):
            return float('-inf')
        weighted_log_likelihoods = [ll + log_p
                for ll, log_p in zip(log_likelihoods, self.log_distribution)]
        return scipy.maxentropy.logsumexp(weighted_log_likelihoods)


class UniformMixture:
    """
    This allows sampling and likelihood calculations for a mixture model.
    Each component of the mixture is equally likely.
    This class can act as a HMM hidden state.
    """

    def __init__(self, states):
        """
        @param states: a sequence of hidden states
        """
        self.states = states

    def sample_observation(self):
        return random.choice(self.states).sample_observation()

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))

    def get_log_likelihood(self, observation):
        log_likelihoods = [state.get_log_likelihood(observation) for state in self.states]
        if all(ll==float('-inf') for ll in log_likelihoods):
            return float('-inf')
        log_likelihood = scipy.maxentropy.logsumexp(log_likelihoods) - math.log(len(self.states))
        return log_likelihood


class SinglePatternState:
    """
    This is for when a single pattern is expected.
    Useful states are mixtures of these states.
    """

    def __init__(self, distribution, expected_coverage):
        """
        @param distribution: the expected nucleotide distribution
        @param expected_coverage: the read coverage at a position is poisson distributed with this expectation
        """
        self.distribution = distribution
        self.expected_coverage = expected_coverage

    def sample_observation(self):
        """
        @return: a sample of counts in each A, C, G, T state
        """
        n = np.random.poisson(self.expected_coverage)
        return np.random.multinomial(n, self.distribution)

    def get_log_likelihood(self, observation):
        if len(observation) != 4:
            raise ValueError('expected the observation to be a vector of four integers')
        n = sum(observation)
        accum = 0
        accum += xgstats.poisson_log_pmf(n, self.expected_coverage)
        accum += xgstats.multinomial_log_pmf(self.distribution, observation)
        return accum

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class FlatState:
    """
    This is supposed to be a somewhat flat distribution.
    Each of the counts is sampled independently
    according to a geometric distribution.
    The distribution of the sum of counts is negative binomial.
    """

    def __init__(self, nstates, expected_coverage):
        """
        Each state has a geometrically distributed count.
        @param nstates: the number of different possible states
        @param expected_coverage: expected read coverage at a position
        """
        # store the arguments
        self.nstates = nstates
        self.expected_coverage = expected_coverage
        # precalculate part of the log likelihood
        self.mu = self.expected_coverage / float(self.nstates)
        self.pr = 1/(self.mu+1)
        if self.pr != 0.0:
            self.log_pr = math.log(self.pr)
        if self.pr != 1.0:
            self.log_not_pr = math.log(1.0 - self.pr)

    def sample_observation(self):
        """
        @return: a sample of counts in each state
        """
        return tuple(scipy.stats.geom.rvs(self.pr, loc=-1, size=self.nstates))

    def get_log_likelihood(self, observation):
        if len(observation) != self.nstates:
            raise ValueError('expected a vector of %d integers' % self.nstates)
        if self.pr == 0.0:
            return float('-inf')
        if self.pr == 1.0:
            if any(observation):
                return float('-inf')
            else:
                return 0
        return sum(observation) * self.log_not_pr + self.nstates * self.log_pr

    def get_likelihood(self, observation):
        return math.exp(self.get_log_likelihood(observation))


class Homozygous(UniformMixture):
    def __init__(self, randomization_rate, expected_coverage):
        distributions = get_homozygous_distributions(randomization_rate)
        states = [SinglePatternState(d, expected_coverage) for d in distributions]
        UniformMixture.__init__(self, states)

class Heterozygous(UniformMixture):
    def __init__(self, randomization_rate, expected_coverage):
        distributions = get_heterozygous_distributions(randomization_rate)
        states = [SinglePatternState(d, expected_coverage) for d in distributions]
        UniformMixture.__init__(self, states)

class Overcovered(UniformMixture):
    def __init__(self, randomization_rate, expected_coverage):
        distributions = get_homozygous_distributions(randomization_rate) + get_heterozygous_distributions(randomization_rate)
        states = [FlatState(4, expected_coverage)] + [SinglePatternState(d, expected_coverage) for d in distributions]
        UniformMixture.__init__(self, states)

*/
