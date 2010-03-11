#include <cmath>

using namespace std;


// Observation interface.
class Obs
{
  virtual ~Obs() = 0;
};

// An observation which is a vector of integers.
class ObsVector: public Obs
{
  public:
    ObsVector();
}

class Mixture: public Model
{
  public:
    get_likelihood(const Obs *pobs);
    get_log_likelihood(const Obs *pobs);
    get_posterior_distribution(const Obs *pobs);
};

class UniformMixture: public Mixture
{
  public:
};

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




/*
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