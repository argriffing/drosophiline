#include <iostream>
#include <vector>
#include <limits>

#include "xgstats.H"

using namespace std;

class WrappedPoisson: public Model<int>
{
  public:
    WrappedPoisson() : Model<int>(), expectation_(0) {}
    WrappedPoisson(int expectation) :
      Model<int>(), expectation_(expectation) {}
    double get_lik(const int &obs) {return exp(poisson_log_pmf(obs));}
    double get_log_lik(const int &obs) {return poisson_log_pmf(obs);}
  private:
    double expectation_;
};

class GoodMultiCoverage: public UniformMixture<int>
{
  public:
    GoodMultiCoverage() : UniformMixture() {}
    GoodMultiCoverage(int nom_coverage, int kmulticoverages) : UniformMixture()
    {
      // create the models
      int i;
      for (i=0; i<kmulticoverages; i++)
      {
        int e = nom_coverage * (i + 1);
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
    SinglePatternState() : Mixture<vector<int> >(),
      pcoverage_model_(NULL) {}
    SinglePatternState(const vector<double> &distn,
      const Model<int> *pcoverage_model) : Model<vector<int> >(),
      distn_(distn), pcoverage_model_(pcoverage_model) {}
    double get_lik(const vector<int> &obs) {return exp(get_log_lik(obs));}
    double get_log_lik(const vector<int> &obs);
  private:
    vector<double> distn_;
    const Model<int> *pcoverage_model_;
};

double SinglePatternState::get_log_lik(const vector<int> &obs)
{
  if (distn_.empty())
    return numeric_limits<double>::signaling_NaN();
  int n = accumulate(obs.begin(), obs.end(), 0);
  double accum = 0;
  accum += pcoverage_model->get_log_lik(n);
  accum += multinomial_log_pmf(distn_.begin(), distn_.end(),
      obs.begin(), obs.end());
  return accum;
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


def gen_RR_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    d = [r/4.0]*4
    d[0] = 1 - 3*r/4.0
    yield d

def gen_RA_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    for i in range(1, 4):
        d = [r/4.0]*4
        d[0] = .5 - r/4.0
        d[i] = .5 - r/4.0
        yield d

def gen_AA_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    for i in range(1, 4):
        d = [r/4.0]*4
        d[i] = 1 - 3*r/4.0
        yield d

def gen_AB_distns(r):
    """
    Yield distributions over the elements of an observation.
    @param r: sequencing randomization rate
    """
    for i in range(1, 4):
        for j in range(i+1, 4):
            d = [r/4.0]*4
            d[i] = .5 - r/4.0
            d[j] = .5 - r/4.0
            yield d

class GoodState(ReadCoverage.Mixture):
    def __init__(self, dref, dchild, seqerr, nomcoverage, kmulticoverages):
        """
        @param dref: a branch length parameter
        @param dchild: a branch length parameter
        @param seqerr: probability of sequence randomization
        @param nomcoverage: nominal coverage
        @param kmulticoverages: allowed multiples of nominal coverage
        """
        mcov = GoodMultiCoverage(nomcoverage, kmulticoverages)
        # define the states
        r = seqerr
        RR_states = [SinglePatternState(d, mcov) for d in gen_RR_distns(r)]
        RA_states = [SinglePatternState(d, mcov) for d in gen_RA_distns(r)]
        AA_states = [SinglePatternState(d, mcov) for d in gen_AA_distns(r)]
        AB_states = [SinglePatternState(d, mcov) for d in gen_AB_distns(r)]
        # define the distributions
        RR = ReadCoverage.UniformMixture(RR_states)
        RA = ReadCoverage.UniformMixture(RA_states)
        AA = ReadCoverage.UniformMixture(AA_states)
        AB = ReadCoverage.UniformMixture(AB_states)
        states = (RR, RA, AA, AB)
        zygo_distn = dgrp.get_zygosity_distribution(dref, dchild)
        ReadCoverage.Mixture.__init__(self, states, zygo_distn)

class HMMRecent(GoodState):
    """
    A predominantly homozygous region.
    """
    def __init__(self, x, y, z, seqerr, nomcoverage, kmulticoverages):
        GoodState.__init__(self, x+y, z, seqerr, nomcoverage, kmulticoverages)

class HMMAncient(GoodState):
    """
    This region has many heterozygous states.
    """
    def __init__(self, x, y, z, seqerr, nomcoverage, kmulticoverages):
        GoodState.__init__(self, x, y+z, seqerr, nomcoverage, kmulticoverages)

class HMMGarbage(ReadCoverage.UniformMixture):
    """
    This region has states with ill-defined zygosity.
    """
    def __init__(self, low, med, high):
        states = [ReadCoverage.FlatState(4, x) for x in (low, med, high)]
        ReadCoverage.UniformMixture.__init__(self, states)


class TestReadCoverageRef:

    def test_gen_xx_distns(self):
        for r in (.001, .01, .1, .9):
            RR = list(gen_RR_distns(r))
            RA = list(gen_RA_distns(r))
            AA = list(gen_AA_distns(r))
            AB = list(gen_AB_distns(r))
            self.assertEqual(len(RR), 1)
            self.assertEqual(len(RA), 3)
            self.assertEqual(len(AA), 3)
            self.assertEqual(len(AB), 6)
            all_distns = RR + RA + AA + AB
            for d in all_distns:
                self.assertAlmostEqual(sum(d), 1.0)


if __name__ == '__main__':
    unittest.main()


/*
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
*/


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


def _gen_observations(n):
    """
    This is a helper function for testing.
    @param n: the total number of counts
    """
    for i in range(n+1):
        for j in range(n+1-i):
            for k in range(n+1-i-j):
                yield i, j, k, n-i-j-k


class TestReadCoverage(unittest.TestCase):

    def test_gen_observations(self):
        """
        Validate the sequences of possible observations.
        """
        # Assert that the right number of sequences is generated.
        # See OEIS A000292.
        for i in range(5):
            expected = ((i+1)*(i+2)*(i+3))/6
            observed = len(list(_gen_observations(i)))
            self.assertEqual(expected, observed)
        # Get a sequence.
        observations = list(_gen_observations(5))
        # Assert that each observation in the sequence is unique.
        self.assertEqual(len(observations), len(set(observations)))
        # Assert that each observation sums to the right number.
        for obs in observations:
            self.assertEqual(sum(obs), 5)

    def test_homozygous_sample(self):
        state = Homozygous(.2, 10)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)

    def test_homozygous_likelihood(self):
        state = Homozygous(.2, 5)
        for n in range(10):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            p_expected = scipy.stats.poisson.pmf(n, 5)
            self.assertTrue(np.allclose(p_observed, p_expected), (n, p_observed, p_expected))

    def test_heterozygous_sample(self):
        state = Heterozygous(.2, 10)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)

    def test_heterozygous_likelihood(self):
        state = Heterozygous(.2, 5)
        for n in range(10):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            p_expected = scipy.stats.poisson.pmf(n, 5)
            self.assertTrue(np.allclose(p_observed, p_expected), (n, p_observed, p_expected))

    def test_flat_sample(self):
        state = FlatState(4, 10)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)
        for x in observation:
            self.assertTrue(0 <= x)

    def test_flat_likelihood(self):
        expected_coverage = 5
        state = FlatState(4, expected_coverage)
        for n in range(10):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            pr = 1.0 / (1.0 + expected_coverage / 4.0)
            p_expected = scipy.stats.nbinom.pmf(n, 4, pr)
            self.assertTrue(np.allclose(p_observed, p_expected), (n, p_observed, p_expected))

    def test_overcovered_sample(self):
        state = Overcovered(.2, 5)
        observation = state.sample_observation()
        self.assertEqual(len(observation), 4)
        for x in observation:
            self.assertTrue(0 <= x)

    def test_overcovered_likelihood(self):
        state = Overcovered(.2, 4)
        for n in range(8):
            p_observed = sum(state.get_likelihood(observation) for observation in _gen_observations(n))
            self.assertTrue(0 <= p_observed <= 1)

    def test_likelihood_ratios(self):
        """
        Assert that the distributions do what we want them to do.
        That is, assert that the objective classification of observations
        matches our subjective classification.
        """
        # define parameters
        good_coverage = 20
        bad_coverage = 60
        error_rate = .1
        # define hidden states
        homozygous = Homozygous(error_rate, good_coverage)
        heterozygous = Heterozygous(error_rate, good_coverage)
        overcovered = Overcovered(error_rate, bad_coverage)
        # define observations
        obs_a = [19, 1, 0, 0]
        obs_b = [10, 7, 0, 0]
        obs_c = [150, 12, 2, 1]
        # assert that the inference we want has the best log likelihood
        self.assertTrue(homozygous.get_log_likelihood(obs_a) > heterozygous.get_log_likelihood(obs_a))
        self.assertTrue(homozygous.get_log_likelihood(obs_a) > overcovered.get_log_likelihood(obs_a))
        self.assertTrue(heterozygous.get_log_likelihood(obs_b) > homozygous.get_log_likelihood(obs_b))
        self.assertTrue(heterozygous.get_log_likelihood(obs_b) > overcovered.get_log_likelihood(obs_b))
        self.assertTrue(overcovered.get_log_likelihood(obs_c) > homozygous.get_log_likelihood(obs_c))
        self.assertTrue(overcovered.get_log_likelihood(obs_c) > heterozygous.get_log_likelihood(obs_c))

    def test_large_observation_sizes(self):
        """
        Assert that the stats module can handle large input.
        """
        # define parameters
        good_coverage = 20
        bad_coverage = 60
        error_rate = .1
        # define hidden states
        homozygous = Homozygous(error_rate, good_coverage)
        heterozygous = Heterozygous(error_rate, good_coverage)
        overcovered = Overcovered(error_rate, bad_coverage)
        models = (homozygous, heterozygous, overcovered)
        # define the large observation
        obs = (0, 0, 605, 2)
        # attempt to get a log likelihood for each model
        for model in models:
            log_likelihood = model.get_log_likelihood(obs)
            self.failIf(math.isnan(log_likelihood))

    def test_myopic_model(self):
        """
        Test a mixture with a degenerately single-minded component.
        """
        # define a hidden state that is a mixture
        components = (FlatState(4, 10), SinglePatternState((.9, .1, 0, 0), 10))
        model = UniformMixture(components)
        # define an observation that is completely incompatible with the second component
        obs = (0, 0, 10, 2)
        # assert that the log likelihood is reasonable
        log_likelihood = model.get_log_likelihood(obs)
        self.failIf(math.isnan(log_likelihood))

    def test_failed_mixture(self):
        """
        Test a mixture where no component can explain the data at all.
        """
        components = (
            SinglePatternState((0.9, 0.1, 0.0, 0.0), 10),
            SinglePatternState((0.8, 0.1, 0.1, 0.0), 10))
        model = UniformMixture(components)
        obs = (0, 0, 10, 2)
        log_likelihood = model.get_log_likelihood(obs)
        self.failIf(math.isnan(log_likelihood))

    def test_weighted_mixture_model_compatibility(self):
        """
        Create a weighted mixture that is the same as a uniform mixture.
        The likelihoods should be the same because it is the same model.
        """
        states = [FlatState(4, 1), FlatState(4, 10)]
        model_a = Mixture(states, [0.5, 0.5])
        model_b = UniformMixture(states)
        observation = (0,1,2,3)
        likelihood_a = model_a.get_likelihood(observation)
        likelihood_b = model_b.get_likelihood(observation)
        self.assertAlmostEqual(likelihood_a, likelihood_b)
        log_likelihood_a = model_a.get_log_likelihood(observation)
        log_likelihood_b = model_b.get_log_likelihood(observation)
        self.assertAlmostEqual(log_likelihood_a, log_likelihood_b)

    def test_weighted_mixture_model_inequality(self):
        """
        Create mixture models that differ only in their mixing proportions.
        The likelihoods should differ in predictable ways.
        """
        states = [FlatState(4, 1), FlatState(4, 10)]
        mixture_a = Mixture(states, [0.5, 0.5])
        mixture_b = Mixture(states, [0.4, 0.6])
        observation = (1,2,3,4)
        likelihood_a = mixture_a.get_likelihood(observation)
        likelihood_b = mixture_b.get_likelihood(observation)
        self.assertTrue(likelihood_a < likelihood_b)

if __name__ == '__main__':
    unittest.main()

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


int main(int argc, const char *argv[])
{
  cout << log(0.0) << endl;
  test_mixture();
  test_finite_distn();
  test_mixture_b();
  test_uniform_mixture();
  demo_poisson_log_pmf();
  demo_geometric_log_pmf();
  return 0;
}

