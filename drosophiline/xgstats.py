
import random
import unittest
from math import log

import scipy.stats
from scipy.special import gammaln

def binomial_log_pmf(observed_n, max_n, p_success):
    #TODO special cases
    accum = 0
    accum += gammaln(max_n + 1)
    accum -= gammaln(observed_n + 1)
    accum -= gammaln((max_n - observed_n) + 1)
    accum += observed_n * log(p_success)
    accum += (max_n - observed_n) * log(1.0 - p_success)
    return accum

def geometric_log_pmf(observed_n, pr):
    """
    @param observed_n: the number of completed events
    @param pr: the probability of quitting
    """
    if pr == 0.0:
        return float('-inf')
    if pr == 1.0:
        if observed_n:
            return float('-inf')
        else:
            return log(pr)
    return observed_n * log(1.0 - pr) + log(pr)

def poisson_log_pmf(observed_n, expected_n):
    if not expected_n:
        if observed_n:
            return float('-inf')
        else:
            return 0.0
    accum = 0
    accum += observed_n * log(expected_n)
    accum -= expected_n
    accum -= gammaln(observed_n+1)
    return accum

def multinomial_log_pmf(distribution, counts):
    """
    This should be in scipy.stats but it isn't.
    @param distribution: the distribution over classes
    @param counts: the observed counts over classes
    """
    # check for a degeneracy
    for d, c in zip(distribution, counts):
        if c and not d:
            return float('-inf')
    n = sum(counts)
    # initialize the log probability mass
    accum = 0
    # add the contribution of n to the multinomial coefficient
    if n > 1:
        accum += gammaln(n+1)
    # add the contribution of the counts to the multinomial coefficient
    accum -= sum(gammaln(count+1) for count in counts if count > 1)
    # add the contribution of probabilities
    for p, count in zip(distribution, counts):
        if count:
            accum += count * log(p)
    return accum

def weights_to_distribution(weights):
    for weight in weights:
        assert weight >= 0
    weight_list = list(weights)
    total = sum(weight_list)
    assert total > 0
    return [weight / total for weight in weight_list]

def random_weighted_int(weights):
    """
    @param weights: an ordered sequence of nonnegative weights summing to one
    """
    x = random.random()
    accum = 0
    for i, w in enumerate(weights):
        accum += w
        if x < accum:
            return i

def weighted_choice(weight_state_pairs):
    if not weight_state_pairs:
        raise ValueError('no choices available')
    if len(weight_state_pairs) == 1:
        weight, state = weight_state_pairs[0]
        return state
    total_weight = sum(weight for weight, state in weight_state_pairs)
    assert total_weight > 0
    r = random.uniform(0, total_weight)
    cumulative_weight = 0
    for weight, state in weight_state_pairs:
        if weight < 0:
            raise ValueError('weights must be non-negative')
        cumulative_weight += weight
        if r < cumulative_weight:
            return state
    raise ValueError('no choice was made')


class TestXgcode(unittest.TestCase):

    def test_poisson_log_pmf(self):
        observed_n = 60
        expected_n = 20
        likelihood = scipy.stats.poisson.pmf(observed_n, expected_n)
        expected = log(likelihood)
        observed = poisson_log_pmf(observed_n, expected_n)
        self.assertAlmostEqual(expected, observed)

    def test_geometric_log_pmf(self):
        obs = 5
        pr = 0.1
        scipy_result = log(scipy.stats.geom.pmf(obs, pr, loc=-1))
        util_result = geometric_log_pmf(obs, pr)
        self.assertAlmostEqual(scipy_result, util_result)


if __name__ == '__main__':
    unittest.main()
