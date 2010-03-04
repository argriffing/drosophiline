"""Functions dealing with the JC69 nucleotide rate matrix.

This is the continuous time Markov model defined by Jukes and Cantor in 1969.
"""

import unittest
import math

def distance_to_probability(d):
    """
    @param d: the expected number of substitutions along a branch
    @return: the probability that the endpoints are different
    """
    return (3.0 / 4.0) * (1.0 - math.exp(-(4.0 / 3.0) * d))

def probability_to_distance(p):
    """
    Note that p can also be an observed proportion of nucleotide differences.
    In this case the returned value is the maximum likelihood estimate
    of the expected number of substitutions along the branch.
    @param p: the probability that the endpoints of a branch are different
    @return: the expected number of substitutions along the branch
    """
    if p >= (3.0 / 4.0):
        return float('inf')
    return -(3.0 / 4.0) * math.log(1.0 - (4.0 / 3.0) * p)

def get_ml_distance(sa, sb):
    """
    Use a closed form maximum likelihood estimator.
    @param sa: a nucleotide sequence
    @param sb: a nucleotide sequence
    @return: the estimated expected number of changes per position
    """
    assert len(sa) == len(sb)
    assert set(sa+sb) <= set('ACGT')
    n = len(sa)
    n2 = sum(1 for x, y in zip(sa, sb) if x != y)
    n1 = n - n2
    if n2 >= 3 * n1:
        return float('inf')
    mle = -(3.0/4.0) * math.log((3.0*n1 - n2) / (3.0*n1 + 3.0*n2))
    return mle


class TestJC69(unittest.TestCase):

    def test_conversion_special_cases(self):
        d = probability_to_distance(0.0)
        self.assertAlmostEqual(d, 0)
        d = probability_to_distance(3.0 / 4.0)
        self.assertEqual(d, float('inf'))
    
    def test_conversion(self):
        probs = (0.0, 0.2, 0.4, 0.6)
        for p_expected in probs:
            d = probability_to_distance(p_expected)
            p_observed = distance_to_probability(d)
            self.assertAlmostEqual(p_expected, p_observed)

    def test_ml(self):
        sa = 'AAAAA'
        sbs = ('AAAAA', 'AAAAC', 'AAACC', 'AACCC')
        probs = (0.0, 0.2, 0.4, 0.6)
        for p_expected, sb in zip(probs, sbs):
            d_expected = probability_to_distance(p_expected)
            d_observed = get_ml_distance(sa, sb)
            self.assertAlmostEqual(d_expected, d_observed)


if __name__ == '__main__':
    unittest.main()
