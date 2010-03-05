
import unittest
import math
from StringIO import StringIO

from scipy.maxentropy import logsumexp

import ReadCoverageRef
import paramhelper

g_param_rows = [
        ('w', '0.5', float,
            'misalignment branch length'),
        ('x', '0.1', float,
            'reference branch length'),
        ('y', '0.01', float,
            'line branch length'),
        ('z', '0.0001', float,
            'mutant branch length'),
        ('low', '2', int,
            'low random coverage per base'),
        ('med', '20', int,
            'medium random coverage per base'),
        ('high', '1000', int,
            'high random coverage per base'),
        ('seqerr', '0.1', float,
            'sequencing error'),
        ('nomcoverage', '20', int,
            'nominal coverage per position'),
        ('kmulticoverages', '4', int,
            'allowed multiples of nominal coverage')]

g_param_definitions = [
        paramhelper.ParamDefinition(*defn) for defn in g_param_rows]

def get_default_string():
    """
    @return: a multiline string with default config
    """
    out = StringIO()
    for defn in g_param_definitions:
        print >> out, '#', defn.description
        print >> out, defn.name, defn.default_string
        print >> out
    return out.getvalue()


class Model:

    def __init__(self):
        pass

    def from_lines(self, lines):
        # extract the parameters
        rows = [line.split() for line in lines]
        rows = [row for row in rows if row]
        rows = [row for row in rows if not row[0].startswith('#')]
        if any(len(r)!=2 for r in rows):
            raise Exception('parameter syntax error')
        paramhelper.read_params(g_param_definitions, rows, self)
        self.states = [
                self.get_recent_state(),
                self.get_ancient_state(),
                self.get_garbage_state(),
                self.get_misaligned_state()]
        self.validate()

    def validate(self):
        # validate the parameter ranges
        if not (1 <= self.kmulticoverages <= 100):
            raise Exception('kmulticoverages is out of range')
        if not (1 <= self.low <= self.med <= self.high):
            raise Exception('expected 1 <= low <= med <= high')

    def get_recent_state(self):
        return ReadCoverageRef.HMMRecent(self.x, self.y, self.z,
                self.seqerr, self.nomcoverage, self.kmulticoverages)
    def get_ancient_state(self):
        return ReadCoverageRef.HMMAncient(self.x, self.y, self.z,
                self.seqerr, self.nomcoverage, self.kmulticoverages)
    def get_garbage_state(self):
        return ReadCoverageRef.HMMGarbage(
                self.low, self.med, self.high)
    def get_misaligned_state(self):
        return ReadCoverageRef.HMMRecent(self.w + self.x, self.y, self.z,
                self.seqerr, self.nomcoverage, self.kmulticoverages)

    def get_nstates(self):
        return len(self.states)

    def get_likelihoods(obs):
        return [s.get_likelihood(obs) for s in self.states]

    def call_polymorphism(self, obs):
        """Get the polymorphism probability.
        This is the posterior probability that the strain is homozygous
        for the non-reference base with the highest count at this position.
        @param obs: one ref base count and three non-ref base counts
        @return: the polymorphism probability
        """
        # get the prior probability of polymorphism conditional on state
        p_recent_AA = self.states[0].get_posterior_distribution(obs)[2]
        p_ancient_AA = self.states[1].get_posterior_distribution(obs)[2]
        # compute the posterior probability of a polymorphism
        posterior_polymorphism = 0
        posterior_polymorphism += p_recent * p_recent_AA
        posterior_polymorphism += p_ancient * p_ancient_AA
        # Given that a polymorphism occurred,
        # get the probability distribution over the
        # three non-reference nucleotides.
        r = self.seqerr
        log_Pr = math.log(r/4.0)
        log_PA = math.log(1 - 3*r/4.0)
        logs = [
                obs[1]*log_PA + obs[2]*log_Pr + obs[3]*log_Pr,
                obs[1]*log_Pr + obs[2]*log_PA + obs[3]*log_Pr,
                obs[1]*log_Pr + obs[2]*log_Pr + obs[3]*log_PA]
        condmaxpost = math.exp(max(logs) - logsumexp(logs))
        # get the posterior probability distribution
        maxpost = posterior_polymorphism * condmaxpost
        return maxpost


class TestModelA(unittest.TestCase):

    def test_default_string(self):
        s = get_default_string()


if __name__ == '__main__':
    unittest.main()
