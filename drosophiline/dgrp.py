"""
Functions to explore the Drosophila genetic reference panel.

Chromosome names and position ranges from flybase 5.13 are included.
Positions in these chromosomes are 1-indexed.
The filtered pileup format is as follows.
Each output file has information from a single chromosome.
Input columns are
{chromosome name, chromosome position, reference base,
calls for the two alleles, coverage, literal A, A count, literal C,
C count, literal G, G count, literal T, T count, first quality score,
second quality score, third quality score}.
While the reads can be reconstructed from a pileup file,
after this filtering has been done the reads can no longer be reconstructed.
"""

import unittest

import jc69
import ambignt
import ReadCoverageRef

class DGRPError(Exception): pass

# from flybase 5.13
g_chromosome_length_pairs = (
        ('YHet', 347038),
        ('2L', 23011544),
        ('X', 22422827),
        ('3L', 24543557),
        ('4', 1351857),
        ('2R', 21146708),
        ('3R', 27905053),
        ('Uextra', 29004656),
        ('2RHet', 3288761),
        ('2LHet', 368872),
        ('3LHet', 2555491),
        ('3RHet', 2517507),
        ('U', 10049037),
        ('XHet', 204112),
        ('dmel_mitochondrion_genome', 19517))

class ChromoSkimmer:
    """
    Skim a data file.
    Each line of the file is a row of whitespace separated values,
    and the first column is the name of a chromosome.
    Remember the number of rows and the names of the chromosomes.
    All rows with the same chromosome name should be in a contiguous block.
    Yield unique chromosome names as they are encountered.
    """

    def __init__(self):
        self.name_set = set()
        self.last_name = None
        self.name_list = []
        self.linecount = 0

    def skim(self, rows):
        """
        @param rows: sequences defining rows of input values
        @raise Exception: when chromosomes are not contiguous blocks
        """
        for row in rows:
            self.linecount += 1
            name = row[0]
            if name != self.last_name:
                if name in self.name_set:
                    msg = 'each chromsome should be a contiguous block'
                    raise Exception(msg)
                yield name
                self.name_set.add(name)
                self.name_list.append(name)
                self.last_name = name

def check_chromo_monotonicity(rows):
    """
    Yield row indices.
    The first value of each input row should be the chromosome name.
    The second value of each input row should be an integer position.
    @param rows: sequences defining rows of typed input values
    @raise Exception: when chromosome positions do not monotonically increase
    """
    # scan the input file for formatting
    name_to_last_pos = {}
    for i, row in enumerate(rows):
        msg = 'the first two values of each row should be name and position'
        if len(row) < 2:
            raise Exception(msg)
        name, pos = row[0], row[1]
        if type(pos) is not int:
            raise Exception('the position should be an integer')
        last_pos = name_to_last_pos.get(name, None)
        msg = 'expected strictly increasing positions per chromosome'
        if last_pos is not None:
            if last_pos >= pos:
                raise Exception(msg)
        name_to_last_pos[name] = pos
        yield i

def filtered_pileup_row_to_typed(values):
    """
    Parse a line and do error checking.
    @param srow: a sequence of strings
    @return: a tuple with proper types
    """
    if len(values) != 16:
        raise DGRPError('expected 16 values per line')
    if values[2] not in ambignt.g_resolve_nt:
        msg = 'the reference allele should be a nucleotide code: ' + values[2]
        raise DGRPError(msg)
    msg = 'literal A, C, G, T letters were not found where expected'
    if values[5] != 'A' or values[7] != 'C':
        raise DGRPError(msg)
    if values[9] != 'G' or values[11] != 'T':
        raise DGRPError(msg)
    typed_values = [
            values[0], int(values[1]), values[2],
            values[3], int(values[4]),
            values[5], int(values[6]), values[7], int(values[8]),
            values[9], int(values[10]), values[11], int(values[12]),
            int(values[13]), int(values[14]), int(values[15])]
    return typed_values

def filtered_pileup_typed_to_obs(row):
    """
    Return a flat tuple consisting of the reference and non-reference counts.
    The reference allele count is the first element of the tuple.
    The remaining allele counts are sorted in decreasing count order.
    @param row: a sequence of values in an expected format
    @return: a sufficient statistic
    """
    name, pos, ref = row[:3]
    A, C, G, T = row[6], row[8], row[10], row[12]
    acgt_counts = (A, C, G, T)
    nt_to_count = dict(zip('ACGT', acgt_counts))
    # hack the reference allele if it is ambiguous
    if ref not in list('ACGT'):
        nts = ambignt.g_resolve_nt[ref]
        count_nt_pairs = [(nt_to_count[nt], nt) for nt in nts]
        ref_count, ref = max(count_nt_pairs)
    # get the count of the reference allele followed by decreasing counts
    R = nt_to_count[ref]
    non_ref_counts = [nt_to_count[c] for c in 'ACGT' if c != ref]
    obs = [R] + list(reversed(sorted(non_ref_counts)))
    return tuple(obs)

def get_zygosity_distribution(ref_length, child_length):
    """
    This is based on the Jukes-Cantor model on a three taxon tree.
    @param ref_length: length of the reference taxon branch
    @param child_length: length of each child taxon branch
    @return: the distribution (RR, RA, AA, AB)
    """
    p_ref_change = jc69.distance_to_probability(ref_length)
    p_child_change = jc69.distance_to_probability(child_length)
    # For now sum over all possibilities of non-reference nodes.
    # This could be done more efficiently using Felsenstein pruning,
    # but I am ignoring this for now.
    p_RR = 0.0
    p_RA = 0.0
    p_AA = 0.0
    p_AB = 0.0
    ref = 0
    for c12 in range(4):
        if c12 == ref:
            p12 = 1.0 - p_ref_change
        else:
            p12 = p_ref_change / 3.0
        for c1 in range(4):
            if c1 == c12:
                p1 = p12 * (1.0 - p_child_change)
            else:
                p1 = p12 * (p_child_change / 3.0)
            for c2 in range(4):
                if c2 == c12:
                    p2 = p1 * (1.0 - p_child_change)
                else:
                    p2 = p1 * (p_child_change / 3.0)
                # Classify the joint distribution
                # and add weight to the appropriate state.
                if c1 == ref and c2 == ref:
                    p_RR += p2
                elif c1 == ref or c2 == ref:
                    p_RA += p2
                elif c1 == c2:
                    p_AA += p2
                else:
                    p_AB += p2
    v = (p_RR, p_RA, p_AA, p_AB)
    total = sum(v)
    if abs(total - 1) > 1e-7:
        raise DGRPError('probabilities do not sum to one')
    return v

def get_foo(p_ancient, p_recent, AA_ancient, AA_recent, x2, x3, x4, r):
    """
    @param p_ancient: posterior probability of ancient state
    @param p_recent: posterior probability of recent state
    @param x2: count of a non-reference nucleotide
    @param x3: count of a non-reference nucleotide
    @param x4: count of a non-reference nucleotide
    @param r: probability of random nucleotide sequencing call
    @return: posterior polymorphism probability
    """
    pass

class Model:
    def __init__(self):
        self.expected_param_names = ['x', 'y', 'z', 'low', 'med', 'high',
                'seqerr', 'nomcoverage', 'kmulticoverages']
    def from_lines(self, lines):
        # extract the parameters
        rows = [line.split() for line in lines]
        if any(len(r)!=2 for r in rows):
            raise Exception('parameter syntax error')
        param_to_value = dict(rows)
        # check the parameters
        if set(self.expected_param_names) != set(param_to_value.keys()):
            raise Exception('invalid or missing parameter')
        # get the typed parameters
        self.x = float(param_to_value['x'])
        self.y = float(param_to_value['y'])
        self.z = float(param_to_value['z'])
        self.low = int(param_to_value['low'])
        self.med = int(param_to_value['med'])
        self.high = int(param_to_value['high'])
        self.seqerr = float(param_to_value['seqerr'])
        self.nomcoverage = int(param_to_value['nomcoverage'])
        self.kmulticoverages = int(param_to_value['kmulticoverages'])
        self.validate()
    def from_fieldstorage(self, fs):
        for name in self.expected_param_names:
            setattr(self, name, getattr(fs, name))
        self.validate()
    def validate(self):
        # validate the parameter ranges
        for name in self.expected_param_names:
            if getattr(self, name) < 0:
                raise Exception('parameters should be nonnegative')
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
    def get_misaligned_state(self, w):
        return ReadCoverageRef.HMMRecent(self.x + w, self.y, self.z,
                self.seqerr, self.nomcoverage, self.kmulticoverages)

class TestDGRP(unittest.TestCase):

    def test_foo(self):
        pass


if __name__ == '__main__':
    unittest.main()
