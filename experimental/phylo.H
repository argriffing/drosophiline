#ifndef PHYLO_H
#define PHYLO_H

#include <vector>

/*
 * @param d: the expected number of substitutions along a branch
 * @return: the probability that the endpoints are different
 */
double jc69_distance_to_probability(double d);

/*
 * This is based on the Jukes-Cantor model on a three taxon tree.
 * @param ref_length: length of the reference taxon branch
 * @param child_length: length of each child taxon branch
 * @return: the distribution (RR, RA, AA, AB)
 */
std::vector<double> get_zygosity_distn(double ref_len, double child_len);

#endif
