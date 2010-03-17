#include <vector>
#include <cmath>

#include "phylo.H"

double jc69_distance_to_probability(double d)
{
  return (3.0 / 4.0) * (1.0 - exp(-(4.0 / 3.0) * d));
}

std::vector<double> get_zygosity_distn(double ref_len, double child_len)
{
  double p_ref_change = jc69_distance_to_probability(ref_len);
  double p_child_change = jc69_distance_to_probability(child_len);
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
  std::vector<double> distn;
  distn.push_back(p_RR);
  distn.push_back(p_RA);
  distn.push_back(p_AA);
  distn.push_back(p_AB);
  return distn;
}
