#include "snppack.h"

int main() {
	double pmm[4][4], pf[4][4][4][4];
	struct read r1;
	double v = 0.01, u = 1e-7, f = 0.1, e = 0.001;
	cout << "SNPPack" << endl;
	snppack_init(pmm, pf, v, u, f);
	//print_matrix44(pmm);
	//print_matrix44(pf[0][1]);
	r1.data[0] = 20; r1.data[1] = 2; r1.data[2] = 1; r1.data[3] = 0;
	snppack_prob(&r1, pmm, pf, e);
	cout << "Log-Weight " << r1.logweight << endl;
	for (int i = 0; i < 4; ++ i) {
		for (int j = 0; j < 4; ++ j) {
			if (i == j)
				continue;
			double sum = 0.0;
			for (int a = 0; a < 4; ++ a) {
				for (int b = a; b < 4; ++ b) {
					sum += r1.prob[i][j][a][b];
				}
			}
			cout << i << "," << j << "  " << sum << endl;
			print_matrix44(r1.prob[i][j]);
			cout << endl;
		}
	}
	return 0;
}

bool snppack_init(double pmm[4][4], double pf[4][4][4][4],
				  double v, double u, double f)
// prepare pmm and pf based on v, u, f
// pmm: P(major, minor|pi, r)
// pf:  P(major, minor, genotype|v, u, f)
//      after mutation and imbreeding
// v:   or theta, diversity
// u:   mutation rate
// f:   imbreeding coefficient
{
	double pg[4][4][4][4], ph[4][4][4][4];
	prep_pmm(pmm);
	prep_pg(pg, v);
	if (!sum_up_to_one(pg)) cout << "Error" << endl;
	//print_matrix44(pg[0][1]);
	prep_ph(ph, pg, u);
	if (!sum_up_to_one(ph)) cout << "Error" << endl;
	//print_matrix44(ph[0][1]);
	prep_pf(pf, ph, f);
	if (!sum_up_to_one(pf)) cout << "Error" << endl;
	return true;
}

bool snppack_prob(struct read *r, double pmm[4][4],
				  double pf[4][4][4][4], double e)
// calculate logweight, prob
// from read data, pmm, pf and e
// r:   read data(numbers, logweight, probabilities)
// pmm: P(major, minor|pi, r)
// pf:  P(major, minor, genotype|v, u, f)
//      after mutation and imbreeding
// e:   sequencing error rate
{
    double prg[4][4], pr[4][4][4][4];
	double pr_sum;
	r->logweight = prep_prg(prg, r->data, e);
	//print_matrix44(prg);
	pr_sum = prep_pr(pr, prg, pmm, pf);
	//print_matrix44(pr[0][1]);
	//cout << pr_sum << endl;
	calc_prob(r->prob, pr, pr_sum);
	r->logweight += log(pr_sum);
	//print_matrix44(r->prob[0][1]);
	if (!sum_up_to_one(r->prob)) cout << "Error" << endl;
	return true;
}

void prep_pmm(double pmm[4][4])
// prepare pmm based on pi and r
// pmm: P(major, minor|pi, r)
// pi:  mutation rate matrix
// r:   background distribution
{
	int i, j;
	double sum = 0;
	for (i = 0; i < 4; i ++)
		for (j = 0; j < 4; j ++)
			if (i != j) sum += pi[i] * r[i][j];
	for (i = 0; i < 4; i ++)
		for (j = 0; j < 4; j ++)
			if (i != j)
				pmm[i][j] = pi[i] * r[i][j] / sum;
			else
				pmm[i][j] = 0.0;
}

void prep_pg(double pg[4][4][4][4], double v)
// prepare pg based on v
// pg: P(major, minor, genotype|v) before mutation
// v:  or theta, diversity
// genotype AC is treated different than CA
//   simplifies prep_ph function
{
	int major, minor, a1, a2;
	for (major = 0; major < 4; major ++) {
		for (minor = 0; minor < 4; minor ++) {
			for (a1 = 0; a1 < 4; a1 ++) {
				for (a2 = 0; a2 < 4; a2 ++) {
					if (major == minor)
						pg[major][minor][a1][a2] = 0.0;
					else if (a1 == major && a2 == major)
						pg[major][minor][a1][a2] = 1.0 - v * 1.5;
					else if (a1 == major && a2 == minor
						|| a1 == minor && a2 == major
						|| a1 == minor && a2 == minor)
						pg[major][minor][a1][a2] = v * 0.5;
					else
						pg[major][minor][a1][a2] = 0.0;
				}
			}
		}
	}
}

void prep_ph(double ph[4][4][4][4], double pg[4][4][4][4], double u)
// prepare ph from pg based on u
// ph: P(major, minor, genotype|v, u)
//     after mutation, before imbreeding
// u:  mutation rate
// genotype AC is treated different than CA
//   simplifies prep_ph function
{
	int major, minor, ga1, ga2, ha1, ha2;
	double item, sum;
	for (major = 0; major < 4; major ++) {
		for (minor = 0; minor < 4; minor ++) {
			for (ha1 = 0; ha1 < 4; ha1 ++) {
				for (ha2 = 0; ha2 < 4; ha2 ++) {
					if (major == minor)
						ph[major][minor][ha1][ha2] = 0.0;
					else {
						sum = 0.0;
						// TODO: Can be simplified to inspect
						// only the major and minor alleles
						for (ga1 = 0; ga1 < 4; ga1 ++) {
							for (ga2 = 0; ga2 < 4; ga2 ++) {
								item = pg[major][minor][ga1][ga2];
								if (ha1 == ga1)
									item *= 1 + u * r[ga1][ga1];
								else
									item *= u * r[ga1][ha1];
								if (ha2 == ga2)
									item *= 1 + u * r[ga2][ga2];
								else
									item *= u * r[ga2][ha2];
								sum += item;
							}
						}
						ph[major][minor][ha1][ha2] = sum;
					}
				}
			}
		}
	}
}

void prep_pf(double pf[4][4][4][4], double ph[4][4][4][4], double f)
// prepare pf from ph based on f
// pf: P(major, minor, genotype|u, v, f) after imbreeding
// f:  imbreeding coefficient
{
	int major, minor, ha1, ha2, fa1, fa2;
	double item, sum;
	for (major = 0; major < 4; major ++) {
		for (minor = 0; minor < 4; minor ++) {
			for (fa1 = 0; fa1 < 4; fa1 ++) {
				for (fa2 = 0; fa2 < 4; fa2 ++) {
					// TODO: Can be unrolled so that you don't need the inner
					//       4x4 loop.  Item not needed.
					if (major == minor)
						pf[major][minor][fa1][fa2] = 0.0;
					else if (fa1 > fa2)
						pf[major][minor][fa1][fa2] = 0.0;
					else if (fa1 < fa2)
						pf[major][minor][fa1][fa2] = (1 - f) *
							(ph[major][minor][fa1][fa2] + ph[major][minor][fa2][fa1]);
					else {
						sum = ph[major][minor][fa1][fa1];
						for (ha1 = 0; ha1 < 4; ha1 ++) {
							for (ha2 = 0; ha2 < 4; ha2 ++) {
								if (fa1 == ha1 && fa1 != ha2 ||	fa1 != ha1 && fa1 == ha2) {
									item = ph[major][minor][ha1][ha2] * f * 0.5;
									sum += item;
								}
							}
						}
						pf[major][minor][fa1][fa2] = sum;
					}
				}
			}
		}
	}
}

double prep_prg(double prg[4][4], int data[4], double e)
// prepare prg from read date based on e
// prg:  P(read|genotype)
// data: read data (# of A, C, G, T)
// e:    sequencing error rate
{
	// TODO: These logs should be set in the init function.
	// TODO: Can be vectorized to get rid of the if(k != i) test
	double log_homo = log(1 - e), log_err = log(e / 3.0), log_het = log(0.5 - e / 3.0);
	int i, j, k;
	double log_max = -DBL_MAX;
	for (i = 0; i < 4; i ++) {
		for (j = 0; j < 4; j ++) {
			if (i > j) {
				prg[i][j] = 0.0;
				continue;
			} else if (i == j) {
				prg[i][j] = log_homo * data[i];
				for (k = 0; k < 4; k ++)
					if (k != i)
						prg[i][j] += log_err * data[k];
			}
			else {
				prg[i][j] = log_het * (data[i] + data[j]);
				for (k = 0; k < 4; k ++)
					if (k != i && k != j)
						prg[i][j] += log_err * data[k];
			}
			if(prg[i][j] > log_max) {
				log_max = prg[i][j];
			}
		}
	}
	for(i = 0; i < 4; i ++) {
		for(j = i; j < 4; j ++) {
			prg[i][j] = exp(prg[i][j] - log_max);
		} 
	}
	return log_max;
}

double prep_pr(double pr[4][4][4][4], double prg[4][4],
			 double pmm[4][4], double pf[4][4][4][4])
// prepare pr and pr_sum from prg, pmm and pf
// pr: P(read, major, minor, genotype|u, v, f, e)
// pr = prg * pf * pmm
{
	int major, minor, i, j;
	double pr_sum = 0.0;
	for (major = 0; major < 4; major ++) {
		for (minor = 0; minor < 4; minor ++) {
			for (i = 0; i < 4; i ++) {
				for (j = 0; j < 4; j ++) {
					if (major == minor)
						pr[major][minor][i][j] = 0.0;
					else {
						pr[major][minor][i][j] =
							prg[i][j] *	pf[major][minor][i][j] *
							pmm[major][minor];
						pr_sum += pr[major][minor][i][j];
					}
				}
			}
		}
	}
	return pr_sum;
}

void calc_prob(double prob[4][4][4][4], double pr[4][4][4][4], double pr_sum)
// calculate read prob from pr and pr_sum
// prob: P(major, minor, genotype|read)
// prob = pr / pr_sum
{
	int major, minor, i, j;
	for (major = 0; major < 4; major ++) {
		for (minor = 0; minor < 4; minor ++) {
			for (i = 0; i < 4; i ++) {
				for (j = 0; j < 4; j ++) {
					if (major == minor)
						prob[major][minor][i][j] = 0.0;
					else if (i > j)
						prob[major][minor][i][j] = 0.0;
					else
						prob[major][minor][i][j] =
							pr[major][minor][i][j] / pr_sum;
				}
			}
		}
	}
}

void print_matrix44(double matrix[4][4]) {
	int i, j;
	for (i = 0; i < 4; i ++) {
		for (j = 0; j < 4; j ++) {
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}

bool sum_up_to_one(double matrix[4][4][4][4]) {
	int major, minor, i, j;
	double sum;
	for (major = 0; major < 4; major ++) {
		for (minor = 0; minor < 4; minor ++) {
			if (major == minor) continue;
			sum = 0.0;
			for (i = 0; i < 4; i ++) {
				for (j = 0; j < 4; j ++) {
					sum += matrix[major][minor][i][j];
				}
			}
			if (fabs(sum - 1) > 1e-6) return false;
		}
	}
	return true;
}
