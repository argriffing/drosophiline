#pragma once
#include <iostream>
#include <math.h>
#include <float.h>

using namespace std;

int main();
bool snppack_init(double pmm[4][4], double pf[4][4][4][4],
				  double v, double u, double f);
bool snppack_prob(struct read *r, double pmm[4][4],
				  double pf[4][4][4][4], double e);

void prep_pmm(double pmm[4][4]);
void prep_pg(double pg[4][4][4][4], double v);
void prep_ph(double ph[4][4][4][4], double pg[4][4][4][4], double u);
void prep_pf(double pf[4][4][4][4], double ph[4][4][4][4], double f);
double prep_prg(double prg[4][4], int data[4], double e);
double prep_pr(double pr[4][4][4][4], double prg[4][4],
			 double pmm[4][4], double pf[4][4][4][4]);
void calc_prob(double prob[4][4][4][4], double pr[4][4][4][4],
			   double pr_sum);

void print_matrix44(double matrix[4][4]);
bool sum_up_to_one(double matrix[4][4][4][4]);

const double pi[4] = {0.285500, 0.214500, 0.214500, 0.285500};
const double r[4][4] =
	{{-0.917070, 0.198534, 0.467593, 0.250944},
	{0.264249, -1.111169, 0.223312, 0.623608},
	{0.622367, 0.223312, -1.109139, 0.263459},
	{0.250944, 0.468525, 0.197941, -0.917409}};

struct read {
	int data[4];
	double logweight;
	double prob[4][4][4][4];
};
