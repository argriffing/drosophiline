#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cstdio>

#include "xgstats.H"
#include "dgrpstats.H"

#define RABC 4

using namespace std;

int main(int argc, const char *argv[])
{
  // define some model parameters
  int low = 2;
  int med = 20;
  int high = 1000;
  double x = 0.1;
  double y = 0.01;
  double z = .0001;
  double seqerr = .01;
  int nomcoverage = 10;
  int kmulticoverages = 8;
  // create the models
  HMMGarbage garbage(low, med, high);
  HMMRecent recent(x, y, z, seqerr, nomcoverage, kmulticoverages);
  HMMAncient ancient(x, y, z, seqerr, nomcoverage, kmulticoverages);
  // let the models be accessed polymorphically in a vector
  vector<Model<vector<int> > *> models;
  models.push_back(&garbage);
  models.push_back(&recent);
  models.push_back(&ancient);
  // convert observations to lilkelihoods
  vector<int32_t> rabc_buffer(RABC);
  vector<double> likelihoods(models.size());
  int incount;
  while (1)
  {
    incount = fread(&rabc_buffer[0], sizeof(int32_t), RABC, stdin);
    if (!incount)
    {
      /* if we read nothing then we are at the end */
      break;
    }
    if (incount < RABC)
    {
      fprintf(stderr, "read a partial binary chunk\n");
      return EXIT_FAILURE;
    }
    /* get the likelihoods */
    int i;
    for (i=0; i<models.size(); i++)
    {
      likelihoods[i] = models[i]->get_lik(rabc_buffer);
    }
    /* write the likelihoods */
    fwrite(&likelihoods[0], sizeof(double), likelihoods.size(), stdout);
  }
  // return success
  return EXIT_SUCCESS;
}
