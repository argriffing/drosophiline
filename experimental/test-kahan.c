#include <stdio.h>
#include <stdlib.h>

/*
 * GCC presumably considers overoptimization of kahan summation to be a bug,
 * because GCC uses Kahan's PARANOIA program for regression testing.
 */

float naive_sum()
{
  float accum = 0;
  int n = 100000000;
  int i;
  for (i=0; i<n; i++)
  {
    accum += 1.0;
  }
  return accum;
}

float kahan_sum()
{
  float accum = 0;
  int n = 100000000;
  int i;
  float c = 0;
  for (i=0; i<n; i++)
  {
    float y = 1.0 - c;
    float t = accum + y;
    c = (t - accum) - y;
    accum = t;
  }
  return accum;
}

int main(int argc, const char **argv)
{
  printf("naive sum: %f\n", naive_sum());
  printf("kahan sum: %f\n", kahan_sum());
  return EXIT_SUCCESS;
}
