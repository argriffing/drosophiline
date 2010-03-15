#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <iostream>

using namespace std;

template <typename T>
double sum_of_sqrts(T begin, T end)
{
  double accum = 0;
  T it;
  for (it=begin; it!=end; it++)
    accum += *it;
  return accum;
}

int main(int argc, const char **argv)
{
  double arr[3] = {1.0, 2.0, 3.5};
  vector<double> myvect(arr, arr+3);
  list<double> mylist(arr, arr+3);
  cout << sum_of_sqrts(myvect.begin(), myvect.end()) << endl;
  cout << sum_of_sqrts(mylist.begin(), mylist.end()) << endl;
  return EXIT_SUCCESS;
}
