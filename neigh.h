#ifndef NEIGH_LIST_H
#define NEIGH_LIST_H

#include "memory.h"
#include <vector>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "crystal.h"

#define MAXLINE 5120
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

using namespace std;

class Neigh {
public:
  Neigh(Crystal *);
  ~Neigh();

  int flag; // 0, read successfully; otherwise, failed
  int natom;
  int *attyp;
  double **x;
  int **neig_list;
  int NumNei, N1stN;

  Memory *memory;
  Crystal *crystal;
 
  virtual double get_alat(int){return 1.;}

  char oneline[MAXLINE];

protected:
  std::vector<int> local_list;
  std::vector<double> local_prop;
  void sort_list(int);
  int count_words(const char *);
};
#endif
