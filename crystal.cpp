#include "crystal.h"

/* -----------------------------------------------------------------------------
 * Constructor, actually does nothing but set NULL to pointers
 * ---------------------------------------------------------------------------*/
Crystal::Crystal()
{
  memory = new Memory();

  pos_ideal = NULL;
}

/* -----------------------------------------------------------------------------
 * deconstructor, free memory
 * ---------------------------------------------------------------------------*/
Crystal::~Crystal()
{
  if (pos_ideal) memory->destroy(pos_ideal);

  delete memory;

return;
}

/* ------------------------------------------------------------------
 * Private method to get the norm of the vector of size 3
 * ------------------------------------------------------------------ */
double Crystal::norm2(double *a)
{
  double la = 0.;
  for (int i=0; i<3; i++) la += a[i]*a[i];

return la;
}

/* ------------------------------------------------------------------
 * Private method to get cosine between two vectors of size 3
 * ------------------------------------------------------------------ */
double Crystal::angle(double *a, double *b)
{
  double la = norm2(a);
  double lb = norm2(b);
  if (la < 1.e-8 || lb < 1.e-8) return 0.;

  double ab = 0.;
  for (int i=0; i<3; i++) ab += a[i]*b[i];

return ab/sqrt(la*lb);
}

/* ------------------------------------------------------------------
 * Private method to get the cross product of two vectors of size 3.
 * a = b x c
 * ------------------------------------------------------------------ */
void Crystal::cross(double *a, double *b, double *c)
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];

return;
}

/* ------------------------------------------------------------------
 * Private method to normalize a vector of size 3
 * ------------------------------------------------------------------ */
void Crystal::normalize(double *a)
{
  double inv_len = 1./sqrt(norm2(a));
  for (int i=0; i<3; i++) a[i] *= inv_len;
}
