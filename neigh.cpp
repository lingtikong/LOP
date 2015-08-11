#include "neigh.h"

/* ------------------------------------------------------------------
 * Constructor, initialize some variables
 * ------------------------------------------------------------------ */
Neigh::Neigh(Crystal *crys)
{
  flag = 1;
  natom = 0;

  x = NULL;
  attyp = NULL;
  neig_list = NULL;

  memory = new Memory();
  crystal = crys;

  NumNei = crystal->NumNei;
  N1stN  = crystal->N1stN;

  local_list.clear(); local_prop.clear();

return;
}

/* ------------------------------------------------------------------
 * free memory
 * ------------------------------------------------------------------ */
Neigh::~Neigh()
{
  if (x) memory->destroy(x);
  if (attyp) memory->destroy(attyp);
  if (neig_list) memory->destroy(neig_list);

  local_list.clear(); local_prop.clear();

  crystal = NULL;
  delete memory;
}

/* ------------------------------------------------------------------
 * Private method to sort the local neighbor list: neighbor with
 * local_prop go first.
 * ------------------------------------------------------------------ */
void Neigh::sort_list(int id)
{
  int nlocal = local_list.size();
  for (int i=0; i<nlocal; i++)
  for (int j=i+1; j<nlocal; j++){
    if (local_prop[j] > local_prop[i]){
      double ftmp = local_prop[i];
      local_prop[i] = local_prop[j];
      local_prop[j] = ftmp;

      int itmp = local_list[i];
      local_list[i] = local_list[j];
      local_list[j] = itmp;
    }
  }

  int nmax = MIN(nlocal, NumNei);
  for (int i=0; i<nmax; i++) neig_list[id][i] = local_list[i];
  for (int i=nmax; i<NumNei; i++) neig_list[id][i] = id;

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int Neigh::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy, n, "copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}
/*----------------------------------------------------------------------------*/
