#include "orderpm.h"

int main(int narg, char **arg)
{

  OrderPM *opm = new OrderPM(narg, arg);
  delete opm;

return 0;
}
