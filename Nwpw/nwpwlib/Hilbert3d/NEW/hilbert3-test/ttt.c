#include <math.h>
#include <stdio.h>

greypath(int i) {
  int j, n, m, tn;
  if (i < 2)
    return i;
  else {
    n = (int)floor(log(1.0 * i) / log(2.0));
    tn = 1;
    for (j = 0; j < n; ++j)
      tn *= 2;
    m = 2 * tn - i - 1;
    return tn + greypath(m);
  }
}

main() {
  int i;
  for (i = 0; i < 16; ++i)
    printf("i = %d path = %d \n", i, greypath(i));
}
