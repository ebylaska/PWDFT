#include <math.h>
#include <stdio.h>

greypath2(int i) {
  int j, n, m, tn;
  if (i < 2)
    return i;
  else {
    n = (int)floor(log(1.0 * i) / log(2.0));
    tn = 1;
    for (j = 0; j < n; ++j)
      tn *= 2;
    m = 2 * tn - i - 1;
    return tn + greypath2(m);
  }
}

main() {
  int i;
  for (i = 0; i < 16; ++i)
    printf("i = %d path = %d \n", i, greypath2(i));
}
