#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXN 10000000

char bool_s(bool a) {
  if(a)
    return 'T';
   
  return 'F';
}

int main(int argc, char** argv) {
  double *a = (double*) malloc(sizeof(double) * MAXN),
      *b = (double*) malloc(sizeof(double) * MAXN),
      *c = (double*) malloc(sizeof(double) * MAXN), 
      *sumb = (double*) malloc(sizeof(double) * MAXN), 
      *sumc = (double*) malloc(sizeof(double) * MAXN);
  long n;
  FILE *f, *g;

  if(argc <= 1) {
    printf("function must have arguments..\n");
    return 2;
  }

  g = fopen("E_result.dat", "w");

  for(long i = 0; i < MAXN; i++ ) {
    sumc[i] = 0; 
    sumb[i] = 0;
    a[i] = 0;
    b[i] = 0;
    c[i] = 0;
  }

  printf("opening files\n");

  for(int i = 1; ; i++) {
    if(i >= argc) break;

    printf("%s\n", argv[i]);
    f = fopen(argv[i], "r");
    if(f == NULL) {
       printf("error opening file\n");
       return 1;
    }

    n = 0;
    while(!feof(f)) {
      if(n >= MAXN)
        continue;

      fscanf(f,"%lf", a + n);
      fscanf(f,"%lf", b + n);
      fscanf(f,"%lf", c + n);

      sumc[n] += c[n];
      sumb[n] += b[n];

      n++;
    }

  }

  // quadratic interpolation of the first point
  sumb[0] = sumb[1] + (-a[1] + a[0])*((sumb[1] - sumb[2])/(a[1] - a[2]) + (((-sumb[1] + sumb[2])/(a[1] - a[2]) + (sumb[2] - sumb[3])/(a[2] - a[3]))*(-a[2] + a[0]))/(-a[1] + a[3]));
  sumc[0] = sumc[1] + (-a[1] + a[0])*((sumc[1] - sumc[2])/(a[1] - a[2]) + (((-sumc[1] + sumc[2])/(a[1] - a[2]) + (sumc[2] - sumc[3])/(a[2] - a[3]))*(-a[2] + a[0]))/(-a[1] + a[3]));

  for(long i = 0; i < n - 1; i++) {
    fprintf(g,"%e %e %e\n",a[i], sumb[i]/sumb[0], -sumc[i]/sumb[0]);
  }


  fclose(g);
  return 0;
}
