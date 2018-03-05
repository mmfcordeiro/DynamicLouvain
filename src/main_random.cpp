#include <stdlib.h>
#include <iostream>

# if defined(WIN32)
#include <stdio.h>
#include <process.h>
#include <stdio.h>
# endif

#include "graph_binary.h"
#include "community.h"


using namespace std;

char *outfile  = NULL;


int
main(int argc, char **argv) {
# if defined(WIN32)
    srand(time(NULL)+_getpid());
# else
    srand(time(NULL)+getpid());
# endif

  int n = atoi(argv[1]);
  int degree = atoi(argv[2]);

  vector<int> v (n);
  for(unsigned int i = 0; i < n; i++) {
        v[i] = rand();
  }

  for (unsigned int i=0 ; i<n ; i++) {
    for (unsigned int j=0 ; j<n ; j++) {
      int r  = rand()%n;
      if (r<degree)
        cout << v[i] << " " << v[j] << endl;
    }
  }
}
