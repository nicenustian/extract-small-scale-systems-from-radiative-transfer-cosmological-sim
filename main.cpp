//INPUT: base directory, redshift, base file name for output

#include <fstream>
#include <iostream>
#include <omp.h>
#include "sys_3d.h"
using namespace std;

int main()
{
	int no_models = 5;
	string base_dir2 = "/";
    string base_dir;
    string model[no_models]= {"zreion12/G3e-13/","zreion8/G3e-13/"};
    omp_set_num_threads(N_THREADS);
  
    FILE *infile;
    string file;

	//dir, redshift fiducial, redshift, NGRID, NDOM, xHI threshold
	systems(base_dir, fiducial_redshift[i][j], redshift[i][j], thresholds[i][j], 1024, 32);

  return 0;
}
