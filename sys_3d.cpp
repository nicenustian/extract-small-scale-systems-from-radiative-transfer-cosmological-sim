#include "sys_3d.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <time.h>

using namespace std;

/*The main function which performs these functions 

1)  Reads the data cubes from Radiative Transfer Hydrodynamical simulation 
    the data cubes are Dark matter, Gas mass, temperature, 
    numbers densities of HI and HII and ne
    
2) Find systems using 4-pixels connectivity in the data cube  

3) Calculate several quantites from these systems

	Position of maximum in HI field
	Total Mass
	Total dark matter mass in solar
	Size in [kpc/h]^3
	Gas temp [K]
	Radiation field Gamma
	HI enutral fraction HI_delta
	HI column densities NHI
*/

void systems(std::string base_dir, std::string fiducial_redshift, std::string redshift,  
double threshold, const long NGRID, const long NDOM)
{

  FILE   *infile, *outfile;
  float  *gammaHI_cube, *delta_cube, *rhog_cube, *temp_cube, *nHI_cube, 
  *nHII_cube, *ne_cube, *xHI_cube, *ratio_cube, *dm_cube;
  int    *labels_cube;
  int    comp;
  double z, z_fid, dx_proper, dx_comov, dV_proper, dV_comov;
  double *size, *mass, *dm_mass;

  //number of skewers for each system paralell to  box axes
  int sk= 6;
  int sk_random = 32;
  int *xi, *yi, *zi;
  double *ratio_max, *nHI_max;
  double *temp_w, *gamma_w, *xHI_delta_w, *delta_w, *NHI;

  //files to read
  string base_dir_in = base_dir;

  const long NSLICE  = NGRID;
  const long NGRID3D = NGRID*NGRID*NGRID;
  const long NGRIDR  = NGRID * NGRID * NSLICE;
  const double dxh   = BOXSIZE/double(NGRID); //in [Mpc/h]
  const long dm_particle_mass = 105; //DM particle mass in solar
  
  //convert redshift string to numerical value and compute units
  //Note that conversion from code units to proper requires the "Redshift"
  z = atof(redshift.c_str());
  z_fid = atof(fiducial_redshift.c_str());
  
  //The dx in proper units [cm]
  dx_proper = BOXSIZE/NGRID/HUBBLE_h/(1+z)*MPCTOCM;
  //The dx in [Mpc/h]
  dx_comov = dx_proper*HUBBLE_h*(1+z_fid)/MPCTOCM;

  //in units of cm-3
  dV_proper = dx_proper * dx_proper * dx_proper;
  dV_comov = dx_comov * dx_comov * dx_comov;


  cout<<"Fiducial z = "<< z_fid  << "  Redshift = " << z <<endl<<
	     "BOXSIZE="<<BOXSIZE<< "  NGRID = "<<NGRID<<" NSLICE ="<<NSLICE<<" NGRIDR="<<NGRIDR<<endl
	     <<" Threshold in xHI = "<< threshold<<endl<<
	     "dx = "<<dx_proper<<"[cm],   "<<dx_comov<<" [Mpc/h]"<<endl;
  
  // Create an output string stream
  std::ostringstream streamObj;
  //Add double to stream
  streamObj << threshold;
  // Get string from output string stream
  std::string strObj = streamObj.str(); 

  //the input filenames for grid files
  string gasdatagridfilename  = base_dir +"/"+ "gas_data.grid.z="+fiducial_redshift;
  string gammagridfilename    = base_dir +"/"+ "gammaHI.grid.z="+fiducial_redshift;
  string dmfilename           = base_dir +"/density_dm_1024.z=07.8890";

  //Output file names
  string sysfilename          = base_dir +"/"+ "sys_3D_vthres_z="+fiducial_redshift;
  string labelsfilename       = base_dir +"/"+ "labels_vthres_z="+fiducial_redshift;
  string pdffilename          = base_dir +"/"+ "pdf_vthres_z="+fiducial_redshift;
  string profilesfilename     = base_dir +"/"+ "profiles_vthres_z="+fiducial_redshift; 

  //cubes to store the varabiles from simulation
  delta_cube = new float[NGRIDR];
  rhog_cube = new float[NGRIDR];
  gammaHI_cube = new float[NGRIDR];
  temp_cube  = new float[NGRIDR];
  nHI_cube   = new float[NGRIDR];
  nHII_cube  = new float[NGRIDR];
  ne_cube    = new float[NGRIDR];
  xHI_cube   = new float[NGRIDR];
  ratio_cube = new float[NGRIDR];
  labels_cube = new int[NGRIDR];

  infile  = fopen(gammagridfilename.c_str(),"rb");
  cout<<"Reading file "<<gammagridfilename.c_str()<< endl;
  fread(gammaHI_cube, sizeof(float), NGRIDR, infile);
  cout << gammaHI_cube[0] <<" "<<gammaHI_cube[NGRIDR-1] << endl;
  fclose(infile);


  infile  = fopen(dmfilename.c_str(),"rb");
  cout<<"Reading file "<<dmfilename.c_str()<< endl;
  fread(dm_cube, sizeof(float), NGRIDR, infile);
  cout << dm_cube[0] <<" "<<dm_cube[NGRIDR-1] << endl;
  fclose(infile);


  infile  = fopen(gasdatagridfilename.c_str(),"rb");
  cout<<"Reading file "<<gasdatagridfilename.c_str()<< endl;

  fread(delta_cube, sizeof(float), NGRIDR, infile);
  fseek ( infile,   NGRID3D*4, SEEK_SET );
  fread(rhog_cube, sizeof(float), NGRIDR, infile);
  fseek ( infile,   2*NGRID3D*4, SEEK_SET );
  fread(temp_cube, sizeof(float), NGRIDR, infile);
  fseek ( infile, 3*NGRID3D*4, SEEK_SET );
  fread(nHI_cube, sizeof(float), NGRIDR, infile);
  fseek ( infile, 4*NGRID3D*4, SEEK_SET );
  fread(nHII_cube, sizeof(float), NGRIDR, infile);
  fseek ( infile, 5*NGRID3D*4, SEEK_SET );
  fread(ne_cube, sizeof(float), NGRIDR, infile);  
  fclose(infile);
  
   ////////////////////////////////////////////////////////////////////////////////////

  //find others drived quantities
  #pragma omp parallel for  
  for (int i=0;i<NGRIDR;i++)
  {
	xHI_cube[i]   = nHI_cube[i] / (nHI_cube[i] + nHII_cube[i]);
	ratio_cube[i] = ( gammaHI_cube[i] * nHI_cube[i] ) / (alphaB(temp_cube[i]) * ne_cube[i] * nHII_cube[i]);
  }

  cout << "Done with drived quantities" << endl;
  cout <<delta_cube[0]<<"  "<< rhog_cube[0] <<" "<< temp_cube[0] <<" "<<nHI_cube[0] <<" "<< nHII_cube[0] <<" "<< ne_cube[0] <<" "<< xHI_cube[0] <<" "<< ratio_cube[0] <<" "   << endl;
  cout << delta_cube[NGRIDR-1] << " " << rhog_cube[NGRIDR-1] <<" "<< temp_cube[NGRIDR-1] <<" "<<nHI_cube[NGRIDR-1] <<" "<< nHII_cube[NGRIDR-1] <<" "<< ne_cube[NGRIDR-1] <<" "<<xHI_cube[NGRIDR-1]<<" "<< ratio_cube[NGRIDR-1] <<" "<< endl;

  ////////////////////////////////////////////////////////////////////////////////////
  //Find connected regions through xHI threshold  

  comp = find_systems(xHI_cube, labels_cube, threshold, NGRID, NGRIDR, NSLICE);

  cout << endl << "DIM = " << NGRID<< " x "<<NGRID<<" x "<<NSLICE<<endl;
  cout << "Done with labels to systems found "<<comp << endl;

  outfile = fopen(labelsfilename.c_str(),"wb");
  cout << "Writing the labels file "<<labelsfilename.c_str() << endl;
  fwrite(labels_cube, sizeof(int), NGRIDR, outfile);
  fclose(outfile);
  
  ////////////////////////////////////////////////////////////////////////////////////
  mass  = new double[comp];
  size = new double[comp];

  xi = new int[comp];
  yi = new int[comp];
  zi = new int[comp];

  temp_w = new double[comp];
  gamma_w = new double[comp];
  xHI_delta_w = new double[comp];
  delta_w  = new double[comp];
  NHI  = new double[comp];

  ratio_max = new double[comp];
  nHI_max = new double[comp];

  for (int i=0;i<comp;i++){
	mass[i] = 0;
	size[i] = 0;
	temp_w[i] = 0;
	gamma_w[i] = 0;
	xHI_delta_w[i] = 0;
	delta_w[i] = 0;
	NHI[i] = 0;
	ratio_max[i] = 0;
    nHI_max[i] = 0;
  }

  //find the location of maximum in NHI which serves as the center of a system
  find_argmax(nHI_cube, labels_cube, nHI_max, comp,  xi, yi, zi, NGRID,  NGRIDR,  NSLICE);

  //find the gas mass for each system in solar mass
  find_total( rhog_cube, labels_cube, mass, comp, dV_proper/SOLARTOGRAMS, NGRIDR);

  //find the dark matter mass for each system in solar mass
  find_total( dm_cube, labels_cube, dm_mass, comp, dm_particle_mass, NGRIDR);

  //the sizes of systems in pixels. dV is in units of [kpc/h]^3 comoving
  cout<<" DV factor=" <<dV_comov*1e9<<endl;
  find_size(labels_cube, size, comp, dV_comov*1e9, NGRIDR);

  //The nHI weighted quantities and column densitites
  skewers(delta_cube, temp_cube, nHI_cube, xHI_cube, gammaHI_cube,
             delta_w, temp_w, NHI, xHI_delta_w, gamma_w,
             labels_cube, comp, sk, xi, yi, zi, dx_proper, NGRID, NSLICE);

  find_max(ratio_cube, labels_cube, ratio_max, comp, NGRIDR);

  for (int i = 0; i < comp; i++){
 	 if (i<100){
		cout<<"x,y,z="<<xi[i]<<","<<yi[i]<<","<<zi[i]<<"\t"
		     <<"log Mass[solar]="<<log10(mass[i])<<"log Mass_dm[solar]="<<log10(dm_mass[i])<<"\t"<<"ratio="<<mass[i]/dm_mass[i];
		     <<"Size [kpc/h]^3="<<size[i]<<"\t"
		     <<"temp ="<<temp_w[i]<<"\t"
		     <<"gamma ="<<gamma_w[i]<<"\t"
		     <<"xHI_delta ="<<xHI_delta_w[i]<<"\t"
		     <<"delta="<<delta_w[i]<<"\t"
		     <<"NHI="<<NHI[i]<<"\t"
		     <<"P/R="<<ratio_max[i]<<endl;
		}
  }

 ////////////////////////////////////////////////////////////////////////////
 
 outfile = fopen(sysfilename.c_str(),"wb");
 cout << "Writing the systems data file "<<sysfilename.c_str() << endl;

 for (int i = 0; i < comp; i++)
 {
	//the position is the index for max. in nHI
	fwrite(&xi[i], sizeof(int),1, outfile);
	fwrite(&yi[i], sizeof(int),1, outfile);
	fwrite(&zi[i], sizeof(int),1, outfile);

	//the mass and size in kpc/h  
	fwrite(&mass[i], sizeof(double),1, outfile);
	fwrite(&size[i], sizeof(double),1, outfile);
	
	fwrite(&NHI[i], sizeof(double),1, outfile);	
	
	fwrite(&delta_w[i], sizeof(double), 1, outfile);
	fwrite(&temp_w[i], sizeof(double),1, outfile);
	fwrite(&xHI_delta_w[i], sizeof(double),1, outfile);
	fwrite(&gamma_w[i], sizeof(double),1, outfile);
	fwrite(&ratio_max[i], sizeof(double),1, outfile);
 }

 fclose(outfile);

 ////////////////////////////////////////////////////////////////////////////
 delete [] xi;
 delete [] yi;
 delete [] zi;
 delete [] mass;
 delete [] size;
 delete [] NHI;
 delete [] temp_w;
 delete [] gamma_w;
 delete [] xHI_delta_w;
 delete [] delta_w;
 delete [] ratio_max;

 delete [] labels_cube;
 delete [] temp_cube;
 delete [] nHI_cube;
 delete [] nHII_cube;
 delete [] xHI_cube;
 delete [] ne_cube;
 delete [] gammaHI_cube;
 delete [] delta_cube;
 delete [] rhog_cube;
 delete [] ratio_cube;

 
 return;
}


/* integral of a quantity, such as mass  if integral of gas density is taken with fact as dV*/
void find_total(float *data, int *labels, double *total, int comp, double fact, const long NGRIDR)
{
  for (int i=0;i<comp;i++) total[i] = 0;

  for (int i=0;i<NGRIDR;i++)
          if ( (labels[i]!=0)) total[labels[i]-1] += data[i];
  
  for (int i=0;i<comp;i++) total[i] *= fact;
}

/*size in pixels but factor convert it to units of kpc/h**3 in this particular*/
void find_size(int *labels, double *size, int comp, double fact, const long NGRIDR)
{
	for (int i=0;i<comp;i++) size[i] = 0;
  
	for (int i=0;i<NGRIDR;i++)
		if ( (labels[i]!=0)) size[labels[i]-1]++;

	for (int i=0;i<comp;i++) size[i] *= fact;

}


/*Find maximum value in the datacube*/
void find_max(float *data, int *labels, double *max, int comp, const long NGRIDR)
{
        //the max array should contain all negatives. As all our values would be positive or zero. 
        for (int i=0;i<comp;i++) max[i] = -1e200;
        	
	for (int i=0;i<NGRIDR;i++){
        		if ( (labels[i]!=0) && (data[i] > max[labels[i]-1]) && isfinite(data[i]) ) 
				max[labels[i]-1] = data[i];
        
        }
}



/*Find the index of maximum values in a data cube*/
void find_argmax(float *data, int *labels, double *max, int comp, int *x, int *y, int *z, const long NGRID, const long NGRIDR, const long NSLICE)
{

  //the max array should contain all negatives. As all our values would be positive or zero.
  for (int i=0;i<comp;i++) 
  {  
	max[i] = -1e200; 
        x[i] = y[i] = z[i] = -1;
  }

  for (int i=0;i<NGRIDR;i++)
  {

	int ix, iy, iz;
        if ( (labels[i]!=0) && (data[i] > max[labels[i]-1]) && isfinite(data[i]) ) 
	{
		max[labels[i]-1] = data[i];

		//corresponding (x,y,z) index from index, column-major
		iz = i/NGRID/NSLICE;
	        iy = (i - iz*NGRID*NSLICE)/NSLICE;
        	ix = i - iz*NGRID*NSLICE - iy*NSLICE;

		x[labels[i]-1] = ix;
		y[labels[i]-1] = iy;
		z[labels[i]-1] = iz;


	}
  }
}






/*Find minimum value in the datacube*/
void find_min(float *data, int *labels, double *min, int comp, const long NGRIDR)
{
	//the max array should contain all negatives. As all our values would be positive or zero. 
	for (int i=0;i<comp;i++) min[i] = 1e200;  

	for (int i=0;i<NGRIDR;i++)
        	if ( (labels[i]!=0) && (data[i] < min[labels[i]-1]) && isfinite(data[i]) ) min[labels[i]-1] = data[i];

}


/*get weighted average quantity over the skewers  for each systems, 
such as HI clumn density and cosmic density*/
void skewers (float *delta, float *temp, float *nHI, float *xHI, float *gamma, 
             double *delta_w, double *temp_w, double *NHI, double *xHI_delta_w, double *gamma_w,
	     int *labels, int comp, int sk, int *x, int *y, int *z, double dx, const long NGRID, const long NSLICE)
{

	int index, index_x, index_y, index_z, j;
    double *theta, *phi, *vecx, *vecy, *vecz;
	double delta_sum, temp_sum, gamma_sum, xHI_delta_sum,  nHI_sum;
	bool flag;

    theta = new double [sk];
	phi   = new double [sk];
	vecx = new double [sk];
	vecy = new double [sk];
	vecz = new double [sk];

	//set the seed to keep same results every time
	srand(1234);
        for (int i=0;i<sk;i++)
        {
                theta[i] =  ((double) rand() / RAND_MAX)*M_PI;
                phi[i]   =  ((double) rand() / RAND_MAX)*2*M_PI;

                //unit vector in random direction
                vecx[i] = cos(phi[i])*sin(theta[i]);
                vecy[i] = sin(phi[i])*sin(theta[i]);
                vecz[i] = cos(theta[i]);
        }

         //loop over the number of objects
        for (int ob=0; ob<comp; ob++ )
        {
                //loop over N skewers starting from peak given by (x, y, z) of each object
                for (int i=0; i<sk; i++)
                {
                	j = 0; 
			flag = true;
			
			delta_sum = 0;
			temp_sum=0;
			gamma_sum=0;
			nHI_sum=0;
			xHI_delta_sum=0;

			//set the peak in absorbers as a starting point
                        index_x = x[ob];
			index_y = y[ob]; 
			index_z = z[ob];

            index = index_z*NGRID*NSLICE + index_y*NSLICE + index_x;

			//do the positive and then negative directions of unit vector to complete a skewer spreading across boundary through the center.
			//positive unit vector direction
			//cout <<" positive" << endl;
            while ( flag  && 
				index_x>=0 && 
				index_x<NGRID && 
				index_y>=0 && 
				index_y<NGRID &&  
				index_z>=0 && 
				index_z<NGRID )
                        {
				if (labels[index]>0)
				{
					nHI_sum += nHI[index];

					//nHI weighted quantities
					delta_sum += (nHI[index] * delta[index]);
					temp_sum += (nHI[index] * temp[index]);
					gamma_sum += (nHI[index] * gamma[index]);
					xHI_delta_sum += (nHI[index] * xHI[index]);///delta[index]);

  
					++j;
                	                //the next index of current skewer
                        	        index_x = ( vecx[i]*j  + x[ob] );  
					index_y = ( vecy[i]*j  + y[ob] );  
					index_z = ( vecz[i]*j  + z[ob] );
                                	index = index_z*NGRID*NSLICE + index_y*NSLICE + index_x;
				}

				else flag=false;
                        }

			j=1; 
			flag=true;
                        index_x = ( -vecx[i]*j  + x[ob] );  
			index_y = ( -vecy[i]*j  + y[ob] );  
			index_z = ( -vecz[i]*j  + z[ob] );
                        index = index_z*NGRID*NSLICE + index_y*NSLICE + index_x;

			//cout <<" negative" << endl;
       			//negative unit direction
			while (flag && 
				index_x>=0 && 
				index_x<NGRID && 
				index_y>0 && 
				index_y<NGRID && 
				index_z>=0 && 
				index_z<NGRID)
                        {
				if (labels[index]>0)
				{

					nHI_sum += nHI[index];

					delta_sum += (nHI[index] * delta[index]);
					temp_sum += (nHI[index] * temp[index]);
					gamma_sum += (nHI[index] * gamma[index]);
					xHI_delta_sum += (nHI[index] * xHI[index]);///delta[index]);

        	                        ++j;
                	                //the staring index of current skewer
                        	        index_x = ( -vecx[i]*j  + x[ob] );  
					index_y = ( -vecy[i]*j  + y[ob] );  
					index_z = ( -vecz[i]*j  + z[ob] );
                                	index = index_z*NGRID*NSLICE + index_y*NSLICE + index_x;
				}

				else flag=false;
                        }

			delta_w[ob] += (delta_sum/nHI_sum);
			temp_w[ob] += (temp_sum/nHI_sum);
			gamma_w[ob] += (gamma_sum/nHI_sum);
			xHI_delta_w[ob] += (xHI_delta_sum/nHI_sum);
			NHI[ob] += (nHI_sum * dx);


		}//skewer loop
		//the average over sk skewers 
		 delta_w[ob] /= sk;
                 temp_w[ob] /= sk;
                 gamma_w[ob] /= sk;
                 xHI_delta_w[ob] /= sk;
                 NHI[ob] /= sk;

	}//object loop

}



/*Find all pixels which are connected using 4-pixel connectiveity 
in a data  cube by picking a certain threshold*/
int find_systems(float *data, int *labels, double thres_v, 
const long NGRID, const long NGRIDR, const long NSLICE)
{
  bool *thres;
  int comp=0, index; 
  thres= new bool[NGRIDR];
  
  // threshold data
  for (int i=0;i<NGRIDR;i++)
  {
    if (data[i] > thres_v) thres[i] = true;
  	else thres[i] = false;
  
  	//also set all the labels to zero. zero means not labelled.
    labels[i] = 0;
  }
 
 
  for (int i = 0; i < NGRID; ++i)
  	for (int j = 0; j < NGRID; ++j)
  		for (int k = 0; k < NSLICE; ++k)
  		{
			index = i*NGRID*NSLICE+j*NSLICE+k;
			//column major index for array
			if (!labels[index] && thres[index])
				dfs(i, j, k, ++comp, thres, labels, NGRID, NSLICE);
		}

  return comp;
}


/*recurvsilvey labelling the pixels based on 4 pixels connectivity*/
void dfs(int i, int j, int k, int current_label, bool *data, int *labels, const long NGRID, const long NSLICE)  
{
   int index = i*NGRID*NSLICE+j*NSLICE+k;
   if (i < 0 || i == NGRID) return; // out of bounds
   if (j < 0 || j == NGRID) return; // out of bounds
   if (k < 0 || k == NSLICE) return; // out of bounds

   if (labels[index] || !data[index]) return; // already labeled or not marked with 1 

   // mark the current cell
   labels[index] = current_label;

   // recursively mark the neighbors, 4 pixel connectivity
   for (int dir = 0; dir < 6; ++dir)
       dfs(i + dx[dir], j + dy[dir],  k + dz[dir], current_label, data, labels, NGRID, NSLICE);

}


//Case-A/B recombination coefficient.  From Hui+Gnedin 1997.
//Input temperature in K, return in cm^3/s
double alphaB(double T)
{
	double x = 2*157807/T;
	//return 1.269e-13*pow(x, 1.503) / pow( 1.0 + pow(x/0.522, 0.470) , 1.923 ); 
	return 2.753e-14*pow(x,1.500)/pow((1.0+pow((x/2.740),0.407)),2.2420);
}



