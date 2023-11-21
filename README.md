# extract-small-scale-systems-from-radiative-transfer-cosmological-sim

The main function which performs these functions 

1)  Reads the data cubes from Radiative Transfer Hydrodynamical simulation 
    the data cubes are Dark matter, Gas mass, temperature, 
    numbers densities of HI and HII and ne


  ```cpp
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
```


    
    
3) Find systems using 4-pixels connectivity in the data cube

```cpp
  comp = find_systems(xHI_cube, labels_cube, threshold, NGRID, NGRIDR, NSLICE);
```

```cpp
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
```


4) Calculate several quantites from these systems

	Position of maximum in HI field
	Total Mass
	Total dark matter mass in solar
	Size in [kpc/h]^3
	Gas temp [K]
	Radiation field Gamma
	HI enutral fraction HI_delta
	HI column densities NHI


```cpp

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

```


Derived quantites over skewers through the cube.

```cpp

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


```
