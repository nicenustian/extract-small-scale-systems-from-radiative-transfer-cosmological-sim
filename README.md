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
