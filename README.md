# extract-small-scale-systems-from-radiative-transfer-cosmological-sim

The main function which performs these functions 

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
