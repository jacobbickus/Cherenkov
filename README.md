# Cherenkov
Cherenkov Radiation Number of Photons produced Calculation 
Input Description:
Script takes up to 4 arguments when called. Inputs In order: 
1. inFile Filename of the .dat file that user would like to read electron energies from. 
2. Min Detector Photon Wavelength(nm) used for dN/dx integration. Default set to 350. 
3. Max Detector Photon Wavelength(nm) used for dN/dx integration Default set to 550.
4. Name of Detector Medium. Default set to Water. 

Script also requires ElementProperities.txt to be in the same directory. ElementProperities.txt goes into dE/dx calculation and provides the correct number of electrons, density, and index of refraction

Output Description:
Script Outputs:
1. NPhotons.txt which is a lookup table for the script to read the electron energy and number of photons emitted at that energy. Not really necessary but was a good check to make sure number of photons emitted per energy is reasonable.
2. Solution.dat which is the table with the electron energies found in .dat file from grasshopper and the number of photons that energy emitted. 
3. Prints to Screen total number of photons emitted

ISSUES:

1. Update stopping power to include additional factor (shell correction etc...) The stopping power formulas are valid at high energies as long as the ym/M <<1 where y = 1/sqrt(1-beta^2). Additionally at higher energies the accuracy of the stopping power formula is decreased because the other physical factors play a larger role. Additionally at lower energies the stopping power would technically become negative however since this is not physically possible the code sets the stopping power = 0. This is not important however since this is a factor below 25 keV which is not of interest for cherenkov radiation.
2. Some Elements require a user input for index of refraction.
3. Ionization energies are approximated with equations from James Turner "Atoms, Radiation and Radiation Protection". Additionally Compounds not included will require input from user.
4. Compounds not included in Element Properties.txt will ask for user input. 
