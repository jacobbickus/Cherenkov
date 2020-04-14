#include "CherenkovPhotons.hh"

// -------------------------------------------------------------------------------------------------------------------------------------- //
int main(int argc, char* argv[]){
    
// --------------------------- Defaults ---------------------------//
    
    if(argc < 2){
        std::cout << "ERROR: Must Input Dat File to Read Electron Energies!"<<std::endl;
        return EXIT_FAILURE;
    }
    else if(argc < 3){
        filename = argv[1];
        lambda1 = 350;
        lambda2 = 550;
        Element.ElementName = "Water";
        
    }
    else if(argc == 3){
        std::cout << "ERROR: Missing INPUT"<<std::endl;
        return EXIT_FAILURE;
    }
    else if(argc < 5){
        filename = argv[1];
        std::string lambda1s = argv[2];
        std::string lambda2s = argv[3];
        lambda1 = std::stod(lambda1s,nullptr);
        lambda2 = std::stod(lambda2s,nullptr);
        Element.ElementName = "Water";
    }
    else{
        filename = argv[1];
        std::string lambda1s = argv[2];
        std::string lambda2s = argv[3];
        Element.ElementName = argv[4];
        lambda1 = std::stod(lambda1s,nullptr);
        lambda2 = std::stod(lambda2s,nullptr);
        
    }
    // Convert lambda to cm
    lambda1 = lambda1*1*pow(10,-7);
    lambda2 = lambda2*1*pow(10,-7);
    
// --------------------------- Read Grasshopper dat File ---------------------------//
    
    inFile.open(filename);
    checkinFile(filename);
    while(inFile && !inFile.eof()){ // read columns
        inFile.ignore(100,'\n'); // ignore headers
        inFile >> I.BeamEnergy >> I.Energy >> I.EventID >> I.ParticleName >> I.CreatorProcessName >> I.Time;
        beamtoArray(I);
    }
    inFile.close(); // close inFile
    
// ---------------------------Read Element Properties File ---------------------------//
    
    filename = "ElementProperties.txt";
    inFile.open(filename); // Open Element Properties text
    checkinFile(filename);

    while(inFile && !inFile.eof()){
        inFile >> Element.Ztot >> Element.Ele >> Element.Molar >> Element.density >> Element.n;
        toArray(Element);
    }
    inFile.close();
    Element.INDEX = findElement(Element);
    toVariable(Element);
    std::cout << "Element " << Element.ElementName << " found."<<std::endl;
    
    if(Element.n == -1){
        std::cout << "Index of Refraction not found in ElementProperties.txt. Please input index of refraction" <<std::endl;
        std::cin >> Element.n;
    }
    
// --------------------------- Write Energy/Photons Table to NPhotons.txt File ---------------------------//
    
    std::ofstream outdata;
    outdata.open("NPhotons.txt");
    MaxE = std::round(*max_element(I.EnergyV.begin(),I.EnergyV.end()))*1000;
    for(int i=0;i<MaxE;i++){
        EnergyList.push_back(i);
        EnergyList[i] = EnergyList[i]/1000; // starts at 1 keV
        Beta = sqrt(1 - pow(1/((EnergyList[i]/.511) + 1),2));
        double  dE_dx_total = stoppingPower(EnergyList[i],Beta,Element.Ion,NA,Element.Ztot);
        double dN_dx = PhotonsProducedPerDistance(EnergyList[i],lambda1,lambda2,Element.n,Beta);
            N_PHOTONS.push_back(dN_dx/dE_dx_total);
        if(std::isinf(N_PHOTONS[i])){
            N_PHOTONS[i] = 0;
        }
    }
    for(int i=0;i<N_PHOTONS.size();i++){
        Answer.push_back(0);
    }
    add = 0;
    for(int i=2;i<N_PHOTONS.size();i++){
        add = N_PHOTONS[i - 1]*(EnergyList[i] - EnergyList[i - 1]);
        Answer[i] = add + Answer[i - 1];
        outdata << EnergyList[i] << "\t"<< Answer[i] << std::endl;
    }
    outdata.close();
    std::cout << std::endl << "Successfully Wrote Number of Photons Per Electron Energy into Table found in NPhotons.txt"<<std::endl;
    // File created with Number of Photons
    
// --------------------------- Write Electron Energies from Grasshopper dat file and Number of Photons into Solution.dat ---------------------------//
    outdata.open("Solution.dat");

    int n = EnergyList.size();
    double TOTALPHOTONS = 0;
    int TOTALELECTRONS = 0;
    // Integrate dN/dx/dE/dx with respect to energy
    for(int i=0;i<I.EnergyV.size();i++){
        if(I.ParticleNameV[i] == "e-"){
            double value = findClosest(EnergyList,n,I.EnergyV[i]);
            int index = find(EnergyList,n,value);
            TOTALPHOTONS = TOTALPHOTONS + Answer[index];
            TOTALELECTRONS = TOTALELECTRONS + 1;
            outdata << I.EnergyV[i] << "\t" << "\t" << Answer[index] <<std::endl;
        }
    }
    outdata.close();
    std::cout << std::endl << "Successfully Wrote Electron Energies and Number of Photons Emitted in Solution.dat" <<std::endl;
    std::cout <<std::endl << "The total number of electrons found was " << TOTALELECTRONS <<std::endl;
    if(!std::isnan(TOTALPHOTONS)){
        std::cout << "The total number of Photons Emitted from Cherenkov was " <<TOTALPHOTONS<<std::endl;
        
    }
    else{
        std::cout << "Cherenkov Radiation Not Prevalent in Material" <<std::endl;
    }
    return 0;
}
