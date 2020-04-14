#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib> // for exit function
#include <cmath>
#include <algorithm>
#include <string>


// Declare variables here

std::string filename;
double Beta;
std::vector<double> N_PHOTONS;
double add;
double lambda1;
double lambda2;
double NA;
std::vector<double> Answer;
std::vector<double> EnergyList;
int MaxE;
bool found;

// Declare Structures here

struct Properties{
    std::string ElementName;
    std::vector<int> ZtotV;//
    std::vector<double> Molar_MassV;
    std::vector<double> NV; // index of refraction
    std::vector<double> densityV;
    std::vector<std::string> Element;
    std::string Ele;
    double density;
    double Molar;
    double n;
    int Ztot;
    double Ion;
    int INDEX;
};
Properties Element;

struct Input{
    double BeamEnergy;
    std::vector<double> BeamEnergyV;
    double Energy;
    std::vector<double> EnergyV;
    int EventID;
    std::vector<int> EventIDV;
    std::string ParticleName;
    std::vector<std::string> ParticleNameV;
    std::string CreatorProcessName;
    std::vector<std::string> CreatorProcessNameV;
    double Time;
    std::vector<double> TimeV;
    
};
Input I;

std::ifstream inFile;
// Define Sub Functions here

void beamtoArray(Input &I){
    I.BeamEnergyV.push_back(I.BeamEnergy);
    I.EnergyV.push_back(I.Energy);
    I.EventIDV.push_back(I.EventID);
    I.ParticleNameV.push_back(I.ParticleName);
    I.CreatorProcessNameV.push_back(I.CreatorProcessName);
    I.TimeV.push_back(I.Time);
    return;
    
}
double stoppingPower(double Energy, double Beta, double Ionization, double NA, double Ztot){
    double Tau = Energy/0.511;
    double F_Beta = ((1 - pow(Beta,2))/2)*(1 + (pow(Tau,2)/8) - (2*Tau + 1)*log(2));
    double G_Beta = log(3.61*pow(10,5)*Tau*sqrt(Tau + 2)) + F_Beta;
    double dE_dx_col = 5.08*pow(10,-31)*NA/pow(Beta,2)*G_Beta - log(Ionization); // Collision stopping power
    //std::cout << dE_dx_col <<std::endl;
    long double F_rad = (pow(1.60*pow(10,-19),4)*pow(8.99*pow(10,9),2)*pow(Ztot,2)*NA*(Energy + .511))/(137*pow(.511,2));
    double ratio = Ztot*(Energy + .511)/800;
    double dE_dx_rad = dE_dx_col*ratio;
    //std::cout << dE_dx_rad <<std::endl;
    double dE_dx_total = dE_dx_col + dE_dx_rad;
    if(dE_dx_total < 0){
        dE_dx_total = 0;
    }
    //std::cout << dE_dx_total << std::endl;
    return dE_dx_total;
}
double PhotonsProducedPerDistance(double Energy, double lambda1, double lambda2, double n, double Beta){
    double Crit_Angle = acos(1/(Beta*n));
    if(std::isnan(Crit_Angle)){
        Crit_Angle = 0;
    }
    double dN_dx = (2*acos(-1)/137)*(-1/lambda2 + 1/lambda1)*pow(sin(Crit_Angle),2);
    return dN_dx;
}
double getClosest(double, double, double);

double findClosest(std::vector<double> EnergyList, int n, double Energy){
    // Corner cases
    if (Energy <= EnergyList[0])
        return EnergyList[0];
    if (Energy >= EnergyList[n - 1])
        return EnergyList[n - 1];
  
    // Doing binary search
    int i = 0, j = n, mid = 0;
    while (i < j) {
        mid = (i + j) / 2;
  
        if (EnergyList[mid] == Energy)
            return EnergyList[mid];
  
        /* If Energy is less than EnergyList element,
            then search in left */
        if (Energy < EnergyList[mid]) {
  
            // If Energy is greater than previous
            // to mid, return closest of two
            if (mid > 0 && Energy > EnergyList[mid - 1])
                return getClosest(EnergyList[mid - 1],
                                  EnergyList[mid], Energy);
  
            /* Repeat for left half */
            j = mid;
        }
  
        // If Energy is greater than mid
        else {
            if (mid < n - 1 && Energy < EnergyList[mid + 1])
                return getClosest(EnergyList[mid],
                                  EnergyList[mid + 1], Energy);
            // update i
            i = mid + 1;
        }
    }
  
    // Only single element left after search
    return EnergyList[mid];
}
  
// Method to compare which one is the more close.
// We find the closest by taking the difference
// between the Energy and both values. It assumes
// that val2 is greater than val1 and Energy lies
// between these two.
double getClosest(double val1, double val2,
               double Energy)
{
    if (Energy - val1 >= val2 - Energy)
        return val2;
    else
        return val1;
}

int find(std::vector<double> EnergyList, int n, double value){ // find index
  int index = -1;

    for(int i=0; i<n; i++){
       if(EnergyList[i]==value){
         index=i;
         break;
       }
    }
   return index;
 }

int findElement(Properties Element){
    for(int i=0;i<Element.Element.size();i++){
        if(Element.ElementName == Element.Element[i] || Element.Element[i] == "end"){
            Element.INDEX = i;
            break;
        }
    }
    return Element.INDEX;
}

double Ionization(Properties Element){
    if(Element.ElementName == "Water"){
        Element.Ion = 4.312; // units eV
    }
    else if(Element.ElementName == "CS2"){
        Element.Ion = 93.8;
    }
    else if(Element.ElementName == "Benzene"){
        Element.Ion = 12.1;
    }
    else if(Element.ElementName == "Sugar"){
        Element.Ion = 7.39;
    }
    else if(Element.ElementName == "GaAs"){
        Element.Ion = 331.8;
    }
    else if(Element.ElementName == "SiO2"){
        Element.Ion = 122.3;
    }
    else if(Element.ElementName == "SiC"){
        Element.Ion = 162.0;
    }
    else if(Element.Ztot == 1){
        Element.Ion = 19;
    }
    else if(Element.Ztot >= 2 && Element.Ztot <= 13){
        Element.Ion = 11.2 + 11.7*Element.Ztot;
    }
    else{
        Element.Ion = 52.8 + 8.71*Element.Ztot;
    }
    return Element.Ion;
}

void toArray(Properties &Element){
    Element.Molar_MassV.push_back(Element.Molar);
    Element.Element.push_back(Element.Ele);
    Element.NV.push_back(Element.n);
    Element.ZtotV.push_back(Element.Ztot);
    Element.densityV.push_back(Element.density);
    return;
}

void toVariable(Properties &Element){
    if(Element.INDEX != 100){
        Element.Ztot = Element.ZtotV[Element.INDEX];
        Element.n = Element.NV[Element.INDEX];
        Element.Molar = Element.Molar_MassV[Element.INDEX];
        Element.density = Element.densityV[Element.INDEX];
        NA = (6.02*pow(10,23)*Element.density*pow(10,6)/Element.Molar)*Element.Ztot;
        Element.Ion = Ionization(Element);
        
    }
    else{
        std::cout << "ERROR Element/Compound Not Found. Please Enter number of electrons, Molar Mass, density, index of refraction, and Mean Excitation Energy" <<std::endl;
        std::cin >> Element.Ztot >> Element.Molar >> Element.density >> Element.n >> Element.Ion;
        NA = (6.02*pow(10,23)*Element.density*pow(10,6)/Element.Molar)*Element.Ztot;
    }
}
void checkinFile(std::string filename){
    if(!inFile){//check to insure file opens
        std::cout<<"Unable to open file "<<filename<<std::endl;
        exit(1);
    }
    else{
        std::cout << "Successfully read " << filename << std::endl;
        
    }
}
