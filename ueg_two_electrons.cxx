/*  This is a piece of code evaluating the correlation energy of the
    uniform electron gas with only two electrons (of opposite spin).
    The implemented equations are motivated by the coupled cluster
    ansatz, viz ccd. For a two electron system the ccsd is able to
    yield the exact correlation energy. (for further details see
    arXiv:2103.06788 and references therein).
    The two parameters are the value for rs ( which defines the
    electron density of this system) and the number of virtual hartree
    fock states.

    For this model system with two electrons, only 5 different types
    of diagram exist, named Mp2, ppl, madelung, ring, and quadratic.
    The logic of the indivual energy contributions follows PRL123, 156401

*/


#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <math.h>
#include <fstream>

using namespace std;

double madelung, volume;

bool meshsort(array<double,3> i, array<double,3> j){
  return (i[0]*i[0]+i[1]*i[1]+i[2]*i[2]<j[0]*j[0]+j[1]*j[1]+j[2]*j[2]);
}


class Mesh{
  public:
    Mesh(int Nv_){
      Nv = Nv_;
      grid.resize(Nv);
      sqrLength.resize(Nv);
      length.resize(Nv);
      hfEnergy.resize(Nv);
    }
    ~Mesh(){
    }
    int Nv;
    vector<array<double,3> > grid;
    vector<double> sqrLength;
    vector<double> length;
    vector<double> hfEnergy;
    void initialize(double boxLength);
    void eigenEnergies();
};

class Container{
  public:
   //Constructor
   Container(int length_){
      length = length_;
	  	data.resize(length);
      for (size_t i(0); i < length; ++i) data[i] = 0.0;
   }
   // Destructor
   ~Container(){
   }
   vector<double> data;
   int length;
   void intitializeMp2Amplitudes(vector<double> sqrLength, vector<double> hfEnergy);
   inline double integral(double sqrLength);
   double getEnergy(vector<double> sqrLength);
   void structureFactorRing(
     vector<double> sqrLength, vector<double> hfEnergy, vector<double> amplitudes
   );
   void structureFactorQuadratic(
     double oldEnergy, vector<double> hfEnergy, vector<double> amplitudes
   );
   void structureFactorPpl(
     vector<array<double,3>> grid, vector<double> hfEnergy, vector<double> amplitudes
   );
   void structureFactorPplSplit(
     bool icut, bool jcut, int NvSmall, vector<array<double,3>> grid,
     vector<double> hfEnergy, vector<double> amplitudes
   );
   void add(vector<double> data);
};


void Container::add(vector<double> data_){
  for (int i(0); i < data.size(); i++)
    data[i] += data_[i];
}


inline double Container::integral(double sqrLength){
  return 4.0*M_PI/sqrLength/volume;
}

void Container::intitializeMp2Amplitudes(vector<double> sqrLength, vector<double> hfEnergy){
  for ( int i(0); i < length; i++)
    data[i] = integral(sqrLength[i])/(-2.0*madelung - 2.0*hfEnergy[i]);
}

double Container::getEnergy(vector<double> sqrLength){
  double energy(0.0);
  for ( int i(0); i < length; i++)
    energy += integral(sqrLength[i])*data[i];
  return energy;
}

void Container::structureFactorRing(
  vector<double> sqrLength, vector<double> hfEnergy, vector<double> amplitudes
){
  for( int i(0); i < length; i++)
    data[i] = 2.0*integral(sqrLength[i])*amplitudes[i]/(-2.0*madelung - 2.0*hfEnergy[i]);
}

void Container::structureFactorQuadratic(
  double oldEnergy, vector<double> hfEnergy, vector<double> amplitudes
){
  for (int i(0); i < length; i++)
    data[i] = amplitudes[i]*oldEnergy/(-2.0*madelung - 2.0*hfEnergy[i]);
}

void Container::structureFactorPpl(
  vector<array<double,3>> grid, vector<double> hfEnergy, vector<double> amplitudes
  ){
  for (int i(0); i < length; i++){
    double value(0.0);
    double ix(grid[i][0]), iy(grid[i][1]), iz(grid[i][2]);
    for (int j(0); j < length; j++){
      if ( i == j ) continue;
      double ijx(ix-grid[j][0]), ijy(iy-grid[j][1]), ijz(iz-grid[j][2]);
      double convolutedKernel = integral(ijx*ijx + ijy*ijy + ijz*ijz);
      value += amplitudes[j]*convolutedKernel;
    }
    value += madelung*amplitudes[i];  // the i=j contribution
    data[i] = value/(-2.0*madelung - 2.0*hfEnergy[i]);
  }
}


void Mesh::initialize(double boxLength){
  int iMax(0);
  int elements(0);
  double reziprocalVector(2.0*M_PI/boxLength);
  vector<int> shell;
  int iMaxCube(pow(Nv,1.0/3.0));
  iMaxCube++;

  while(elements< 8*iMaxCube*iMaxCube*iMaxCube){
    iMax++;
    for (int i(-iMax); i <= iMax; i++)
      for( int j(-iMax); j <= iMax; j++)
        for (int k(-iMax); k <= iMax; k++){
          if ( (abs(i) < iMax) && (abs(j) < iMax) && (abs(k) < iMax) ) continue;
          int absVal(i*i+j*j+k*k);
          if(find(shell.begin(), shell.end(), absVal) == shell.end()) {
            shell.push_back(absVal);
          }
          elements++;
     }
  }
  sort(shell.begin(), shell.end());
  vector<int> shellCounter;
  for (int v(0); v < shell.size(); v++){
    elements = 0;
    for(int w(1); w <= iMax; w++){
      for (int i(-w); i <= w; i++)
        for( int j(-w); j <= w; j++)
          for (int k(-w); k <= w; k++){
            if ( (abs(i) < w) && (abs(j) < w) && (abs(k) < w) ) continue;
            int absVal(i*i+j*j+k*k);
            if (absVal <= shell[v]) elements++;
       }
    }
    shellCounter.push_back(elements);
  }
  int threshold(-1);
  for (int i(0); i < shell.size(); i++){
    if ( shellCounter[i] == Nv ){
      threshold = shell[i];
    }
  }
  if ( threshold < 0 ){
    for ( int i(0); i < shell.size(); i++){
      if ( shellCounter[i] > Nv){
        cout << "Possible values: " << shellCounter[i-1] << " " << shellCounter[i] << endl;
        break;
      }
    }
    throw invalid_argument("wrong number of virtuals");
  }
  elements = 0;
  for(int w(1); w <= iMax; w++){
    for (int i(-w); i <= w; i++)
      for (int j(-w); j <= w; j++)
        for (int k(-w); k <= w; k++){
          if ( (abs(i) < w) && (abs(j) < w) && (abs(k) < w) ) continue;
          double absval(i*i+j*j+k*k);
          if ( absval > threshold) continue;
          grid[elements][0] = i*reziprocalVector;
          grid[elements][1] = j*reziprocalVector;
          grid[elements][2] = k*reziprocalVector;
          sqrLength[elements] = absval*reziprocalVector*reziprocalVector;
          length[elements] = sqrt(absval)*reziprocalVector;
          elements++;
    }
  }
  sort(grid.begin(), grid.end(), meshsort);
  sort(length.begin(), length.end());
  sort(sqrLength.begin(), sqrLength.end());
}

void Mesh::eigenEnergies(){
  for ( int i(0); i < Nv; i++){
    hfEnergy[i] = 0.5*sqrLength[i] - 4.0*M_PI/volume/sqrLength[i];
  }
}

int main(int argc, char *argv[]){
  cout.precision(12);
  if (argc < 3 ){
    throw invalid_argument("Usage: give me rs +  number of virtual");
  }
  double rs(atof(argv[1]));
  int argNv(atoi(argv[2]));

  cout << "This two electron UEG result is for rs: " << rs << endl;
  int iterations(30);
  int Nv(argNv);
  Mesh grid(Nv);

  // initialize system
  volume = 8.0/3.0*M_PI*pow(rs,3.);
  double boxLength(pow(volume, 1.0/3.0));
  // Note: the madelung number is not unabiguously defined/used in the field
  madelung = 2.8372/boxLength;
  cout << "With volume " << volume << endl;
  cout << "madelung " << madelung << endl;
  grid.initialize(boxLength);
  grid.eigenEnergies();

  // eval Mp2 Energy
  Container amplitudes(Nv);
  amplitudes.intitializeMp2Amplitudes(grid.sqrLength, grid.hfEnergy);
  double mp2Energy(amplitudes.getEnergy(grid.sqrLength));

  double energy(mp2Energy);
  cout << "Iteration 0" << " Energy: " << energy << endl;

  // ccsd iterations
  Container newAmplitudes(Nv), container(Nv);
  for (int iter(1); iter < iterations; iter++){
    newAmplitudes.intitializeMp2Amplitudes(grid.sqrLength, grid.hfEnergy);
    container.structureFactorRing(grid.sqrLength, grid.hfEnergy, amplitudes.data);
    newAmplitudes.add(container.data);
    container.structureFactorQuadratic(-energy, grid.hfEnergy, amplitudes.data);
    newAmplitudes.add(container.data);
    container.structureFactorPpl(grid.grid, grid.hfEnergy, amplitudes.data);
    newAmplitudes.add(container.data);
    container.structureFactorQuadratic(-madelung*3.0, grid.hfEnergy, amplitudes.data);
    newAmplitudes.add(container.data);
    energy = newAmplitudes.getEnergy(grid.sqrLength);
    for (int i(0); i < amplitudes.data.size(); i++)
      amplitudes.data[i] = newAmplitudes.data[i];
    cout << "Iteration " << iter << " Energy: " << energy << endl;
  }


  // evalute energy contribution from different channels using converged amplitutes
  Container structureFactor(Nv);

  structureFactor.structureFactorRing(grid.sqrLength, grid.hfEnergy, amplitudes.data);
  double ringEnergy(structureFactor.getEnergy(grid.sqrLength));

  structureFactor.structureFactorQuadratic(-energy, grid.hfEnergy, amplitudes.data);
  double quadraticEnergy(structureFactor.getEnergy(grid.sqrLength));

  structureFactor.structureFactorPpl(grid.grid, grid.hfEnergy, amplitudes.data);
  double pplEnergy(structureFactor.getEnergy(grid.sqrLength));

  structureFactor.structureFactorQuadratic(-madelung*3.0, grid.hfEnergy, amplitudes.data);
  double madelungEnergy(structureFactor.getEnergy(grid.sqrLength));


  cout << "Number of iterations: " << iterations << endl;
  cout << "MP2 Energy for " << Nv  << " virtuals:       " <<
           mp2Energy << "\n";
  cout << "ring Energy for " << Nv  << " virtuals:      " <<
           ringEnergy << "\n";
  cout << "quadratic Energy for " << Nv  << " virtuals: " <<
           quadraticEnergy << "\n";
  cout << "ppl Energy for " << Nv  << " virtuals:       " <<
           pplEnergy << "\n";
  cout << "Madelung Energy for " << Nv  << " virtuals:  " <<
           madelungEnergy << "\n";
  cout << "Total energy for   " << Nv << " virtuals:    "
       << (mp2Energy+ringEnergy+quadraticEnergy+pplEnergy+madelungEnergy)
       << "\n";

  cout << "---------------------------------------------" << endl;

  return 0;
}
