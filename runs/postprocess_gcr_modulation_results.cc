#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Spectrum;

const double J0 = 1.0;
const double mu = 1.8;
const double T0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT;
const double specie = SPECIES_PROTON_BEAM;

// Unmodulated spectrum
//inline double unmod_spectrum(double T) {return J0 * pow(T / T0, -mu);};
//inline double unmod_spectrum(double T) {return 15.0 * pow(T / T0, -2.6)/ (1.0 + 4.3*pow(T / T0, -0.8) + 0.14*pow(T / T0, -2.85) + 1.0*pow(T / T0, -2.0));};
inline double unmod_spectrum(double T) {return 12.0 * pow(T / T0, -2.6)/ (1.0 + 5.3*pow(T / T0, -1.22) + 1.3*pow(T / T0, -2.8) + 0.0087*pow(T / T0, -4.32));};

int main(int argc, char** argv)
{
   std::string infilename1 = "output/modulated_gcr_spectrum_SSSW_fixednewKdriftsNLGC_10au.dat";
   std::string outfilename = "output/modulated_gcr_spectrum_SSSW_fixednewKdriftsNLGC_10au_PP.dat";
   std::string line;
   int i, N = 100;
   int sum_c1[N];
   double energy1[N], distro1[N], sum_w1[N];
   double scal, v, p2, energy0 = SPC_CONST_CGSM_MEGA_ELECTRON_VOLT;

// Open input distro file
   std::ifstream input_spectrum_file1(infilename1);

// Read first two lines of distro file
   std::getline(input_spectrum_file1, line);
   std::getline(input_spectrum_file1, line);

// Read data
   for(i = 0; i < N; i++) {
      input_spectrum_file1 >> energy1[i];
      input_spectrum_file1 >> distro1[i];
      input_spectrum_file1 >> sum_w1[i];
      input_spectrum_file1 >> sum_c1[i];
   };

// Close input cartesian distro file
   input_spectrum_file1.close();

// Open output distro file
   std::ofstream output_spectrum_file(outfilename);

   std::cout << unmod_spectrum(T0) << std::endl;

// Output data
   output_spectrum_file << std::setprecision(8);
   scal = Vel(Mom(T0 / unit_energy_particle, specie), specie) * unmod_spectrum(T0);
   for(i = 0; i < N; i++) {
      p2 = Sqr(Mom(energy1[i] / unit_energy_particle, specie));
      v = Vel(Mom(energy1[i] / unit_energy_particle, specie), specie);
      output_spectrum_file << std::setw(20) << energy1[i] / energy0;
      output_spectrum_file << std::setw(20) << unmod_spectrum(energy1[i]);
      output_spectrum_file << std::setw(20) << p2 * distro1[i];
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

   return 0;
};
