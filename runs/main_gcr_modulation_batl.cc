#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_solarwind.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include <iostream>
#include <iomanip>

using namespace Spectrum;

int main(int argc, char** argv)
{

   DataContainer container;


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create a simulation object
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<SimulationWorker> simulation;
   simulation = CreateSimulation(argc, argv);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle type
//----------------------------------------------------------------------------------------------------------------------------------------------------

    int specie = SPECIES_PROTON_BEAM;
    simulation->SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Background
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Initial time
   double t0 = 0.0;
   container.Insert(t0);

// Origin
   container.Insert(gv_zeros);

// Velocity
   container.Insert(gv_zeros);

// Magnetic field
   container.Insert(gv_zeros);

// Effective "mesh" resolution
   double dmax = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(dmax);

   simulation->AddBackground(BackgroundServerBATL(), container, "3d__var_3_t30000330_n00456411");

//----------------------------------------------------------------------------------------------------------------------------------------------------

    container.Clear();

// Initial time
    double init_t = 0.0;
    container.Insert(init_t);

    simulation->AddInitial(InitialTimeFixed(), container);

    container.Clear();

// Initial position
    GeoVector init_pos(3 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid, 0.0, 0.0);
    container.Insert(init_pos);

    simulation->AddInitial(InitialSpaceFixed(), container);

    container.Clear();

// Lower bound for momentum
    double momentum1 = Mom(10.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
    container.Insert(momentum1);

// Upper bound for momentum
    double momentum2 = Mom(5000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
    container.Insert(momentum2);

// Log bias
    bool log_bias = true;
    container.Insert(log_bias);

    simulation->AddInitial(InitialMomentumThickShell(), container);

    container.Clear();

    //// Inner boundary ----------------------------------------------------

    // Max crossings
    int max_crossings_inner = 1;
    container.Insert(max_crossings_inner);

    // Action
    std::vector<int> actions_inner;
    actions_inner.push_back(-1);
    actions_inner.push_back(-1);
    container.Insert(actions_inner);

    // Origin
    container.Insert(gv_zeros);

    // Radius
    double inner_boundary = 2 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
    container.Insert(inner_boundary);

    simulation->AddBoundary(BoundarySphereReflect(), container);

    container.Clear();

    //// Outer boundary ----------------------------------------------------
    
    // Max crossings
    int max_crossings_outer = 1;
    container.Insert(max_crossings_outer);

    // Action
    std::vector<int> actions_outer;
    actions_outer.push_back(0);
    actions_outer.push_back(0);
    container.Insert(actions_outer);

    // Region Index
    int region_index = 1;
    container.Insert(region_index);

    // Region Value
    double region_value = 0.99;
    container.Insert(region_value);

    simulation->AddBoundary(BoundaryRegionAbsorb(), container);

    container.Clear();

    //// Time limit ----------------------------------------------------

    // Not needed because this class sets the value to -1
    int max_crossings_time = 1;
    container.Insert(max_crossings_time);

    // Action
    std::vector<int> actions_time;
    actions_time.push_back(-1);
    actions_time.push_back(-1);
    container.Insert(actions_time);

    // Max duration of the trajectory
    double maxtime = 3 * 60.0 * 60.0 * 24.0 * 365.0 / unit_time_fluid;
    container.Insert(maxtime);

    simulation->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion
//----------------------------------------------------------------------------------------------------------------------------------------------------

    container.Clear();

// Diffusion coefficient normalization factor
   double kap0 = 1.5e21 / unit_diffusion_fluid;
   container.Insert(kap0);

// Rigidity normalization factor
   double T0_d = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0_d);

// Magnetic field normalization factor
   double r0 = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(r0);

// Power law slope for rigidity
   double pow_law_rigidity = 1.0;
   container.Insert(pow_law_rigidity);

// Power law slope for radius field
   double pow_law_r = 1.2;
   container.Insert(pow_law_r);

// Ratio of kappa_perp to kappa_para
   double kap_rat = 0.05;
   container.Insert(kap_rat);

// Flow dependency flag
   int stream_dep_idx = 0;
   container.Insert(stream_dep_idx);

// Upstream velocity before shock
   double u_up = 4.0e7 / unit_velocity_fluid;
   container.Insert(u_up);

// Shock width
   double w_sh = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(w_sh);

// Shock strength
   double s_sh = 4.0;
   container.Insert(s_sh);

// Pass ownership of "diffusion" to simulation
   simulation->AddDiffusion(DiffusionKineticEnergyRadialDistancePowerLaw(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 1 (kinetic energy)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Number of bins
   MultiIndex n_bins1(100, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(EnrKin(momentum1, specie), 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1(EnrKin(momentum2, specie), 0.0, 0.0);
   container.Insert(maxval1);

// Linear or logarithmic bins
   MultiIndex log_bins1(1, 0, 0);
   container.Insert(log_bins1);

// Add outlying events to the end bins
   MultiIndex bin_outside1(0, 0, 0);
   container.Insert(bin_outside1);

// Physical units of the distro variable
   double unit_distro1 = 1.0 / (Sqr(unit_length_fluid) * unit_time_fluid * M_4PI * unit_energy_particle);
   container.Insert(unit_distro1);

// Physical units of the bin variable
   GeoVector unit_val1 = {unit_energy_particle, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

// Normalization for the "hot" boundary
   double J0 = 1.0 / unit_distro1;
   container.Insert(J0);

// Characteristic energy
   double T0 = 1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0);

// Spectral power law
   double pow_law_T = -1.8;
   container.Insert(pow_law_T);

// Constant value for the "cold" condition
   double val_cold1 = 0.0;
   container.Insert(val_cold1);

   simulation->AddDistribution(DistributionSpectrumKineticEnergyPowerLaw(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 2 (exit time)
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins2(100, 0, 0);
   container.Insert(n_bins2);
   
// Smallest value
   GeoVector minval2(0.0, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(maxtime, 0.0, 0.0);
   container.Insert(maxval2);

// Linear or logarithmic bins
   MultiIndex log_bins2(0, 0, 0);
   container.Insert(log_bins2);

// Add outlying events to the end bins
   MultiIndex bin_outside2(1, 0, 0);
   container.Insert(bin_outside2);

// Physical units of the distro variable
   double unit_distro2 = 1.0;
   container.Insert(unit_distro2);

// Physical units of the bin variable
   GeoVector unit_val2 = {unit_time_fluid, 1.0, 1.0};
   container.Insert(unit_val2);

// Keep records
   bool keep_records2 = false;
   container.Insert(keep_records2);

// Value for the "hot" condition
   double val_hot2 = 1.0;
   container.Insert(val_hot2);

// Value for the "cold" condition
   double val_cold2 = 0.0;
   container.Insert(val_cold2);

// Value for which time to bin
   int val_time2 = 1;
   container.Insert(val_time2);

   simulation->AddDistribution(DistributionTimeUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------

    int n_traj;
    int batch_size;

    batch_size = n_traj = 1;
    if(argc > 1) {
        n_traj = atoi(argv[1]);
    }
    if(argc > 2) {
        batch_size = atoi(argv[2]);
    }

    simulation->SetTasks(n_traj, batch_size);

    simulation->MainLoop();

    std::string simulation_files_prefix = "single_particle_batl_test_distro_";
    simulation->DistroFileName(simulation_files_prefix);

    simulation->PrintDistro1D(0, 0, "single_particle_batl_test_distro_0.dat", true);

    return 0;
};
