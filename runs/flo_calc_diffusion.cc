#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace Spectrum;

inline GeoVector drift_numer(double r_L, double vel, SpatialData spdata, int specie)
{
   GeoVector drift = (r_L * vel / 3.0) * (spdata.curlB() - 2.0 * (spdata.gradBmag ^ spdata.bhat)) / spdata.Bmag;
   // Correct magnitude if necessary
   if (drift.Norm() > 0.5 * vel)
   {
      drift.Normalize();
      drift *= 0.5 * vel;
   };
   return drift;
};

int main(int argc, char **argv)
{
   int active_local_workers, workers_stopped;
   BackgroundServerBATL background;

   SpatialData spdata;
   double t = 0.0;
   int i, j, k;
   GeoVector pos = gv_zeros, vel = gv_zeros, mom = gv_zeros;
   int specie = SPECIES_PROTON_BEAM;

   std::ofstream data_file;

   //--------------------------------------------------------------------------------------------------
   // Server
   //--------------------------------------------------------------------------------------------------

   std::shared_ptr<MPI_Config> mpi_config = std::make_shared<MPI_Config>(argc, argv);
   std::shared_ptr<ServerBaseBack> server_back = nullptr;

   std::string fname_pattern = "3d__var_4_t30000330_n00456411";

   if (mpi_config->is_boss)
   {
      server_back = std::make_unique<ServerBackType>(fname_pattern);
      active_local_workers = mpi_config->workers_in_node;
      server_back->ServerStart();
   };

   DataContainer container;

   //--------------------------------------------------------------------------------------------------
   // Background
   //--------------------------------------------------------------------------------------------------

   container.Clear();

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

   background.SetupObject(container);

   //--------------------------------------------------------------------------------------------------
   if (mpi_config->is_boss)
   {
      while (active_local_workers)
      {
         workers_stopped = server_back->ServerFunctions();
         active_local_workers -= workers_stopped;
      };
      server_back->ServerFinish();
   }
   else if (mpi_config->is_worker)
   {
      int i, j, k, N = 401;
      double Rs = 6.957e+10 / unit_length_fluid;
      double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
      double x_min = -200.0 * one_au;
      double y_min = 0 * one_au;
      double z_min = -400.0 * one_au;
      double dx = 800.0 * one_au / (N - 1);
      double dy = 400.0 * one_au / (N - 1);
      double dz = 800.0 * one_au / (N - 1);
      double kappa_lzp_perp;
      double kappa_lzp_par;
      double kappa_nlgc_perp;
      double kappa_nlgc_par;
      double lambda_lzp_perp;
      double lambda_nlgc_perp;
      double lambda_lzp_par;
      double lambda_nlgc_par;
      double lam0 = 0.1 * one_au;
      double R0 = 3.33e6 / unit_rigidity_particle;
      double B0 = 5.0e-5 / unit_magnetic_fluid;
      int l_perp_index = 2;
      int del_b_index = 3;
      int in_HS = 0;
      double delb_lzp;
      double delb_nlgc;
      spdata._mask = BACKGROUND_ALL | BACKGROUND_gradB;

      pos[0] = one_au * 40.0;
      mom[0] = Mom(1000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
      vel[0] = Vel(mom[0], specie);

      //--------------------------------------------------------------------------------------------------
      std::cout << "3D plots..." << std::endl;
      data_file.open("flo_diffusion_data.txt");

      pos = gv_zeros;
      for (i = 0; i < N; i++)
      {
         for (j = 0; j < 1; j++)
         {
            for (k = 0; k < N; k++)
            {
               pos[0] = x_min + i * dx;
               pos[1] = y_min + j * dy;
               pos[2] = z_min + k * dz;

               double vmag = Vel(mom[0]);

               if (pos.Norm() < 1.0 * one_au)
               {
                  spdata.n_dens = 0.0;
                  spdata.Uvec = gv_zeros;
                  spdata.Bmag = 0.0;
                  spdata.dmax = 0.0;
                  kappa_lzp_par = 0.0;
                  kappa_nlgc_par = 0.0;
                  kappa_lzp_perp = 0.0;
                  kappa_nlgc_perp = 0.0;
               }
               else
               {
                    background.GetFields(t, pos, mom, spdata);
                    delb_lzp = 3e-12 * spdata.n_dens / (0.06); //*sin(M_PI/2.0 - acos(abs(pos[2])/pos.Norm())/1.67));
                    delb_nlgc = 5e-12 * spdata.n_dens / (0.06); //*sin(M_PI/2.0 - acos(abs(pos[2])/pos.Norm())/1.67));

                    kappa_lzp_par = (27.0/35.0) * pow(LarmorRadius(mom[0], spdata.Bmag, specie)*unit_length_fluid * 0.05 * 0.05 * unit_length_fluid * unit_length_fluid, 0.333) * vmag * unit_velocity_fluid * spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid * (10/delb_lzp)
                                * ((7.0/27.0) * pow(LarmorRadius(mom[0], spdata.Bmag, specie)/0.05, 1.67) + 1.0);
                    kappa_lzp_perp = 0.03 * delb_lzp * kappa_lzp_par / (spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid);

                    if (spdata.Uvec.Norm() > 3.0) {
                        kappa_nlgc_par = (27.0/35.0) * pow(LarmorRadius(mom[0], spdata.Bmag, specie)*unit_length_fluid * 0.1 * 0.1 * unit_length_fluid * unit_length_fluid, 0.333) * vmag * unit_velocity_fluid * spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid * (10/delb_nlgc)
                                * ((7.0/27.0) * pow(LarmorRadius(mom[0], spdata.Bmag, specie)/0.1, 1.67) + 1.0);
                        kappa_nlgc_perp = 1.1 * pow((delb_nlgc * 0.1 * unit_length_fluid * vmag * unit_velocity_fluid)/(spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid * 3.0), 0.667) * pow(kappa_nlgc_par, 0.333);
                    }
                    else{
                        kappa_nlgc_par = (27.0/35.0) * pow(LarmorRadius(mom[0], spdata.Bmag, specie)*unit_length_fluid * 0.033 * 0.033 * unit_length_fluid * unit_length_fluid, 0.333) * vmag * unit_velocity_fluid * spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid * (10/delb_nlgc)
                                * ((7.0/27.0) * pow(LarmorRadius(mom[0], spdata.Bmag, specie)/0.033, 1.67) + 1.0);
                        kappa_nlgc_perp = 1.1 * pow((delb_nlgc * 0.033 * unit_length_fluid * vmag * unit_velocity_fluid)/(spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid * 3.0), 0.667) * pow(kappa_nlgc_par, 0.333);
                    }
               };

                    

               data_file << std::setw(15) << pos[0] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT
                         << std::setw(15) << pos[1] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT
                         << std::setw(15) << pos[2] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT
                         << std::setw(15) << delb_lzp/5e-12
                         << std::setw(15) << 3*kappa_lzp_par/vmag/unit_diffusion_fluid
                         << std::setw(15) << 3*kappa_nlgc_par/vmag/unit_diffusion_fluid
                         << std::setw(15) << 3*kappa_lzp_perp/vmag/unit_diffusion_fluid
                         << std::setw(15) << 3*kappa_nlgc_perp/vmag/unit_diffusion_fluid
                         << std::endl;
            };
         };
      };
      data_file.close();

      background.StopServerFront();
      std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "Node " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};
