#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include "src/diffusion_other.hh"
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
   DiffusionSammy diffusion_model;

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

   std::string fname_pattern = "3d__var_3_t30000330_n00456411";

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

   container.Clear();

   int l_perp_index = 2;
   container.Insert(l_perp_index);

   int del_b_index = 3;
   container.Insert(del_b_index);

   double kap_rat = 0.05;
   container.Insert(kap_rat);

   double R_max = Rigidity(Mom(5000 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT), specie);
   container.Insert(R_max);

   double A = 1.0;
   container.Insert(A);

   double mu = 0.0;
   container.Insert(mu);

   diffusion_model.SetSpecie(specie);
   diffusion_model.SetupObject(container);

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
      int i, j, k, N = 101;
      double Rs = 6.957e+10 / unit_length_fluid;
      double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
      double x_min = -200.0 * one_au;
      double y_min = 0 * one_au;
      double z_min = -400.0 * one_au;
      double dx = 800.0 * one_au / (N - 1);
      double dy = 400.0 * one_au / (N - 1);
      double dz = 800.0 * one_au / (N - 1);
      double kappa_para, kappa_perp;
      double kappa_old;
      double lam0 = 0.1 * one_au;
      double R0 = 3.33e6 / unit_rigidity_particle;
      double B0 = 5.0e-5 / unit_magnetic_fluid;
      int l_perp_index = 2;
      int del_b_index = 3;
      int in_HS = 0;
      double l_perp;
      double del_b;
      double e_idx;
      double d_eidx = 3.0 / (N-1);
      spdata._mask = BACKGROUND_ALL | BACKGROUND_gradB;
      double Omega;

      pos[0] = one_au * 40.0;
      mom[0] = Mom(1000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
      vel[0] = Vel(mom[0], specie);

      //--------------------------------------------------------------------------------------------------
      std::cout << "3D plots..." << std::endl;
      data_file.open("diffusion_data_rigidity3au.txt");

      pos = gv_zeros;
      for (i = 0; i < N; i++)
      {
               e_idx = 1.0 + i*d_eidx;
               mom[0] = Mom(pow(10.0, e_idx) * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
               vel[0] = Vel(mom[0], specie);
               pos[0] = 3.0;
               pos[1] = 0.0;
               pos[2] = 0.0;

               if (pos.Norm() < 1.0 * one_au)
               {
                  spdata.n_dens = 0.0;
                  spdata.Uvec = gv_zeros;
                  spdata.Bmag = 0.0;
                  spdata.dmax = 0.0;
                  kappa_para = 0.0;
                  kappa_perp = 0.0;
                  kappa_old = 0.0;
               }
               else
               {
                    background.GetFields(t, pos, mom, spdata);
                    kappa_para = diffusion_model.GetComponent(1, t, pos, mom, spdata);
                    kappa_perp = diffusion_model.GetComponent(0, t, pos, mom, spdata);

                    /* l_perp = spdata.region[l_perp_index] * 0.128;
                    del_b = (sqrt(spdata.region[del_b_index])/21.81) * 1e-5 / unit_magnetic_fluid;
                    Omega = vel[0] / LarmorRadius(mom[0], spdata.Bmag, specie);

                    if (spdata.Uvec.Norm() > 3.0) {
                        in_HS = 0;
                        kappa_para = ((3.0 * vel[0] * Sqr(Rigidity(mom[0], specie)))/(20.0 * l_perp * Sqr(del_b) * sin(3.0 * M_1_PI/5.0)))
                                    * (1.0 + (72.0/7.0)*pow(l_perp / LarmorRadius(mom[0], spdata.Bmag, specie) , 1.667));
                        // kappa_para = ((3.0 * Cube(vel[0]) * Sqr(spdata.Bmag))/(20.0 * l_perp * Sqr(Omega) * Sqr(del_b) * sin(3.0 * M_1_PI/5.0)))
                        //            * (1.0 + (72.0/7.0)*pow(l_perp / LarmorRadius(mom[0], spdata.Bmag, specie) , 1.667));


                        std::cout << "Vel  " << vel[0] << std::endl;
                        std::cout << "Bmag  " << spdata.Bmag << std::endl;
                        std::cout << "del_b  " << del_b << std::endl;
                        std::cout << "Larmor  " << LarmorRadius(mom[0], spdata.Bmag, specie) << std::endl;
                        std::cout << "Part 1" << ((3.0 * Cube(vel[0]) * Sqr(spdata.Bmag))/(20.0 * l_perp * Sqr(Omega) * Sqr(del_b) * sin(3.0 * M_1_PI/5.0))) << std::endl;
                        std::cout << "Part 2" << (1.0 + (72.0/7.0)*pow(l_perp / LarmorRadius(mom[0], spdata.Bmag, specie) , 1.667)) << std::endl;
                        kappa_perp = 1.1 * pow((del_b * del_b * 1e-5 * 1e-5 * l_perp * unit_length_fluid * vel[0] * unit_velocity_fluid)/(3.0 * spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid), 0.667) * pow(kappa_para, 0.333);
                    }
                    else {
                        in_HS = 1;
                        kappa_para = (
                           (1 / 9.7705) * pow(GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_MASS_PROTON / SPC_CONST_CGSM_ELECTRON_CHARGE, 0.333)
                           * pow(RelFactor1(mom[0]), 0.333) * vel[0] * vel[0] * unit_velocity_fluid * unit_velocity_fluid
                           * pow(l_perp * unit_length_fluid, 0.667)
                           * 0.144
                           / (0.03 * unit_magnetic_server)
                        )/3;
                        kappa_perp = 1.1 * pow((0.03 * 0.03 * 1e-5 * 1e-5 * l_perp * unit_length_fluid * vel[0] * unit_velocity_fluid)/(3.0 * spdata.Bmag * spdata.Bmag * unit_magnetic_fluid * unit_magnetic_fluid), 0.667) * pow(kappa_para, 0.333);

                    };
                    kappa_old = (lam0 * vel[0] * unit_length_fluid * unit_velocity_fluid / 3.0) * pow(Rigidity(mom[0], specie) / R0, 0.5) * pow(spdata.Bmag / B0, -1.0); */
               };

               data_file << std::setw(15) << Rigidity(mom[0], specie)
                         << std::setw(15) << pow(10.0, e_idx)
//                         << std::setw(15) << spdata.Bmag
//                         << std::setw(15) << spdata.region[0]
//                         << std::setw(15) << spdata.region[1]
//                         << std::setw(15) << l_perp
//                         << std::setw(15) << del_b
                         << std::setw(15) << kappa_para * unit_diffusion_fluid
                         << std::setw(15) << 3.0 * kappa_para / vel[0]
                         << std::setw(15) << 3.0 * kappa_perp / vel[0]
                         << std::endl;

      };
      data_file.close();

      background.StopServerFront();
      std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "Node " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;

};