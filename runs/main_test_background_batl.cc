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
      int i, j, k, N = 40;
      double Rs = 6.957e+10 / unit_length_fluid;
      double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
      double x_min = -200.0 * one_au;
      double y_min = -200.0 * one_au;
      double z_min = -200.0 * one_au;
      double dx = 400.0 * one_au / (N - 1);
      double dy = 400.0 * one_au / (N - 1);
      double dz = 400.0 * one_au / (N - 1);
      double polarity, r_L;
      GeoVector drift_vel;
      spdata._mask = BACKGROUND_ALL | BACKGROUND_gradB;

      pos[0] = one_au * 40.0;
      mom[0] = Mom(1000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
      vel[0] = Vel(mom[0], specie);
      background.GetFields(t, pos, mom, spdata);
      r_L = LarmorRadius(mom[0], spdata.Bmag, specie);
      std::cout << "|B| @ 1au = "
                << std::setw(18) << spdata.Bmag * unit_magnetic_fluid * 1.0e5 << " nT"
                << std::endl;
      std::cout << "rL (1 GeV) = "
                << std::setw(18) << r_L << " au"
                << std::endl;

      //--------------------------------------------------------------------------------------------------
      std::cout << "3D plots..." << std::endl;
      data_file.open("output/test_read_grid_.dat");

      pos = gv_zeros;
      for (i = 0; i < N; i++)
      {
         for (j = 0; j < N; j++)
         {
            for (k = 0; k < N; k++)
            {
               pos[0] = x_min + i * dx;
               pos[1] = y_min + j * dy;
               pos[2] = z_min + k * dz;

               if (pos.Norm() < 1.0 * one_au)
               {
                  spdata.n_dens = 0.0;
                  spdata.Uvec = gv_zeros;
                  spdata.Bmag = 0.0;
                  spdata.region = SimpleArray<double, 4>();
               }
               else
               {
                  background.GetFields(t, pos, mom, spdata);
               };

               data_file << std::setw(18) << pos[0] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT
                         << std::setw(18) << pos[1] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT
                         << std::setw(18) << pos[2] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT
                         << std::setw(18) << spdata.n_dens * unit_number_density_fluid
                         << std::setw(18) << spdata.Uvec[0] * unit_velocity_fluid
                         << std::setw(18) << spdata.Uvec[1] * unit_velocity_fluid
                         << std::setw(18) << spdata.Uvec[2] * unit_velocity_fluid
                         << std::setw(18) << spdata.Bmag * unit_magnetic_fluid
                         << std::setw(18) << spdata.region[0]
                         << std::setw(18) << spdata.region[1]
                         << std::setw(18) << spdata.region[2]
                         << std::setw(18) << spdata.region[3]
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
