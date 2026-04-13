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
   double lat, lon;
   GeoVector pos = gv_zeros, vel = gv_zeros, mom = gv_zeros;
   int specie = SPECIES_PROTON_BEAM;
   std::ofstream CurrentSheet_file;

   //    if (argc > 1) {
   //       cir_date = argv[1];
   //       std::cout << "CIR date: " << cir_date << std::endl;
   //    } else {
   //       std::cout << "ERROR: No CIR date provided." << std::endl;
   //       exit(1);
   //    };

   //--------------------------------------------------------------------------------------------------
   // Server
   //--------------------------------------------------------------------------------------------------

   std::shared_ptr<MPI_Config> mpi_config = std::make_shared<MPI_Config>(argc, argv);
   std::shared_ptr<ServerBaseBack> server_back = nullptr;

   std::string fname_pattern = "/data001/cosmicrays_vf/Sammy/SPECTRUM/runs/3d__var_3_t30000330_n00456411";

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

      int lat_i, lon_i, time_i, N = 108;
      int M = 1 + (N/2);
      double Rs = 6.957e+10 / unit_length_fluid;
      double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
      double lat_max = M_PI/6;
      double lon_min = 0;
      double d_lat = (M_PI/3) / (M-1);
      double d_lon = 2 * M_PI / N;
      int max_time = N * 10;
      double dt = 86400*(27.0/N);
      std::cout << "dt is..." << dt << std::endl;
      spdata._mask = BACKGROUND_ALL | BACKGROUND_gradB;
      int lat_selector = 0;
      pos[0] = one_au * 40.0;
      mom[0] = Mom(1000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
      vel[0] = Vel(mom[0], specie);
      background.GetFields(t, pos, mom, spdata);

      //--------------------------------------------------------------------------------------------------
      std::cout << "3D plots..." << std::endl;
      CurrentSheet_file.open("current_sheet_points_subset.txt");

      pos = gv_zeros;
      for (lon_i = 0; lon_i < N; lon_i++)
      {
	 std::cout << "Longitude " << lon_i << " of " << N << "." << std::endl;
         lon = lon_min + lon_i * d_lon;
         for (lat_i = 0; lat_i < M; lat_i++)
         {
	    lat_selector = abs((lat_i)%(N) - (N/2));
            lat = lat_max * cos(lat_i * M_PI / num_lats);
            pos[0] = one_au * cos(lat) * cos(lon);
            pos[1] = one_au * cos(lat) * sin(lon);
            pos[2] = one_au * sin(lat);
            background.GetFields(t, pos, mom, spdata);
	    if (abs((lon_i)%(N) - (N/2)) == lat_selector)
		{	    
            		CurrentSheet_file << std::setprecision(6) << pos[0] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT << " "
                        << pos[1] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT << " "
                     	<< pos[2] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT //<< " "
//                     << spdata.Uvec[0] * unit_velocity_fluid << " "
//                     << spdata.Uvec[1] * unit_velocity_fluid << " "
//                     << spdata.Uvec[2] * unit_velocity_fluid << " "
//                     << spdata.Bvec[0] * unit_magnetic_fluid << " "
//                     << spdata.Bvec[1] * unit_magnetic_fluid << " "
//                     << spdata.Bvec[2] * unit_magnetic_fluid << " "
//                     << spdata.region[0] << " "
//                     << spdata.region[1]
                     << std::endl;
		};
            for (time_i = 0; time_i < max_time-1; time_i++)
            {  
               background.GetFields(t, pos, mom, spdata);

               pos[0] = pos[0] + (spdata.Uvec[0]*unit_velocity_fluid*dt/unit_length_fluid);
               pos[1] = pos[1] + (spdata.Uvec[1]*unit_velocity_fluid*dt/unit_length_fluid);
               pos[2] = pos[2] + (spdata.Uvec[2]*unit_velocity_fluid*dt/unit_length_fluid);
               if (abs((time_i+lon_i)%(N) - (N/2)) == lat_selector)
	       {
		       CurrentSheet_file << std::setprecision(6) << pos[0] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT << " "
                         << pos[1] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT << " "
                         << pos[2] * unit_length_fluid / GSL_CONST_CGSM_ASTRONOMICAL_UNIT // << " "
//                         << spdata.Uvec[0] * unit_velocity_fluid << " "
//                         << spdata.Uvec[1] * unit_velocity_fluid << " "
//                         << spdata.Uvec[2] * unit_velocity_fluid << " "
//                         << spdata.Bvec[0] * unit_magnetic_fluid << " "
//                         << spdata.Bvec[1] * unit_magnetic_fluid << " "
//                         << spdata.Bvec[2] * unit_magnetic_fluid << " "
//                         << spdata.region[0] << " "
//                         << spdata.region[1]
                         << std::endl;
	       };
            };
	 };
      };
      CurrentSheet_file.close();

      background.StopServerFront();
      std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "Node " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};
