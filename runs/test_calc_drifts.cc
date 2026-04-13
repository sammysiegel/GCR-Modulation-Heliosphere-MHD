#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

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

std::vector<int> i_to_streamline(int idx, int num_lats=108, int num_steps=500) {
    int lon_i = idx/(num_steps*num_lats);
    int lat_i = (idx%(num_steps*num_lats))/num_steps;
    int step_i = idx%num_steps;
    return {lon_i, lat_i, step_i};
};

int streamline_to_i(int lon_i, int lat_i, int step_i, int num_lons=108, int num_lats=108, int num_steps=500) {
    int N = num_lats*num_lons*num_steps;
    return (lon_i*num_lats*num_steps + lat_i*num_steps + step_i + N)%N;
};

int get_lon_next(int idx, int num_lons=108) {
    int extra_step = 0;
    int lon_i, lat_i, step_i;
    std::vector<int> streamline = i_to_streamline(idx);
    lon_i = streamline[0];
    lat_i = streamline[1];
    step_i = streamline[2];
    if ((lon_i+lat_i)%num_lons==0) {
        extra_step = -1;
    }
    return(streamline_to_i((lon_i+1)%num_lons, lat_i, step_i+extra_step));
};

int get_lon_prev(int idx, int num_lons=108) {
    int extra_step = 0;
    int lon_i, lat_i, step_i;
    std::vector<int> streamline = i_to_streamline(idx);
    lon_i = streamline[0];
    lat_i = streamline[1];
    step_i = streamline[2];
    if ((lon_i+lat_i-1)%num_lons==0) {
        extra_step = 1;
    }
    return(streamline_to_i((lon_i-1+num_lons)%num_lons, lat_i, step_i+extra_step));
};

int get_lat_next(int idx, int num_lons=108) {
    int extra_step = 0;
    int lon_i, lat_i, step_i;
    std::vector<int> streamline = i_to_streamline(idx);
    lon_i = streamline[0];
    lat_i = streamline[1];
    step_i = streamline[2];
    if ((lon_i+lat_i)%num_lons==0) {
        extra_step = -1;
    }
    return(streamline_to_i(lon_i, (lat_i+1)%num_lons, step_i+extra_step));
};

int get_lat_prev(int idx, int num_lons=108) {
    int extra_step = 0;
    int lon_i, lat_i, step_i;
    std::vector<int> streamline = i_to_streamline(idx);
    lon_i = streamline[0];
    lat_i = streamline[1];
    step_i = streamline[2];
    if ((lon_i+lat_i-1)%num_lons==0) {
        extra_step = 1;
    }
    return(streamline_to_i(lon_i, (lat_i-1+num_lons)%num_lons, step_i+extra_step));
};

double get_distance_idx(std::vector<double> p1, std::vector<double> p0) {
    double dx, dy, dz;
    dx = p0[0] - p1[0];
    dy = p0[1] - p1[1];
    dz = p0[2] - p1[2];
    return (dx*dx + dy*dy + dz*dz);
};

double get_cs_drift_speed(double d_rg) {
    return 0.457 - 0.412*d_rg + 0.0915*d_rg*d_rg;
};

int main(int argc, char **argv)
{
   int active_local_workers, workers_stopped;
   BackgroundServerBATL background;

   SpatialData spdata;
   double t = 0.0;
   GeoVector pos = gv_zeros, vel = gv_zeros, mom = gv_zeros;
   int specie = SPECIES_PROTON_BEAM;

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

        double Rs = 6.957e+10 / unit_length_fluid;
        double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
        spdata._mask = BACKGROUND_ALL | BACKGROUND_gradB;
        pos[0] = one_au * 40.0;
        mom[0] = Mom(1000.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
        vel[0] = Vel(mom[0], specie);
        background.GetFields(t, pos, mom, spdata);

        //--------------------------------------------------------------------------------------------------

        std::cout << "3D plots..." << std::endl;

        pos = gv_zeros;
      
        int num_lats = 108;
        int num_lons = 108;
        int num_steps = 500;
        double polarity = 1.0;
        int i;
        int x_i, z_i;
        int N_x = 200;
        int N_z = 400;
        double x_start = -100.0;
        double z_start = 0.0;
        double delta_x = 0.025;
        double delta_z = 0.025;
        int N = num_lats*num_lons*num_steps;
        std::cout << "N is " << N << std::endl;
        std::vector<double> X_vals, Y_vals, Z_vals;

        std::string infilename = "current_sheet_points_subset.txt";
        std::ifstream input_file(infilename);
        std::ofstream drift_results_file;
        double x, y, z;

        for (i=0; i<N; i++) {
            input_file >> x;
            input_file >> y;
            input_file >> z;
            X_vals.push_back(x);
            Y_vals.push_back(y);
            Z_vals.push_back(z);
        }

        std::vector<double> test_point = {52.8, 54.0, 20.0};
        double distance; 
        double distance_nearest;
        int idx_0_0;
        double dx, dy, dz;
        std::vector<double> p_a, p_b, p_c, p_ab, p_ac, p_a0, cross_product, final_drift_vector;
        double final_distance, direction_mag, r_L, drift_speed;
        int idx_1_0, idx_m1_0, idx_0_1, idx_0_m1, idx_1_1, idx_m1_1, idx_1_m1, idx_m1_m1;
        double distance_1_0, distance_m1_0, distance_0_1, distance_0_m1, distance_1_1, distance_m1_1, distance_1_m1, distance_m1_m1;
	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	double side_of_sheet;
        drift_results_file.open("drift_data.txt");
	
        for (x_i=0; x_i<N_x; x_i++) {
            for (z_i=0; z_i<N_z; z_i++) {
		test_point = {x_start+x_i*delta_x, 0.0, z_start+z_i*delta_z};
                if (test_point[0]*test_point[0]+test_point[2]*test_point[2] < 1){
		    continue;
		}
		pos[0] = test_point[0]*one_au;
                pos[1] = test_point[1]*one_au;
                pos[2] = test_point[2]*one_au;
                background.GetFields(t, pos, mom, spdata);
		if (spdata.region[0] < -0.1) {
                    continue;
                };
                r_L = LarmorRadius(mom[0], spdata.Bmag, specie)/one_au;
                distance_nearest = 1000000.0;
                idx_0_0 = -1;
		if (z_i==0) {
                    start = std::chrono::high_resolution_clock::now();
		};
                for (i=0; i<N; i++) {
                    dx = test_point[0]-X_vals[i];
                    dy = test_point[1]-Y_vals[i];
                    dz = test_point[2]-Z_vals[i];
                    distance = dx*dx + dy*dy + dz*dz;
                    if(distance<distance_nearest) {
                        distance_nearest = distance;
                        idx_0_0 = i;
                    };
                };

		if(z_i==0) {
                    stop = std::chrono::high_resolution_clock::now();
                    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
                    std::cout << duration.count() << " microseconds" << std::endl;
                };

		p_a = {X_vals[idx_0_0], Y_vals[idx_0_0], Z_vals[idx_0_0]};
                if(idx_0_0%num_steps==0) {
                    final_distance = sqrt(distance_nearest);
		    side_of_sheet = (test_point[2]-p_a[2])/abs(test_point[2]-p_a[2]);
                    if (final_distance/r_L < 2.0) {
                        drift_speed = get_cs_drift_speed(final_distance/r_L);
                        idx_0_m1 = get_lat_prev(idx_0_0);
                        final_drift_vector = {X_vals[idx_0_0]-X_vals[idx_0_m1], Y_vals[idx_0_0]-Y_vals[idx_0_m1], Z_vals[idx_0_0]-Z_vals[idx_0_m1]};
                        direction_mag = sqrt(final_drift_vector[0]*final_drift_vector[0]+final_drift_vector[1]*final_drift_vector[1]+final_drift_vector[2]*final_drift_vector[2]);
                        final_drift_vector[0] = polarity*final_drift_vector[0]*drift_speed/direction_mag;
                        final_drift_vector[1] = polarity*final_drift_vector[1]*drift_speed/direction_mag;
                        final_drift_vector[2] = polarity*final_drift_vector[2]*drift_speed/direction_mag;
                    }
                    else {
			drift_speed = 0.0;
                        final_drift_vector = {0.0, 0.0, 0.0};
                    };
                }
                else {
                    idx_1_0 = get_lon_next(idx_0_0);
                    distance_1_0 = get_distance_idx({X_vals[idx_1_0], Y_vals[idx_1_0], Z_vals[idx_1_0]}, test_point);
                    idx_m1_0 = get_lon_prev(idx_0_0);
                    distance_m1_0 = get_distance_idx({X_vals[idx_m1_0], Y_vals[idx_m1_0], Z_vals[idx_m1_0]}, test_point);
                    idx_0_1 = get_lat_next(idx_0_0);
                    distance_0_1 = get_distance_idx({X_vals[idx_0_1], Y_vals[idx_0_1], Z_vals[idx_0_1]}, test_point);
                    idx_0_m1 = get_lat_prev(idx_0_0);
                    distance_0_m1 = get_distance_idx({X_vals[idx_0_m1], Y_vals[idx_0_m1], Z_vals[idx_0_m1]}, test_point);
                    idx_1_1 = get_lat_next(idx_1_0);
                    distance_1_1 = get_distance_idx({X_vals[idx_1_1], Y_vals[idx_1_1], Z_vals[idx_1_1]}, test_point);
                    idx_m1_1 = get_lat_next(idx_m1_0);
                    distance_m1_1 = get_distance_idx({X_vals[idx_m1_1], Y_vals[idx_m1_1], Z_vals[idx_m1_1]}, test_point);
                    idx_1_m1 = get_lat_prev(idx_1_0);
                    distance_1_m1 = get_distance_idx({X_vals[idx_1_m1], Y_vals[idx_1_m1], Z_vals[idx_1_m1]}, test_point);
                    idx_m1_m1 = get_lat_prev(idx_m1_m1);
                    distance_m1_m1 = get_distance_idx({X_vals[idx_m1_m1], Y_vals[idx_m1_m1], Z_vals[idx_m1_m1]}, test_point);

                    if ((distance_1_0 < distance_m1_0) and (distance_1_0 < distance_0_1) and (distance_1_0 < distance_0_m1)) {
                        if ((distance_1_1 < distance_1_m1) and (distance_1_1 < distance_0_1) and (distance_1_1 < distance_0_m1)) {
                            p_b = {X_vals[idx_1_0], Y_vals[idx_1_0], Z_vals[idx_1_0]};
                            p_c = {X_vals[idx_1_1], Y_vals[idx_1_1], Z_vals[idx_1_1]};
                        }
                        else if ((distance_1_m1 < distance_1_1) and (distance_1_m1 < distance_0_1) and (distance_1_m1 < distance_0_m1)) {
                            p_b = {X_vals[idx_1_m1], Y_vals[idx_1_m1], Z_vals[idx_1_m1]};
                            p_c = {X_vals[idx_1_0], Y_vals[idx_1_0], Z_vals[idx_1_0]};
                        }
                        else if ((distance_0_1 < distance_1_1) and (distance_0_1 < distance_1_m1) and (distance_0_1 < distance_0_m1)) {
                            p_b = {X_vals[idx_1_0], Y_vals[idx_1_0], Z_vals[idx_1_0]};
                            p_c = {X_vals[idx_0_1], Y_vals[idx_0_1], Z_vals[idx_0_1]};
                        }
                        else if ((distance_0_m1 < distance_1_1) and (distance_0_m1 < distance_1_m1) and (distance_0_m1 < distance_0_1)) {
                            p_b = {X_vals[idx_0_m1], Y_vals[idx_0_m1], Z_vals[idx_0_m1]};
                            p_c = {X_vals[idx_1_0], Y_vals[idx_1_0], Z_vals[idx_1_0]};
                        };
                    }
                    else if ((distance_m1_0 < distance_1_0) and (distance_m1_0 < distance_0_1) and (distance_m1_0 < distance_0_m1)) {
                        if ((distance_0_1 < distance_0_m1) and (distance_0_1 < distance_m1_1) and (distance_0_1 < distance_m1_m1)) {
                            p_b = {X_vals[idx_0_1], Y_vals[idx_0_1], Z_vals[idx_0_1]};
                            p_c = {X_vals[idx_m1_0], Y_vals[idx_m1_0], Z_vals[idx_m1_0]};
                        }
                        else if ((distance_0_m1 < distance_0_1) and (distance_0_m1 < distance_m1_1) and (distance_0_m1 < distance_m1_m1)) {
                            p_b = {X_vals[idx_m1_0], Y_vals[idx_m1_0], Z_vals[idx_m1_0]};
                            p_c = {X_vals[idx_0_m1], Y_vals[idx_0_m1], Z_vals[idx_0_m1]};
                        }
                        else if ((distance_m1_1 < distance_0_1) and (distance_m1_1 < distance_0_m1) and (distance_m1_1 < distance_m1_m1)) {
                            p_b = {X_vals[idx_m1_1], Y_vals[idx_m1_1], Z_vals[idx_m1_1]};
                            p_c = {X_vals[idx_m1_0], Y_vals[idx_m1_0], Z_vals[idx_m1_0]};
                        }
                        else if ((distance_m1_m1 < distance_0_1) and (distance_m1_m1 < distance_0_m1) and (distance_m1_m1 < distance_m1_1)) {
                            p_b = {X_vals[idx_m1_0], Y_vals[idx_m1_0], Z_vals[idx_m1_0]};
                            p_c = {X_vals[idx_m1_m1], Y_vals[idx_m1_m1], Z_vals[idx_m1_m1]};
                        };
                    }
                    else if ((distance_0_1 < distance_1_0) and (distance_0_1 < distance_m1_0) and (distance_0_1 < distance_0_m1)) {
                        if ((distance_m1_1 < distance_1_1) and (distance_m1_1 < distance_m1_0) and (distance_m1_1 < distance_1_0)) {
                            p_b = {X_vals[idx_0_1], Y_vals[idx_0_1], Z_vals[idx_0_1]};
                            p_c = {X_vals[idx_m1_1], Y_vals[idx_m1_1], Z_vals[idx_m1_1]};
                        }
                        else if ((distance_1_1 < distance_m1_1) and (distance_1_1 < distance_m1_0) and (distance_1_1 < distance_1_0)) {
                            p_b = {X_vals[idx_1_1], Y_vals[idx_1_1], Z_vals[idx_1_1]};
                            p_c = {X_vals[idx_0_1], Y_vals[idx_0_1], Z_vals[idx_0_1]};
                        }
                        else if ((distance_m1_0 < distance_m1_1) and (distance_m1_0 < distance_1_1) and (distance_m1_0 < distance_1_0)) {
                            p_b = {X_vals[idx_0_1], Y_vals[idx_0_1], Z_vals[idx_0_1]};
                            p_c = {X_vals[idx_m1_1], Y_vals[idx_m1_1], Z_vals[idx_m1_1]};
                        }
                        else if ((distance_1_0 < distance_m1_1) and (distance_1_0 < distance_1_1) and (distance_1_0 < distance_m1_0)) {
                            p_b = {X_vals[idx_1_0], Y_vals[idx_1_0], Z_vals[idx_1_0]};
                            p_c = {X_vals[idx_0_1], Y_vals[idx_0_1], Z_vals[idx_0_1]};
                        };
                    }
                    else if ((distance_0_m1 < distance_1_0) and (distance_0_m1 < distance_m1_0) and (distance_0_m1 < distance_0_1)) {
                        if ((distance_m1_0 < distance_1_0) and (distance_m1_0 < distance_m1_m1) and (distance_m1_0 < distance_1_m1)) {
                            p_b = {X_vals[idx_m1_0], Y_vals[idx_m1_0], Z_vals[idx_m1_0]};
                            p_c = {X_vals[idx_0_m1], Y_vals[idx_0_m1], Z_vals[idx_0_m1]};
                        }
                        else if ((distance_1_0 < distance_m1_0) and (distance_1_0 < distance_m1_m1) and (distance_1_0 < distance_1_m1)) {
                            p_b = {X_vals[idx_0_m1], Y_vals[idx_0_m1], Z_vals[idx_0_m1]};
                            p_c = {X_vals[idx_1_0], Y_vals[idx_1_0], Z_vals[idx_1_0]};
                        }
                        else if ((distance_m1_m1 < distance_m1_0) and (distance_m1_m1 < distance_1_0) and (distance_m1_m1 < distance_1_m1)) {
                            p_b = {X_vals[idx_m1_m1], Y_vals[idx_m1_m1], Z_vals[idx_m1_m1]};
                            p_c = {X_vals[idx_0_m1], Y_vals[idx_0_m1], Z_vals[idx_0_m1]};
                        }
                        else if ((distance_1_m1 < distance_m1_0) and (distance_1_m1 < distance_1_0) and (distance_1_m1 < distance_m1_m1)) {
                            p_b = {X_vals[idx_0_m1], Y_vals[idx_0_m1], Z_vals[idx_0_m1]};
                            p_c = {X_vals[idx_1_m1], Y_vals[idx_1_m1], Z_vals[idx_1_m1]};
                        };
                    };

                    p_ab = {p_b[0]-p_a[0], p_b[1]-p_a[1], p_b[2]-p_a[2]};
                    p_ac = {p_c[0]-p_a[0], p_c[1]-p_a[1], p_c[2]-p_a[2]};
                    p_a0 = {test_point[0]-p_a[0], test_point[1]-p_a[1], test_point[2]-p_a[2]};

                    cross_product = {
                        p_ab[1]*p_ac[2] - p_ac[1]*p_ab[2],
                        p_ac[0]*p_ab[2] - p_ab[0]*p_ac[2],
                        p_ab[0]*p_ac[1] - p_ac[0]*p_ab[1]
                    };

                    final_distance = sqrt(cross_product[0]*cross_product[0] + cross_product[1]*cross_product[1] + cross_product[2]*cross_product[2]);
                    cross_product = {cross_product[0]/final_distance, cross_product[1]/final_distance, cross_product[2]/final_distance};
                    final_distance = (p_a0[0]*cross_product[0] + p_a0[1]*cross_product[1] + p_a0[2]*cross_product[2]);
		    side_of_sheet = final_distance/abs(final_distance);
		    final_distance = final_distance/side_of_sheet;
                    
		    if (final_distance/r_L < 2.0) {
                        drift_speed = get_cs_drift_speed(final_distance/r_L);
                        final_drift_vector = {X_vals[idx_0_0]-X_vals[idx_0_m1], Y_vals[idx_0_0]-Y_vals[idx_0_m1], Z_vals[idx_0_0]-Z_vals[idx_0_m1]};
                        direction_mag = sqrt(final_drift_vector[0]*final_drift_vector[0]+final_drift_vector[1]*final_drift_vector[1]+final_drift_vector[2]*final_drift_vector[2]);
                        final_drift_vector[0] = polarity*final_drift_vector[0]*drift_speed/direction_mag;
                        final_drift_vector[1] = polarity*final_drift_vector[1]*drift_speed/direction_mag;
                        final_drift_vector[2] = polarity*final_drift_vector[2]*drift_speed/direction_mag;
                    }
                    else {
			drift_speed = 0.0;
                        final_drift_vector = {0.0, 0.0, 0.0};
                    };
                };

                drift_results_file << std::setprecision(6) 
                << test_point[0] << ' '
                << test_point[1] << ' '
                << test_point[2] << ' '
                << final_drift_vector[0] << ' '
                << final_drift_vector[1] << ' '
                << final_drift_vector[2] << ' '
		<< drift_speed << ' '
		<< final_distance << ' '
		<< r_L << ' '
		<< side_of_sheet << ' '
		<< idx_0_0 << std::endl;
            };
        };
        drift_results_file.close();
        background.StopServerFront();
        std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "Node " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};
