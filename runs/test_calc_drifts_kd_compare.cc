#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <ANN/ANN.h>

using namespace Spectrum;

int				k				= 1;			// number of nearest neighbors
int				dim				= 3;			// dimension
double			eps				= 0;			// error bound
int				N			= 1000;			// maximum number of data points

std::istream*		dataIn			= NULL;			// input for data points

bool readPt(std::istream &in, ANNpoint p)			// read point (false on EOF)
{
	for (int i = 0; i < dim; i++) {
		if(!(in >> p[i])) return false;
	}
	return true;
}

void printPt(std::ostream &out, ANNpoint p)			// print point
{
	out << "(" << p[0];
	for (int i = 1; i < dim; i++) {
		out << ", " << p[i];
	}
	out << ")\n";
}

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

   std::string fname_pattern = "/data001/cosmicrays_vf/Sammy/SPECTRUM_devel/runs/3d__var_3_t30000330_n00456411";

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
        int N_x = 2000;
        int N_z = 4000;
        double x_start = -100.0;
        double z_start = 0.0;
        double delta_x = 0.025;
        double delta_z = 0.025;
        int N = num_lats*num_lons*num_steps;
        std::cout << "N is " << N << std::endl;

        std::string infilename = "current_sheet_points_subset.txt";
        std::ifstream input_file;
        input_file.open(infilename, std::ios::in);
        dataIn = &input_file;
        std::ofstream drift_results_file;

        int					nPts;					// actual number of data points
        ANNpointArray		dataPts;				// data points
        ANNpoint			queryPt;				// query point
        ANNidxArray			nnIdx;					// near neighbor indices
        ANNdistArray		dists;					// near neighbor distances
        ANNkd_tree*			kdTree;					// search structure

        // getArgs(argc, argv);						// read command-line arguments

        queryPt = annAllocPt(dim);					// allocate query point
        dataPts = annAllocPts(N, dim);			// allocate data points
        nnIdx = new ANNidx[k];						// allocate near neigh indices
        dists = new ANNdist[k];						// allocate near neighbor dists

        nPts = 0;									// read data points

        std::cout << "Reading Data Points...\n";
        while (nPts < N && readPt(*dataIn, dataPts[nPts])) {
            nPts++;
        }

        kdTree = new ANNkd_tree(					// build search structure
                    dataPts,					// the data points
                    nPts,						// number of points
                    dim);						// dimension of space

        std::vector<double> test_point;
        double distance_nearest;
        int idx_0_0;
        std::vector<double> p_a, p_b, p_c, p_ab, p_ac, p_a0, cross_product, final_drift_vector;
        double final_distance, direction_mag, r_L, drift_speed;
        int idx_1_0, idx_m1_0, idx_0_1, idx_0_m1, idx_1_1, idx_m1_1, idx_1_m1, idx_m1_m1;
        double distance_1_0, distance_m1_0, distance_0_1, distance_0_m1, distance_1_1, distance_m1_1, distance_1_m1, distance_m1_m1;
        double side_of_sheet;
        drift_results_file.open("drift_data.txt");

        int num_query_points;

        std::vector<double> X_vals, Y_vals, Z_vals;
        std::string queryfilename = "query_points.txt";
        std::ifstream query_file(queryfilename);
        double x, y, z, dummy1, dummy2, dummy3;

        for (i=0; i<578210; i++) {
            query_file >> x;
            query_file >> y;
            query_file >> z;
	    query_file >> dummy1;
	    query_file >> dummy2;
	    query_file >> dummy3;
            X_vals.push_back(x);
            Y_vals.push_back(y);
            Z_vals.push_back(z);
        }

        num_query_points = X_vals.size();

        for (x_i=0; x_i<num_query_points; x_i++) {
                test_point = {X_vals[x_i], Y_vals[x_i], Z_vals[x_i]};
                pos[0] = test_point[0]*one_au;
                pos[1] = test_point[1]*one_au;
                pos[2] = test_point[2]*one_au;
                queryPt[0] = test_point[0];
                queryPt[1] = test_point[1];
                queryPt[2] = test_point[2];
                background.GetFields(t, pos, mom, spdata);
                r_L = LarmorRadius(mom[0], spdata.Bmag, specie)/one_au;
                kdTree->annkSearch(						// search
                        queryPt,						// query point
                        k,								// number of near neighbors
                        nnIdx,							// nearest neighbors (returned)
                        dists,							// distance (returned)
                        eps);							// error bound
                idx_0_0 = nnIdx[0];
                distance_nearest = sqrt(dists[0]);
                p_a = {dataPts[idx_0_0][0], dataPts[idx_0_0][1], dataPts[idx_0_0][2]};

                if(idx_0_0%num_steps==0) {
                    final_distance = sqrt(distance_nearest);
                    side_of_sheet = (test_point[2]-p_a[2])/abs(test_point[2]-p_a[2]);
                    if (final_distance/r_L < 2.0) {
                        drift_speed = get_cs_drift_speed(final_distance/r_L);
                        idx_0_m1 = get_lat_prev(idx_0_0);
                        final_drift_vector = {dataPts[idx_0_0][0]-dataPts[idx_0_m1][0], dataPts[idx_0_0][1]-dataPts[idx_0_m1][1], dataPts[idx_0_0][2]-dataPts[idx_0_m1][2]};
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
                    distance_1_0 = get_distance_idx({dataPts[idx_1_0][0], dataPts[idx_1_0][1], dataPts[idx_1_0][2]}, test_point);
                    idx_m1_0 = get_lon_prev(idx_0_0);
                    distance_m1_0 = get_distance_idx({dataPts[idx_m1_0][0], dataPts[idx_m1_0][1], dataPts[idx_m1_0][2]}, test_point);
                    idx_0_1 = get_lat_next(idx_0_0);
                    distance_0_1 = get_distance_idx({dataPts[idx_0_1][0], dataPts[idx_0_1][1], dataPts[idx_0_1][2]}, test_point);
                    idx_0_m1 = get_lat_prev(idx_0_0);
                    distance_0_m1 = get_distance_idx({dataPts[idx_0_m1][0], dataPts[idx_0_m1][1], dataPts[idx_0_m1][2]}, test_point);
                    idx_1_1 = get_lat_next(idx_1_0);
                    distance_1_1 = get_distance_idx({dataPts[idx_1_1][0], dataPts[idx_1_1][1], dataPts[idx_1_1][2]}, test_point);
                    idx_m1_1 = get_lat_next(idx_m1_0);
                    distance_m1_1 = get_distance_idx({dataPts[idx_m1_1][0], dataPts[idx_m1_1][1], dataPts[idx_m1_1][2]}, test_point);
                    idx_1_m1 = get_lat_prev(idx_1_0);
                    distance_1_m1 = get_distance_idx({dataPts[idx_1_m1][0], dataPts[idx_1_m1][1], dataPts[idx_1_m1][2]}, test_point);
                    idx_m1_m1 = get_lat_prev(idx_m1_m1);
                    distance_m1_m1 = get_distance_idx({dataPts[idx_m1_m1][0], dataPts[idx_m1_m1][1], dataPts[idx_m1_m1][2]}, test_point);

                    if ((distance_1_0 < distance_m1_0) and (distance_1_0 < distance_0_1) and (distance_1_0 < distance_0_m1)) {
                        if ((distance_1_1 < distance_1_m1) and (distance_1_1 < distance_0_1) and (distance_1_1 < distance_0_m1)) {
                            p_b = {dataPts[idx_1_0][0], dataPts[idx_1_0][1], dataPts[idx_1_0][2]};
                            p_c = {dataPts[idx_1_1][0], dataPts[idx_1_1][1], dataPts[idx_1_1][2]};
                        }
                        else if ((distance_1_m1 < distance_1_1) and (distance_1_m1 < distance_0_1) and (distance_1_m1 < distance_0_m1)) {
                            p_b = {dataPts[idx_1_m1][0], dataPts[idx_1_m1][1], dataPts[idx_1_m1][2]};
                            p_c = {dataPts[idx_1_0][0], dataPts[idx_1_0][1], dataPts[idx_1_0][2]};
                        }
                        else if ((distance_0_1 < distance_1_1) and (distance_0_1 < distance_1_m1) and (distance_0_1 < distance_0_m1)) {
                            p_b = {dataPts[idx_1_0][0], dataPts[idx_1_0][1], dataPts[idx_1_0][2]};
                            p_c = {dataPts[idx_0_1][0], dataPts[idx_0_1][1], dataPts[idx_0_1][2]};
                        }
                        else if ((distance_0_m1 < distance_1_1) and (distance_0_m1 < distance_1_m1) and (distance_0_m1 < distance_0_1)) {
                            p_b = {dataPts[idx_0_m1][0], dataPts[idx_0_m1][1], dataPts[idx_0_m1][2]};
                            p_c = {dataPts[idx_1_0][0], dataPts[idx_1_0][1], dataPts[idx_1_0][2]};
                        }
                    }
                    else if ((distance_m1_0 < distance_1_0) and (distance_m1_0 < distance_0_1) and (distance_m1_0 < distance_0_m1)) {
                        if ((distance_0_1 < distance_0_m1) and (distance_0_1 < distance_m1_1) and (distance_0_1 < distance_m1_m1)) {
                            p_b = {dataPts[idx_0_1][0], dataPts[idx_0_1][1], dataPts[idx_0_1][2]};
                            p_c = {dataPts[idx_m1_0][0], dataPts[idx_m1_0][1], dataPts[idx_m1_0][2]};
                        }
                        else if ((distance_0_m1 < distance_0_1) and (distance_0_m1 < distance_m1_1) and (distance_0_m1 < distance_m1_m1)) {
                            p_b = {dataPts[idx_m1_0][0], dataPts[idx_m1_0][1], dataPts[idx_m1_0][2]};
                            p_c = {dataPts[idx_0_m1][0], dataPts[idx_0_m1][1], dataPts[idx_0_m1][2]};
                        }
                        else if ((distance_m1_1 < distance_0_1) and (distance_m1_1 < distance_0_m1) and (distance_m1_1 < distance_m1_m1)) {
                            p_b = {dataPts[idx_m1_1][0], dataPts[idx_m1_1][1], dataPts[idx_m1_1][2]};
                            p_c = {dataPts[idx_m1_0][0], dataPts[idx_m1_0][1], dataPts[idx_m1_0][2]};
                        }
                        else if ((distance_m1_m1 < distance_0_1) and (distance_m1_m1 < distance_0_m1) and (distance_m1_m1 < distance_m1_1)) {
                            p_b = {dataPts[idx_m1_0][0], dataPts[idx_m1_0][1], dataPts[idx_m1_0][2]};
                            p_c = {dataPts[idx_m1_m1][0], dataPts[idx_m1_m1][1], dataPts[idx_m1_m1][2]};
                        }
                    }
                    else if ((distance_0_1 < distance_1_0) and (distance_0_1 < distance_m1_0) and (distance_0_1 < distance_0_m1)) {
                        if ((distance_m1_1 < distance_1_1) and (distance_m1_1 < distance_m1_0) and (distance_m1_1 < distance_1_0)) {
                            p_b = {dataPts[idx_0_1][0], dataPts[idx_0_1][1], dataPts[idx_0_1][2]};
                            p_c = {dataPts[idx_m1_1][0], dataPts[idx_m1_1][1], dataPts[idx_m1_1][2]};
                        }
                        else if ((distance_1_1 < distance_m1_1) and (distance_1_1 < distance_m1_0) and (distance_1_1 < distance_1_0)) {
                            p_b = {dataPts[idx_1_1][0], dataPts[idx_1_1][1], dataPts[idx_1_1][2]};
                            p_c = {dataPts[idx_0_1][0], dataPts[idx_0_1][1], dataPts[idx_0_1][2]};
                        }
                        else if ((distance_m1_0 < distance_m1_1) and (distance_m1_0 < distance_1_1) and (distance_m1_0 < distance_1_0)) {
                            p_b = {dataPts[idx_0_1][0], dataPts[idx_0_1][1], dataPts[idx_0_1][2]};
                            p_c = {dataPts[idx_m1_1][0], dataPts[idx_m1_1][1], dataPts[idx_m1_1][2]};
                        }
                        else if ((distance_1_0 < distance_m1_1) and (distance_1_0 < distance_1_1) and (distance_1_0 < distance_m1_0)) {
                            p_b = {dataPts[idx_1_0][0], dataPts[idx_1_0][1], dataPts[idx_1_0][2]};
                            p_c = {dataPts[idx_0_1][0], dataPts[idx_0_1][1], dataPts[idx_0_1][2]};
                        }
                    }
                    else if ((distance_0_m1 < distance_1_0) and (distance_0_m1 < distance_m1_0) and (distance_0_m1 < distance_0_1)) {
                        if ((distance_m1_0 < distance_1_0) and (distance_m1_0 < distance_m1_m1) and (distance_m1_0 < distance_1_m1)) {
                            p_b = {dataPts[idx_m1_0][0], dataPts[idx_m1_0][1], dataPts[idx_m1_0][2]};
                            p_c = {dataPts[idx_0_m1][0], dataPts[idx_0_m1][1], dataPts[idx_0_m1][2]};
                        }
                        else if ((distance_1_0 < distance_m1_0) and (distance_1_0 < distance_m1_m1) and (distance_1_0 < distance_1_m1)) {
                            p_b = {dataPts[idx_0_m1][0], dataPts[idx_0_m1][1], dataPts[idx_0_m1][2]};
                            p_c = {dataPts[idx_1_0][0], dataPts[idx_1_0][1], dataPts[idx_1_0][2]};
                        }
                        else if ((distance_m1_m1 < distance_m1_0) and (distance_m1_m1 < distance_1_0) and (distance_m1_m1 < distance_1_m1)) {
                            p_b = {dataPts[idx_m1_m1][0], dataPts[idx_m1_m1][1], dataPts[idx_m1_m1][2]};
                            p_c = {dataPts[idx_0_m1][0], dataPts[idx_0_m1][1], dataPts[idx_0_m1][2]};
                        }
                        else if ((distance_1_m1 < distance_m1_0) and (distance_1_m1 < distance_1_0) and (distance_1_m1 < distance_m1_m1)) {
                            p_b = {dataPts[idx_0_m1][0], dataPts[idx_0_m1][1], dataPts[idx_0_m1][2]};
                            p_c = {dataPts[idx_1_m1][0], dataPts[idx_1_m1][1], dataPts[idx_1_m1][2]};
                        }
                    }

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
                        final_drift_vector = {dataPts[idx_0_0][0]-dataPts[idx_0_m1][0], dataPts[idx_0_0][1]-dataPts[idx_0_m1][1], dataPts[idx_0_0][2]-dataPts[idx_0_m1][2]};
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

                drift_results_file << std::setprecision(15) 
                << test_point[0] << ' '
                << test_point[1] << ' '
                << test_point[2] << ' '
                << final_drift_vector[0] << ' '
                << final_drift_vector[1] << ' '
                << final_drift_vector[2] << std::endl;
            };
        drift_results_file.close();
        background.StopServerFront();
        std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "Node " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};
