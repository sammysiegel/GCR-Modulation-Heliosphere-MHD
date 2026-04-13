#include "src/server_config.hh"
#include "src/background_server_batl.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <ANN/ANN.h>
#include "src/drift_current_sheet.hh"


using namespace Spectrum;

int				k				= 1;			// number of nearest neighbors
int				dim				= 3;			// dimension
double			eps				= 0;			// error bound
int				N			= 1000;			// maximum number of data points

std::istream*		dataIn			= NULL;			// input for data points

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
      
        int num_lats = 216;
        int num_lons = 216;
        int num_steps = 100;
        double polarity = 1.0;
        int i;
        int x_i, z_i;
        int N_x = 4000;
        int N_z = 4000;
        double x_start = -100.0;
        double z_start = -100.0;
        double delta_x = 0.05;
        double delta_z = 0.05;
        int N = num_lats*num_lons*num_steps;
        std::cout << "N is " << N << std::endl;

        std::string infilename = "current_sheet_points_subset_hires.txt";
        std::ifstream input_file;
        input_file.open(infilename, std::ios::in);
        dataIn = &input_file;
        std::ofstream drift_results_file;

        int					nPts;					// actual number of data points
        ANNpointArray		dataPts;				// data points
        ANNkd_tree*			kdTree;					// search structure

        // getArgs(argc, argv);						// read command-line arguments

        dataPts = annAllocPts(N, dim);			// allocate data points

        nPts = 0;									// read data points

        std::cout << "Reading Data Points...\n";
        while (nPts < N && readPt(*dataIn, dataPts[nPts])) {
            nPts++;
        }

        kdTree = new ANNkd_tree(					// build search structure
                    dataPts,					// the data points
                    nPts,						// number of points
                    dim);						// dimension of space

        GeoVector test_point;
        auto start = std::chrono::high_resolution_clock::now();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        double side_of_sheet;
        double r_L;
        GeoVector cs_drift_vel_unit;
        GeoVector gc_drift_vel;
        GeoVector drift_vel;
        drift_results_file.open("drift_data.txt");

        for (x_i=0; x_i<N_x; x_i++) {
            for (z_i=0; z_i<N_z; z_i++) {
                test_point = {x_start+x_i*delta_x, 0.0, z_start+z_i*delta_z};
                if (test_point[0]*test_point[0]+test_point[2]*test_point[2] < 1){
                    cs_drift_vel_unit = gv_zeros;
                    gc_drift_vel = gv_zeros;
                    drift_vel = gv_zeros;
                    side_of_sheet = 0;
                }
                else{
                    pos[0] = test_point[0]*one_au;
                    pos[1] = test_point[1]*one_au;
                    pos[2] = test_point[2]*one_au;
                    background.GetFields(t, pos, mom, spdata);
                    if (spdata.region[0] < -99.9) {
                        continue;
                    };
                    r_L = LarmorRadius(mom[0], spdata.Bmag, specie)/one_au;
                    if (z_i==0) {
                        start = std::chrono::high_resolution_clock::now();
                    };

                    CalculateCurrentSheetDrift(pos, r_L, &cs_drift_vel_unit, &side_of_sheet, &dataPts, kdTree);

                    if (z_i==0) {
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
                        std::cout << duration.count() << " microseconds" << std::endl;
                    };

                    gc_drift_vel = (spdata.curlB() - 2.0 * (spdata.gradBmag ^ spdata.bhat)) * side_of_sheet / spdata.Bmag;
                    gc_drift_vel *= r_L * vel[0] / 3.0;
                    drift_vel = gc_drift_vel + cs_drift_vel_unit * vel[0];
                };

                drift_results_file << std::setprecision(6) 
                << test_point[0] << ' '
                << test_point[1] << ' '
                << test_point[2] << ' '
                << cs_drift_vel_unit[0] << ' '
                << cs_drift_vel_unit[1] << ' '
                << cs_drift_vel_unit[2] << ' '
                << gc_drift_vel[0] << ' '
                << gc_drift_vel[1] << ' '
                << gc_drift_vel[2] << ' '
                << drift_vel[0] << ' '
                << drift_vel[1] << ' '
                << drift_vel[2] << ' '
                << side_of_sheet << std::endl;
            };
        };
        drift_results_file.close();
        background.StopServerFront();
        std::cout << "Background samples outputted to file." << std::endl;
   };

   std::cout << "Node " << mpi_config->glob_comm_rank << " exited." << std::endl;
   return 0;
};
