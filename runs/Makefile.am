## Process this file with automake to produce Makefile.in
ANN_HOME = /data001/cosmicrays_vf/Sammy/SPECTRUM_devel/ann_1.1.2
AM_CXXFLAGS = $(MPI_CFLAGS) -I$(ANN_HOME)/include
AM_LDFLAGS = $(MPI_LIBS)
AM_FCFLAGS = -J$(BATL_HOME)/share/include -I$(ANN_HOME)/include

bin_PROGRAMS = gcr_modulation_example postprocess_gcr_modulation_results main_test_background_batl main_current_sheet_3d main_gcr_modulation_batl test_calc_diffusion test_calc_drifts_kd_header

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src
BATL_HOME = /data001/cosmicrays_vf/Sammy/BATL

gcr_modulation_example_SOURCES = gcr_modulation_test.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/trajectory_parker.cc \
   $(SPBL_SOURCE_DIR)/trajectory_parker.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_time.cc \
   $(SPBL_SOURCE_DIR)/initial_time.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/source_base.cc \
   $(SPBL_SOURCE_DIR)/source_base.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_SOURCE_DIR)/server_config.hh \
   $(SPBL_COMMON_DIR)/workload_manager.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

gcr_modulation_example_LDADD = $(MPI_LIBS) $(GSL_LIBS)

postprocess_gcr_modulation_results_SOURCES = postprocess_gcr_modulation_results.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

postprocess_gcr_modulation_results_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_test_background_batl_SOURCES = main_test_background_batl.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.hh \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_server.cc \
   $(SPBL_SOURCE_DIR)/background_server.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_test_background_batl_LDADD = -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh -lmpi

main_current_sheet_3d_SOURCES = main_current_sheet_3d.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.hh \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_server.cc \
   $(SPBL_SOURCE_DIR)/background_server.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

main_current_sheet_3d_LDADD = -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh -lmpi

main_gcr_modulation_batl_SOURCES = gcr_modulation_test_batl.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/trajectory_parker.cc \
   $(SPBL_SOURCE_DIR)/trajectory_parker.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_time.cc \
   $(SPBL_SOURCE_DIR)/initial_time.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/drift_current_sheet.hh \
   $(SPBL_SOURCE_DIR)/source_base.cc \
   $(SPBL_SOURCE_DIR)/source_base.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_SOURCE_DIR)/server_config.hh \
   $(SPBL_COMMON_DIR)/workload_manager.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_SOURCE_DIR)/background_server_batl.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.hh \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_server.cc \
   $(SPBL_SOURCE_DIR)/background_server.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(ANN_HOME)/include/ANN/ANN.h

main_gcr_modulation_batl_LDADD = -L$(ANN_HOME)/lib $(GSL_LIBS) -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lmpi_mpifh -lmpi -lgfortran -lANN

test_calc_diffusion_SOURCES = test_calc_diffusion_rigidity.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.hh \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_server.cc \
   $(SPBL_SOURCE_DIR)/background_server.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh

test_calc_diffusion_LDADD = -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh -lmpi

test_calc_drifts_kd_header_SOURCES = test_calc_drifts_kd_header.cc \
   $(SPBL_SOURCE_DIR)/drift_current_sheet.hh \
   $(SPBL_SOURCE_DIR)/background_server_batl.cc \
   $(SPBL_SOURCE_DIR)/background_server_batl.hh \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/background_server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_server.cc \
   $(SPBL_SOURCE_DIR)/background_server.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_batl.cc \
   $(SPBL_SOURCE_DIR)/server_batl.hh \
   $(SPBL_SOURCE_DIR)/server_cartesian.cc \
   $(SPBL_SOURCE_DIR)/server_cartesian.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_batl.cc \
   $(SPBL_SOURCE_DIR)/block_batl.hh \
   $(SPBL_SOURCE_DIR)/block_cartesian.cc \
   $(SPBL_SOURCE_DIR)/block_cartesian.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/spectrum_interface.f90 \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(ANN_HOME)/include/ANN/ANN.h

test_calc_drifts_kd_header_LDADD = -L$(ANN_HOME)/lib -L$(BATL_HOME)/lib -lREADAMR -lTIMING -lSHARE -lgfortran -lmpi_mpifh -lmpi -lANN
