
/*!
\file drift_current_sheet.hh
\brief Defines functions for calculating current sheet drift velocity
\author Sammy Siegel

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef DRIFT_CURRENT_SHEET_HH
#define DRIFT_CURRENT_SHEET_HH

#include <ANN/ANN.h>
#include "common/print_warn.hh"
#include "common/vectors.hh"
#include <vector>

namespace Spectrum {

std::vector<int> i_to_streamline(int idx, int num_lats=216, int num_steps=30) {
    int lon_i = idx/(num_steps*num_lats);
    int lat_i = (idx%(num_steps*num_lats))/num_steps;
    int step_i = idx%num_steps;
    return {lon_i, lat_i, step_i};
};

int streamline_to_i(int lon_i, int lat_i, int step_i, int num_lons=216, int num_lats=216, int num_steps=30) {
    int N = num_lats*num_lons*num_steps;
    return (lon_i*num_lats*num_steps + lat_i*num_steps + step_i + N)%N;
};

int get_lon_next(int idx, int num_lons=216) {
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

int get_lon_prev(int idx, int num_lons=216) {
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

int get_lat_next(int idx, int num_lons=216) {
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

int get_lat_prev(int idx, int num_lons=216) {
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

double get_distance_idx(GeoVector p1, GeoVector p0) {
    double dx, dy, dz;
    dx = p0[0] - p1[0];
    dy = p0[1] - p1[1];
    dz = p0[2] - p1[2];
    return (dx*dx + dy*dy + dz*dz);
};

double get_cs_drift_speed(double d_rg) {
    return 0.457 - 0.412*d_rg + 0.0915*d_rg*d_rg;
};

bool readPt(std::istream &in, ANNpoint p)			// read point (false on EOF)
{
	for (int i = 0; i < 3; i++) {
		if(!(in >> p[i])) return false;
	}
	return true;
}

void printPt(std::ostream &out, ANNpoint p)			// print point
{
	out << "(" << p[0];
	for (int i = 1; i < 3; i++) {
		out << ", " << p[i];
	}
	out << ")\n";
}

void CalculateCurrentSheetDrift(GeoVector test_point, double r_L, GeoVector* cs_drift_vel_unit, double* side_of_sheet, ANNpointArray* dataPts, ANNkd_tree* kdTree) {
    ANNpoint			queryPt;				// query point
    ANNidxArray			nnIdx;					// near neighbor indices
    ANNdistArray		dists;                  // nearest distances
    int				    k				= 1;	// number of nearest neighbors
    int				    dim				= 3;	// dimension
    double			    eps				= 0;	// error bound
    int                 num_lons        = 216;
    int                 num_lats        = 216;
    int                 num_steps       = 30;
    double              polarity        = 1.0;

    double                  distance_nearest, final_distance, direction_mag, drift_speed;
    double                  distance_1_0, distance_m1_0, distance_0_1, distance_0_m1, distance_1_1, distance_m1_1, distance_1_m1, distance_m1_m1;
    int                     idx_0_0, idx_1_0, idx_m1_0, idx_0_1, idx_0_m1, idx_1_1, idx_m1_1, idx_1_m1, idx_m1_m1;
    std::vector<double>     p_a, p_b, p_c, p_ab, p_ac, p_a0, cross_product, final_drift_vector;   
    
    queryPt = annAllocPt(dim);                  // allocate query point
    nnIdx = new ANNidx[k];						// allocate near neighbor indices
    dists = new ANNdist[k];						// allocate near neighbor dists

    queryPt[0] = test_point[0];
    queryPt[1] = test_point[1];
    queryPt[2] = test_point[2];

    kdTree->annkSearch(						// search
            queryPt,						// query point
            k,								// number of near neighbors
            nnIdx,							// nearest neighbors (returned)
            dists,							// distance (returned)
            eps);							// error bound

    idx_0_0 = nnIdx[0];
    distance_nearest = sqrt(dists[0]);
    p_a = {(*dataPts)[idx_0_0][0], (*dataPts)[idx_0_0][1], (*dataPts)[idx_0_0][2]};

    if(idx_0_0%num_steps==0) {
        final_distance = sqrt(distance_nearest);
        *side_of_sheet = (test_point[2]-p_a[2])/abs(test_point[2]-p_a[2]);
        if (final_distance/r_L < 2.0) {
            drift_speed = get_cs_drift_speed(final_distance/r_L);
            idx_0_m1 = get_lat_prev(idx_0_0);
            final_drift_vector = {(*dataPts)[idx_0_0][0]-(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_0][1]-(*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_0][2]-(*dataPts)[idx_0_m1][2]};
            direction_mag = sqrt(final_drift_vector[0]*final_drift_vector[0]+final_drift_vector[1]*final_drift_vector[1]+final_drift_vector[2]*final_drift_vector[2]);
            final_drift_vector[0] = polarity*final_drift_vector[0]*drift_speed/direction_mag;
            final_drift_vector[1] = polarity*final_drift_vector[1]*drift_speed/direction_mag;
            final_drift_vector[2] = polarity*final_drift_vector[2]*drift_speed/direction_mag;
            *cs_drift_vel_unit = GeoVector(final_drift_vector[0], final_drift_vector[1], final_drift_vector[2]);
        }
        else {
            drift_speed = 0.0;
            final_drift_vector = {0.0, 0.0, 0.0};
            *cs_drift_vel_unit = gv_zeros;
        };
    }
    else {
        idx_1_0 = get_lon_next(idx_0_0);
        distance_1_0 = get_distance_idx({(*dataPts)[idx_1_0][0], (*dataPts)[idx_1_0][1], (*dataPts)[idx_1_0][2]}, test_point);
        idx_m1_0 = get_lon_prev(idx_0_0);
        distance_m1_0 = get_distance_idx({(*dataPts)[idx_m1_0][0], (*dataPts)[idx_m1_0][1], (*dataPts)[idx_m1_0][2]}, test_point);
        idx_0_1 = get_lat_next(idx_0_0);
        distance_0_1 = get_distance_idx({(*dataPts)[idx_0_1][0], (*dataPts)[idx_0_1][1], (*dataPts)[idx_0_1][2]}, test_point);
        idx_0_m1 = get_lat_prev(idx_0_0);
        distance_0_m1 = get_distance_idx({(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_m1][2]}, test_point);
        idx_1_1 = get_lat_next(idx_1_0);
        distance_1_1 = get_distance_idx({(*dataPts)[idx_1_1][0], (*dataPts)[idx_1_1][1], (*dataPts)[idx_1_1][2]}, test_point);
        idx_m1_1 = get_lat_next(idx_m1_0);
        distance_m1_1 = get_distance_idx({(*dataPts)[idx_m1_1][0], (*dataPts)[idx_m1_1][1], (*dataPts)[idx_m1_1][2]}, test_point);
        idx_1_m1 = get_lat_prev(idx_1_0);
        distance_1_m1 = get_distance_idx({(*dataPts)[idx_1_m1][0], (*dataPts)[idx_1_m1][1], (*dataPts)[idx_1_m1][2]}, test_point);
        idx_m1_m1 = get_lat_prev(idx_m1_m1);
        distance_m1_m1 = get_distance_idx({(*dataPts)[idx_m1_m1][0], (*dataPts)[idx_m1_m1][1], (*dataPts)[idx_m1_m1][2]}, test_point);

        if ((distance_1_0 < distance_m1_0) and (distance_1_0 < distance_0_1) and (distance_1_0 < distance_0_m1)) {
            if ((distance_1_1 < distance_1_m1) and (distance_1_1 < distance_0_1) and (distance_1_1 < distance_0_m1)) {
                p_b = {(*dataPts)[idx_1_0][0], (*dataPts)[idx_1_0][1], (*dataPts)[idx_1_0][2]};
                p_c = {(*dataPts)[idx_1_1][0], (*dataPts)[idx_1_1][1], (*dataPts)[idx_1_1][2]};
            }
            else if ((distance_1_m1 < distance_1_1) and (distance_1_m1 < distance_0_1) and (distance_1_m1 < distance_0_m1)) {
                p_b = {(*dataPts)[idx_1_m1][0], (*dataPts)[idx_1_m1][1], (*dataPts)[idx_1_m1][2]};
                p_c = {(*dataPts)[idx_1_0][0], (*dataPts)[idx_1_0][1], (*dataPts)[idx_1_0][2]};
            }
            else if ((distance_0_1 < distance_1_1) and (distance_0_1 < distance_1_m1) and (distance_0_1 < distance_0_m1)) {
                p_b = {(*dataPts)[idx_1_0][0], (*dataPts)[idx_1_0][1], (*dataPts)[idx_1_0][2]};
                p_c = {(*dataPts)[idx_0_1][0], (*dataPts)[idx_0_1][1], (*dataPts)[idx_0_1][2]};
            }
            else if ((distance_0_m1 < distance_1_1) and (distance_0_m1 < distance_1_m1) and (distance_0_m1 < distance_0_1)) {
                p_b = {(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_m1][2]};
                p_c = {(*dataPts)[idx_1_0][0], (*dataPts)[idx_1_0][1], (*dataPts)[idx_1_0][2]};
            }
        }
        else if ((distance_m1_0 < distance_1_0) and (distance_m1_0 < distance_0_1) and (distance_m1_0 < distance_0_m1)) {
            if ((distance_0_1 < distance_0_m1) and (distance_0_1 < distance_m1_1) and (distance_0_1 < distance_m1_m1)) {
                p_b = {(*dataPts)[idx_0_1][0], (*dataPts)[idx_0_1][1], (*dataPts)[idx_0_1][2]};
                p_c = {(*dataPts)[idx_m1_0][0], (*dataPts)[idx_m1_0][1], (*dataPts)[idx_m1_0][2]};
            }
            else if ((distance_0_m1 < distance_0_1) and (distance_0_m1 < distance_m1_1) and (distance_0_m1 < distance_m1_m1)) {
                p_b = {(*dataPts)[idx_m1_0][0], (*dataPts)[idx_m1_0][1], (*dataPts)[idx_m1_0][2]};
                p_c = {(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_m1][2]};
            }
            else if ((distance_m1_1 < distance_0_1) and (distance_m1_1 < distance_0_m1) and (distance_m1_1 < distance_m1_m1)) {
                p_b = {(*dataPts)[idx_m1_1][0], (*dataPts)[idx_m1_1][1], (*dataPts)[idx_m1_1][2]};
                p_c = {(*dataPts)[idx_m1_0][0], (*dataPts)[idx_m1_0][1], (*dataPts)[idx_m1_0][2]};
            }
            else if ((distance_m1_m1 < distance_0_1) and (distance_m1_m1 < distance_0_m1) and (distance_m1_m1 < distance_m1_1)) {
                p_b = {(*dataPts)[idx_m1_0][0], (*dataPts)[idx_m1_0][1], (*dataPts)[idx_m1_0][2]};
                p_c = {(*dataPts)[idx_m1_m1][0], (*dataPts)[idx_m1_m1][1], (*dataPts)[idx_m1_m1][2]};
            }
        }
        else if ((distance_0_1 < distance_1_0) and (distance_0_1 < distance_m1_0) and (distance_0_1 < distance_0_m1)) {
            if ((distance_m1_1 < distance_1_1) and (distance_m1_1 < distance_m1_0) and (distance_m1_1 < distance_1_0)) {
                p_b = {(*dataPts)[idx_0_1][0], (*dataPts)[idx_0_1][1], (*dataPts)[idx_0_1][2]};
                p_c = {(*dataPts)[idx_m1_1][0], (*dataPts)[idx_m1_1][1], (*dataPts)[idx_m1_1][2]};
            }
            else if ((distance_1_1 < distance_m1_1) and (distance_1_1 < distance_m1_0) and (distance_1_1 < distance_1_0)) {
                p_b = {(*dataPts)[idx_1_1][0], (*dataPts)[idx_1_1][1], (*dataPts)[idx_1_1][2]};
                p_c = {(*dataPts)[idx_0_1][0], (*dataPts)[idx_0_1][1], (*dataPts)[idx_0_1][2]};
            }
            else if ((distance_m1_0 < distance_m1_1) and (distance_m1_0 < distance_1_1) and (distance_m1_0 < distance_1_0)) {
                p_b = {(*dataPts)[idx_0_1][0], (*dataPts)[idx_0_1][1], (*dataPts)[idx_0_1][2]};
                p_c = {(*dataPts)[idx_m1_1][0], (*dataPts)[idx_m1_1][1], (*dataPts)[idx_m1_1][2]};
            }
            else if ((distance_1_0 < distance_m1_1) and (distance_1_0 < distance_1_1) and (distance_1_0 < distance_m1_0)) {
                p_b = {(*dataPts)[idx_1_0][0], (*dataPts)[idx_1_0][1], (*dataPts)[idx_1_0][2]};
                p_c = {(*dataPts)[idx_0_1][0], (*dataPts)[idx_0_1][1], (*dataPts)[idx_0_1][2]};
            }
        }
        else if ((distance_0_m1 < distance_1_0) and (distance_0_m1 < distance_m1_0) and (distance_0_m1 < distance_0_1)) {
            if ((distance_m1_0 < distance_1_0) and (distance_m1_0 < distance_m1_m1) and (distance_m1_0 < distance_1_m1)) {
                p_b = {(*dataPts)[idx_m1_0][0], (*dataPts)[idx_m1_0][1], (*dataPts)[idx_m1_0][2]};
                p_c = {(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_m1][2]};
            }
            else if ((distance_1_0 < distance_m1_0) and (distance_1_0 < distance_m1_m1) and (distance_1_0 < distance_1_m1)) {
                p_b = {(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_m1][2]};
                p_c = {(*dataPts)[idx_1_0][0], (*dataPts)[idx_1_0][1], (*dataPts)[idx_1_0][2]};
            }
            else if ((distance_m1_m1 < distance_m1_0) and (distance_m1_m1 < distance_1_0) and (distance_m1_m1 < distance_1_m1)) {
                p_b = {(*dataPts)[idx_m1_m1][0], (*dataPts)[idx_m1_m1][1], (*dataPts)[idx_m1_m1][2]};
                p_c = {(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_m1][2]};
            }
            else if ((distance_1_m1 < distance_m1_0) and (distance_1_m1 < distance_1_0) and (distance_1_m1 < distance_m1_m1)) {
                p_b = {(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_m1][2]};
                p_c = {(*dataPts)[idx_1_m1][0], (*dataPts)[idx_1_m1][1], (*dataPts)[idx_1_m1][2]};
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
        *side_of_sheet = final_distance / abs(final_distance);
        final_distance = final_distance / *side_of_sheet;
        
        if (final_distance/r_L < 2.0) {
            drift_speed = get_cs_drift_speed(final_distance/r_L);
            final_drift_vector = {(*dataPts)[idx_0_0][0]-(*dataPts)[idx_0_m1][0], (*dataPts)[idx_0_0][1]-(*dataPts)[idx_0_m1][1], (*dataPts)[idx_0_0][2]-(*dataPts)[idx_0_m1][2]};
            direction_mag = sqrt(final_drift_vector[0]*final_drift_vector[0]+final_drift_vector[1]*final_drift_vector[1]+final_drift_vector[2]*final_drift_vector[2]);
            final_drift_vector[0] = polarity*final_drift_vector[0]*drift_speed/direction_mag;
            final_drift_vector[1] = polarity*final_drift_vector[1]*drift_speed/direction_mag;
            final_drift_vector[2] = polarity*final_drift_vector[2]*drift_speed/direction_mag;
            *cs_drift_vel_unit = GeoVector(final_drift_vector[0], final_drift_vector[1], final_drift_vector[2]);
        }
        else {
            drift_speed = 0.0;
            final_drift_vector = {0.0, 0.0, 0.0};
            *cs_drift_vel_unit = gv_zeros;
        };
    };
};

};

#endif


