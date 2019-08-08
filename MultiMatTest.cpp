/*
 * Copyright (c) 2017-2019, Triad National Security, LLC.
 * All rights Reserved.
 *
 * This is the code released under LANL Copyright Disclosure C17041/LA-CC-17-041
 * Copyright 2017-2019.  Triad National Security, LLC. This material was
 * produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 * National Laboratory (LANL), which is operated by Triad National Security, LLC
 * for the U.S. Department of Energy. See LICENSE file for details.
 *
 * Released under the New BSD License
 *
 * Bob Robey brobey@lanl.gov and Rao Garimella rao@lanl.gov
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "genmalloc.h"
#include "input.h"
#include "timer.h"

#define CLOCK_RATE 2.7e9    // GHz laptop
#define STREAM_RATE 13375;  // MB/sec

#define MININT(i, j) (((i) < (j)) ? (i) : (j))

void get_neighbors(int ncells, int *&num_nbrs, int **&nbrs);
void get_centroids(double (*&cen)[2]);
void get_vol_frac_matrix_rand(double **&Volfrac, float &filled_percentage);
void get_vol_frac_matrix_file(double **&Volfrac, float &filled_percentage);
void get_vol_frac_matrix(int method, double **&Volfrac,
                         float &filled_percentage);
void setup_cell_dominant_data_structure(int method, double *&Vol,
                                        double *&Density, double *&Temperature,
                                        double *&Pressure, double **&Volfrac,
                                        double **&Densityfrac,
                                        double **&Temperaturefrac,
                                        double **&Pressurefrac,
                                        float &filled_percentage);
void setup_material_dominant_data_structure(
    int method, double *&Vol, double *&Density, double *&Temperature,
    double *&Pressure, double **&Volfrac, double **&Densityfrac,
    double **&Temperaturefrac, double **&Pressurefrac,
    float &filled_percentage);
void setup_cell_dominant_compact_data_structure(
    int method, int *&imaterial, int *&nmaterials, double *&Vol,
    double *&Density, double *&Temperature, double *&Pressure,
    int *&imaterialfrac, int *&nextfrac, int *&frac2cell, double *&Volfrac,
    double *&Densityfrac, double *&Temperaturefrac, double *&Pressurefrac,
    float &filled_percentage);
void setup_mat_dominant_compact_data_structure(
    int method, int **&subset2mesh, int **&mesh2subset, int *&nmatscell,
    int *&matids, int *&ncellsmat, double *&Vol, double *&Density,
    double **&Volfrac, double **&Densityfrac, double **&Temperaturefrac,
    double **&Pressurefrac, float &filled_percentage);

void convert_compact_material_2_compact_cell(
    int ncells, int **subset2mesh, int **mesh2subset, int *nmatscell,
    int *matids, int *ncellsmat, double **Volfrac, double **Densityfrac,
    double **Temperaturefrac, double **Pressurefrac, int *&Cimaterial,
    int *&Cnmaterials, int *&Cimaterialfrac, int *&Cnextfrac, int *&Cfrac2cell,
    double *&CVolfrac, double *&CDensityfrac, double *&CTemperaturefrac,
    double *&CPressurefrac);
void convert_compact_cell_2_compact_material(
    int ncells, int nmats, int *Cimaterial, int *Cnmaterials,
    int *Cimaterialfrac, double *CVol, double *CDensity, double *CTemperature,
    double *CVolfrac, double *CDensityfrac, double *CTemperaturefrac,
    int **&subset2mesh, int **&mesh2subset, int *&nmatscell, int *&matids,
    int *&ncellsmat, double **&Volfrac, double **&Densityfrac,
    double **&Temperaturefrac, double **&Pressurefrac);
void print_performance_estimates(float est_perf, int64_t memops8byte,
                                 int64_t memops4byte, int64_t flops,
                                 float penalty_msecs);

bool verbose = false;
bool memory_verbose = true;
int itermax = 100;
int ncells = 1000000;
int nmats = 50;
float est_perf, act_perf, model_error;

int main(int argc, char **argv) {
  struct timespec tstart_cpu;
  double density_ave;
  double time_sum;
  double VolTotal;
  double *Density_average;
  double nmatconst = 5.0;
  double nmatconsts[nmats];
  for (int i = 0; i < nmats; i++) {
    nmatconsts[i] = 5.0;
  }
  int64_t memops, memops4byte, memops8byte;
  int64_t flops;
  float penalty_msecs;
  float filled_percentage, filled_fraction;
  int method = 0;  // VF initialization: 0 - random, 1 - read volfrac.dat

  if (argc > 1) sscanf(argv[1], "%d", &method);

  printf("Run stream benchmark for your system\n");
  printf(
      "L3 Cache on Macbook Pro is 6MB so problem size is just bigger at 16MB "
      "min\n");
  printf(
      "First test should give Stream Benchmark or problem size is too small\n");
  printf("Second problem should give about twice the first\n");

  // Some globals

  float L_f = method ? 0.5 : 1.0;  // ave frac of nbrs containing material
  int nnbrs_ave = 8;  // nearly so; 4000 boundary cells in 1 million cells
  //                  // in 3D, nnbrs_ave would be 26

  // Build up list of neighbors for each cell
  // Assuming a 2D structured mesh, each cell will have a maximum of 8 nbrs

  int *nnbrs = (int *)genvector("Num_Neighbors", ncells, sizeof(int));
  int **nbrs = (int **)genmatrix("Neighbors", ncells, nmats, sizeof(int));

  get_neighbors(ncells, nnbrs, nbrs);

  // Compute centroids of cells
  double(*cen)[2] =
      (double(*)[2])genvector("Centroids", ncells, sizeof(double[2]));
  get_centroids(cen);

  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Starting Single Material Data Structure
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("=======================================\n");
    printf("Starting Single Material Data Structure\n");
    printf("=======================================\n\n");

    double *Density = (double *)genvector("Density", ncells, sizeof(double));
    double *Temperature =
        (double *)genvector("Temperature", ncells, sizeof(double));
    double *Pressure = (double *)genvector("Pressure", ncells, sizeof(double));
    double *Vol = (double *)genvector("Volume", ncells, sizeof(double));

    VolTotal = 0.0;
    for (int ic = 0; ic < ncells; ic++) {
      Density[ic] = 2.0;
      Temperature[ic] = 0.5;
      Vol[ic] = 1.0;
      VolTotal += Vol[ic];
    }

    if (memory_verbose) {
      genmalloc_MB_memory_report();
    }
    genmalloc_MB_memory_total();
    printf("\n");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with cell densities (pure cells)
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      density_ave = 0.0;
      for (int ic = 0; ic < ncells; ic++) {
        density_ave += Density[ic] * Vol[ic];  // Nc stores, Nc*3 loads
      }
      density_ave /= VolTotal;

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf("Average Density of pure cells    %lf, compute time is %lf msecs\n",
           density_ave, act_perf);
    printf("Updated prediction: %d\n", ncells + ncells*3);
    memops = 2 * ncells;
    flops = 2 * ncells;
    penalty_msecs = 0.0;
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Calculate pressure using ideal gas law
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        //Nc*4 loads, Nc stores
        Pressure[ic] = (nmatconst * Density[ic] * Temperature[ic]) / Vol[ic];
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Pressure Calculation for cell,             compute time is %lf "
        "msecs\n",
        act_perf);
    printf("Updated prediction: %d\n", ncells*5);

    memops = 4 * ncells;
    flops = 3 * ncells;
    penalty_msecs = 0.0;
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    genvectorfree(Vol);
    genvectorfree(Density);
    genvectorfree(Temperature);
    genvectorfree(Pressure);
  }

  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Starting Cell-Dominant Full Matrix Data Structure
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n");
    printf("=================================================\n");
    printf("Starting Cell-Dominant Full Matrix Data Structure\n");
    printf("=================================================\n\n");

    double *Vol, *Density, *Temperature, *Pressure;
    double **Densityfrac, **Temperaturefrac, **Pressurefrac, **Volfrac;

    setup_cell_dominant_data_structure(
        method, Vol, Density, Temperature, Pressure, Volfrac, Densityfrac,
        Temperaturefrac, Pressurefrac, filled_percentage);

    filled_fraction = filled_percentage / 100.0;

    if (memory_verbose) {
      genmalloc_MB_memory_report();
    }
    genmalloc_MB_memory_total();
    printf("\n");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Density_average =
        (double *)genvector("Density_average", ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        density_ave = 0.0; // Nc stores
        for (int m = 0; m < nmats; m++) {
          density_ave += Densityfrac[ic][m] * Volfrac[ic][m]; // Nc*Nm*4 + Nc*2 memops
        }
        Density_average[ic] = density_ave / Vol[ic]; // Nc*2 loads, Nc stores
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of mixed material cells    compute time is %lf "
        "msecs\n",
        act_perf);
    printf("Updated prediction: %d\n", ncells * nmats * 4 + ncells * 6);

    genvectorfree(Density_average);

    memops = 2 * ncells * nmats;  // line 4 loads
    memops += 2 * ncells;         // line 6 stores
    flops = 2 * ncells * nmats;   // line 4 flops
    penalty_msecs = 0.0;
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities with if test
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Density_average =
        (double *)genvector("Density_average", ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        density_ave = 0.0; // Nc stores
        for (int m = 0; m < nmats; m++) {
          if (Volfrac[ic][m] > 0.0) { // Nc*Nm + Nc loads
            // (Nc*Nm*3 + Nc)*frac_filled memops
            density_ave += Densityfrac[ic][m] * Volfrac[ic][m];
          }
        }
        Density_average[ic] = density_ave / Vol[ic]; // Nc*2 loads, Nc stores
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of frac with if            compute time is %lf "
        "msecs\n",
        act_perf);
    printf("Updated prediction: %d\n", ncells + ncells * nmats + ncells * nmats * 6 * 0.13 + ncells * 4);

    genvectorfree(Density_average);

    float cache_miss_freq = method ? 0.7 : 1.0;
    memops = ncells * nmats;  // line 4 loads
    memops += (int64_t)(filled_fraction *
                        (float)(2 * ncells * nmats));  // line 5 loads
    memops += 2 * ncells;  // line 8 stores and loads
    flops = (int64_t)(filled_fraction *
                      (float)(2 * ncells * nmats));  // line 5 flops
    flops += ncells;                                 // line 8 flops
    float branch_wait =
        1.0 / CLOCK_RATE * 16;  // Estimate a 16 cycle wait for branch
                                // misprediction for a 2.7 GHz processor
    float cache_wait =
        1.0 / CLOCK_RATE * 7 * 16;  // Estimate a 7*16 or 112 cycle wait for
                                    // missing prefetch for a 2.7 GHz processor
    penalty_msecs = 1000.0 * cache_miss_freq * (branch_wait + cache_wait) *
                    (filled_fraction * (float)(ncells * nmats));  // line 4 if
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //   Calculate pressure using ideal gas law
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        for (int m = 0; m < nmats; m++) {
          if (Volfrac[ic][m] > 0.) { // Nc*Nm + Nc loads
            // (Nc*Nm*4 + Nc*3) * filled_frac memops
            Pressurefrac[ic][m] =
                (nmatconsts[m] * Densityfrac[ic][m] * Temperaturefrac[ic][m]) /
                (Volfrac[ic][m]);
          } else {
            // (Nc*Nm + Nc) * 1-filled_frac memops
            Pressurefrac[ic][m] = 0.0;
          }
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Pressure Calculation of mixed material cells with if compute time is "
        "%lf msecs\n",
        act_perf);

    float sparsity_fraction = 1.0 - filled_fraction;
    memops = ncells * nmats;  // line 3 loads
    memops +=
        (int64_t)(filled_fraction * (float)(ncells * nmats));  // line 5 stores
    memops += (int64_t)(filled_fraction *
                        (float)(3 * ncells * nmats));  // line 6 loads
    memops += (int64_t)(sparsity_fraction *
                        (float)(ncells * nmats));  // line 8 stores
    flops = (int64_t)(filled_fraction *
                      (float)(3 * ncells * nmats));  // line 6 flops
    branch_wait =
        1.0 / CLOCK_RATE * 16;  // Estimate a 16 cycle wait for branch
                                // misprediction for a 2.7 GHz processor
    cache_wait =
        1.0 / CLOCK_RATE * 7 * 16;  // Estimate a 7*16 or 112 cycle wait for
                                    // missing prefetch for a 2.7 GHz processor
    penalty_msecs = 1000.0 * cache_miss_freq * (branch_wait + cache_wait) *
                    (filled_fraction * (float)(ncells * nmats));  // line 3 if
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average material density over neighborhood of each cell
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double **MatDensity_average = (double **)genmatrix(
        "MatDensity_average", ncells, nmats, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        double xc[2];
        xc[0] = cen[ic][0];
        xc[1] = cen[ic][1];
        int nn = nnbrs[ic];
        int cnbrs[8];
        double dsqr[8];
        for (int n = 0; n < nn; n++) cnbrs[n] = nbrs[ic][n];
        for (int n = 0; n < nn; n++) {
          dsqr[n] = 0.0;
          for (int d = 0; d < 1; d++) {
            double ddist = (xc[d] - cen[cnbrs[n]][d]);
            dsqr[n] += ddist * ddist;
          }
        }
        for (int m = 0; m < nmats; m++) {
          if (Volfrac[ic][m] > 0.0) {
            int nnm = 0;  // number of nbrs with this material
            for (int n = 0; n < nn; n++) {
              int jc = cnbrs[n];
              if (Volfrac[jc][m] > 0.0) {
                MatDensity_average[ic][m] += Densityfrac[ic][m] / dsqr[n];
                nnm++;
              }
            }
            MatDensity_average[ic][m] /= nnm;
          } else {
            MatDensity_average[ic][m] = 0.0;
          }
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf("Average Material Density            compute time is %lf msecs\n",
           act_perf);

    genmatrixfree((void **)MatDensity_average);

    memops = (2 + 2 * nmats + (0.5 + 16) * nnbrs_ave) * ncells;
    //        // Formula differs from paper because it is 2D here
    memops +=
        (int64_t)(filled_fraction * 8 * (1 + L_f) * ncells * nmats * nnbrs_ave);
    flops = 6 * ncells * nnbrs_ave;
    flops += (int64_t)(filled_fraction * 3 * ncells * nmats * nnbrs_ave * L_f);
    branch_wait =
        1.0 / CLOCK_RATE * 16;  // Estimate a 16 cycle wait for branch
                                // misprediction for a 2.7 GHz processor
    cache_wait =
        1.0 / CLOCK_RATE * 7 * 16;  // Estimate a 7*16 or 112 cycle wait for
                                    // missing prefetch for a 2.7 GHz processor
    penalty_msecs = 1000.0 * cache_miss_freq * (branch_wait + cache_wait) *
                    (filled_fraction * (float)(ncells * nmats));
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    genvectorfree((void *)Vol);
    genvectorfree((void *)Density);
    genvectorfree((void *)Temperature);
    genvectorfree((void *)Pressure);
    genmatrixfree((void **)Volfrac);
    genmatrixfree((void **)Densityfrac);
    genmatrixfree((void **)Temperaturefrac);
    genmatrixfree((void **)Pressurefrac);
  }

  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Starting Material-Dominant Full Matrix Data Structure
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n");
    printf("===================================================\n");
    printf("Starting Material-Dominant Full Matrix Data Structure\n");
    printf("===================================================\n");

    double *Vol, *Density, *Temperature, *Pressure;
    double **Volfrac, **Densityfrac, **Temperaturefrac, **Pressurefrac;

    setup_material_dominant_data_structure(
        method, Vol, Density, Temperature, Pressure, Volfrac, Densityfrac,
        Temperaturefrac, Pressurefrac, filled_percentage);

    filled_fraction = filled_percentage / 100.0;

    if (memory_verbose) {
      genmalloc_MB_memory_report();
    }
    genmalloc_MB_memory_total();
    printf("\n");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Density_average =
        (double *)genvector("Density_average", ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        Density_average[ic] = 0.0;
      }
      for (int m = 0; m < nmats; m++) {
        for (int ic = 0; ic < ncells; ic++) {
          Density_average[ic] += Densityfrac[m][ic] * Volfrac[m][ic];
        }
      }
      for (int ic = 0; ic < ncells; ic++) {
        Density_average[ic] /= Vol[ic];
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of mixed material cells    compute time is %lf "
        "msecs\n",
        act_perf);

    genvectorfree(Density_average);

    memops = ncells;               // line 3 loads
    memops += ncells * nmats;      // line 6 stores
    memops += 3 * ncells * nmats;  // line 7 loads
    flops = 2 * ncells * nmats;    // line 7 flops
    memops += 3 * ncells;          // line 11 loads/stores
    flops += ncells;               // line 11 flops
    penalty_msecs = 0.0;
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities with if test
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Density_average =
        (double *)genvector("Density_average", ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        Density[ic] = 0.0;
      }
      for (int m = 0; m < nmats; m++) {
        for (int ic = 0; ic < ncells; ic++) {
          if (Volfrac[m][ic] > 0.0) {
            Density[ic] += Densityfrac[m][ic] * Volfrac[m][ic];
          }
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of mixed material cells with if compute time is %lf "
        "msecs\n",
        act_perf);

    genvectorfree(Density_average);

    float cache_miss_freq = method ? 0.2 : 1.0;
    memops = ncells;           // line 2 loads
    memops += ncells * nmats;  // line 6 loads
    memops +=
        (int64_t)(filled_fraction * (float)(ncells * nmats));  // line 7 stores
    memops +=
        (int64_t)(filled_fraction * (float)(ncells * nmats));  // line 8 loads
    flops = (int64_t)(filled_fraction *
                      (float)(2 * ncells * nmats));  // line 8 flops
    memops += 2 * ncells;                            // line 11 stores and loads
    flops += ncells;                                 // line 11 flops
    float branch_wait =
        1.0 / CLOCK_RATE * 16;  // Estimate a 16 cycle wait for branch
                                // misprediction for a 2.7 GHz processor
    float cache_wait =
        1.0 / CLOCK_RATE * 7 * 16;  // Estimate a 7*16 or 112 cycle wait for
                                    // missing prefetch for a 2.7 GHz processor
    penalty_msecs = 1000.0 * cache_miss_freq * (branch_wait + cache_wait) *
                    (filled_fraction * (float)(ncells * nmats));  // line 6 if
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Calculate pressure using ideal gas law
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int m = 0; m < nmats; m++) {
        for (int ic = 0; ic < ncells; ic++) {
          if (Volfrac[m][ic] > 0.0) {
            Pressurefrac[m][ic] =
                (nmatconsts[m] * Densityfrac[m][ic] * Temperaturefrac[m][ic]) /
                Volfrac[m][ic];
          } else {
            Pressurefrac[m][ic] = 0.0;
          }
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Pressure Calculation of frac with if       compute time is %lf "
        "msecs\n",
        act_perf);

    float sparsity_fraction = 1.0 - filled_fraction;
    memops = nmats;            // line 2 loads
    memops += ncells * nmats;  // line 4 loads
    memops +=
        (int64_t)(filled_fraction * (float)(ncells * nmats));  // line 5 stores
    memops += (int64_t)(filled_fraction *
                        (float)(2 * ncells * nmats));  // line 6 loads
    flops = (int64_t)(filled_fraction *
                      (float)(3 * ncells * nmats));  // line 6 flops
    memops += (int64_t)(sparsity_fraction *
                        (float)(ncells * nmats));  // line 8 stores
    branch_wait =
        1.0 / CLOCK_RATE * 16;  // Estimate a 16 cycle wait for branch
                                // misprediction for a 2.7 GHz processor
    cache_wait =
        1.0 / CLOCK_RATE * 7 * 16;  // Estimate a 7*16 or 112 cycle wait for
                                    // missing prefetch for a 2.7 GHz processor
    penalty_msecs = 1000.0 * cache_miss_freq * (branch_wait + cache_wait) *
                    (filled_fraction * (float)(ncells * nmats));  // line 6 if
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average material density over neighborhood of each cell
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double **MatDensity_average = (double **)genmatrix(
        "MatDensity_average", nmats, ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int m = 0; m < nmats; m++) {
        for (int ic = 0; ic < ncells; ic++) {
          if (Volfrac[m][ic] > 0.0) {
            double xc[2];
            xc[0] = cen[ic][0];
            xc[1] = cen[ic][1];
            int nn = nnbrs[ic];
            int cnbrs[8];
            double dsqr[8];
            for (int n = 0; n < nn; n++) cnbrs[n] = nbrs[ic][n];
            for (int n = 0; n < nn; n++) {
              dsqr[n] = 0.0;
              for (int d = 0; d < 1; d++) {
                double ddist = (xc[d] - cen[cnbrs[n]][d]);
                dsqr[n] += ddist * ddist;
              }
            }

            int nnm = 0;  // number of nbrs with this material
            for (int n = 0; n < nn; n++) {
              int jc = cnbrs[n];
              if (Volfrac[m][jc] > 0.0) {
                MatDensity_average[m][ic] += Densityfrac[m][ic] / dsqr[n];
                nnm++;
              }
            }
            MatDensity_average[m][ic] /= nnm;
          } else {
            MatDensity_average[m][ic] = 0.0;
          }
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf("Average Material Density            compute time is %lf msecs\n",
           act_perf);

    genmatrixfree((void **)MatDensity_average);

    memops = 2 * ncells * nmats;
    //        // Formula differs from paper because it is 2D here
    memops += (int64_t)(2 * filled_fraction * ncells * nmats);
    memops += (int64_t)(8.5 * filled_fraction * ncells * nmats * nnbrs_ave);
    memops +=
        (int64_t)(24 * filled_fraction * L_f * ncells * nmats * nnbrs_ave);
    flops = (int64_t)(filled_fraction * ncells * nmats);
    flops += (int64_t)(9 * filled_fraction * ncells * nmats * nnbrs_ave * L_f);
    branch_wait =
        1.0 / CLOCK_RATE * 16;  // Estimate a 16 cycle wait for branch
                                // misprediction for a 2.7 GHz processor
    cache_wait =
        1.0 / CLOCK_RATE * 7 * 16;  // Estimate a 7*16 or 112 cycle wait for
                                    // missing prefetch for a 2.7 GHz processor
    penalty_msecs = 1000.0 * cache_miss_freq * (branch_wait + cache_wait) *
                    (filled_fraction * ncells * nmats);
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    genvectorfree((void *)Vol);
    genvectorfree((void *)Density);
    genvectorfree((void *)Temperature);
    genvectorfree((void *)Pressure);
    genmatrixfree((void **)Volfrac);
    genmatrixfree((void **)Densityfrac);
    genmatrixfree((void **)Temperaturefrac);
    genmatrixfree((void **)Pressurefrac);
  }

  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Starting Cell-Dominant Compact Data Structure
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n");
    printf("===================================================\n");
    printf("Starting Cell-Dominant Compact Data Structure\n");
    printf("===================================================\n");

    int *imaterial, *nmaterials, *imaterialfrac, *nextfrac, *frac2cell;
    double *Vol, *Density, *Temperature, *Pressure, *Volfrac, *Densityfrac,
        *Temperaturefrac, *Pressurefrac;

    setup_cell_dominant_compact_data_structure(
        method, imaterial, nmaterials, Vol, Density, Temperature, Pressure,
        imaterialfrac, nextfrac, frac2cell, Volfrac, Densityfrac,
        Temperaturefrac, Pressurefrac, filled_percentage);

    filled_fraction = filled_percentage / 100.0;

    int nmixlength = 0;
    int pure_cell_count = 0;
    int mixed_cell_count = 0;
    for (int ic = 0; ic < ncells; ic++) {
      int ix = imaterial[ic];
      if (ix <= 0) {
        for (ix = -ix; ix >= 0; ix = nextfrac[ix]) nmixlength++;
        mixed_cell_count++;
      } else {
        pure_cell_count++;
      }
    }
    float mixed_cell_fraction = (float)mixed_cell_count / (float)ncells;
    float pure_cell_fraction = (float)pure_cell_count / (float)ncells;
    int nmats_ave =
        ((float)pure_cell_count + (float)nmixlength) / (float)ncells;

    if (memory_verbose) {
      genmalloc_MB_memory_report();
    }
    genmalloc_MB_memory_total();
    printf("\n");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        density_ave = 0.0;
        int ix = imaterial[ic];
        if (ix <= 0) {  // material numbers for clean cells start at 1
          for (ix = -ix; ix >= 0; ix = nextfrac[ix]) {
            density_ave += Densityfrac[ix] * Volfrac[ix];
          }
          Density[ic] = density_ave / Vol[ic];
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of mixed material cells    compute time is %lf "
        "msecs\n",
        act_perf);

    float cache_miss_freq = method ? 0.1 : 1.0;
    memops4byte = ncells;                         // line 3 loads
    memops4byte += nmixlength;                    // line 5 loads
    memops8byte = 2 * nmixlength;                 // line 6 loads
    flops = 2 * nmixlength;                       // line 6 flops
    memops8byte += mixed_cell_fraction * ncells;  // line 8 stores
    memops8byte += mixed_cell_fraction * ncells;  // line 8 loads
    flops += mixed_cell_fraction * ncells;        // line 8 flops
    float loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs = 1000.0 * loop_overhead * mixed_cell_fraction *
                    (float)ncells;  // line 5 for
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities using nmats
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        density_ave = 0.0;
        int mstart = imaterial[ic];
        if (mstart <= 0) {  // material numbers for clean cells start at 1
          mstart = -mstart;
          for (int ix = 0; ix < nmaterials[ic]; ix++) {
            density_ave += Densityfrac[mstart + ix] * Volfrac[mstart + ix];
          }
          Density[ic] = density_ave / Vol[ic];
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of mixed material cells with nmats  compute time is "
        "%lf msecs\n",
        act_perf);

    memops4byte = ncells;                         // line 3 loads
    memops4byte += mixed_cell_fraction * ncells;  // line 5 loads
    memops8byte = 2 * nmixlength;                 // line 6 loads
    flops = 2 * nmixlength;                       // line 6 flops
    memops8byte += mixed_cell_fraction * ncells;  // line 8 stores
    memops8byte += mixed_cell_fraction * ncells;  // line 8 loads
    flops += mixed_cell_fraction * ncells;        // line 8 flops
    loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs = 1000.0 * cache_miss_freq * loop_overhead *
                    mixed_cell_fraction * (float)ncells;  // line 5 for
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Density_average =
        (double *)genvector("Density_average", ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        density_ave = 0.0;
        int ix = imaterial[ic];
        if (ix <= 0) {  // material numbers for clean cells start at 1
          for (ix = -ix; ix >= 0; ix = nextfrac[ix]) {
            density_ave += Densityfrac[ix] * Volfrac[ix];
          }
          Density_average[ic] = density_ave / Vol[ic];
        } else {
          Density_average[ic] = Density[ic];
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of mixed material cells    compute time is %lf "
        "msecs\n",
        act_perf);

    genvectorfree(Density_average);

    memops4byte = ncells;                         // line 3 loads
    memops4byte += nmixlength;                    // line 5 loads
    memops8byte = 2 * nmixlength;                 // line 6 loads
    flops = 2 * nmixlength;                       // line 6 flops
    memops8byte += mixed_cell_fraction * ncells;  // line 8 stores
    memops8byte += mixed_cell_fraction * ncells;  // line 8 loads
    flops += mixed_cell_fraction * ncells;        // line 8 flops
    memops8byte += pure_cell_fraction * ncells;   // line 10 stores
    memops8byte += pure_cell_fraction * ncells;   // line 10 loads
    loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs = 1000.0 * cache_miss_freq * loop_overhead *
                    mixed_cell_fraction * (float)ncells;  // line 5 for
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities with pure calculation filler
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Density_average =
        (double *)genvector("Density_average", ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        int ix = imaterial[ic];
        if (ix <= 0) {  // material numbers for clean cells start at 1
          density_ave = 0.0;
          for (ix = -ix; ix >= 0; ix = nextfrac[ix]) {
            density_ave += Densityfrac[ix] * Volfrac[ix];
          }
        } else {  // Pure cell
          density_ave = Density[ic] * Vol[ic];
        }
        Density_average[ic] = density_ave / Vol[ic];
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of mixed materials cells pure filler  compute time is "
        "%lf msecs\n",
        act_perf);

    genvectorfree(Density_average);

    memops4byte = ncells;                                     // line 3 loads
    memops4byte += nmixlength;                                // line 5 loads
    memops8byte = 2 * nmixlength;                             // line 6 loads
    flops = 2 * nmixlength;                                   // line 6 flops
    memops8byte += (int64_t)mixed_cell_fraction * ncells;     // line 8 stores
    memops8byte += (int64_t)mixed_cell_fraction * ncells;     // line 8 loads
    flops += (int64_t)mixed_cell_fraction * ncells;           // line 8 flops
    memops8byte += (int64_t)pure_cell_fraction * ncells;      // line 10 stores
    memops8byte += (int64_t)2 * pure_cell_fraction * ncells;  // line 10 loads
    flops += (int64_t)pure_cell_fraction * ncells;            // line 8 flops
    loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs = 1000.0 * cache_miss_freq * loop_overhead *
                    mixed_cell_fraction * (float)ncells;  // line 5 for
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //   Calculate pressure using ideal gas law
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        int ix = imaterial[ic];
        if (ix <= 0) {  // material numbers for clean cells start at 1
          for (ix = -ix; ix >= 0; ix = nextfrac[ix]) {
            int m = imaterialfrac[ix];
            Pressurefrac[ix] =
                (nmatconsts[m] * Densityfrac[ix] * Temperaturefrac[ix]) /
                Volfrac[ix];
          }
        } else {
          Pressure[ic] =
              nmatconsts[ix] * Density[ic] * Temperature[ic] / Vol[ic];
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Pressure Calculation of mixed material cells compute time is %lf "
        "msecs\n",
        act_perf);

    memops4byte = ncells;                                     // line 2 loads
    memops4byte += nmixlength;                                // line 4 loads
    memops4byte += nmixlength;                                // line 5 loads
    memops8byte = 3 * nmixlength;                             // line 6 loads
    memops8byte += nmixlength;                                // line 6 stores
    flops = 3 * nmixlength;                                   // line 6 flops
    memops8byte += (int64_t)pure_cell_fraction * ncells;      // line 9 stores
    memops8byte += (int64_t)4 * pure_cell_fraction * ncells;  // line 9 loads
    flops += (int64_t)3 * pure_cell_fraction * ncells;        // line 9 flops
    loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs = 1000.0 * cache_miss_freq * loop_overhead *
                    mixed_cell_fraction * (float)ncells;  // line 5 for
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average material density over neighborhood of each cell
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double **MatDensity_average = (double **)genmatrix(
        "MatDensity_average", ncells, nmats, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int ic = 0; ic < ncells; ic++) {
        for (int m = 0; m < nmats; m++) MatDensity_average[ic][m] = 0.0;

        double xc[2];
        xc[0] = cen[ic][0];
        xc[1] = cen[ic][1];
        int nn = nnbrs[ic];
        int cnbrs[8];
        double dsqr[8];
        for (int n = 0; n < nn; n++) cnbrs[n] = nbrs[ic][n];
        for (int n = 0; n < nn; n++) {
          dsqr[n] = 0.0;
          for (int d = 0; d < 1; d++) {
            double ddist = (xc[d] - cen[cnbrs[n]][d]);
            dsqr[n] += ddist * ddist;
          }
        }

        int ix = imaterial[ic];
        if (ix <= 0) {
          for (ix = -ix; ix >= 0; ix = nextfrac[ix]) {
            int m = imaterialfrac[ix];

            int nnm = 0;  // number of nbrs with this material
            for (int n = 0; n < nn; n++) {
              int jc = cnbrs[n];

              int jx = imaterial[jc];
              if (jx <= 0) {
                for (jx = -jx; jx >= 0; jx = nextfrac[jx]) {
                  if (imaterialfrac[jx] == m) {
                    MatDensity_average[ic][m] += Densityfrac[jx] / dsqr[n];
                    nnm++;
                    break;
                  }
                }
              } else {
                if (imaterialfrac[jx] == m) {
                  MatDensity_average[ic][m] += Densityfrac[jx] / dsqr[n];
                  nnm++;
                }
              }
            }
            MatDensity_average[ic][m] /= nnm;
          }
        } else {
          int m = imaterialfrac[ix];

          int nnm = 0;  // number of nbrs with this material
          for (int n = 0; n < nn; n++) {
            int jc = cnbrs[n];

            int jx = imaterial[jc];
            if (jx <= 0) {
              for (jx = -jx; jx >= 0; jx = nextfrac[jx]) {
                if (imaterialfrac[jx] == m) {
                  MatDensity_average[ic][m] += Densityfrac[jx] / dsqr[n];
                  nnm++;
                  break;
                }
              }
            } else {
              if (imaterialfrac[jx] == m) {
                MatDensity_average[ic][m] += Densityfrac[jx] / dsqr[n];
                nnm++;
              }
            }
          }
          MatDensity_average[ic][m] /= nnm;
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf("Average Material Density            compute time is %lf msecs\n",
           act_perf);

    genmatrixfree((void **)MatDensity_average);

    filled_fraction = filled_percentage / 100.0;
    // Formula differs a bit from paper because it is 2D here
    memops = (int64_t)(2.5 * ncells * (1 + nnbrs_ave) + 0.5 * nmixlength);
    memops += (int64_t)(ncells * nmats * (1 + 1.5 * filled_fraction));
    memops += (int64_t)(4 * filled_fraction * ncells * nmats * nnbrs_ave);
    memops +=
        (int64_t)(8 * filled_fraction * ncells * nmats * nnbrs_ave * nmats_ave);
    memops += (int64_t)(8 * filled_fraction * ncells * nmats * nnbrs_ave * L_f);
    flops = 6 * ncells * nnbrs_ave;
    flops += (int64_t)(3 * filled_fraction * ncells * nmats * nnbrs_ave * L_f);
    flops += (int64_t)(filled_fraction * ncells * nmats);
    float branch_wait =
        1.0 / CLOCK_RATE * 16;  // Estimate a 16 cycle wait for branch
                                // misprediction for a 2.7 GHz processor
    float cache_wait =
        1.0 / CLOCK_RATE * 7 * 16;  // Estimate a 7*16 or 112 cycle wait for
                                    // missing prefetch for a 2.7 GHz processor
    penalty_msecs = 1000.0 * cache_miss_freq * (branch_wait + cache_wait) *
                    mixed_cell_fraction * ncells;
    print_performance_estimates(act_perf, memops, 0, flops, penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    genvectorfree(imaterial);
    genvectorfree(nmaterials);
    genvectorfree(Vol);
    genvectorfree(Density);
    genvectorfree(Temperature);
    genvectorfree(Pressure);
    genvectorfree(imaterialfrac);
    genvectorfree(nextfrac);
    genvectorfree(frac2cell);
    genvectorfree(Volfrac);
    genvectorfree(Densityfrac);
    genvectorfree(Temperaturefrac);
    genvectorfree(Pressurefrac);
  }

  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Starting Material-Centric Compact Data Structure
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n");
    printf("===================================================\n");
    printf("Starting Material-Centric Compact Data Structure\n");
    printf("===================================================\n");

    int *nmatscell, *matids, *ncellsmat;
    int **subset2mesh, **mesh2subset;
    double *Vol, *Density;
    double **Volfrac, **Densityfrac, **Temperaturefrac, **Pressurefrac;
    float filled_percentage;
    float cache_miss_freq = method ? 0.2 : 1.0;

    setup_mat_dominant_compact_data_structure(
        method, subset2mesh, mesh2subset, nmatscell, matids, ncellsmat, Vol,
        Density, Volfrac, Densityfrac, Temperaturefrac, Pressurefrac,
        filled_percentage);

    filled_fraction = filled_percentage / 100.0;

    if (memory_verbose) {
      genmalloc_MB_memory_report();
    }
    genmalloc_MB_memory_total();
    printf("\n");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities - MAT-DOMINANT LOOP
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int C = 0; C < ncells; C++) Density[C] = 0.0;

      for (int m = 0; m < nmats; m++) {
        for (int c = 0; c < ncellsmat[m]; c++) {  // Note that this is c not C
          int C = subset2mesh[m][c];
          Density[C] += Densityfrac[m][c] * Volfrac[m][c];
        }
      }

      for (int C = 0; C < ncells; C++) Density[C] /= Vol[C];

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of cells - Mat Dominant loop  -  compute time is %lf "
        "msecs\n",
        act_perf);

    float ninner = 0.0;  // number of times inner loop executed
    for (int m = 0; m < nmats; m++)
      for (int c = 0; c < ncellsmat[m]; c++) ninner++;
    // ninner = F_f*N_m*N_c
    memops8byte = ncells;  // Initialization of Density

    memops4byte = nmats;   // load ncmat
    memops8byte += nmats;  // load subset2mesh

    memops8byte += 4 * ninner;  // load subset, Density, Densityfrac, Volfrac
    memops8byte += ninner;      // store Density

    memops8byte += 2 * ncells;  // load density, Vol
    memops8byte += ncells;      // Store Density

    flops = 2 * ninner;  // multiply and add
    flops += ncells;     // divide cell density by cell volume

    // memops8byte  += (int64_t) (8*cache_miss_freq +
    // (1-cache_miss_freq))*ninner ;    // load Density (cache miss, reload 8
    // doubles)
    float loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs =
        0.0;  // Only if we account for removed materials with an if-check
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average density with fractional densities - CELL-DOMINANT LOOP
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int C = 0; C < ncells; C++) {
        double density_ave = 0.0;
        for (int im = 0; im < nmatscell[C]; im++) {
          int m = matids[4 * C + im];
          int c = mesh2subset[m][C];
          density_ave += Densityfrac[m][c] * Volfrac[m][c];
        }
        Density[C] = density_ave / Vol[C];
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Average Density of cells - Cell Dominant loop  -  compute time is %lf "
        "msecs\n",
        act_perf);

    ninner = 0;
    for (int C = 0; C < ncells; C++)
      for (int im = 0; im < nmatscell[C]; im++) ninner++;

    memops4byte = ncells;   // load nmatscells
    memops4byte += ninner;  // load matids
    memops4byte += (int64_t)(16 * cache_miss_freq + (1 - cache_miss_freq)) *
                   ninner;  // load mesh2subset (cache miss, reload 16 integers)
    memops8byte =
        (int64_t)(8 * cache_miss_freq + (1 - cache_miss_freq)) * 2 *
        ninner;  // load Densityfrac, Volfrac (cache miss, reload 8 doubles)
    memops8byte += ncells;  // load of cell volume Vol
    memops8byte += ncells;  // Store of Density
    flops = 2 * ninner;     // multiply and add
    flops += ncells;        // divide density_ave by Vol
    loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs =
        0.0;  // Only if we account for removed materials with an if-check
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //   Calculate pressure using ideal gas law - MAT-CENTRIC COMPACT STRUCTURE
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int m = 0; m < nmats; m++) {
        double matconst = nmatconsts[m];
        for (int c = 0; c < ncellsmat[m]; c++) {
          Pressurefrac[m][c] =
              (matconst * Densityfrac[m][c] * Temperaturefrac[m][c]) /
              Volfrac[m][c];
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf(
        "Pressure Calculation of cells - Mat-Centric -  compute time is %lf "
        "msecs\n",
        act_perf);

    ninner = 0;
    for (int m = 0; m < nmats; m++)
      for (int c = 0; c < ncellsmat[m]; c++) ninner++;
    memops4byte = 0;
    memops8byte = nmats;        // load of nmatsconsts
    memops8byte += 3 * ninner;  // load DensityFrac, TemperatureFrac, VolFrac
    memops8byte += ninner;      // store PressureFrac
    flops = 3 * ninner;
    loop_overhead =
        1.0 / CLOCK_RATE * 20;  // Estimate a 20 cycle loop exit overhead
    penalty_msecs =
        0.0;  // Only if we account for removed materials with an if-check
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    Average material density over neighborhood of each cell
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double **MatDensity_average = (double **)genmatrix(
        "MatDensity_average", nmats, ncells, sizeof(double));

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      cpu_timer_start(&tstart_cpu);

      for (int m = 0; m < nmats; m++) {
        for (int C = 0; C < ncells; C++) MatDensity_average[m][C] = 0.0;

        for (int c = 0; c < ncellsmat[m]; c++) {  // Note that this is c not C
          int C = subset2mesh[m][c];
          double xc[2];
          xc[0] = cen[C][0];
          xc[1] = cen[C][1];
          int nn = nnbrs[C];
          int cnbrs[9];
          double dsqr[8];
          for (int n = 0; n < nn; n++) cnbrs[n] = nbrs[C][n];
          for (int n = 0; n < nn; n++) {
            dsqr[n] = 0.0;
            for (int d = 0; d < 1; d++) {
              double ddist = (xc[d] - cen[cnbrs[n]][d]);
              dsqr[n] += ddist * ddist;
            }
          }

          int nnm = 0;  // number of nbrs with this material
          for (int n = 0; n < nn; n++) {
            int C_j = cnbrs[n];
            int c_j = mesh2subset[m][C_j];
            if (c_j >= 0) {
              MatDensity_average[m][C] += Densityfrac[m][c_j] / dsqr[n];
              nnm++;
            }
          }
          MatDensity_average[m][C] /= nnm;
        }
      }

      time_sum += cpu_timer_stop(tstart_cpu);
    }
    act_perf = time_sum * 1000.0 / itermax;
    printf("Average Material Density  -  compute time is %lf msecs\n",
           act_perf);

    memops8byte = ncells * nmats;
    memops8byte += (int64_t)24.5 * filled_fraction * ncells * nmats;
    memops8byte += (int64_t)8 * filled_fraction * ncells * nmats * nnbrs_ave;
    memops8byte +=
        (int64_t)17 * filled_fraction * ncells * nmats * nnbrs_ave * L_f;
    flops = (int64_t)8 * filled_fraction * ncells * nmats * nnbrs_ave * L_f;
    flops += (int64_t)filled_fraction * ncells * nmats;
    penalty_msecs = 0.0;
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Convert from MATERIAL-CENTRIC COMPACT DATA STRUCTURE to CELL_CENTRIC
//   COMPACT DATA STRUCTURE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#define CONVERSION_CHECK 1
#ifdef CONVERSION_CHECK
    int *CCimaterial, *CCnmaterials, *CCimaterialfrac, *CCnextfrac,
        *CCfrac2cell;
    double *CCVol, *CCDensity, *CCTemperature, *CCPressure, *CCVolfrac,
        *CCDensityfrac, *CCTemperaturefrac, *CCPressurefrac;

    setup_cell_dominant_compact_data_structure(
        method, CCimaterial, CCnmaterials, CCVol, CCDensity, CCTemperature,
        CCPressure, CCimaterialfrac, CCnextfrac, CCfrac2cell, CCVolfrac,
        CCDensityfrac, CCTemperaturefrac, CCPressurefrac, filled_percentage);
#endif

    time_sum = 0;
    for (int iter = 0; iter < itermax; iter++) {
      int *Cimaterial;
      int *Cnmaterials;
      int *Cimaterialfrac;
      int *Cnextfrac;
      int *Cfrac2cell;
      double *CVolfrac;
      double *CDensityfrac;
      double *CTemperaturefrac;
      double *CPressurefrac;

      cpu_timer_start(&tstart_cpu);
      convert_compact_material_2_compact_cell(
          ncells, subset2mesh, mesh2subset, nmatscell, matids, ncellsmat,
          Volfrac, Densityfrac, Temperaturefrac, Pressurefrac, Cimaterial,
          Cnmaterials, Cimaterialfrac, Cnextfrac, Cfrac2cell, CVolfrac,
          CDensityfrac, CTemperaturefrac, CPressurefrac);
      time_sum += cpu_timer_stop(tstart_cpu);

#ifdef CONVERSION_CHECK
      int ix = 0;
      for (int ic = 0; ic < ncells; ic++) {
        int CNmats = nmatscell[ic];
        if (CNmats == 1) {
          if (CCimaterial[ic] != Cimaterial[ic] && ic < 100) {
            printf(
                "DEBUG line %d ic %d CNmats %d CCimaterial %d Cimaterial %d\n",
                __LINE__, ic, CNmats, CCimaterial[ic], Cimaterial[ic]);
          }
          if (matids[ic * 4] + 1 != Cimaterial[ic]) {
            printf("DEBUG line %d ic %d CNmats %d matids %d Cimaterial %d\n",
                   __LINE__, ic, CNmats, matids[ic * 4], Cimaterial[ic]);
            exit(1);
          }
        } else {
          ix = abs(Cimaterial[ic]);
          int m = 0;
          while (ix > 0) {
            if (matids[ic * 4 + m] + 1 != Cimaterialfrac[ix]) {
              printf(
                  "DEBUG CNmats %d ix %d mixed material %d Cimaterialfrac %d\n",
                  CNmats, ix, matids[ic * 4 + m] + 1, Cimaterialfrac[ix]);
            }
            if (CCimaterialfrac[ix] != Cimaterialfrac[ix]) {
              printf(
                  "DEBUG CNmats %d ix %d CCimaterialfrac %d Cimaterialfrac "
                  "%d\n",
                  CNmats, ix, CCimaterialfrac[ix], Cimaterialfrac[ix]);
            }
            ix = Cnextfrac[ix];
            m++;
          }
        }
      }
      exit(0);
#endif

      genvectorfree(Cimaterial);
      genvectorfree(Cimaterialfrac);
      genvectorfree(Cnextfrac);
      genvectorfree(Cfrac2cell);
      genvectorfree(CVolfrac);
      genvectorfree(CDensityfrac);
      genvectorfree(CTemperaturefrac);
      genvectorfree(CPressurefrac);
    }

#ifdef CONVERSION_CHECK
    genvectorfree(CCimaterial);
    genvectorfree(CCimaterialfrac);
    genvectorfree(CCnextfrac);
    genvectorfree(CCfrac2cell);
    genvectorfree(CCVolfrac);
    genvectorfree(CCDensityfrac);
    genvectorfree(CCTemperaturefrac);
    genvectorfree(CCPressurefrac);
#endif

    printf(
        "Conversion from compact material data structure to compact cell data "
        "structure\n");
    act_perf = time_sum * 1000.0 / itermax;
    // 4 arrays read from and stored to. Add 8x penalty for non-contiguous reads
    memops8byte = (4 * 8 + 4) * 5.0e5;
    // reads with 8x penalty for non-contiguous reads
    memops4byte = 3 * 8 * ncells + 2 * 8 * 5.0e5 + 3 * 4 * (ncells + 5.0e5);
    flops = 0.1;
    print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                                penalty_msecs);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    genvectorfree((void *)nmatscell);
    genvectorfree((void *)matids);
    genvectorfree((void *)Vol);
    genvectorfree((void *)Density);

    for (int m = 0; m < nmats; m++) {
      genvectorfree((void *)Volfrac[m]);
      genvectorfree((void *)Densityfrac[m]);
      genvectorfree((void *)Temperaturefrac[m]);
      genvectorfree((void *)Pressurefrac[m]);
      genvectorfree((void *)subset2mesh[m]);
      genvectorfree((void *)mesh2subset[m]);
    }
    genvectorfree((void *)Volfrac);
    genvectorfree((void *)Densityfrac);
    genvectorfree((void *)Temperaturefrac);
    genvectorfree((void *)Pressurefrac);
    genvectorfree((void *)subset2mesh);
    genmatrixfree((void **)mesh2subset);
  }

  int *Cimaterial, *Cnmaterials, *Cimaterialfrac, *Cnextfrac, *Cfrac2cell;
  double *CVol, *CDensity, *CTemperature, *CPressure, *CVolfrac, *CDensityfrac,
      *CTemperaturefrac, *CPressurefrac;

  setup_cell_dominant_compact_data_structure(
      method, Cimaterial, Cnmaterials, CVol, CDensity, CTemperature, CPressure,
      Cimaterialfrac, Cnextfrac, Cfrac2cell, CVolfrac, CDensityfrac,
      CTemperaturefrac, CPressurefrac, filled_percentage);
  for (int ic = 0; ic < ncells; ic++) {
    if (Cnmaterials[ic] < 1 || Cnmaterials[ic] > nmats) {
      printf("DEBUG -- ic %d Cnmaterials %d\n", ic, Cnmaterials[ic]);
    }
  }

  time_sum = 0;
  for (int iter = 0; iter < itermax; iter++) {
    int **subset2mesh;
    int **mesh2subset;
    int *nmatscell;
    int *matids;
    int *ncellsmat;
    double **Volfrac;
    double **Densityfrac;
    double **Temperaturefrac;
    double **Pressurefrac;

    cpu_timer_start(&tstart_cpu);

    convert_compact_cell_2_compact_material(
        ncells, nmats, Cimaterial, Cnmaterials, Cimaterialfrac, CVol, CDensity,
        CTemperature, CVolfrac, CDensityfrac, CTemperaturefrac, subset2mesh,
        mesh2subset, nmatscell, matids, ncellsmat, Volfrac, Densityfrac,
        Temperaturefrac, Pressurefrac);

    time_sum += cpu_timer_stop(tstart_cpu);

    for (int m = 0; m < nmats; m++) {
      genvectorfree((void *)subset2mesh[m]);
      genvectorfree((void *)Volfrac[m]);
      genvectorfree((void *)Densityfrac[m]);
      genvectorfree((void *)Temperaturefrac[m]);
      genvectorfree((void *)Pressurefrac[m]);
    }
    genvectorfree((void **)subset2mesh);
    genvectorfree((void **)Volfrac);
    genvectorfree((void **)Densityfrac);
    genvectorfree((void **)Temperaturefrac);
    genvectorfree((void **)Pressurefrac);

    genvectorfree((void *)matids);
    genvectorfree((void *)ncellsmat);

    genmatrixfree((void **)mesh2subset);
  }

  printf(
      "Conversion from compact cell data structure to compact material data "
      "structure\n");
  act_perf = time_sum * 1000.0 / itermax;
  // 4 arrays read from and stored to. Add 8x penalty for non-contiguous reads
  memops8byte = (3 * 8 + 3) * (ncells + 5.0e5);
  // reads with 8x penalty for non-contiguous reads
  memops4byte = 3 * 8 * ncells + 2 * 8 * 5.0e5 + 3 * 4 * (ncells + 5.0e5);
  memops4byte += 1 * (ncells * nmats);
  flops = 0.1;
  print_performance_estimates(act_perf, memops8byte, memops4byte, flops,
                              penalty_msecs);

  genvectorfree(Cimaterial);
  genvectorfree(Cimaterialfrac);
  genvectorfree(Cnextfrac);
  genvectorfree(Cfrac2cell);
  genvectorfree(CVolfrac);
  genvectorfree(CDensityfrac);
  genvectorfree(CTemperaturefrac);
  genvectorfree(CPressurefrac);
}

void setup_cell_dominant_data_structure(int method, double *&Vol,
                                        double *&Density, double *&Temperature,
                                        double *&Pressure, double **&Volfrac,
                                        double **&Densityfrac,
                                        double **&Temperaturefrac,
                                        double **&Pressurefrac,
                                        float &filled_percentage) {
  Vol = (double *)genvector("Volume", ncells, sizeof(double));
  Density = (double *)genvector("Density", ncells, sizeof(double));
  Temperature = (double *)genvector("Temperature", ncells, sizeof(double));
  Pressure = (double *)genvector("Pressure", ncells, sizeof(double));
  Densityfrac =
      (double **)genmatrix("DensityFrac", ncells, nmats, sizeof(double));
  Temperaturefrac =
      (double **)genmatrix("TemperatureFrac", ncells, nmats, sizeof(double));
  Pressurefrac =
      (double **)genmatrix("PressureFrac", ncells, nmats, sizeof(double));

  get_vol_frac_matrix(method, Volfrac, filled_percentage);

  for (int ic = 0; ic < ncells; ic++) {
    Vol[ic] = 0.0;
    for (int m = 0; m < nmats; m++) {
      if (Volfrac[ic][m] > 0.0) {
        Densityfrac[ic][m] = 2.0;
        Temperaturefrac[ic][m] = 0.5;
        Vol[ic] += Volfrac[ic][m];
      } else {
        Densityfrac[ic][m] = 0.0;
        Temperaturefrac[ic][m] = 0.0;
      }
      Pressurefrac[ic][m] = 0.0;
    }
  }
}

void setup_material_dominant_data_structure(
    int method, double *&Vol, double *&Density, double *&Temperature,
    double *&Pressure, double **&Volfrac, double **&Densityfrac,
    double **&Temperaturefrac, double **&Pressurefrac,
    float &filled_percentage) {
  Vol = (double *)genvector("Volume", ncells, sizeof(double));
  Density = (double *)genvector("Density", ncells, sizeof(double));
  Temperature = (double *)genvector("Temperature", ncells, sizeof(double));
  Pressure = (double *)genvector("Pressure", ncells, sizeof(double));
  Volfrac = (double **)genmatrix("VolumeFrac", nmats, ncells, sizeof(double));
  Densityfrac =
      (double **)genmatrix("DensityFrac", nmats, ncells, sizeof(double));
  Temperaturefrac =
      (double **)genmatrix("TemperatureFrac", nmats, ncells, sizeof(double));
  Pressurefrac =
      (double **)genmatrix("PressureFrac", nmats, ncells, sizeof(double));

  double **Volfrac_fullcc;  // cell centric full matrix of fractional volumes
  get_vol_frac_matrix(method, Volfrac_fullcc, filled_percentage);

  for (int m = 0; m < nmats; m++) {
    for (int ic = 0; ic < ncells; ic++) {
      if (Volfrac_fullcc[ic][m] > 0.0) {
        Volfrac[m][ic] = Volfrac_fullcc[ic][m];
        Densityfrac[m][ic] = 2.0;
        Temperaturefrac[m][ic] = 0.5;
      } else {
        Volfrac[m][ic] = Volfrac_fullcc[ic][m];
        Densityfrac[m][ic] = 0.0;
        Temperaturefrac[m][ic] = 0.0;
      }
      Pressurefrac[m][ic] = 0.0;
    }
  }

  // Now free the full data structures
  genmatrixfree((void **)Volfrac_fullcc);
}

void setup_cell_dominant_compact_data_structure(
    int method, int *&imaterial, int *&nmaterials, double *&Vol,
    double *&Density, double *&Temperature, double *&Pressure,
    int *&imaterialfrac, int *&nextfrac, int *&frac2cell, double *&Volfrac,
    double *&Densityfrac, double *&Temperaturefrac, double *&Pressurefrac,
    float &filled_percentage) {
  imaterial = (int *)genvector("imaterial", ncells, sizeof(int));
  nmaterials = (int *)genvector("nmaterials", ncells, sizeof(int));
  Vol = (double *)genvector("Vol", ncells, sizeof(double));
  Density = (double *)genvector("Density", ncells, sizeof(double));
  Temperature = (double *)genvector("Temperature", ncells, sizeof(double));
  Pressure = (double *)genvector("Pressure", ncells, sizeof(double));

  double **Volfrac_fullcc;  // full cell-centric matrix of fractional volumes
  get_vol_frac_matrix(method, Volfrac_fullcc, filled_percentage);

  int ix = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int nnz = 0;
    for (int im = 0; im < nmats; im++) {
      if (Volfrac_fullcc[ic][im] > 0.0) nnz++;
    }
    if (nnz > 1) ix += nnz;
  }
  int mxsize = ix;

  imaterialfrac = (int *)genvector("imaterialfrac", mxsize, sizeof(int));
  nextfrac = (int *)genvector("nextfrac", mxsize, sizeof(int));
  frac2cell = (int *)genvector("frac2cell", mxsize, sizeof(int));
  Volfrac = (double *)genvector("Volfrac", mxsize, sizeof(double));
  Densityfrac = (double *)genvector("Densityfrac", mxsize, sizeof(double));
  Temperaturefrac =
      (double *)genvector("Temperaturefrac", mxsize, sizeof(double));
  Pressurefrac = (double *)genvector("Pressurefrac", mxsize, sizeof(double));

  ix = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int m1, m2, m3, m4;
    int nnz = 0;
    for (int im = 0; im < nmats; im++) {
      if (Volfrac_fullcc[ic][im] > 0.0) {
        nnz++;
        if (nnz == 1)
          m1 = im;
        else if (nnz == 2)
          m2 = im;
        else if (nnz == 3)
          m3 = im;
        else if (nnz == 4)
          m4 = im;
      }
    }
    if (nnz == 4) {
      imaterial[ic] = -ix;
      nmaterials[ic] = 4;
      imaterialfrac[ix] = (m1 + 1);
      imaterialfrac[ix + 1] = (m2 + 1);
      imaterialfrac[ix + 2] = (m3 + 1);
      imaterialfrac[ix + 3] = (m4 + 1);
      Volfrac[ix] = 0.4;
      Volfrac[ix + 1] = 0.3;
      Volfrac[ix + 2] = 0.2;
      Volfrac[ix + 3] = 0.1;
      Densityfrac[ix] = 2.0;
      Densityfrac[ix + 1] = 2.0;
      Densityfrac[ix + 2] = 2.0;
      Densityfrac[ix + 3] = 2.0;
      Temperaturefrac[ix] = 0.5;
      Temperaturefrac[ix + 1] = 0.5;
      Temperaturefrac[ix + 2] = 0.5;
      Temperaturefrac[ix + 3] = 0.5;
      nextfrac[ix] = ix + 1;
      nextfrac[ix + 1] = ix + 2;
      nextfrac[ix + 2] = ix + 3;
      nextfrac[ix + 3] = -1;
      frac2cell[ix] = ic;
      frac2cell[ix + 1] = ic;
      frac2cell[ix + 2] = ic;
      frac2cell[ix + 3] = ic;
      ix += 4;
    } else if (nnz == 3) {
      imaterial[ic] = -ix;
      nmaterials[ic] = 3;
      imaterialfrac[ix] = (m1 + 1);
      imaterialfrac[ix + 1] = (m2 + 1);
      imaterialfrac[ix + 2] = (m3 + 1);
      Volfrac[ix] = 0.5;
      Volfrac[ix + 1] = 0.3;
      Volfrac[ix + 2] = 0.2;
      Densityfrac[ix] = 2.0;
      Densityfrac[ix + 1] = 2.0;
      Densityfrac[ix + 2] = 2.0;
      Temperaturefrac[ix] = 0.5;
      Temperaturefrac[ix + 1] = 0.5;
      Temperaturefrac[ix + 2] = 0.5;
      nextfrac[ix] = ix + 1;
      nextfrac[ix + 1] = ix + 2;
      nextfrac[ix + 2] = -1;
      frac2cell[ix] = ic;
      frac2cell[ix + 1] = ic;
      frac2cell[ix + 2] = ic;
      ix += 3;
    } else if (nnz == 2) {
      imaterial[ic] = -ix;
      nmaterials[ic] = 2;
      imaterialfrac[ix] = (m1 + 1);
      imaterialfrac[ix + 1] = (m2 + 1);
      Volfrac[ix] = 0.5;
      Volfrac[ix + 1] = 0.5;
      Densityfrac[ix] = 2.0;
      Densityfrac[ix + 1] = 2.0;
      Temperaturefrac[ix] = 0.5;
      Temperaturefrac[ix + 1] = 0.5;
      nextfrac[ix] = ix + 1;
      nextfrac[ix + 1] = -1;
      frac2cell[ix] = ic;
      frac2cell[ix + 1] = ic;
      ix += 2;
    } else {
      imaterial[ic] = (m1 + 1);
      nmaterials[ic] = 1;
    }
  }

  int filled_count = 0;
  int pure_cell_count = 0;
  int mixed_cell_count = 0;
  int mixed_frac_count = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int ix = imaterial[ic];
    if (ix > 0) {  // material numbers for clean cells start at 1
                   // clean cells
      Vol[ic] = 1.0;
      Density[ic] = 2.0;
      Temperature[ic] = 0.5;
      pure_cell_count++;
      filled_count++;
    } else {
      // multimaterial cells
      mixed_cell_count++;
      for (ix = -ix; ix >= 0; ix = nextfrac[ix]) {
        Vol[ic] += Volfrac[ix];
        filled_count++;
        mixed_frac_count++;
      }
    }
  }

#ifdef XXX
  for (int ic = 0; ic < ncells; ic++) {
    int ix = imaterial[ic];
    // clean cells
    printf("DEBUG cell is %3d Density %8.2lf material is %d\n", ic, Density[ic],
           ix);
  }

  for (int ix = 0; ix < mxsize; ix++) {
    printf("DEBUG mx is %d Dfrac %lf irf %d nxtf %d ixf %d\n", ix,
           Densityfrac[ix], imaterialfrac[ix], nextfrac[ix], frac2cell[ix]);
  }

  for (int ic = 0; ic < ncells; ic++) {
    int ix = imaterial[ic];
    if (ix > 0) {  // material numbers for clean cells start at 1
      // clean cells
      printf("DEBUG cell is %3d Density %8.2lf material is %d\n", ic,
             Density[ic], ix);
    } else {
      // multimaterial cells
      for (ix = -ix; ix >= 0; ix = nextfrac[ix]) {
        printf(
            "DEBUG cell is %3d mx is %3d Dfrac %8.2lf irf %d nxtf %4d ixf "
            "%3d\n",
            ic, ix, Densityfrac[ix], imaterialfrac[ix], nextfrac[ix],
            frac2cell[ix]);
        if (frac2cell[ix] != ic) {
          printf(
              "DEBUG -- error!!! mix item %d points to wrong cell %d, should "
              "be ic %d\n",
              ic, frac2cell[ix], ic);
          break;
        }
      }
    }
  }
#endif

  // Now free the full data structures
  genmatrixfree((void **)Volfrac_fullcc);
}

void setup_mat_dominant_compact_data_structure(
    int method, int **&subset2mesh, int **&mesh2subset, int *&nmatscell,
    int *&matids, int *&ncellsmat, double *&Vol, double *&Density,
    double **&Volfrac, double **&Densityfrac, double **&Temperaturefrac,
    double **&Pressurefrac, float &filled_percentage) {
  Vol = (double *)genvector("Volume", ncells, sizeof(double));
  Density = (double *)genvector("Density", ncells, sizeof(double));

  mesh2subset = (int **)genmatrix("mesh2subset", nmats, ncells, sizeof(int));
  for (int m = 0; m < nmats; m++)
    for (int C = 0; C < ncells; C++) mesh2subset[m][C] = -1;
  nmatscell = (int *)genvector("nmatscell", ncells, sizeof(int));
  for (int C = 0; C < ncells; C++) nmatscell[C] = 0;
  matids = (int *)genvector("matids", 4 * ncells, sizeof(int));
  for (int C = 0; C < 4 * ncells; C++) matids[C] = -1;
  ncellsmat = (int *)genvector("ncellsmat", nmats, sizeof(int));
  for (int m = 0; m < nmats; m++) ncellsmat[m] = 0;

  subset2mesh = (int **)genvector("subset2mesh", nmats, sizeof(int *));
  Volfrac = (double **)genvector("VolumeFrac", nmats, sizeof(double *));
  Densityfrac = (double **)genvector("DensityFrac", nmats, sizeof(double *));
  Temperaturefrac =
      (double **)genvector("TemperatureFrac", nmats, sizeof(double *));
  Pressurefrac = (double **)genvector("PressureFrac", nmats, sizeof(double *));

  double **Volfrac_fullcc;
  get_vol_frac_matrix(method, Volfrac_fullcc, filled_percentage);

  for (int ic = 0; ic < ncells; ic++) {
    nmatscell[ic] = 0;
    for (int m = 0; m < nmats; m++) {
      if (Volfrac_fullcc[ic][m] > 0.0) {
        matids[4 * ic + nmatscell[ic]] = m;
        nmatscell[ic]++;
        ncellsmat[m]++;
      }
    }
  }

  // Allocate compact data structures

  for (int m = 0; m < nmats; m++) {
    subset2mesh[m] =
        (int *)genvector("subset2mesh_m", ncellsmat[m], sizeof(int));
    Volfrac[m] = (double *)genvector("VolFrac_m", ncellsmat[m], sizeof(double));
    Densityfrac[m] =
        (double *)genvector("DensityFrac_m", ncellsmat[m], sizeof(double));
    Temperaturefrac[m] =
        (double *)genvector("TemperatureFrac_m", ncellsmat[m], sizeof(double));
    Pressurefrac[m] =
        (double *)genvector("PressureFrac_m", ncellsmat[m], sizeof(double));
  }

  // Now populate the compact data structures
  for (int m = 0; m < nmats; m++) ncellsmat[m] = 0;
  for (int C = 0; C < ncells; C++) {
    for (int im = 0; im < nmatscell[C]; im++) {
      int m = matids[4 * C + im];
      int c = ncellsmat[m];
      subset2mesh[m][c] = C;
      mesh2subset[m][C] = c;
      Volfrac[m][c] = Volfrac_fullcc[C][m];
      Densityfrac[m][c] = 2.0;
      Temperaturefrac[m][c] = 0.5;
      (ncellsmat[m])++;
    }
  }

  // Now free the full data structures
  genmatrixfree((void **)Volfrac_fullcc);
}

void convert_compact_material_2_compact_cell(
    int ncells, int **subset2mesh, int **mesh2subset, int *nmatscell,
    int *matids, int *ncellsmat, double **Volfrac, double **Densityfrac,
    double **Temperaturefrac, double **Pressurefrac, int *&Cimaterial,
    int *&Cnmaterials, int *&Cimaterialfrac, int *&Cnextfrac, int *&Cfrac2cell,
    double *&CVolfrac, double *&CDensityfrac, double *&CTemperaturefrac,
    double *&CPressurefrac) {
  Cimaterial = (int *)genvector("Cimaterial", ncells, sizeof(int));
  Cnmaterials = nmatscell;  // This is the same in both, but we assign pointer
                            // to maintain the naming convention

  // First we count up the number of mixed cells for the mix cell data
  // structure. We skip single material cells since they do not get put into the
  // compact structure.
  int mxsize = 0;
  for (int ic = 0; ic < ncells; ic++) {
    if (nmatscell[ic] > 1) mxsize += nmatscell[ic];
  }

  Cimaterialfrac = (int *)genvector("Cimaterialfrac", mxsize, sizeof(int));
  Cnextfrac = (int *)genvector("Cnextfrac", mxsize, sizeof(int));
  Cfrac2cell = (int *)genvector("Cfrac2cell", mxsize, sizeof(int));
  CVolfrac = (double *)genvector("CVolfrac", mxsize, sizeof(double));
  CDensityfrac = (double *)genvector("CDensityfrac", mxsize, sizeof(double));
  CTemperaturefrac =
      (double *)genvector("CTemperaturefrac", mxsize, sizeof(double));
  CPressurefrac = (double *)genvector("CPressurefrac", mxsize, sizeof(double));

  // This assumes that Vol, Density, Temperature, Pressure hold the correct
  // values for single material values. If not, they should be copied
  // from CVolfrac, CDensityfrac, CTemperaturefrac, and CPressurefrac

  int ix = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int nmats = nmatscell[ic];
    if (nmats == 1) {
      int m = matids[4 * ic];
      Cimaterial[ic] = m + 1;
    } else {  // nmatscell > 1
      Cimaterial[ic] = -ix;
      for (int im = 0; im < nmats; im++) {
        int m = matids[4 * ic + im];
        int c = mesh2subset[m][ic];
        Cimaterialfrac[ix] = m + 1;
        Cnextfrac[ix] = ix + 1;
        Cfrac2cell[ix] = ic;
        CVolfrac[ix] = Volfrac[m][c];
        CDensityfrac[ix] = Densityfrac[m][c];
        CTemperaturefrac[ix] = Temperaturefrac[m][c];
        CPressurefrac[ix] = Pressurefrac[m][c];
        ix++;
      }
      Cnextfrac[ix - 1] = -1;
    }  // nmatscell > 1
  }
}

void convert_compact_cell_2_compact_material(
    int ncells, int nmats, int *Cimaterial, int *Cnmaterials,
    int *Cimaterialfrac, double *CVol, double *CDensity, double *CTemperature,
    double *CVolfrac, double *CDensityfrac, double *CTemperaturefrac,
    int **&subset2mesh, int **&mesh2subset, int *&nmatscell, int *&matids,
    int *&ncellsmat, double **&Volfrac, double **&Densityfrac,
    double **&Temperaturefrac, double **&Pressurefrac) {
  // Already setup and just needs a name change
  nmatscell = Cnmaterials;

  mesh2subset = (int **)genmatrix("mesh2subset", nmats, ncells, sizeof(int));
  for (int m = 0; m < nmats; m++)
    for (int C = 0; C < ncells; C++) mesh2subset[m][C] = -1;

  matids = (int *)genvector("matids", 4 * ncells, sizeof(int));
  for (int C = 0; C < 4 * ncells; C++) matids[C] = -1;
  ncellsmat = (int *)genvector("ncellsmat", nmats, sizeof(int));
  for (int m = 0; m < nmats; m++) ncellsmat[m] = 0;

  // We need ncellsmat for each material
  int ix = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int CNmats = Cnmaterials[ic];
    if (CNmats == 1) {
      int m = Cimaterial[ic] - 1;
      (ncellsmat[m])++;
    } else {
      for (int im = 0; im < CNmats; im++) {
        int m = Cimaterialfrac[ix] - 1;
        (ncellsmat[m])++;
        ix++;
      }
    }
  }

  subset2mesh = (int **)genvector("subset2mesh", nmats, sizeof(int *));
  Volfrac = (double **)genvector("VolumeFrac", nmats, sizeof(double *));
  Densityfrac = (double **)genvector("DensityFrac", nmats, sizeof(double *));
  Temperaturefrac =
      (double **)genvector("TemperatureFrac", nmats, sizeof(double *));
  Pressurefrac = (double **)genvector("PressureFrac", nmats, sizeof(double *));

  // Allocate compact data structures

  for (int m = 0; m < nmats; m++) {
    subset2mesh[m] =
        (int *)genvector("subset2mesh_m", ncellsmat[m], sizeof(int));
    Volfrac[m] = (double *)genvector("VolFrac_m", ncellsmat[m], sizeof(double));
    Densityfrac[m] =
        (double *)genvector("DensityFrac_m", ncellsmat[m], sizeof(double));
    Temperaturefrac[m] =
        (double *)genvector("TemperatureFrac_m", ncellsmat[m], sizeof(double));
    Pressurefrac[m] =
        (double *)genvector("PressureFrac_m", ncellsmat[m], sizeof(double));
  }

  // Now populate the compact data structures
  for (int m = 0; m < nmats; m++) ncellsmat[m] = 0;
  ix = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int CNmats = Cnmaterials[ic];
    if (CNmats == 1) {
      int m = Cimaterial[ic] - 1;
      int c = ncellsmat[m];
      matids[4 * ic] = m;
      mesh2subset[m][ic] = c;
      subset2mesh[m][c] = ic;
      Volfrac[m][c] = CVol[ic];
      Densityfrac[m][c] = CDensity[ic];
      Temperaturefrac[m][c] = CTemperature[ic];
      (ncellsmat[m])++;
    } else {
      for (int im = 0; im < CNmats; im++) {
        int m = Cimaterialfrac[ix] - 1;
        int c = ncellsmat[m];
        matids[4 * ic + im] = m;
        mesh2subset[m][ic] = c;
        subset2mesh[m][c] = ic;
        Volfrac[m][c] = CVolfrac[ix];
        Densityfrac[m][c] = CDensityfrac[ix];
        Temperaturefrac[m][c] = CTemperaturefrac[ix];
        (ncellsmat[m])++;
        ix++;
      }
    }
  }
}

void print_performance_estimates(float est_perf, int64_t memops8byte,
                                 int64_t memops4byte, int64_t flops,
                                 float penalty_msecs) {
  float Megamemops, Megabytes, Megaflops;
  int STREAM = STREAM_RATE;

  Megamemops = (float)(memops8byte + memops4byte) / 1000000.0;

  // First divide by 1000000 and then multiply by bytes to avoid overflow
  // Using floats to make sure this works on GPUs as well?
  Megabytes = 8 * (((float)memops8byte) / 1000000) +
              4 * (((float)memops4byte) / 1000000.);
  Megaflops = (float)flops / 1000000.;
  printf(
      "Memory Operations  are %.1f M memops, %.1f Mbytes, %.1f Mflops, %0.2f:1 "
      "memops:flops\n",
      Megamemops, Megabytes, Megaflops, (float)Megamemops / (float)Megaflops);
  est_perf = (float)Megabytes / (float)STREAM * 1000.0 + penalty_msecs;
  model_error = (est_perf - act_perf) / act_perf * 100.0;
  printf(
      "Estimated performance %.2f msec, actual %.2f msec, model error %f "
      "%%\n\n",
      est_perf, act_perf, model_error);
}

void get_vol_frac_matrix(int method, double **&Volfrac,
                         float &filled_percentage) {
  if (method == 0)
    get_vol_frac_matrix_rand(Volfrac, filled_percentage);
  else if (method == 1)
    get_vol_frac_matrix_file(Volfrac, filled_percentage);
}

void get_vol_frac_matrix_rand(double **&Volfrac, float &filled_percentage) {
  Volfrac = (double **)genmatrix("VolumeBase", ncells, nmats, sizeof(double));
  int *mf_rand = (int *)genvector("mf_rand", ncells, sizeof(int));

  srand(0);
  for (int ic = 0; ic < ncells; ic++) {
    mf_rand[ic] =
        (int)((float)rand() * 1000.0 / (float)((long long)RAND_MAX + 1));
  }

  for (int ic = 0; ic < ncells; ic++)
    for (int m = 0; m < nmats; m++) Volfrac[ic][m] = 0.0;

  double VolTotal = 0.0;
  int filled_count = 0;
  int mixed_cell_count = 0;
  int mixed_frac_count = 0;
  int pure_frac_count = 0;
  int pure_cell_count = 0;
  int onematcell = 0;
  int twomatcell = 0;
  int threematcell = 0;
  int fourmatcell = 0;
  int fiveplusmatcell = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int m1 =
        (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
    m1 = MININT(m1, nmats - 1);
    Volfrac[ic][m1] = 1.0;
    int mf = mf_rand[ic];
    if (mf < 25) {
      int m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while (m2 == m1) {
        m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      int m3 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while (m3 == m2 || m3 == m1) {
        m3 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      int m4 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while (m4 == m3 || m4 == m2 || m4 == m1) {
        m4 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      m2 = MININT(m2, nmats - 1);
      m3 = MININT(m3, nmats - 1);
      m4 = MININT(m4, nmats - 1);
      Volfrac[ic][m1] = 0.4;
      Volfrac[ic][m2] = 0.3;
      Volfrac[ic][m3] = 0.2;
      Volfrac[ic][m4] = 0.1;
    } else if (mf < 75) {
      int m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while (m2 == m1) {
        m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      int m3 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while (m3 == m2 || m3 == m1) {
        m3 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      m2 = MININT(m2, nmats - 1);
      m3 = MININT(m3, nmats - 1);
      Volfrac[ic][m1] = 0.5;
      Volfrac[ic][m2] = 0.3;
      Volfrac[ic][m3] = 0.2;
    } else if (mf < 200) {
      int m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while (m2 == m1) {
        m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      m2 = MININT(m2, nmats - 1);
      Volfrac[ic][m1] = 0.5;
      Volfrac[ic][m2] = 0.5;
    }
    int mat_count = 0;
    for (int m = 0; m < nmats; m++) {
      if (Volfrac[ic][m] > 0.0) {
        filled_count++;
        mat_count++;
      }
      VolTotal += Volfrac[ic][m];
    }
    if (mat_count >= 2) {
      mixed_cell_count++;
      mixed_frac_count += mat_count;
    } else {
      pure_frac_count++;
    }
    if (mat_count == 1) pure_cell_count++;
    if (mat_count == 1) onematcell++;
    if (mat_count == 2) twomatcell++;
    if (mat_count == 3) threematcell++;
    if (mat_count == 4) fourmatcell++;
    if (mat_count >= 5) fiveplusmatcell++;
  }

  genvectorfree(mf_rand);

  printf("Ratios to Full Data Structure\n");
  filled_percentage = (float)filled_count * 100.0 / (float)(ncells * nmats);
  float sparsity_percentage =
      (float)(ncells * nmats - filled_count) * 100.0 / (float)(ncells * nmats);
  printf("Sparsity %lf percent/Filled %lf percent\n\n", sparsity_percentage,
         filled_percentage);

  printf("Ratios to Number of Cells\n");
  float pure_cell_percentage = (float)pure_cell_count * 100.0 / (float)ncells;
  float mixed_cell_percentage = (float)mixed_cell_count * 100.0 / (float)ncells;
  printf("Pure cell %lf percentage/Mixed material %lf percentage\n\n",
         pure_cell_percentage, mixed_cell_percentage);

  printf("Ratios to Mixed Material Data Structure\n");
  float mixed_material_sparsity_percentage =
      (float)mixed_frac_count * 100.0 / (float)(mixed_cell_count * nmats);
  float mixed_material_filled_percentage =
      (float)(mixed_cell_count * nmats - mixed_frac_count) * 100.0 /
      (float)(mixed_cell_count * nmats);
  printf(
      "Mixed material Sparsity %lf percent/Mixed material Filled %lf "
      "percent\n\n",
      mixed_material_sparsity_percentage, mixed_material_filled_percentage);

  // printf("Vol Total %lf\n",VolTotal);
  // printf("%f percent of the cells are
  // filled\n",(float)filled_count*100.0/(float)(ncells*nmats)); printf("%f
  // percent of the cells are
  // mixed\n",(float)mixed_cell_count*100.0/(float)ncells); printf("%f percent of
  // the total are mixed\n",(float)mixed_frac_count*100.0/(float)(ncells*nmats));
  // printf("%f percent of the frac are
  // mixed\n",(float)mixed_frac_count*100.0/(float)(mixed_cell_count*nmats));
  // printf("%f percent
  // sparsity\n",(float)(ncells*nmats-mixed_frac_count)*100.0/(float)(ncells*nmats));
  // printf("%f percent of the frac are
  // pure\n",(float)pure_frac_count*100.0/(float)ncells);
  printf("1 matcell %d 2 matcell %d 3 matcell %d 4 matcell %d 5 matcell %d\n\n",
         onematcell, twomatcell, threematcell, fourmatcell, fiveplusmatcell);
  // printf("Total cells %d\n\n",
  // onematcell+2*twomatcell+3*threematcell+4*fourmatcell+5*fiveplusmatcell);
}

void get_vol_frac_matrix_file(double **&Volfrac, float &filled_percentage) {
  int status;
  FILE *fp;
  fp = fopen("volfrac.dat", "r");
  if (!fp) {
    fprintf(stderr, "unable to read volume fractions from file \"%s\"\n",
            "volfrac.dat");
    exit(-1);
  }

  status = fscanf(fp, "%d", &nmats);
  if (status < 0) {
    printf("error in read at line %d\n", __LINE__);
    exit(1);
  }
  Volfrac = (double **)genmatrix("VolumeBase", ncells, nmats, sizeof(double));

  for (int ic = 0; ic < ncells; ic++)
    for (int m = 0; m < nmats; m++) Volfrac[ic][m] = 0.0;

  char matname[256];
  for (int m = 0; m < nmats; m++) {
    status = fscanf(fp, "%s", matname);  // read and discard
    if (status < 0) {
      printf("error in read at line %d\n", __LINE__);
      exit(1);
    }
  }

  double VolTotal = 0.0;
  int filled_count = 0;
  int mixed_cell_count = 0;
  int mixed_frac_count = 0;
  int pure_frac_count = 0;
  int pure_cell_count = 0;
  int onematcell = 0;
  int twomatcell = 0;
  int threematcell = 0;
  int fourmatcell = 0;
  int fiveplusmatcell = 0;
  for (int ic = 0; ic < ncells; ic++) {
    int mat_count = 0;
    for (int m = 0; m < nmats; m++) {
      status = fscanf(fp, "%lf", &(Volfrac[ic][m]));
      if (status < 0) {
        printf("error in read at line %d\n", __LINE__);
        exit(1);
      }
      if (Volfrac[ic][m] > 0.0) {
        filled_count++;
        mat_count++;
      }
      VolTotal += Volfrac[ic][m];
    }
    if (mat_count >= 2) {
      mixed_cell_count++;
      mixed_frac_count += mat_count;
    } else {
      pure_frac_count++;
    }
    if (mat_count == 1) pure_cell_count++;
    if (mat_count == 1) onematcell++;
    if (mat_count == 2) twomatcell++;
    if (mat_count == 3) threematcell++;
    if (mat_count == 4) fourmatcell++;
    if (mat_count >= 5) fiveplusmatcell++;
  }
  fclose(fp);

  printf("Ratios to Full Data Structure\n");
  filled_percentage = (float)filled_count * 100.0 / (float)(ncells * nmats);
  float sparsity_percentage =
      (float)(ncells * nmats - filled_count) * 100.0 / (float)(ncells * nmats);
  printf("Sparsity %lf percent/Filled %lf percent\n\n", sparsity_percentage,
         filled_percentage);

  printf("Ratios to Number of Cells\n");
  float pure_cell_percentage = (float)pure_cell_count * 100.0 / (float)ncells;
  float mixed_cell_percentage = (float)mixed_cell_count * 100.0 / (float)ncells;
  printf("Pure cell %lf percentage/Mixed material %lf percentage\n\n",
         pure_cell_percentage, mixed_cell_percentage);

  printf("Ratios to Mixed Material Data Structure\n");
  float mixed_material_sparsity_percentage =
      (float)mixed_frac_count * 100.0 / (float)(mixed_cell_count * nmats);
  float mixed_material_filled_percentage =
      (float)(mixed_cell_count * nmats - mixed_frac_count) * 100.0 /
      (float)(mixed_cell_count * nmats);
  printf(
      "Mixed material Sparsity %lf percent/Mixed material Filled %lf "
      "percent\n\n",
      mixed_material_sparsity_percentage, mixed_material_filled_percentage);

  printf("Vol Total %lf\n", VolTotal);
  printf("%f percent of the cells are filled\n",
         (float)filled_count * 100.0 / (float)(ncells * nmats));
  printf("%f percent of the cells are mixed\n",
         (float)mixed_cell_count * 100.0 / (float)ncells);
  printf("%f percent of the total are mixed\n",
         (float)mixed_frac_count * 100.0 / (float)(ncells * nmats));
  printf("%f percent of the frac are mixed\n",
         (float)mixed_frac_count * 100.0 / (float)(mixed_cell_count * nmats));
  printf("%f percent sparsity\n", (float)(ncells * nmats - mixed_frac_count) *
                                      100.0 / (float)(ncells * nmats));
  printf("%f percent of the frac are pure\n",
         (float)pure_frac_count * 100.0 / (float)ncells);
  printf("1 matcell %d 2 matcell %d 3 matcell %d 4 matcell %d 5 matcell %d\n\n",
         onematcell, twomatcell, threematcell, fourmatcell, fiveplusmatcell);
  printf("Total cells %d\n\n", onematcell + 2 * twomatcell + 3 * threematcell +
                                   4 * fourmatcell + 5 * fiveplusmatcell);
}

void get_neighbors(int ncells, int *&num_nbrs, int **&nbrs) {
  int ncells1 = (int)sqrt(ncells);  // assumes ncells is a perfect square
  if (ncells1 * ncells1 != ncells) {
    fprintf(stderr, "Number of cells in mesh is not a perfect square");
    exit(-1);
  }

  for (int i = 0; i < ncells1; i++) {
    for (int j = 0; j < ncells1; j++) {
      int c = i * ncells1 + j;
      int ilo = i == 0 ? i : i - 1;
      int jlo = j == 0 ? j : j - 1;
      int ihi = i == ncells1 - 1 ? i : i + 1;
      int jhi = j == ncells1 - 1 ? j : j + 1;
      int n = 0;
      for (int i1 = ilo; i1 <= ihi; i1++)
        for (int j1 = jlo; j1 <= jhi; j1++) {
          int c2 = i1 * ncells1 + j1;
          if (c2 != c) {
            nbrs[c][n] = i1 * ncells1 + j1;
            n++;
          }
        }
      num_nbrs[c] = n;
    }
  }
}

void get_centroids(double (*&cen)[2]) {
  int ncells1 = (int)sqrt(ncells);  // assumes ncells is a perfect square
  if (ncells1 * ncells1 != ncells) {
    fprintf(stderr, "Number of cells in mesh is not a perfect square");
    exit(-1);
  }

  // Assume domain is a unit square

  double XLO = 0.0, YLO = 0.0, XHI = 1.0, YHI = 1.0;
  double dx = (XHI - XLO) / ncells1, dy = (YHI - YLO) / ncells1;

  for (int i = 0; i < ncells1; i++) {
    for (int j = 0; j < ncells1; j++) {
      int c = i * ncells1 + j;
      cen[c][0] = XLO + i * dx;
      cen[c][1] = YLO + j * dy;
    }
  }
}
