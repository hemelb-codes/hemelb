/* This file contains the globals.
 * Should be removed ASAP
 */

double *f_old = 0, *f_new = 0;
int *f_id = 0;
int *f_recv_iv = 0;

int sites_x, sites_y, sites_z;
int blocks_x, blocks_y, blocks_z;
int blocks;
int block_size;
int shift;
int sites_in_a_block;

double lbm_stress_type;
double lbm_stress_par;

int lbm_terminate_simulation;
