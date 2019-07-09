/**
 * @file    integrator_saba.c
 * @brief   SABA integrator family (Laskar and Robutel 2001).
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the family of symplectic integrators
 * of Laskar and Robutel (2001). 
 * 
 * @section LICENSE
 * Copyright (c) 2019 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "rebound.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_saba.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

// Some coefficients appear multiple times to simplify the for loops below. 
const static double reb_saba_c[4][4] = {
        {0.5, 0., 0., 0.}, // SABA1
        {0.2113248654051871177454256097490212721762, 0.5773502691896257645091487805019574556476, 0., 0.}, // SABA2
        {0.1127016653792583114820734600217600389167, 0.3872983346207416885179265399782399610833,  0.3872983346207416885179265399782399610833, 0.}, // SABA3
        {0.06943184420297371238802675555359524745214, 0.2605776340045981552106403648947824089476, 0.3399810435848562648026657591032446872006, 0.2605776340045981552106403648947824089476}, // SABA4
}; 
const static double reb_saba_d[4][4] = {
        {1., 0., 0., 0.},
        {0.5, 0.5, 0., 0.},
        {0.2777777777777777777777777777777777777778, 0.4444444444444444444444444444444444444444,0.2777777777777777777777777777777777777778, 0.},
        {0.1739274225687269286865319746109997036177, 0.3260725774312730713134680253890002963823, 0.3260725774312730713134680253890002963823, 0.1739274225687269286865319746109997036177},
};
const static double reb_saba_cc[4] = {
        0.08333333333333333333333333333333333333333, // SABA1
        0.01116454968463011276968973577058865137738, // SABA2
        0.005634593363122809402267823769797538671562, // SABA3
        0.003396775048208601331532157783492144, // SABA4
}; 
    

void reb_saba_corrector_step(struct reb_simulation* r, double cc){
    double dt = r->dt;
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* const p_j = ri_whfast->p_jh;
	struct reb_particle* const particles = r->particles;
	const int N = r->N;
	const double G = r->G;

    // Calculate normal kick
    reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N);
    r->gravity_ignore_terms = 1;
    reb_update_acceleration(r);
    reb_transformations_inertial_to_jacobi_acc(r->particles, p_j, r->particles, N);

    // make copy of normal kick, also stores original positions
    struct reb_particle* temp_pj = r->ri_saba.temp_pj;
    memcpy(temp_pj,p_j,r->N*sizeof(struct reb_particle));

    double eta = particles[0].m;
    for (unsigned int i=1;i<N;i++){
        const struct reb_particle pji = p_j[i];
        eta += pji.m;
        const double prefac1 = dt*dt/12.; 
        if (i>1){
            const double rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z);
            const double rji  = sqrt(rj2i);
            const double rj3iM = rji*rj2i*G*eta;
            temp_pj[i].ax += rj3iM*temp_pj[i].x;
            temp_pj[i].ay += rj3iM*temp_pj[i].y;
            temp_pj[i].az += rj3iM*temp_pj[i].z;
        }
        p_j[i].x += prefac1 * temp_pj[i].ax;
        p_j[i].y += prefac1 * temp_pj[i].ay;
        p_j[i].z += prefac1 * temp_pj[i].az;
    }
   
    // recalculate kick 
    reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N);
    reb_update_acceleration(r);
    reb_transformations_inertial_to_jacobi_acc(r->particles, p_j, r->particles, N);

    dt = 12.*dt*cc;
    eta = particles[0].m;
    for (unsigned int i=1;i<N;i++){
        const struct reb_particle pji = p_j[i];
        eta += pji.m;
        if (i>1){
            const double rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z);
            const double rji  = sqrt(rj2i);
            const double rj3iM = rji*rj2i*G*eta;
            p_j[i].ax += rj3iM*pji.x;
            p_j[i].ay += rj3iM*pji.y;
            p_j[i].az += rj3iM*pji.z;
        }
        // commutator is difference between modified and original kick
        p_j[i].vx += dt * (p_j[i].ax - temp_pj[i].ax);
        p_j[i].vy += dt * (p_j[i].ay - temp_pj[i].ay);
        p_j[i].vz += dt * (p_j[i].az - temp_pj[i].az);
        // reset positions
        p_j[i].x = temp_pj[i].x;
        p_j[i].y = temp_pj[i].y;
        p_j[i].z = temp_pj[i].z;
    }

}

void reb_integrator_saba_part1(struct reb_simulation* const r){
    // If the total step consistes of   AB(BA)^(k-1)[A]
    // do only                          A                  here in part1, 
    // and                               B(AB)^(k-1)[A]    in part2.
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    const int k = ri_saba->k;
    const int corrector = ri_saba->corrector;
    const int N = r->N;
    if (reb_integrator_whfast_init(r)){
        // Non recoverable error occured.
        return;
    }
    if (corrector && ri_saba->allocated_N != N){
        ri_saba->allocated_N = N;
        ri_saba->temp_pj = realloc(ri_saba->temp_pj,sizeof(struct reb_particle)*N);
    }
    
    // Only recalculate Jacobi coordinates if needed
    if (ri_saba->safe_mode || ri_whfast->recalculate_coordinates_this_timestep){
        reb_integrator_whfast_from_inertial(r);
        ri_whfast->recalculate_coordinates_this_timestep = 0;
    }
    if (corrector){
        if (ri_saba->is_synchronized){
            reb_saba_corrector_step(r, reb_saba_cc[k-1]);
        }else{
            reb_saba_corrector_step(r, 2.*reb_saba_cc[k-1]);
        }
        // First half DRIFT step
        reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);   
        reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
    }else{
        if (ri_saba->is_synchronized){
            // First half DRIFT step
            reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);   
            reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
        }else{
            // Combined DRIFT step
            reb_whfast_kepler_step(r, 2.*reb_saba_c[k-1][0]*r->dt);   
            reb_whfast_com_step(r, 2.*reb_saba_c[k-1][0]*r->dt);
        }
    }

    reb_integrator_whfast_to_inertial(r);
}

void reb_integrator_saba_synchronize(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    int k = ri_saba->k;
    int corrector = ri_saba->corrector;
    if (ri_saba->is_synchronized == 0){
        const int N = r->N;
        struct reb_particle* sync_pj  = NULL;
        if (ri_whfast->keep_unsynchronized){
            sync_pj = malloc(sizeof(struct reb_particle)*r->N);
            memcpy(sync_pj,r->ri_whfast.p_jh,r->N*sizeof(struct reb_particle));
        }
        if (corrector){
            // DRIFT ALREADY DONE
            reb_saba_corrector_step(r, reb_saba_cc[k-1]);
        }else{
            reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);
            reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
        }
        reb_transformations_jacobi_to_inertial_posvel(r->particles, ri_whfast->p_jh, r->particles, N);
        if (ri_whfast->keep_unsynchronized){
            memcpy(r->ri_whfast.p_jh,sync_pj,r->N*sizeof(struct reb_particle));
            free(sync_pj);
        }else{
            ri_saba->is_synchronized = 1;
        }
    }
}

void reb_integrator_saba_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    struct reb_particle* restrict const particles = r->particles;
    int k = ri_saba->k;
    int corrector = ri_saba->corrector;
    const int N = r->N;
    if (ri_whfast->p_jh==NULL){
        // Non recoverable error occured earlier. 
        // Skipping rest of integration to avoid segmentation fault.
        return;
    }
    
    reb_whfast_interaction_step(r, reb_saba_d[k-1][0]*r->dt);
  
    for(int i=1;i<k;i++){
        reb_whfast_kepler_step(r, reb_saba_c[k-1][i]*r->dt);   
        reb_whfast_com_step(r, reb_saba_c[k-1][i]*r->dt);
        reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N);
        r->gravity_ignore_terms = 1;
        reb_update_acceleration(r);
        reb_whfast_interaction_step(r, reb_saba_d[k-1][i]*r->dt);
    } 

    if (corrector){
        // Always need to do DRIFT step if correctors are turned on
        reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);
        reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
    }

    ri_saba->is_synchronized = 0;
    if (ri_saba->safe_mode){
        reb_integrator_saba_synchronize(r);
    }
    
    r->t+=r->dt;
    r->dt_last_done = r->dt;
}
    
void reb_integrator_saba_reset(struct reb_simulation* const r){
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    ri_saba->k = 1;
    ri_saba->corrector = 0;
    ri_saba->safe_mode = 1;
    ri_saba->is_synchronized = 1;
    reb_integrator_whfast_reset(r);
    if (ri_saba->temp_pj){
        free(ri_saba->temp_pj);
        ri_saba->temp_pj = NULL;
    }
    ri_saba->allocated_N = 0;
}
