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

void reb_saba_corrector_step(struct reb_simulation* r, double dt){
	struct reb_particle* const particles = r->particles;
	const int N = r->N;
	const double G = r->G;
    for (int i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
    for (int i=0; i<N; i++){
        if (reb_sigint) return;
        for (int j=0; j<N; j++){
            if (((j==1 && i==0) || (i==1 && j==0) )) continue;
            if (i==j) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
            const double prefact = -G/(_r*_r*_r)*particles[j].m;
            
            particles[i].ax    += prefact*dx;
            particles[i].ay    += prefact*dy;
            particles[i].az    += prefact*dz;
        }
    }
    const double m0 = particles[0].m;
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* const p_j = ri_whfast->p_jh;
    reb_transformations_inertial_to_jacobi_acc(r->particles, p_j, r->particles, N);
    double eta = m0;
    for (unsigned int i=1;i<N;i++){
        // Eq 132
        const struct reb_particle pji = p_j[i];
        eta += pji.m;
        static double rj2i;
        static double rj3iM;
        static double prefac1;
        p_j[i].vx += dt * pji.ax;
        p_j[i].vy += dt * pji.ay;
        p_j[i].vz += dt * pji.az;
        if (i>1){
            rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z);
            const double rji  = sqrt(rj2i);
            rj3iM = rji*rj2i*G*eta;
            prefac1 = dt*rj3iM;
            p_j[i].vx += prefac1*pji.x;
            p_j[i].vy += prefac1*pji.y;
            p_j[i].vz += prefac1*pji.z;
        }
    }
}

void reb_integrator_saba_part1(struct reb_simulation* const r){
    // If the total step consistes of   AB(BA)^(k-1)[A]
    // do only                          A                  here in part1, 
    // and                               B(AB)^(k-1)[A]    in part2.
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    const int k = ri_saba->k;
    if (reb_integrator_whfast_init(r)){
        // Non recoverable error occured.
        return;
    }
    
    // Only recalculate Jacobi coordinates if needed
    if (ri_whfast->safe_mode || ri_whfast->recalculate_coordinates_this_timestep){
        if (ri_whfast->is_synchronized==0){
            reb_integrator_saba_synchronize(r);
            if (ri_whfast->recalculate_coordinates_but_not_synchronized_warning==0){
                reb_warning(r,"Recalculating coordinates but pos/vel were not synchronized before.");
                ri_whfast->recalculate_coordinates_but_not_synchronized_warning++;
            }
        }
        reb_integrator_whfast_from_inertial(r);
        ri_whfast->recalculate_coordinates_this_timestep = 0;
    }
    if (ri_whfast->is_synchronized){
        // First half DRIFT step
        //if (ri_whfast->corrector){
        //    reb_whfast_apply_corrector(r, 1., ri_whfast->corrector, reb_whfast_corrector_Z);
        //}
        reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);   
        reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
    }else{
        // Combined DRIFT step
        reb_whfast_kepler_step(r, 2.*reb_saba_c[k-1][0]*r->dt);   
        reb_whfast_com_step(r, 2.*reb_saba_c[k-1][0]*r->dt);
    }

    reb_integrator_whfast_to_inertial(r);

}

void reb_integrator_saba_synchronize(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    const int k = ri_saba->k;
    if (ri_whfast->is_synchronized == 0){
        const int N_real = r->N-r->N_var;
        struct reb_particle* sync_pj  = NULL;
        if (ri_whfast->keep_unsynchronized){
            sync_pj = malloc(sizeof(struct reb_particle)*r->N);
            memcpy(sync_pj,r->ri_whfast.p_jh,r->N*sizeof(struct reb_particle));
        }
        reb_whfast_kepler_step(r, reb_saba_c[k-1][0]*r->dt);
        reb_whfast_com_step(r, reb_saba_c[k-1][0]*r->dt);
        //if (ri_whfast->corrector){
        //    reb_whfast_apply_corrector(r, -1., ri_whfast->corrector, reb_whfast_corrector_Z);
        //}
        reb_transformations_jacobi_to_inertial_posvel(r->particles, ri_whfast->p_jh, r->particles, N_real);
        if (ri_whfast->keep_unsynchronized){
            memcpy(r->ri_whfast.p_jh,sync_pj,r->N*sizeof(struct reb_particle));
            free(sync_pj);
        }else{
            ri_whfast->is_synchronized = 1;
        }
    }
}

void reb_integrator_saba_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    struct reb_particle* restrict const particles = r->particles;
    const int k = ri_saba->k;
    const int N_real = r->N-r->N_var;
    if (ri_whfast->p_jh==NULL){
        // Non recoverable error occured earlier. 
        // Skipping rest of integration to avoid segmentation fault.
        return;
    }
    
    reb_whfast_interaction_step(r, reb_saba_d[k-1][0]*r->dt);
  
    for(int i=1;i<k;i++){
        reb_whfast_kepler_step(r, reb_saba_c[k-1][i]*r->dt);   
        reb_whfast_com_step(r, reb_saba_c[k-1][i]*r->dt);
        reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N_real);
        r->gravity_ignore_terms = 1;
        reb_update_acceleration(r);
        reb_whfast_interaction_step(r, reb_saba_d[k-1][i]*r->dt);
    } 
     

    ri_whfast->is_synchronized = 0;
    if (ri_whfast->safe_mode){
        reb_integrator_saba_synchronize(r);
    }
    
    r->t+=r->dt;
    r->dt_last_done = r->dt;

}
    
void reb_integrator_saba_reset(struct reb_simulation* const r){
    struct reb_simulation_integrator_saba* const ri_saba = &(r->ri_saba);
    ri_saba->k = 1;
    reb_integrator_whfast_reset(r);
}
