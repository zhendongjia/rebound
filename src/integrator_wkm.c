/**
 * @file 	integrator_wkm.c
 * @brief   Wisdom Kernel Method 
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the Kernel method of Wisdom, 
 * Holman, and Touma (1996)
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
#include "integrator_wkm.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

void reb_wkm_jerk_step(struct reb_simulation* r){
    // Assume particles.a calculated.
	struct reb_particle* const particles = r->particles;
	const int N = r->N;
    struct reb_particle* jerk = malloc(sizeof(struct reb_particle)*N);
	const double G = r->G;
    for (int i=0; i<N; i++){
        jerk[i].ax = 0; 
        jerk[i].ay = 0; 
        jerk[i].az = 0; 
    }
    /////////////////
    // A+B Term
    for (int j=0; j<N; j++){
        for (int i=0; i<N; i++){
            double Rkm1x = particles[0].x;
            double Rkm1y = particles[0].y;
            double Rkm1z = particles[0].z;
            double Mkm1 = particles[0].m;
            for (int k=1; k<N; k++){
                double Qkx = particles[k].x - Rkm1x;
                double Qky = particles[k].y - Rkm1y;
                double Qkz = particles[k].z - Rkm1z;
                
                const double dr = sqrt(Qkx*Qkx + Qky*Qky + Qkz*Qkz);
                double alphasum = particles[i].ax*Qkx+particles[i].ay*Qky+particles[i].az*Qkz;
                double dQkrj = 1.;
                if (j>k){
                    dQkrj = 0.;
                }
                if (j<k){
                    dQkrj = -particles[j].m/Mkm1;
                }
                double dQkri = 1.;
                if (i>k){
                    dQkri = 0.;
                }
                if (i<k){
                    dQkri = -particles[i].m/Mkm1;
                }
                
                double prefact1 = G*Mkm1*particles[k].m *dQkrj * dQkri /(dr*dr*dr*dr*dr);
                jerk[j].ax    += 6.* prefact1*alphasum*Qkx; 
                jerk[j].ay    += 6.* prefact1*alphasum*Qky;
                jerk[j].az    += 6.* prefact1*alphasum*Qkz; 
                
                double prefact2 = G*Mkm1*particles[k].m *dQkrj * dQkri /(dr*dr*dr);
                jerk[j].ax    -= 2.*prefact2*particles[i].ax;
                jerk[j].ay    -= 2.*prefact2*particles[i].ay;
                jerk[j].az    -= 2.*prefact2*particles[i].az;
                
                Rkm1x = (Rkm1x*Mkm1+particles[k].x*particles[k].m)/(Mkm1+particles[k].m);
                Rkm1y = (Rkm1y*Mkm1+particles[k].y*particles[k].m)/(Mkm1+particles[k].m);
                Rkm1z = (Rkm1z*Mkm1+particles[k].z*particles[k].m)/(Mkm1+particles[k].m);
                Mkm1 += particles[k].m;
            }
        }
    }
    /////////////////
    // C Term
    for (int j=0; j<N; j++){
        for (int i=0; i<N; i++){
            if (i!=j){
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                double prefact1 = G*particles[j].m*particles[i].m /(dr*dr*dr*dr*dr);
                double alphasum = dax*dx+day*dy+daz*dz;
                jerk[j].ax    += 6.*alphasum*dx*prefact1;
                jerk[j].ay    += 6.*alphasum*dy*prefact1;
                jerk[j].az    += 6.*alphasum*dz*prefact1;
                double prefact2 = G*particles[j].m*particles[i].m /(dr*dr*dr);
                jerk[j].ax    -= 2.*dax*prefact2;
                jerk[j].ay    -= 2.*day*prefact2;
                jerk[j].az    -= 2.*daz*prefact2;
            }
        }
    }

    ///////////////////
    for (int i=0; i<N; i++){
        particles[i].ax -= 0.1*r->dt*r->dt/24.*jerk[i].ax/particles[i].m; 
        particles[i].ay -= 0.1*r->dt*r->dt/24.*jerk[i].ay/particles[i].m; 
        particles[i].az -= 0.1*r->dt*r->dt/24.*jerk[i].az/particles[i].m; 
    }
    free(jerk);
}


void reb_wkm_interaction_step(struct reb_simulation* const r, const double _dt){
    const unsigned int N_real = r->N-r->N_var;
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* const p_j = ri_whfast->p_jh;
    reb_transformations_inertial_to_jacobi_acc(r->particles, p_j, r->particles, N_real);
    for (unsigned int i=0;i<N_real;i++){
        const struct reb_particle pji = p_j[i];
        p_j[i].vx += _dt * pji.ax;
        p_j[i].vy += _dt * pji.ay;
        p_j[i].vz += _dt * pji.az;
    }
}

static void reb_wkm_corrector_Z(struct reb_simulation* r, const double a, const double b){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* restrict const particles = r->particles;
    const int N_real = r->N;
    reb_whfast_kepler_step(r, a);
    reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N_real);
    reb_update_acceleration(r);
    reb_wkm_interaction_step(r, -b);
    reb_whfast_kepler_step(r, -2.*a);
    reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N_real);
    reb_update_acceleration(r);
    reb_wkm_interaction_step(r, b);
    reb_whfast_kepler_step(r, a);
}
void reb_wkm_apply_C(struct reb_simulation* const r, double a, double b){
    reb_whfast_kepler_step(r, a);   
    reb_whfast_com_step(r, a);
    
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* restrict const particles = r->particles;
    const int N_real = r->N;
    reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N_real);
    reb_update_acceleration(r);
    reb_wkm_interaction_step(r, b);
    
    reb_whfast_kepler_step(r, -a);   
    reb_whfast_com_step(r, -a);
}

void reb_wkm_apply_Y(struct reb_simulation* const r, double a, double b){
    reb_wkm_apply_C(r, a, b); 
    reb_wkm_apply_C(r, -a, -b); 
}
void reb_wkm_apply_U(struct reb_simulation* const r, double a, double b){
    reb_whfast_kepler_step(r, a);   
    reb_whfast_com_step(r, a);
    reb_wkm_apply_Y(r, a, b); 
    reb_wkm_apply_Y(r, a, -b); 
    reb_whfast_kepler_step(r, -a);   
    reb_whfast_com_step(r, -a);
}
void reb_wkm_apply_corrector2(struct reb_simulation* const r, double dt){
    double a = 0.5 * dt;
    double b = 0.03486083443891981449909050107438281205803 * dt;
    reb_wkm_apply_U(r, a, b); 
    reb_wkm_apply_U(r, -a, b);
}

void reb_integrator_wkm_part1(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_wkm* const ri_wkm = &(r->ri_wkm);
    const int corrector = ri_wkm->corrector%10;
    const int kernel = ri_wkm->corrector/10;
    r->gravity = REB_GRAVITY_JACOBI;
    if (r->var_config_N>0 && ri_whfast->coordinates!=REB_WHFAST_COORDINATES_JACOBI){
        reb_error(r, "Variational particles are not supported in the WKM integrator.");
        return; 
    }
    if (ri_whfast->coordinates!=REB_WHFAST_COORDINATES_JACOBI){
        reb_error(r, "WKM integrator requires ri_whfast.coordinates to be set to Jacobi coordinates.");
        return; 
    }
    if (kernel>2){
        reb_error(r, "WKM Kernel not implemented");
        return;
    }
    if (reb_integrator_whfast_init(r)){
        // Non recoverable error occured.
        return;
    }
    
    // Only recalculate Jacobi coordinates if needed
    if (ri_wkm->safe_mode || ri_whfast->recalculate_coordinates_this_timestep){
        reb_integrator_whfast_from_inertial(r);
        ri_whfast->recalculate_coordinates_this_timestep = 0;
    }
    if (ri_wkm->is_synchronized){
        if (corrector){
            reb_whfast_apply_corrector(r, 1., 11, reb_wkm_corrector_Z);
            if (corrector>=2){
                reb_wkm_apply_corrector2(r, r->dt);
            }
        }
        if (kernel==0){ // composition
            reb_whfast_kepler_step(r, 5./8.*r->dt);   
            reb_whfast_com_step(r, 5./8.*r->dt);
        }else if(kernel==1){ // lazy implementers method
            reb_whfast_kepler_step(r, 1./2.*r->dt);   
            reb_whfast_com_step(r, 1./2.*r->dt);
        }else if(kernel==2){ // jerk method
            reb_whfast_kepler_step(r, 1./2.*r->dt);   
            reb_whfast_com_step(r, 1./2.*r->dt);
        }
    }else{
        reb_whfast_kepler_step(r, r->dt);   
        reb_whfast_com_step(r, r->dt);
    }

    reb_integrator_whfast_to_inertial(r);
}

void reb_integrator_wkm_synchronize(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_wkm* const ri_wkm = &(r->ri_wkm);
    const int corrector = ri_wkm->corrector%10;
    const int kernel = ri_wkm->corrector/10;
    if (ri_wkm->is_synchronized == 0){
        const int N = r->N;
        struct reb_particle* sync_pj  = NULL;
        r->gravity = REB_GRAVITY_JACOBI; // needed here in case of SimulationArchive calls
        if (ri_whfast->keep_unsynchronized){
            sync_pj = malloc(sizeof(struct reb_particle)*N);
            memcpy(sync_pj,r->ri_whfast.p_jh,N*sizeof(struct reb_particle));
        }
        if (kernel==0){ // composition
            reb_whfast_kepler_step(r, 3./8.*r->dt);
            reb_whfast_com_step(r, 3./8.*r->dt);
        }else if (kernel==1){ // lazy implementers method
            reb_whfast_kepler_step(r, 1./2.*r->dt);
            reb_whfast_com_step(r, 1./2.*r->dt);
        }else if (kernel==2){ // jerk method
            reb_whfast_kepler_step(r, 1./2.*r->dt);
            reb_whfast_com_step(r, 1./2.*r->dt);
        }

        if (corrector){
            if (corrector>=2){
                reb_wkm_apply_corrector2(r, -r->dt);
            }
            reb_whfast_apply_corrector(r, -1., 11, reb_wkm_corrector_Z);
        }
        reb_transformations_jacobi_to_inertial_posvel(r->particles, ri_whfast->p_jh, r->particles, N);
        if (ri_whfast->keep_unsynchronized){
            memcpy(r->ri_whfast.p_jh,sync_pj,N*sizeof(struct reb_particle));
            free(sync_pj);
        }else{
            ri_wkm->is_synchronized = 1;
        }
    }
}

void reb_integrator_wkm_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_simulation_integrator_wkm* const ri_wkm = &(r->ri_wkm);
    struct reb_particle* restrict const particles = r->particles;
    const int kernel = ri_wkm->corrector/10;
    const int N = r->N;
    if (ri_whfast->p_jh==NULL){
        // Non recoverable error occured earlier. 
        // Skipping rest of integration to avoid segmentation fault.
        return;
    }
    if (kernel==0){ // composition
        
        // -1/6 B
        reb_wkm_interaction_step(r, -1./6.*r->dt);
      
        { // -1/4 A
            reb_whfast_kepler_step(r, -1./4.*r->dt);   
            reb_whfast_com_step(r, -1./4.*r->dt);
        }
        { // 1/6 B
            reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N);
            reb_update_acceleration(r);
            reb_wkm_interaction_step(r, 1./6.*r->dt);
        }
        { // 1/8 A
            reb_whfast_kepler_step(r, 1./8.*r->dt);   
            reb_whfast_com_step(r, 1./8.*r->dt);
        }
        { // B
            reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N);
            reb_update_acceleration(r);
            reb_wkm_interaction_step(r, r->dt);
        }
        { // -1/8 A
            reb_whfast_kepler_step(r, -1./8.*r->dt);   
            reb_whfast_com_step(r, -1./8.*r->dt);
        }
        { // -1/6 B
            reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N);
            reb_update_acceleration(r);
            reb_wkm_interaction_step(r, -1./6.*r->dt);
        }
        { // 1/4 A
            reb_whfast_kepler_step(r, 1./4.*r->dt);   
            reb_whfast_com_step(r, 1./4.*r->dt);
        }
        { // 1/6 B
            reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N);
            reb_update_acceleration(r);
            reb_wkm_interaction_step(r, 1./6.*r->dt);
        }
    } else if (kernel==1){ //lazy implementers method
        double dt = r->dt;
        struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
        struct reb_particle* const p_j = ri_whfast->p_jh;
        struct reb_particle* const particles = r->particles;
        struct reb_simulation_integrator_wkm* const ri_wkm = &(r->ri_wkm);
        const int N = r->N;
        if (ri_wkm->allocated_N != N){
            ri_wkm->allocated_N = N;
            ri_wkm->temp_pj = realloc(ri_wkm->temp_pj,sizeof(struct reb_particle)*N);
        }
        struct reb_particle* temp_pj = ri_wkm->temp_pj;

        // Calculate normal kick
        // Accelertions were already calculated before part2 gets called
        reb_transformations_inertial_to_jacobi_acc(r->particles, p_j, r->particles, N);

        // make copy of normal kick, also stores original positions
        memcpy(temp_pj,p_j,r->N*sizeof(struct reb_particle));

        // WHT Eq 10.6
        for (unsigned int i=1;i<N;i++){
            const double prefac1 = dt*dt/12.; 
            p_j[i].x += prefac1 * temp_pj[i].ax;
            p_j[i].y += prefac1 * temp_pj[i].ay;
            p_j[i].z += prefac1 * temp_pj[i].az;
        }
       
        // recalculate kick 
        reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N);
        reb_update_acceleration(r);
        reb_wkm_interaction_step(r, dt);

        for (unsigned int i=1;i<N;i++){
            // reset positions
            p_j[i].x = temp_pj[i].x;
            p_j[i].y = temp_pj[i].y;
            p_j[i].z = temp_pj[i].z;
        }
    } else if (kernel==2){ //jerk method
        reb_wkm_jerk_step(r); // adds jerk terms to acceleration
        reb_wkm_interaction_step(r, r->dt);
    }
    
    ri_wkm->is_synchronized = 0;
    if (ri_wkm->safe_mode){
        reb_integrator_wkm_synchronize(r);
    }
    
    r->t+=r->dt;
    r->dt_last_done = r->dt;
}
    
void reb_integrator_wkm_reset(struct reb_simulation* const r){
    struct reb_simulation_integrator_wkm* const ri_wkm = &(r->ri_wkm);
    ri_wkm->corrector = 1;
    ri_wkm->safe_mode = 1;
    ri_wkm->is_synchronized = 1;
    reb_integrator_whfast_reset(r);
    if (ri_wkm->temp_pj){
        free(ri_wkm->temp_pj);
        ri_wkm->temp_pj = NULL;
    }
    ri_wkm->allocated_N = 0;
}
