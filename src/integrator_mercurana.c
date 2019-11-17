/**
 * @file    integrator_mercurana.c
 * @brief   MERCURANA, a modified version of John Chambers' MERCURY algorithm
 *          using the IAS15 integrator and WHFast. It works with planet-planry
 *          collisions, test particles, and additional forces.
 * @author  Hanno Rein, Dan Tamayo
 * 
 * @section LICENSE
 * Copyright (c) 2019 Hanno Rein, Dan Tamayo
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
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_mercurana.h"
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

double reb_integrator_mercurana_L_mercury(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function used by the Mercury integrator.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return 10.*(y*y*y) - 15.*(y*y*y*y) + 6.*(y*y*y*y*y);
    }
}

double reb_integrator_mercurana_dLdr_mercury(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function used by the Mercury integrator.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    double dydr = 1./(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 0.;
    }else{
        return dydr*(20.*(y*y) - 60.*(y*y*y) + 30.*(y*y*y*y));
    }
}

static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}

double reb_integrator_mercurana_L_infinity(const struct reb_simulation* const r, double d, double dcrit){
    // Infinitely differentiable function.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return f(y) /(f(y) + f(1.-y));
    }
}


static void reb_mercurana_encounter_predict(struct reb_simulation* const r, double dt){
    double dts = copysign(1.,dt); 
    dt = fabs(dt);
    struct reb_simulation_integrator_mercurana* rim = &(r->ri_mercurana);
    struct reb_particle* const particles = r->particles;
    const double* const dcrit = rim->dcrit;
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    rim->encounterN = 0;
    for (int i=0; i<N; i++){
        rim->encounter_map[i] = 0;
    }
    for (int i=0; i<N_active; i++){
        for (int j=i+1; j<N; j++){
            const double dx1 = particles[i].x - particles[j].x; // distance at beginning
            const double dy1 = particles[i].y - particles[j].y;
            const double dz1 = particles[i].z - particles[j].z;
            const double r1 = (dx1*dx1 + dy1*dy1 + dz1*dz1);
            const double dvx1 = dts*(particles[i].vx - particles[j].vx); 
            const double dvy1 = dts*(particles[i].vy - particles[j].vy);
            const double dvz1 = dts*(particles[i].vz - particles[j].vz);
            const double dx2 = dx1 +dt*dvx1; // distance at end
            const double dy2 = dy1 +dt*dvy1;
            const double dz2 = dz1 +dt*dvz1;
            const double r2 = (dx2*dx2 + dy2*dy2 + dz2*dz2);
            const double t_closest = (dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);
            const double dx3 = dx1+t_closest*dvx1; // closest approach
            const double dy3 = dy1+t_closest*dvy1;
            const double dz3 = dz1+t_closest*dvz1;
            const double r3 = (dx3*dx3 + dy3*dy3 + dz3*dz3);

            double rmin = MIN(r1,r2);
            if (t_closest/dt>-0.5 && t_closest/dt<1.5){
                rmin = MIN(rmin, r3);
            }

            if (sqrt(rmin)< 2.*MAX(dcrit[i],dcrit[j])){
                if (rim->encounter_map[i]==0){
                    rim->encounter_map[i] = i;
                    rim->encounterN++;
                }
                if (rim->encounter_map[j]==0){
                    rim->encounter_map[j] = j;
                    rim->encounterN++;
                }
            }
        }
    }
}
    
void reb_integrator_mercurana_interaction_step(struct reb_simulation* r, double y, double v, int recalculate){
	const double dt = r->dt;
    if (recalculate){
        reb_update_acceleration(r);
    }
    // Assume particles.a calculated.
	struct reb_particle* const particles = r->particles;
    const int N = r->N;
    if (r->ri_whfast.allocated_N != N){
        r->ri_whfast.allocated_N = N;
        r->ri_whfast.p_jh = realloc(r->ri_whfast.p_jh,sizeof(struct reb_particle)*N);
    }
    struct reb_particle* jerk = r->ri_whfast.p_jh; // Used as a temporary buffer for accelerations
	const double G = r->G;
    for (int j=0; j<N; j++){
        jerk[j].ax = 0; 
        jerk[j].ay = 0; 
        jerk[j].az = 0; 
    }

    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = r->ri_mercurana.L;
    const double* const dcrit = r->ri_mercurana.dcrit;
    if (v!=0.){
        //LAZY
       //     if (r->ri_whfast.allocated_Ntemp != N){
       //         r->ri_whfast.allocated_Ntemp = N;
       //         r->ri_whfast.p_temp = realloc(r->ri_whfast.p_temp,sizeof(struct reb_particle)*N);
       //     }
       //     struct reb_particle* p_temp = r->ri_whfast.p_temp;

       //     // make copy of original positions
       //     memcpy(p_temp,particles,r->N*sizeof(struct reb_particle));

       //     // WHT Eq 10.6
       //     for (unsigned int i=0;i<N;i++){
       //         const double prefac1 = dt*dt*v/y; 
       //         particles[i].x += prefac1 * p_temp[i].ax;
       //         particles[i].y += prefac1 * p_temp[i].ay;
       //         particles[i].z += prefac1 * p_temp[i].az;
       //     }
       //    
       //     // recalculate kick 
       //     reb_update_acceleration(r);

       //     for (unsigned int i=0;i<N;i++){
       //         // reset positions
       //         particles[i].x = p_temp[i].x;
       //         particles[i].y = p_temp[i].y;
       //         particles[i].z = p_temp[i].z;
       //     }
    for (int j=0; j<N; j++){
        for (int i=0; i<j; i++){
            /////////////////
            // Direct Term
            const double dx = particles[j].x - particles[i].x; 
            const double dy = particles[j].y - particles[i].y; 
            const double dz = particles[j].z - particles[i].z; 
            
            const double dax = particles[j].ax - particles[i].ax; 
            const double day = particles[j].ay - particles[i].ay; 
            const double daz = particles[j].az - particles[i].az; 

            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = _L(r,dr,dcritmax);
            const double alphasum = dax*dx+day*dy+daz*dz;
            const double prefact2 = L*G /(dr*dr*dr);
            const double prefact2i = prefact2*particles[i].m;
            const double prefact2j = prefact2*particles[j].m;
            jerk[j].ax    -= dax*prefact2i;
            jerk[j].ay    -= day*prefact2i;
            jerk[j].az    -= daz*prefact2i;
            jerk[i].ax    += dax*prefact2j;
            jerk[i].ay    += day*prefact2j;
            jerk[i].az    += daz*prefact2j;
            const double prefact1 = 3.*alphasum*prefact2 /(dr*dr);
            const double prefact1i = prefact1*particles[i].m;
            const double prefact1j = prefact1*particles[j].m;
            jerk[j].ax    += dx*prefact1i;
            jerk[j].ay    += dy*prefact1i;
            jerk[j].az    += dz*prefact1i;
            jerk[i].ax    -= dx*prefact1j;
            jerk[i].ay    -= dy*prefact1j;
            jerk[i].az    -= dz*prefact1j;
        }
    }
    }

	for (int i=0;i<N;i++){
        const double prefact = dt*dt*dt*v;
		particles[i].vx += dt * y * particles[i].ax + prefact*jerk[i].ax;
		particles[i].vy += dt * y * particles[i].ay + prefact*jerk[i].ay;
		particles[i].vz += dt * y * particles[i].az + prefact*jerk[i].az;
    }
}


static void reb_mercurana_encounter_step(struct reb_simulation* const r, const double _dt){
    // Only particles having a close encounter are integrated by IAS15.
    struct reb_simulation_integrator_mercurana* rim = &(r->ri_mercurana);
    if (rim->encounterN<2){
        return; // If there are no particles (other than the star) having a close encounter, then there is nothing to do.
    }

    int i_enc = 0;
    rim->encounterNactive = 0;
    for (unsigned int i=0; i<r->N; i++){
        if(rim->encounter_map[i]){  
            rim->encounter_map[i_enc] = i;
            i_enc++;
            if (r->N_active==-1 || i<r->N_active){
                rim->encounterNactive++;
            }
        }
    }

    rim->mode = 1;
    
    // run
    const double old_dt = r->dt;
    const double old_t = r->t;
    double t_needed = r->t + _dt; 
        
    reb_integrator_ias15_reset(r);
    
    r->dt = 0.0001*_dt; // start with a small timestep.
    double s = copysign(1.,_dt); 
    
    while(s*r->t < s*t_needed && fabs(r->dt/old_dt)>1e-14 ){
        reb_update_acceleration(r);
        reb_integrator_ias15_part2(r);
        
        if (s*(r->t+r->dt) >  s*t_needed){
            r->dt = t_needed-r->t;
        }

        // Search and resolve collisions
        reb_collision_search(r);

        // Do any additional post_timestep_modifications.
        // Note: post_timestep_modifications is called here but also
        // at the end of the full timestep. The function thus needs
        // to be implemented with care as not to do the same 
        // modification multiple times. To do that, check the value of
        // r->ri_mercurana.mode
        if (r->post_timestep_modifications){
            r->post_timestep_modifications(r);
        }
    }

    // Reset constant for global particles
    r->t = old_t;
    r->dt = old_dt;
    rim->mode = 0;

}

void reb_integrator_mercurana_drift_step(struct reb_simulation* const r, double a){
    const int N = r->N;
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    struct reb_particle* restrict const particles = r->particles;
    const double dt = r->dt;
    reb_mercurana_encounter_predict(r, a*dt);
    for (int i=0;i<N;i++){
        if(rim->encounter_map[i]==0){  // only advance non-encounter particles
            particles[i].x += a*dt*particles[i].vx;
            particles[i].y += a*dt*particles[i].vy;
            particles[i].z += a*dt*particles[i].vz;
        }
    }
    reb_mercurana_encounter_step(r,a*dt);
}
void reb_integrator_mercurana_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurius does not work with variational equations.");
    }
    
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    const int N = r->N;
    
    if (rim->dcrit_allocatedN<N){
        // Need to safe these arrays in SimulationArchive
        rim->dcrit              = realloc(rim->dcrit, sizeof(double)*N);
        rim->dcrit_allocatedN = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        rim->recalculate_dcrit_this_timestep        = 1;
    }
    if (rim->allocatedN<N){
        // These arrays are only used within one timestep. 
        // Can be recreated without loosing bit-wise reproducibility
        rim->encounter_map      = realloc(rim->encounter_map,sizeof(int)*N);
        rim->allocatedN = N;
    }

    if (rim->recalculate_dcrit_this_timestep){
        rim->recalculate_dcrit_this_timestep = 0;
        if (rim->is_synchronized==0){
            reb_integrator_mercurana_synchronize(r);
            reb_warning(r,"MERCURANA: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        rim->dcrit[0] = 2.*r->particles[0].r; // central object only uses physical radius
        const double m0 = r->particles[0].m;
        for (int i=1;i<N;i++){
            const double dx  = r->particles[i].x;  // in dh
            const double dy  = r->particles[i].y;
            const double dz  = r->particles[i].z;
            const double dvx = r->particles[i].vx - r->particles[0].vx; 
            const double dvy = r->particles[i].vy - r->particles[0].vy; 
            const double dvz = r->particles[i].vz - r->particles[0].vz; 
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
            const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;

            const double GM = r->G*(m0+r->particles[i].m);
            const double a = GM*_r / (2.*GM - _r*v2);
            const double vc = sqrt(GM/fabs(a));
            double dcrit = 0;
            // Criteria 1: average velocity
            dcrit = MAX(dcrit, vc*0.4*r->dt);
            // Criteria 2: current velocity
            dcrit = MAX(dcrit, sqrt(v2)*0.4*r->dt);
            // Criteria 3: Hill radius
            dcrit = MAX(dcrit, rim->hillfac*a*pow(r->particles[i].m/(3.*r->particles[0].m),1./3.));
            // Criteria 4: physical radius
            dcrit = MAX(dcrit, 2.*r->particles[i].r);

            rim->dcrit[i] = dcrit;
        }
    }
    
    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurius only works with a direct collision search.");
    }
    
    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_MERCURANA){
        reb_warning(r,"Mercurius has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }
    r->gravity = REB_GRAVITY_MERCURANA;
    rim->mode = 0;
    
    if (rim->L == NULL){
        // Setting default switching function
        rim->L = reb_integrator_mercurana_L_mercury;
    }
    
    // Make copy of particles before the kepler step.
    // Then evolve all particles in kepler step.
    // Result will be used in encounter prediction.
    // Particles having a close encounter will be overwritten 
    // later by encounter step.
    if (rim->is_synchronized){
        // Preprocessor
        double z[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.};
        z[5] = -(z[0]+z[1]+z[2]+z[3]+z[4]);
        double y[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, 0.};
        y[5] = -(y[0]+y[1]+y[2]+y[3]+y[4]);
        double v[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};
        for (int i=0;i<6;i++){
            reb_integrator_mercurana_drift_step(r, z[i]);
            reb_integrator_mercurana_interaction_step(r, y[i], v[i]*2.,1); // recalculates accelerations
        }
        double a1 = -0.0682610383918630;
        reb_integrator_mercurana_drift_step(r,a1);
    }else{
        double a1 = -0.0682610383918630;
        reb_integrator_mercurana_drift_step(r,a1*2.);
    }
}

void reb_integrator_mercurana_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    
    double b1 = 0.2621129352517028;
    reb_integrator_mercurana_interaction_step(r, b1, 0.,0); // does NOT recalculate accelerations
 
    double a1 = -0.0682610383918630;
    double a2 =  0.5-a1;
    reb_integrator_mercurana_drift_step(r, a2);
    
    double b2 = 1.-2.*b1;
    double c2 = 0.0164011128160783;
    reb_integrator_mercurana_interaction_step(r, b2, c2*2.,1); // recalculates accelerations
    
    reb_integrator_mercurana_drift_step(r, a2);
    
    
    reb_integrator_mercurana_interaction_step(r, b1, 0.,1); // recalculates accelerations
        
    rim->is_synchronized = 0;
    if (rim->safe_mode){
        reb_integrator_mercurana_synchronize(r);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

void reb_integrator_mercurana_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    if (rim->is_synchronized == 0){
        r->gravity = REB_GRAVITY_MERCURANA; // needed here again for SimulationArchive
        rim->mode = 0;
        if (rim->L == NULL){
            // Setting default switching function
            rim->L = reb_integrator_mercurana_L_mercury;
        }
        double a1 = -0.0682610383918630;
        reb_integrator_mercurana_drift_step(r,a1);
        
        double z[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.};
        z[5] = -(z[0]+z[1]+z[2]+z[3]+z[4]);
        double y[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, 0.};
        y[5] = -(y[0]+y[1]+y[2]+y[3]+y[4]);
        double v[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};
        for (int i=5;i>=0;i--){
            reb_integrator_mercurana_interaction_step(r, -y[i], -v[i]*2.,1); // recalculates accelerations
            reb_integrator_mercurana_drift_step(r, -z[i]);
        }
        
        rim->is_synchronized = 1;
    }
}

void reb_integrator_mercurana_reset(struct reb_simulation* r){
    r->ri_mercurana.L = NULL;
    r->ri_mercurana.mode = 0;
    r->ri_mercurana.encounterN = 0;
    r->ri_mercurana.is_synchronized = 1;
    r->ri_mercurana.encounterNactive = 0;
    r->ri_mercurana.hillfac = 3;
    // Internal arrays (only used within one timestep)
    free(r->ri_mercurana.particles_backup_additionalforces);
    r->ri_mercurana.particles_backup_additionalforces = NULL;
    free(r->ri_mercurana.encounter_map);
    r->ri_mercurana.encounter_map = NULL;
    r->ri_mercurana.allocatedN = 0;
    r->ri_mercurana.allocatedN_additionalforces = 0;
    // dcrit array
    free(r->ri_mercurana.dcrit);
    r->ri_mercurana.dcrit = NULL;
    r->ri_mercurana.dcrit_allocatedN = 0;
}

