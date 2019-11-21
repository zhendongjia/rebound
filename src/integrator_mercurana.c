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
        return dydr*(30.*(y*y) - 60.*(y*y*y) + 30.*(y*y*y*y));
    }
}

static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}
static double dfdy(double x){
    if (x<0) return 0;
    return exp(-1./x)/(x*x);
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
double reb_integrator_mercurana_dLdr_infinity(const struct reb_simulation* const r, double d, double dcrit){
    // Infinitely differentiable function.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    double dydr = 1./(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 0.;
    }else{
        return dydr*(
                dfdy(y) /(f(y) + f(1.-y))
                -f(y) /(f(y) + f(1.-y))/(f(y) + f(1.-y)) * (dfdy(y) - dfdy(1.-y))
                );
    }
}

const double a_6[2] = {-0.0682610383918630,0.568261038391863038121699}; // a1 a2
const double b_6[2] = {0.2621129352517028, 0.475774129496594366806050}; // b1 b2
const double c_6[2] = {0., 0.0164011128160783}; // c1 c2
const double z_6[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.564398710666239478150885};
const double y_6[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, -0.680252073928462652752103};
const double v_6[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};

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
        rim->hasencounter[i] = 0;
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
                if (rim->hasencounter[i]==0){
                    rim->hasencounter[i] = 1;
                    rim->encounter_map[rim->encounterN] = i;
                    rim->encounterN++;
                }
                if (rim->hasencounter[j]==0){
                    rim->hasencounter[j] = 1;
                    rim->encounter_map[rim->encounterN] = j;
                    rim->encounterN++;
                }
            }
        }
    }
}
    
void reb_integrator_mercurana_drift_step(struct reb_simulation* const r, double a, unsigned int shell){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    const int N = rim.shellN[shell];
    int* map = rim.map[shell];
    struct reb_particle* restrict const particles = r->particles;
    const double dt = r->dt;
    reb_mercurana_encounter_predict(r, a*dt, shell);
    for (int i=0;i<N;i++){  // loop over all particles in shell (includes subshells)
        int mi = map[i]; 
        if(rim->inshell[mi]==shell){  // only advance in-shell particles
            particles[mi].x += a*dt*particles[mi].vx;
            particles[mi].y += a*dt*particles[mi].vy;
            particles[mi].z += a*dt*particles[mi].vz;
        }
    }
    if (shell+1<rim->Nmaxshells){ // does sub-shell exist?
        // advance all sub-shell particles
        //                               a*dt not wrking yet
        dt / Nstepspershell
        reb_integrator_mercurana_preprocessor(r, a*dt, shell+1);
        for (int i=0;i<rim->Nstepspershell;i++){
            reb_integrator_mercurana_step(r, a*dt, shell+1);
        }
        reb_integrator_mercurana_postprocessor(r, a*dt, shell+1);
    }
}

void reb_integrator_mercurana_interaction_step(struct reb_simulation* r, double y, double v, int remove_TODO_recalculate, int shell){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    const int N = rim->shellN[shell];
    const int N_active = rim->shellN_active[shell];
    struct reb_particle* REBOUND_RESTRICT jerk = rim->jerk; 
    struct reb_particle* const particles = r->particles;
    const int testparticle_type   = r->testparticle_type;
    const double G = r->G;
    int* map = rim->map[shell];
    const double* const dcrit_inner = r->ri_mercurana.dcrit[shell];
    const double* const dcrit_outer = NULL;
    if (shell!=0){ // outermost shell has no outer shell
        dcrit_outer = r->ri_mercurana.dcrit[shell-1];
    }

    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = r->ri_mercurana.L;
    double (*_dLdr) (const struct reb_simulation* const r, double d, double dcrit) = r->ri_mercurana.dLdr;
    const double dt = r->dt;
    // Normal force calculation 
    for (int i=0; i<N; i++){
        if (reb_sigint) return;
        int mi = map[i];
        particles[mi].ax = 0; 
        particles[mi].ay = 0; 
        particles[mi].az = 0; 
        for (int j=0; j<N_active; j++){
            int mj = map[j];
            if (mi==mj) continue;
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
            const double dcritmax_inner = MAX(dcrit_inner[mi],dcrit_inner[mj]);
            const double L_inner = _L(r,_r,dcritmax_inner);
            double L_outer = 0.;
            if (shell!=0){
                const double dcritmax_outer = MAX(dcrit_outer[mi],dcrit_outer[mj]);
                L_outer = _L(r,_r,dcritmax_outer);
            }
            const double prefact = -G*particles[mj].m*L_inner*(1.-L_outer)/(_r*_r*_r);
            particles[mi].ax    += prefact*dx;
            particles[mi].ay    += prefact*dy;
            particles[mi].az    += prefact*dz;
        }
    }
    if (_testparticle_type){
        for (int i=0; i<N_active; i++){
            if (reb_sigint) return;
            int mi = map[i];
            for (int j=N_active; j<N; j++){
                int mj = map[j];
                const double dx = particles[mi].x - particles[mj].x;
                const double dy = particles[mi].y - particles[mj].y;
                const double dz = particles[mi].z - particles[mj].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz);
                const double dcritmax_inner = MAX(dcrit_inner[mi],dcrit_inner[mj]);
                const double L_inner = _L(r,_r,dcritmax_inner);
                double L_outer = 0.;
                if (shell!=0){
                    const double dcritmax_outer = MAX(dcrit_outer[mi],dcrit_outer[mj]);
                    L_outer = _L(r,_r,dcritmax_outer);
                }
                const double prefact = -G*particles[mj].m*L_inner*(1.-L_outer)/(_r*_r*_r);
                particles[mi].ax    += prefact*dx;
                particles[mi].ay    += prefact*dy;
                particles[mi].az    += prefact*dz;
            }
        }
    }
    // Jerk calculation
    for (int j=0; j<N; j++){
        jerk[j].ax = 0; 
        jerk[j].ay = 0; 
        jerk[j].az = 0; 
    }
    if (v!=0.){ // is jerk even used?
        for (int j=0; j<N_active; j++){
            int mj = map[j];
            for (int i=0; i<j; i++){
                int mi = map[i];
                const double dx = particles[mj].x - particles[mi].x; 
                const double dy = particles[mj].y - particles[mi].y; 
                const double dz = particles[mj].z - particles[mi].z; 
                
                const double dax = particles[mj].ax - particles[mi].ax; 
                const double day = particles[mj].ay - particles[mi].ay; 
                const double daz = particles[mj].az - particles[mi].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double dcritmax_inner = MAX(dcrit_inner[mi],dcrit_inner[mj]);
                const double L_inner = _L(r,_r,dcritmax_inner);
                const double dLdr_inner = _dLdr(r, dr, dcritmax_inner);
                double L_outer = 0.;
                double dLdr_outer = 0.;
                if (shell!=0){
                    const double dcritmax_outer = MAX(dcrit_outer[mi],dcrit_outer[mj]);
                    L_outer = _L(r,_r,dcritmax_outer);
                }
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2i = L_inner*(1.-L_outer)*prefact2*particles[mi].m;
                const double prefact2j = L_inner*(1.-L_outer)*prefact2*particles[mj].m;
                jerk[j].ax    -= dax*prefact2i;
                jerk[j].ay    -= day*prefact2i;
                jerk[j].az    -= daz*prefact2i;
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2 *(3.*L_inner*(1.-L_outer)/(dr*dr)-(dLdr_inner-dLdr_outer)/dr);
                const double prefact1i = prefact1*particles[mi].m;
                const double prefact1j = prefact1*particles[mj].m;
                jerk[j].ax    += dx*prefact1i;
                jerk[j].ay    += dy*prefact1i;
                jerk[j].az    += dz*prefact1i;
                jerk[i].ax    -= dx*prefact1j;
                jerk[i].ay    -= dy*prefact1j;
                jerk[i].az    -= dz*prefact1j;
            }
        }
        for (int j=0; j<N_active; j++){
            int mj = map[j];
            for (int i=N_active; i<N; i++){
                int mi = map[i];
                const double dx = particles[mj].x - particles[mi].x; 
                const double dy = particles[mj].y - particles[mi].y; 
                const double dz = particles[mj].z - particles[mi].z; 
                
                const double dax = particles[mj].ax - particles[mi].ax; 
                const double day = particles[mj].ay - particles[mi].ay; 
                const double daz = particles[mj].az - particles[mi].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
                const double L_inner = _L(r,_r,dcritmax_inner);
                const double dLdr_inner = _dLdr(r, dr, dcritmax_inner);
                double L_outer = 0.;
                double dLdr_outer = 0.;
                if (shell!=0){
                    const double dcritmax_outer = MAX(dcrit_outer[mi],dcrit_outer[mj]);
                    L_outer = _L(r,_r,dcritmax_outer);
                }
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2j = L_inner*(1.-L_outer)*prefact2*particles[mj].m;
                if (testparticle_type){
                    const double prefact2i = L_inner*(1.-L_outer)*prefact2*particles[mi].m;
                    jerk[j].ax    -= dax*prefact2i;
                    jerk[j].ay    -= day*prefact2i;
                    jerk[j].az    -= daz*prefact2i;
                }
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2 *(3.*L_inner*(1.-L_outer)/(dr*dr)-(dLdr_inner-dLdr_outer)/dr);
                const double prefact1j = prefact1*particles[mj].m;
                if (testparticle_type){
                    const double prefact1i = prefact1*particles[mi].m;
                    jerk[j].ax    += dx*prefact1i;
                    jerk[j].ay    += dy*prefact1i;
                    jerk[j].az    += dz*prefact1i;
                }
                jerk[i].ax    -= dx*prefact1j;
                jerk[i].ay    -= dy*prefact1j;
                jerk[i].az    -= dz*prefact1j;
            }
        }
    }

    for (int i=0;i<N;i++){
        const double prefact = dt*dt*dt*v;
        int mi = map[i];
        particles[mi].vx += dt * y * particles[mi].ax + prefact*jerk[i].ax;
        particles[mi].vy += dt * y * particles[mi].ay + prefact*jerk[i].ay;
        particles[mi].vz += dt * y * particles[mi].az + prefact*jerk[i].az;
    }
}
void reb_integrator_mercurana_preprocessor(struct reb_simulation* const r, int shell, double dt){
    for (int i=0;i<6;i++){
        reb_integrator_mercurana_drift_step(r, dt*z_6[i], shell);
        reb_integrator_mercurana_interaction_step(r, dt*y_6[i], v_6[i]*2., 1, shell);
    }
}
void reb_integrator_mercurana_postprocessor(struct reb_simulation* const r, int shell, double dt){
    for (int i=5;i>=0;i--){
        reb_integrator_mercurana_interaction_step(r, -dt*y_6[i], -v_6[i]*2., 1, shell); TODO need to put dt in front of v_6
        reb_integrator_mercurana_drift_step(r, -dt*z_6[i], shell);
    }
}
void reb_integrator_mercurana_step(struct reb_simulation* const r, int shell, double dt){
    reb_integrator_mercurana_drift_step(r, dt*a_6[0], shell); //TODO combine drift steps
    reb_integrator_mercurana_interaction_step(r, dt*b_6[0], c_6[0]*2., 1, shell);  TODO Need to add dt in front of c_6 
    reb_integrator_mercurana_drift_step(r, dt*a_6[1], shell);
    reb_integrator_mercurana_interaction_step(r, dt*b_6[1], c_6[1]*2., 1, shell); 
    reb_integrator_mercurana_drift_step(r, dt*a_6[1], shell);
    reb_integrator_mercurana_interaction_step(r, dt*b_6[0], c_6[0]*2., 1, shell);
    reb_integrator_mercurana_drift_step(r, dt*a_6[0], shell);
}

void reb_integrator_mercurana_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurius does not work with variational equations.");
    }
    
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    const int N = r->N;
    
    if (rim->allocatedN<N){
        // dcrit
        if (rim->dcrit){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->dcrit[i]);
            }
        }
        rim->dcrit = realloc(rim->dcrit, sizeof(unsigned int*)*(rim->Nmaxshells));
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->dcrit[i] = realloc(rim->dcrit[i], sizeof(unsigned int)*N);
        }
        // map
        if (rim->map){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->map[i]);
            }
        }
        rim->map = realloc(rim->map, sizeof(unsigned int*)*rim->Nmaxshells);
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->map[i] = realloc(rim->map[i], sizeof(unsigned int)*N);
        }
        // inshell
        rim->inshell = realloc(rim->inshell, sizeof(unsigned int)*N);
        // jerk
        rim->jerk = realloc(rim->jerk, sizeof(struct reb_particle)*N);
        // shellN
        rim->shellN = realloc(rim->shellN, sizeof(unsigned int)*rim->Nmaxshells);
        // shellN_active
        rim->shellN_active = realloc(rim->shellN_active, sizeof(unsigned int)*rim->Nmaxshells);

        rim->allocatedN = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        rim->recalculate_dcrit_this_timestep = 1;
    }

    if (rim->recalculate_dcrit_this_timestep){
        rim->recalculate_dcrit_this_timestep = 0;
        if (rim->is_synchronized==0){
            reb_integrator_mercurana_synchronize(r);
            reb_warning(r,"MERCURANA: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        // TODO: need to think about cetral object more
        rim->dcrit[0] = 2.*r->particles[0].r; // central object only uses physical radius
        double critical_timescale = r->dt*rim->dt_frac;
        for (int s=0;s<rim->Nmaxshells-1;s++){ // innermost shell has no dcrit
            for (int i=1;i<N;i++){
                // radius at which escape speed is equal to critical timescale
                // TODO: make machine independent
                double dcrit = pow(critical_timescale*critical_timescale*r->G*r->particles[i].m,1./3.);
                dcrit = MAX(dcrit, 2.*r->particles[i].r);
                rim->dcrit[s][i] = dcrit;
            }
            critical_timescale /= Nstepspershell;
        }
        // Set dcrit for innermost shell to 0 
        for (int i=0;i<N;i++){
            rim->dcrit[rim->Nmaxshells-1][i] = 0;
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
    
    if (rim->L == NULL){
        // Setting default switching function
        rim->L = reb_integrator_mercurana_L_infinity;
        rim->dLdr = reb_integrator_mercurana_dLdr_infinity;
    }
}

void reb_integrator_mercurana_part2(struct reb_simulation* const r){
    if (rim->is_synchronized){
        reb_integrator_mercurana_preprocessor(r, 0, r->dt);
    }
    reb_integrator_mercurana_step(r, 0, r->dt);

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
        if (rim->L == NULL){
            // Setting default switching function
            rim->L = reb_integrator_mercurana_L_infinity;
            rim->dLdr = reb_integrator_mercurana_dLdr_infinity;
        }
        reb_integrator_mercurana_postprocessor(r, 0, r->dt);
        rim->is_synchronized = 1;
    }
}

void reb_integrator_mercurana_reset(struct reb_simulation* r){
    r->ri_mercurana.L = NULL;
    r->ri_mercurana.mode = 0;
    r->ri_mercurana.encounterN = 0;
    r->ri_mercurana.is_synchronized = 1;
    r->ri_mercurana.encounterNactive = 0;
    r->ri_mercurana.Nsubsteps = 10;
    // Internal arrays (only used within one timestep)
    free(r->ri_mercurana.particles_backup_additionalforces);
    r->ri_mercurana.particles_backup_additionalforces = NULL;
    free(r->ri_mercurana.encounter_map);
    r->ri_mercurana.encounter_map = NULL;
    free(r->ri_mercurana.hasencounter);
    r->ri_mercurana.hasencounter = NULL;
    r->ri_mercurana.allocatedN = 0;
    r->ri_mercurana.allocatedN_additionalforces = 0;
    // dcrit array
    free(r->ri_mercurana.dcrit);
    r->ri_mercurana.dcrit = NULL;
    r->ri_mercurana.dcrit_allocatedN = 0;
}

