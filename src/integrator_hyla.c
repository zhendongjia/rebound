/**
 * @file    integrator_hyla.c
 * @brief   HYLA
 * @author  Hanno Rein, Dan Tamayo
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
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_hyla.h"
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

static const double a_6[2] = {-0.0682610383918630,0.568261038391863038121699}; // a1 a2
static const double b_6[2] = {0.2621129352517028, 0.475774129496594366806050}; // b1 b2
static const double c_6[2] = {0., 0.0164011128160783}; // c1 c2
static const double z_6[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.564398710666239478150885};
static const double y_6[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, -0.680252073928462652752103};
static const double v_6[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};

static const double y_4[3] = {0.1859353996846055, 0.0731969797858114, -0.1576624269298081};
static const double z_4[3] = {0.8749306155955435, -0.237106680151022, -0.5363539829039128};

void reb_integrator_hyla_interaction_shell0(struct reb_simulation* r, double y, double v){
    struct reb_simulation_integrator_hyla* const rim = &(r->ri_hyla);
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    const int N_central = rim->Ncentral;
    struct reb_particle* jerk = rim->jerk; 
    struct reb_particle* const particles = r->particles;
    const int testparticle_type   = r->testparticle_type;
    const double G = r->G;
    for (int i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
    
    // Normal force calculation 
    for (int i=N_central; i<N; i++){
        if (reb_sigint) return;
        for (int j=N_central; j<N_active; j++){
            if (i==j) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = -G*particles[j].m/(dr*dr*dr);
            particles[i].ax    += prefact*dx;
            particles[i].ay    += prefact*dy;
            particles[i].az    += prefact*dz;
        }
    }
    if (testparticle_type){
        for (int i=N_central; i<N_active; i++){
            if (reb_sigint) return;
            for (int j=N_active; j<N; j++){
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double dr = sqrt(dx*dx + dy*dy + dz*dz);

                const double prefact = -G*particles[j].m/(dr*dr*dr);
                particles[i].ax    += prefact*dx;
                particles[i].ay    += prefact*dy;
                particles[i].az    += prefact*dz;
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
        for (int j=N_central; j<N_active; j++){
            if (reb_sigint) return;
            for (int i=N_central; i<j; i++){
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2i = prefact2*particles[i].m;
                const double prefact2j = prefact2*particles[j].m;
                jerk[j].ax    -= dax*prefact2i;
                jerk[j].ay    -= day*prefact2i;
                jerk[j].az    -= daz*prefact2i;
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2/dr *3./dr;
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
        for (int j=N_central; j<N_active; j++){
            if (reb_sigint) return;
            for (int i=N_active; i<N; i++){
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2j = prefact2*particles[j].m;
                if (testparticle_type){
                    const double prefact2i = prefact2*particles[i].m;
                    jerk[j].ax    -= dax*prefact2i;
                    jerk[j].ay    -= day*prefact2i;
                    jerk[j].az    -= daz*prefact2i;
                }
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2/dr*3./dr;
                const double prefact1j = prefact1*particles[j].m;
                if (testparticle_type){
                    const double prefact1i = prefact1*particles[i].m;
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
        particles[i].vx += y*particles[i].ax + v*jerk[i].ax;
        particles[i].vy += y*particles[i].ay + v*jerk[i].ay;
        particles[i].vz += y*particles[i].az + v*jerk[i].az;
    }
}
void reb_integrator_hyla_interaction_shell1(struct reb_simulation* r, double y, double v){
    struct reb_simulation_integrator_hyla* const rim = &(r->ri_hyla);
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    const int N_central = rim->Ncentral;
    struct reb_particle* jerk = rim->jerk; 
    struct reb_particle* const particles = r->particles;
    const int testparticle_type   = r->testparticle_type;
    const double G = r->G;
    
    for (int i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }

    // Normal force calculation 
    for (int i=0; i<N_central; i++){
        for (int j=0; j<N_active; j++){
            if (i==j) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = -G*particles[j].m/(dr*dr*dr);
            particles[i].ax    += prefact*dx;
            particles[i].ay    += prefact*dy;
            particles[i].az    += prefact*dz;
        }
    }
    for (int j=0; j<N_central; j++){
        for (int i=0; i<N; i++){
            if (i==j) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = -G*particles[j].m/(dr*dr*dr);
            particles[i].ax    += prefact*dx;
            particles[i].ay    += prefact*dy;
            particles[i].az    += prefact*dz;
        }
    }
    if (testparticle_type){
        // TODO
        printf("not working");
        for (int i=N_central; i<N_active; i++){
            if (reb_sigint) return;
            for (int j=N_active; j<N; j++){
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double dr = sqrt(dx*dx + dy*dy + dz*dz);

                const double prefact = -G*particles[j].m/(dr*dr*dr);
                particles[i].ax    += prefact*dx;
                particles[i].ay    += prefact*dy;
                particles[i].az    += prefact*dz;
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
        for (int j=0; j<N_central; j++){
            if (reb_sigint) return;
            for (int i=0; i<N_active; i++){
                if (i==j) continue;
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2i = prefact2*particles[i].m;
                jerk[j].ax    -= dax*prefact2i;
                jerk[j].ay    -= day*prefact2i;
                jerk[j].az    -= daz*prefact2i;
                const double prefact1 = alphasum*prefact2/dr *3./dr;
                const double prefact1i = prefact1*particles[i].m;
                jerk[j].ax    += dx*prefact1i;
                jerk[j].ay    += dy*prefact1i;
                jerk[j].az    += dz*prefact1i;
            }
        }
        for (int i=0; i<N_central; i++){
            if (reb_sigint) return;
            for (int j=0; j<N_active; j++){
                if (i==j) continue;
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2i = prefact2*particles[i].m;
                jerk[j].ax    -= dax*prefact2i;
                jerk[j].ay    -= day*prefact2i;
                jerk[j].az    -= daz*prefact2i;
                const double prefact1 = alphasum*prefact2/dr *3./dr;
                const double prefact1i = prefact1*particles[i].m;
                jerk[j].ax    += dx*prefact1i;
                jerk[j].ay    += dy*prefact1i;
                jerk[j].az    += dz*prefact1i;
            }
        }
        for (int j=N_central; j<N_active; j++){
            if (reb_sigint) return;
            printf("not working"); //TODO
            for (int i=N_active; i<N; i++){
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2j = prefact2*particles[j].m;
                if (testparticle_type){
                    const double prefact2i = prefact2*particles[i].m;
                    jerk[j].ax    -= dax*prefact2i;
                    jerk[j].ay    -= day*prefact2i;
                    jerk[j].az    -= daz*prefact2i;
                }
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2/dr*3./dr;
                const double prefact1j = prefact1*particles[j].m;
                if (testparticle_type){
                    const double prefact1i = prefact1*particles[i].m;
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
        particles[i].vx += y*particles[i].ax + v*jerk[i].ax;
        particles[i].vy += y*particles[i].ay + v*jerk[i].ay;
        particles[i].vz += y*particles[i].az + v*jerk[i].az;
    }
}
void reb_integrator_hyla_preprocessor(struct reb_simulation* const r, double dt, int order, void (*drift_step)(struct reb_simulation* const r, double a), void (*interaction_step)(struct reb_simulation* const r, double y, double v)){
    switch(order){
        case 6:
            for (int i=0;i<6;i++){
                drift_step(r, dt*z_6[i]);
                interaction_step(r, dt*y_6[i], dt*dt*dt*v_6[i]*2.);
            }
            break;
        case 4:
            for (int i=0;i<3;i++){
                interaction_step(r, dt*y_4[i], 0.);
                drift_step(r, dt*z_4[i]);
            }
            break;
        case 2:
        default:
            break;
    }
}
void reb_integrator_hyla_postprocessor(struct reb_simulation* const r, double dt, int order, void (*drift_step)(struct reb_simulation* const r, double a), void (*interaction_step)(struct reb_simulation* const r, double y, double v)){
    switch(order){
        case 6:
            for (int i=5;i>=0;i--){
                interaction_step(r, -dt*y_6[i], -dt*dt*dt*v_6[i]*2.); 
                drift_step(r, -dt*z_6[i]);
             }
            break;
        case 4:
            for (int i=2;i>=0;i--){
                drift_step(r, -dt*z_4[i]);
                interaction_step(r, -dt*y_4[i], 0.); 
             }
            break;
        case 2:
        default:
            break;
    }
}
void reb_integrator_hyla_drift_shell1(struct reb_simulation* const r, double a){
    struct reb_particle* restrict const particles = r->particles;
    unsigned int N = r->N;
    for (int i=0;i<N;i++){  
        particles[i].x += a*particles[i].vx;
        particles[i].y += a*particles[i].vy;
        particles[i].z += a*particles[i].vz;
    } 
}

void reb_integrator_hyla_step(struct reb_simulation* const r, double dt, int order, void (*drift_step)(struct reb_simulation* const r, double a), void (*interaction_step)(struct reb_simulation* const r, double y, double v)){
    switch(order){
        case 6:
            drift_step(r, dt*a_6[0]); //TODO combine drift steps
            interaction_step(r, dt*b_6[0], dt*dt*dt*c_6[0]*2.); 
            drift_step(r, dt*a_6[1]);
            interaction_step(r, dt*b_6[1], dt*dt*dt*c_6[1]*2.); 
            drift_step(r, dt*a_6[1]);
            interaction_step(r, dt*b_6[0], dt*dt*dt*c_6[0]*2.);
            drift_step(r, dt*a_6[0]);
            break;
        case 4:
            drift_step(r, dt*0.5); //TODO combine drift steps
            interaction_step(r, dt, dt*dt*dt/24.*2); 
            drift_step(r, dt*0.5);
            break;
        case 2:
        default:
            drift_step(r, dt*0.5); //TODO combine drift steps
            interaction_step(r, dt, 0.);
            drift_step(r, dt*0.5);
            break;
    }
}
void reb_integrator_hyla_drift_shell0(struct reb_simulation* const r, double a){
    struct reb_simulation_integrator_hyla* const rim = &(r->ri_hyla);
    double as = a/rim->Nstepspershell;
    reb_integrator_hyla_preprocessor(r, as, rim->ordersubsteps, reb_integrator_hyla_drift_shell1, reb_integrator_hyla_interaction_shell1);
    for (int i=0;i<rim->Nstepspershell;i++){
        reb_integrator_hyla_step(r, as, rim->ordersubsteps, reb_integrator_hyla_drift_shell1, reb_integrator_hyla_interaction_shell1);
    }
    reb_integrator_hyla_postprocessor(r, as, rim->ordersubsteps, reb_integrator_hyla_drift_shell1, reb_integrator_hyla_interaction_shell1);
}

void reb_integrator_hyla_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurius does not work with variational equations.");
    }
    
    struct reb_simulation_integrator_hyla* const rim = &(r->ri_hyla);
    const int N = r->N;
    
    if (rim->allocatedN<N){
        // dcrit
        if (rim->dcrit){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->dcrit[i]);
            }
        }
        rim->dcrit = realloc(rim->dcrit, sizeof(double*)*(rim->Nmaxshells));
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->dcrit[i] = malloc(sizeof(double)*N);
        }
        // map
        if (rim->map){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->map[i]);
            }
        }
        rim->map = realloc(rim->map, sizeof(unsigned int*)*rim->Nmaxshells);
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->map[i] = malloc(sizeof(unsigned int)*N);
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
    }

    r->gravity = REB_GRAVITY_NONE;

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurius only works with a direct collision search.");
    }
    
}

void reb_integrator_hyla_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_hyla* const rim = &(r->ri_hyla);
    rim->shellN[0] = r->N;
    rim->shellN_active[0] = r->N_active==-1?r->N:r->N_active;

    if (rim->is_synchronized){
        reb_integrator_hyla_preprocessor(r, r->dt, rim->order, reb_integrator_hyla_drift_shell0, reb_integrator_hyla_interaction_shell0);
    }
    reb_integrator_hyla_step(r, r->dt, rim->order, reb_integrator_hyla_drift_shell0, reb_integrator_hyla_interaction_shell0);

    rim->is_synchronized = 0;
    if (rim->safe_mode){
        reb_integrator_hyla_synchronize(r);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

void reb_integrator_hyla_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_hyla* const rim = &(r->ri_hyla);
    if (rim->is_synchronized == 0){
        r->gravity_ignore_terms = 2;
        reb_integrator_hyla_postprocessor(r, r->dt, rim->order, reb_integrator_hyla_drift_shell0, reb_integrator_hyla_interaction_shell0);
        rim->is_synchronized = 1;
    }
}

void reb_integrator_hyla_reset(struct reb_simulation* r){
    if (r->ri_hyla.allocatedN){
        for (int i=0;i<r->ri_hyla.Nmaxshells;i++){
            free(r->ri_hyla.map[i]);
        }
        free(r->ri_hyla.map);
        free(r->ri_hyla.inshell);
        free(r->ri_hyla.shellN);
        free(r->ri_hyla.shellN_active);
        free(r->ri_hyla.jerk);
    }
    r->ri_hyla.allocatedN = 0;
    r->ri_hyla.map = NULL;
    r->ri_hyla.inshell = NULL;
    r->ri_hyla.shellN = NULL;
    r->ri_hyla.shellN_active = NULL;
    r->ri_hyla.jerk = NULL;
    
    r->ri_hyla.Ncentral = 1;
    r->ri_hyla.order = 2;
    r->ri_hyla.ordersubsteps = 2;
    r->ri_hyla.safe_mode = 1;
    r->ri_hyla.Nmaxshells = 10;
    r->ri_hyla.Nmaxshellused = 1;
    r->ri_hyla.Nstepspershell = 10;
    r->ri_hyla.is_synchronized = 1;
    
}

