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

// Machine independent implementation of pow(*,1./3.) using Newton's method.
// Speed is not an issue. Only used to calculate dcrit.
static double sqrt3(double a){
    double x = 1.;
    for (int k=0; k<200;k++){  // A smaller number should be ok too.
        double x2 = x*x;
        x += (a/x2-x)/3.;
    }
    return x;
}

static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}
static double dfdy(double x){
    if (x<0) return 0;
    return exp(-1./x)/(x*x);
}

double reb_integrator_mercurana_L_infinity(const struct reb_simulation* const r, double d, double ri, double ro){
    // Infinitely differentiable function.
    double y = (d-ri)/(ro-ri);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return f(y) /(f(y) + f(1.-y));
    }
}
double reb_integrator_mercurana_dLdr_infinity(const struct reb_simulation* const r, double d, double ri, double ro){
    // Infinitely differentiable function.
    double y = (d-ri)/(ro-ri);
    double dydr = 1./(ro-ri);
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

void reb_integrator_mercurana_preprocessor(struct reb_simulation* const r, double dt, int shell, int order);
void reb_integrator_mercurana_postprocessor(struct reb_simulation* const r, double dt, int shell, int order);
void reb_integrator_mercurana_step(struct reb_simulation* const r, double dt, int shell, int order);

static const double a_6[2] = {-0.0682610383918630,0.568261038391863038121699}; // a1 a2
static const double b_6[2] = {0.2621129352517028, 0.475774129496594366806050}; // b1 b2
static const double c_6[2] = {0., 0.0164011128160783}; // c1 c2
static const double z_6[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.564398710666239478150885};
static const double y_6[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, -0.680252073928462652752103};
static const double v_6[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};

static const double y_4[3] = {0.1859353996846055, 0.0731969797858114, -0.1576624269298081};
static const double z_4[3] = {0.8749306155955435, -0.237106680151022, -0.5363539829039128};

void reb_mercurana_predict_rmin2(struct reb_particle p1, struct reb_particle p2, double dt, double* rmin2_ab, double* rmin2_abc){ 
    double dts = copysign(1.,dt); 
    dt = fabs(dt);
    const double dx1 = p1.x - p2.x; // distance at beginning
    const double dy1 = p1.y - p2.y;
    const double dz1 = p1.z - p2.z;
    const double r1 = (dx1*dx1 + dy1*dy1 + dz1*dz1);
    const double dvx1 = dts*(p1.vx - p2.vx); 
    const double dvy1 = dts*(p1.vy - p2.vy);
    const double dvz1 = dts*(p1.vz - p2.vz);
    const double dx2 = dx1 +dt*dvx1; // distance at end
    const double dy2 = dy1 +dt*dvy1;
    const double dz2 = dz1 +dt*dvz1;
    const double r2 = (dx2*dx2 + dy2*dy2 + dz2*dz2);
    const double t_closest = (dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);
    const double dx3 = dx1+t_closest*dvx1; // closest approach
    const double dy3 = dy1+t_closest*dvy1;
    const double dz3 = dz1+t_closest*dvz1;
    const double r3 = (dx3*dx3 + dy3*dy3 + dz3*dz3);

    *rmin2_ab = MIN(r1,r2);
    if (t_closest/dt>=0. && t_closest/dt<=1.){
        *rmin2_abc = MIN(*rmin2_ab, r3);
    }else{
        *rmin2_abc = *rmin2_ab;
    }
}

static void reb_mercurana_encounter_predict(struct reb_simulation* const r, double dt, int shell){
    struct reb_simulation_integrator_mercurana* rim = &(r->ri_mercurana);
    struct reb_particle* const particles = r->particles;
    const double* const dcrit = rim->dcrit[shell];
    const int N = rim->shellN[shell];
    const int N_active = rim->shellN_active[shell];
    unsigned int* map = rim->map[shell];

    // Put all particles in current shell by default
    for (int i=0; i<N; i++){
        int mi = map[i]; 
        rim->inshell[mi] = 1;
    }
    
    if (shell+1>=rim->Nmaxshells){ // does sub-shell exist?
        return;
    }
	
    int collisions_N = 0;
    
    // Check if particles are in sub-shell
    rim->shellN[shell+1] = 0;
    rim->shellN_active[shell+1] = 0;

    // Note: Need to find the active particles first
    // TODO make this O(N^2/2)
    for (int i=0; i<N_active; i++){
        int mi = map[i]; 
        for (int j=0; j<N; j++){
            int mj = map[j]; 
            if (i==j) continue;
            double rmin2_ab, rmin2_abc;
            reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt,&rmin2_ab,&rmin2_abc);
            double rsum = particles[mi].r+particles[mj].r;
            double dcritsum = dcrit[mi]+dcrit[mj];
            // Only check if particles are overlapping at the beginning or end of the timestep.
            // TODO: Add flag to also check for overlap during the timestep.
            if (rmin2_ab<rsum*rsum){
                ////printf("col found i=%d j=%d s=%d\n", mi,mj, shell);
                //printf("part i m=%f xyz=%f %f %f \n", r->particles[mi].m, r->particles[mi].x, r->particles[mi].y, r->particles[mi].z );
                //printf("part j m=%f xyz=%f %f %f \n", r->particles[mj].m, r->particles[mj].x, r->particles[mj].y, r->particles[mj].z );
                // Collision occured.
                // Note: No check if particles are approaching each other. 
                // Add particles to collision array.
                if (r->collisions_allocatedN<=collisions_N){
                    // Allocate memory if there is no space in array.
                    // Init to 32 if no space has been allocated yet, otherwise double it.
                    r->collisions_allocatedN = r->collisions_allocatedN ? r->collisions_allocatedN * 2 : 32;
                    r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->collisions_allocatedN);
                }
                r->collisions[collisions_N].p1 = mi;
                r->collisions[collisions_N].p2 = mj;
                //r->collisions[collisions_N].gb = 0;
                collisions_N++;
                // Only mark collisions here. Resolve later.
            }
            if (rmin2_abc< dcritsum*dcritsum){ 
                // j particle will be added later (active particles need to be in array first)
                rim->inshell[mi] = 0;
                rim->map[shell+1][rim->shellN[shell+1]] = mi;
                rim->shellN[shell+1]++;
                rim->shellN_active[shell+1]++;
                break; // only add particle i once
            }

        }
    }
    for (int i=N_active; i<N; i++){
        int mi = map[i]; 
        for (int j=0; j<N_active; j++){
            int mj = map[j]; 
            double rmin2_ab, rmin2_abc;
            reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt,&rmin2_ab,&rmin2_abc);
            double rsum = particles[mi].r+particles[mj].r;
            double dcritsum = dcrit[mi]+dcrit[mj];
            // Only check if particles are overlapping at the beginning or end of the timestep.
            // TODO: Add flag to also check for overlap during the timestep.
            if (rmin2_ab<rsum*rsum){
                // Collision occured.
                // Note: No check if particles are approaching each other. 
                // Add particles to collision array.
                if (r->collisions_allocatedN<=collisions_N){
                    // Allocate memory if there is no space in array.
                    // Init to 32 if no space has been allocated yet, otherwise double it.
                    r->collisions_allocatedN = r->collisions_allocatedN ? r->collisions_allocatedN * 2 : 32;
                    r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->collisions_allocatedN);
                }
                r->collisions[collisions_N].p1 = mi;
                r->collisions[collisions_N].p2 = mj;
                //r->collisions[collisions_N].gb = 0;
                collisions_N++;
                // Only mark collisions here. Resolve later.
            }
            if (rmin2_abc< dcritsum*dcritsum){ 
                rim->inshell[mi] = 0;
                rim->map[shell+1][rim->shellN[shell+1]] = mi;
                rim->shellN[shell+1]++;
                break; // only add particle i once
            }
        }
    }
	// randomize collisions
	for (int i=0;i<collisions_N;i++){
		int new = rand()%collisions_N;
		struct reb_collision c1 = r->collisions[i];
		r->collisions[i] = r->collisions[new];
		r->collisions[new] = c1;
	}
	
	int (*resolve) (struct reb_simulation* const r, struct reb_collision c) = r->collision_resolve;
	if (resolve==NULL){
		// Default is hard sphere
		resolve = reb_collision_resolve_halt;
	}
    int N_orig = r->N;
	for (int i=0;i<collisions_N;i++){
        
        struct reb_collision c = r->collisions[i];
        if (c.p1 != -1 && c.p2 != -1){
            // Resolve collision
            int outcome = resolve(r, c);
            //printf("resolce %d\n",outcome);
            
            // Remove particles
            if (outcome & 1){
                // Remove p1
                if (c.p2==r->N-r->N_var-1 && !(r->tree_root)){
                    // Particles swapped
                    c.p2 = c.p1;
                }
                reb_remove(r,c.p1,r->collision_resolve_keep_sorted);
                // Check for pair
                for (int j=i+1;j<collisions_N;j++){
                    struct reb_collision cp = r->collisions[j];
                    if (cp.p1==c.p1 || cp.p2==c.p1){
                        r->collisions[j].p1 = -1;
                        r->collisions[j].p2 = -1;
                        // Will be skipped.
                    }
                    if (cp.p1==r->N-r->N_var){
                        r->collisions[j].p1 = c.p1;
                    }
                    if (cp.p2==r->N-r->N_var){
                        r->collisions[j].p2 = c.p1;
                    }
                }
            }
            if (outcome & 2){
                // Remove p2
                reb_remove(r,c.p2,r->collision_resolve_keep_sorted);
                // Check for pair
                for (int j=i+1;j<collisions_N;j++){
                    struct reb_collision cp = r->collisions[j];
                    if (cp.p1==c.p2 || cp.p2==c.p2){
                        r->collisions[j].p1 = -1;
                        r->collisions[j].p2 = -1;
                        // Will be skipped.
                    }
                    if (cp.p1==r->N-r->N_var){
                        r->collisions[j].p1 = c.p2;
                    }
                    if (cp.p2==r->N-r->N_var){
                        r->collisions[j].p2 = c.p2;
                    }
                }
            }
        }
	}
    if (N_orig!=r->N){
        // redo prediction. Substeps might not be needed anymore
        reb_mercurana_encounter_predict(r, dt,shell);
    }
}
    
void reb_integrator_mercurana_drift_step(struct reb_simulation* const r, double a, unsigned int shell){
    //printf("drift s=%d\n",shell);
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    struct reb_particle* restrict const particles = r->particles;
    reb_mercurana_encounter_predict(r, a, shell);
    unsigned int* map = rim->map[shell];
    unsigned int N = rim->shellN[shell];
    for (int i=0;i<N;i++){  // loop over all particles in shell (includes subshells)
        int mi = map[i]; 
        if(rim->inshell[mi]){  // only advance in-shell particles
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
        }
    }
    if (shell+1<rim->Nmaxshells){ // does sub-shell exist?
        if (rim->shellN[shell+1]>0){ // are there particles in it?
            rim->Nmaxshellused = MAX(rim->Nmaxshellused, shell+2);
            // advance all sub-shell particles
            double as = a/rim->Nstepspershell;
            reb_integrator_mercurana_preprocessor(r, as, shell+1, rim->ordersubsteps);
            for (int i=0;i<rim->Nstepspershell;i++){
                reb_integrator_mercurana_step(r, as, shell+1, rim->ordersubsteps);
            }
            reb_integrator_mercurana_postprocessor(r, as, shell+1, rim->ordersubsteps);
        }
    }
}

void reb_integrator_mercurana_interaction_step(struct reb_simulation* r, double y, double v, int shell){
    //printf("inter s=%d\n",shell);
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    const int N = rim->shellN[shell];
    const int N_active = rim->shellN_active[shell];
    struct reb_particle* jerk = rim->jerk; 
    struct reb_particle* const particles = r->particles;
    const int testparticle_type   = r->testparticle_type;
    const double G = r->G;
    unsigned int* map = rim->map[shell];
    const double* dcrit_ii = NULL;
    const double* dcrit_i = NULL;
    const double* dcrit_o = NULL;
    if (shell<rim->Nmaxshells-1){
        dcrit_ii = r->ri_mercurana.dcrit[shell+1];
    }
    dcrit_i = r->ri_mercurana.dcrit[shell];
    if (shell!=0){
        dcrit_o = r->ri_mercurana.dcrit[shell-1];
    }

    double (*_L) (const struct reb_simulation* const r, double d, double dcrit, double fracin) = r->ri_mercurana.L;
    double (*_dLdr) (const struct reb_simulation* const r, double d, double dcrit, double fracin) = r->ri_mercurana.dLdr;
    // Normal force calculation 
    for (int i=0; i<N; i++){
        if (reb_sigint) return;
        int mi = map[i];
        particles[mi].ax = 0; 
        particles[mi].ay = 0; 
        particles[mi].az = 0; 
        if (reb_sigint) return;
        for (int j=0; j<N_active; j++){
            int mj = map[j];
            if (mi==mj) continue;
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);
            const double dc_i = dcrit_i[mi]+dcrit_i[mj];
            double Lsum = 0.;
            double dc_o = 0;
            if (dcrit_o){
                dc_o = dcrit_o[mi]+dcrit_o[mj];
                Lsum -= _L(r,dr,dc_i,dc_o);
            }
            double dc_ii = 0;
            if (dcrit_ii){
                dc_ii = dcrit_ii[mi]+dcrit_ii[mj];
                Lsum += _L(r,dr,dc_ii,dc_i);
            }else{
                Lsum += 1; // Innermost
            }

            const double prefact = -G*particles[mj].m*Lsum/(dr*dr*dr);
            particles[mi].ax    += prefact*dx;
            particles[mi].ay    += prefact*dy;
            particles[mi].az    += prefact*dz;
        }
    }
    if (testparticle_type){
        for (int i=0; i<N_active; i++){
            if (reb_sigint) return;
            int mi = map[i];
            for (int j=N_active; j<N; j++){
                int mj = map[j];
                const double dx = particles[mi].x - particles[mj].x;
                const double dy = particles[mi].y - particles[mj].y;
                const double dz = particles[mi].z - particles[mj].z;
                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double dc_i = dcrit_i[mi]+dcrit_i[mj];
                double Lsum = 0.;
                double dc_o = 0;
                if (dcrit_o){
                    dc_o = dcrit_o[mi]+dcrit_o[mj];
                    Lsum -= _L(r,dr,dc_i,dc_o);
                }
                double dc_ii = 0;
                if (dcrit_ii){
                    dc_ii = dcrit_ii[mi]+dcrit_ii[mj];
                    Lsum += _L(r,dr,dc_ii,dc_i);
                }else{
                    Lsum += 1; // Innermost
                }

                const double prefact = -G*particles[mj].m*Lsum/(dr*dr*dr);
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
            if (reb_sigint) return;
            for (int i=0; i<j; i++){
                int mi = map[i];
                const double dx = particles[mj].x - particles[mi].x; 
                const double dy = particles[mj].y - particles[mi].y; 
                const double dz = particles[mj].z - particles[mi].z; 
                
                const double dax = particles[mj].ax - particles[mi].ax; 
                const double day = particles[mj].ay - particles[mi].ay; 
                const double daz = particles[mj].az - particles[mi].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double dc_i = dcrit_i[mi]+dcrit_i[mj];
                double Lsum = 0.;
                double dLdrsum = 0.;
                if (dcrit_o){
                    double dc_o = dcrit_o[mi]+dcrit_o[mj];
                    Lsum    -=    _L(r,dr,dc_i,dc_o);
                    dLdrsum -= _dLdr(r,dr,dc_i,dc_o);
                }
                if (dcrit_ii){
                    double dc_ii = dcrit_ii[mi]+dcrit_ii[mj];
                    Lsum    +=    _L(r,dr,dc_ii,dc_i);
                    dLdrsum += _dLdr(r,dr,dc_ii,dc_i);
                }else{
                    Lsum += 1; // Innermost
                }
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2i = Lsum*prefact2*particles[mi].m;
                const double prefact2j = Lsum*prefact2*particles[mj].m;
                jerk[j].ax    -= dax*prefact2i;
                jerk[j].ay    -= day*prefact2i;
                jerk[j].az    -= daz*prefact2i;
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2/dr *(3.*Lsum/dr-dLdrsum);
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
            if (reb_sigint) return;
            for (int i=N_active; i<N; i++){
                int mi = map[i];
                const double dx = particles[mj].x - particles[mi].x; 
                const double dy = particles[mj].y - particles[mi].y; 
                const double dz = particles[mj].z - particles[mi].z; 
                
                const double dax = particles[mj].ax - particles[mi].ax; 
                const double day = particles[mj].ay - particles[mi].ay; 
                const double daz = particles[mj].az - particles[mi].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double dc_i = dcrit_i[mi]+dcrit_i[mj];
                double Lsum = 0.;
                double dLdrsum = 0.;
                if (dcrit_o){
                    double dc_o = dcrit_o[mi]+dcrit_o[mj];
                    Lsum    -=    _L(r,dr,dc_i,dc_o);
                    dLdrsum -= _dLdr(r,dr,dc_i,dc_o);
                }
                if (dcrit_ii){
                    double dc_ii = dcrit_ii[mi]+dcrit_ii[mj];
                    Lsum    +=    _L(r,dr,dc_ii,dc_i);
                    dLdrsum += _dLdr(r,dr,dc_ii,dc_i);
                }else{
                    Lsum += 1; // Innermost
                }
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2j = Lsum*prefact2*particles[mj].m;
                if (testparticle_type){
                    const double prefact2i = Lsum*prefact2*particles[mi].m;
                    jerk[j].ax    -= dax*prefact2i;
                    jerk[j].ay    -= day*prefact2i;
                    jerk[j].az    -= daz*prefact2i;
                }
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2/dr*(3.*Lsum/dr-dLdrsum);
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
        int mi = map[i];
        particles[mi].vx += y*particles[mi].ax + v*jerk[i].ax;
        particles[mi].vy += y*particles[mi].ay + v*jerk[i].ay;
        particles[mi].vz += y*particles[mi].az + v*jerk[i].az;
    }
}
void reb_integrator_mercurana_preprocessor(struct reb_simulation* const r, double dt, int shell, int order){
    switch(order){
        case 6:
            for (int i=0;i<6;i++){
                reb_integrator_mercurana_drift_step(r, dt*z_6[i], shell);
                reb_integrator_mercurana_interaction_step(r, dt*y_6[i], dt*dt*dt*v_6[i]*2., shell);
            }
            break;
        case 4:
            for (int i=0;i<3;i++){
                reb_integrator_mercurana_interaction_step(r, dt*y_4[i], 0., shell);
                reb_integrator_mercurana_drift_step(r, dt*z_4[i], shell);
            }
            break;
        case 2:
        default:
            break;
    }
}
void reb_integrator_mercurana_postprocessor(struct reb_simulation* const r, double dt, int shell, int order){
    switch(order){
        case 6:
            for (int i=5;i>=0;i--){
                reb_integrator_mercurana_interaction_step(r, -dt*y_6[i], -dt*dt*dt*v_6[i]*2., shell); 
                reb_integrator_mercurana_drift_step(r, -dt*z_6[i], shell);
             }
            break;
        case 4:
            for (int i=2;i>=0;i--){
                reb_integrator_mercurana_drift_step(r, -dt*z_4[i], shell);
                reb_integrator_mercurana_interaction_step(r, -dt*y_4[i], 0., shell); 
             }
            break;
        case 2:
        default:
            break;
    }
}
void reb_integrator_mercurana_step(struct reb_simulation* const r, double dt, int shell, int order){
    switch(order){
        case 6:
            reb_integrator_mercurana_drift_step(r, dt*a_6[0], shell); //TODO combine drift steps
            reb_integrator_mercurana_interaction_step(r, dt*b_6[0], dt*dt*dt*c_6[0]*2., shell); 
            reb_integrator_mercurana_drift_step(r, dt*a_6[1], shell);
            reb_integrator_mercurana_interaction_step(r, dt*b_6[1], dt*dt*dt*c_6[1]*2., shell); 
            reb_integrator_mercurana_drift_step(r, dt*a_6[1], shell);
            reb_integrator_mercurana_interaction_step(r, dt*b_6[0], dt*dt*dt*c_6[0]*2., shell);
            reb_integrator_mercurana_drift_step(r, dt*a_6[0], shell);
            break;
        case 4:
            reb_integrator_mercurana_drift_step(r, dt*0.5, shell); //TODO combine drift steps
            reb_integrator_mercurana_interaction_step(r, dt, dt*dt*dt/24.*2, shell); 
            reb_integrator_mercurana_drift_step(r, dt*0.5, shell);
            break;
        case 2:
        default:
            reb_integrator_mercurana_drift_step(r, dt*0.5, shell); //TODO combine drift steps
            reb_integrator_mercurana_interaction_step(r, dt, 0., shell);
            reb_integrator_mercurana_drift_step(r, dt*0.5, shell);
            break;
    }
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
        rim->recalculate_dcrit_this_timestep = 1;
    }

    if (rim->recalculate_dcrit_this_timestep){
        rim->recalculate_dcrit_this_timestep = 0;
        if (rim->is_synchronized==0){
            reb_integrator_mercurana_synchronize(r);
            reb_warning(r,"MERCURANA: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        double dt_shell = r->dt;
        for (int s=0;s<rim->Nmaxshells;s++){ // innermost shell has no dcrit
            for (int i=0;i<N;i++){
                // distance where dt/dt_frac is equal to dynamical timescale
                // Note: particle radius not needed here.
                double T = dt_shell/(rim->dt_frac*2.*M_PI);
                double dcrit = sqrt3(T*T*r->G*r->particles[i].m);
                rim->dcrit[s][i] = dcrit;
            }
            double longest_drift_step_in_shell = 0.5;                        // 2nd + 4th order
            if ((s==0 && rim->order==6) || (s>0 && rim->ordersubsteps==6)){  // 6th order
                longest_drift_step_in_shell = a_6[1];
            }
            dt_shell *= longest_drift_step_in_shell;
            dt_shell /= rim->Nstepspershell;
            // Initialize shell numbers to zero (not needed, but helps debugging)
            rim->shellN[s] = 0;
            rim->shellN_active[s] = 0;
        }
        for (int i=0;i<N;i++){
            // Set map to identity for outer-most shell
            rim->map[0][i] = i;
        }

    }
    
    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurius only works with a direct collision search.");
    }
    
    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_NONE){
        reb_warning(r,"Mercurius has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }
    r->gravity = REB_GRAVITY_NONE;
    
    if (rim->L == NULL){
        // Setting default switching function
        rim->L = reb_integrator_mercurana_L_infinity;
        rim->dLdr = reb_integrator_mercurana_dLdr_infinity;
    }
}

void reb_integrator_mercurana_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    rim->shellN[0] = r->N;
    rim->shellN_active[0] = r->N_active==-1?r->N:r->N_active;

    if (rim->is_synchronized){
        reb_integrator_mercurana_preprocessor(r, r->dt, 0, rim->order);
    }
    //double alpha1 = 1.35120719195965763404768780897;
    //reb_integrator_mercurana_step(r, alpha1*r->dt, 0, rim->order);
    //reb_integrator_mercurana_step(r, (1.-2.*alpha1)*r->dt, 0, rim->order);
    //reb_integrator_mercurana_step(r, alpha1*r->dt, 0, rim->order);
    reb_integrator_mercurana_step(r, r->dt, 0, rim->order);

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
        r->gravity = REB_GRAVITY_NONE; // needed here again for SimulationArchive
        if (rim->L == NULL){
            // Setting default switching function
            rim->L = reb_integrator_mercurana_L_infinity;
            rim->dLdr = reb_integrator_mercurana_dLdr_infinity;
        }
        reb_integrator_mercurana_postprocessor(r, r->dt, 0, rim->order);
        rim->is_synchronized = 1;
    }
}

void reb_integrator_mercurana_reset(struct reb_simulation* r){
    if (r->ri_mercurana.allocatedN){
        for (int i=0;i<r->ri_mercurana.Nmaxshells;i++){
            free(r->ri_mercurana.map[i]);
            free(r->ri_mercurana.dcrit[i]);
        }
        free(r->ri_mercurana.map);
        free(r->ri_mercurana.dcrit);
        free(r->ri_mercurana.inshell);
        free(r->ri_mercurana.shellN);
        free(r->ri_mercurana.shellN_active);
        free(r->ri_mercurana.jerk);
    }
    r->ri_mercurana.allocatedN = 0;
    r->ri_mercurana.map = NULL;
    r->ri_mercurana.dcrit = NULL;
    r->ri_mercurana.inshell = NULL;
    r->ri_mercurana.shellN = NULL;
    r->ri_mercurana.shellN_active = NULL;
    r->ri_mercurana.jerk = NULL;
    
    r->ri_mercurana.order = 2;
    r->ri_mercurana.ordersubsteps = 2;
    r->ri_mercurana.safe_mode = 1;
    r->ri_mercurana.dt_frac = 0.1;
    r->ri_mercurana.Nmaxshells = 10;
    r->ri_mercurana.Nmaxshellused = 1;
    r->ri_mercurana.Nstepspershell = 10;
    r->ri_mercurana.recalculate_dcrit_this_timestep = 0;
    r->ri_mercurana.is_synchronized = 1;
    r->ri_mercurana.L = NULL;
    r->ri_mercurana.dLdr = NULL;
    
}

