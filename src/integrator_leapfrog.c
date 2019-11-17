/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include "rebound.h"

void reb_leapfrog_apply_C(struct reb_simulation* r, double y, double v, int recalculate){
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

    if (v!=0.){
        //LAZY
    //        // Need temporary array to store old positions
    //        if (r->ri_whfast.allocated_Ntemp != N){
    //            r->ri_whfast.allocated_Ntemp = N;
    //            r->ri_whfast.p_temp = realloc(r->ri_whfast.p_temp,sizeof(struct reb_particle)*N);
    //        }
    //        struct reb_particle* p_temp = r->ri_whfast.p_temp;

    //        // make copy of original positions
    //        memcpy(p_temp,particles,r->N*sizeof(struct reb_particle));

    //        // WHT Eq 10.6
    //        for (unsigned int i=0;i<N;i++){
    //            const double prefac1 = dt*dt*v/y; 
    //            particles[i].x += prefac1 * p_temp[i].ax;
    //            particles[i].y += prefac1 * p_temp[i].ay;
    //            particles[i].z += prefac1 * p_temp[i].az;
    //        }
    //       
    //        // recalculate kick 
    //        reb_update_acceleration(r);

    //        for (unsigned int i=0;i<N;i++){
    //            // reset positions
    //            particles[i].x = p_temp[i].x;
    //            particles[i].y = p_temp[i].y;
    //            particles[i].z = p_temp[i].z;
    //        }
    //}
    //MODIFIED
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
void reb_leapfrog_apply_A(struct reb_simulation* r, double a){
	const int N = r->N;
	struct reb_particle* restrict const particles = r->particles;
	const double dt = r->dt;

    for (int i=0;i<N;i++){
        particles[i].x  += a* dt * particles[i].vx;
        particles[i].y  += a* dt * particles[i].vy;
        particles[i].z  += a* dt * particles[i].vz;
    }
}

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void reb_integrator_leapfrog_part1(struct reb_simulation* r){
    r->gravity_ignore_terms = 0;
	const int N = r->N;
	struct reb_particle* restrict const particles = r->particles;
	const double dt = r->dt;

    if (r->ri_whfast.kernel == REB_WHFAST_KERNEL_6ABA363){
        if (r->ri_whfast.is_synchronized){
            // Preprocessor
            double z[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.};
            z[5] = -(z[0]+z[1]+z[2]+z[3]+z[4]);
            double y[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, 0.};
            y[5] = -(y[0]+y[1]+y[2]+y[3]+y[4]);
            double v[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};
            for (int i=0;i<6;i++){
                reb_leapfrog_apply_A(r, z[i]);
                reb_leapfrog_apply_C(r, y[i], v[i]*2.,1); // recalculates accelerations
            }
            double a1 = -0.0682610383918630;
            reb_leapfrog_apply_A(r, a1);
        }else{
            double a1 = -0.0682610383918630;
            reb_leapfrog_apply_A(r, 2.*a1);
        }
    }else{
        // Traditional leapfrog
        for (int i=0;i<N;i++){
            particles[i].x  += 0.5* dt * particles[i].vx;
            particles[i].y  += 0.5* dt * particles[i].vy;
            particles[i].z  += 0.5* dt * particles[i].vz;
        }
    }
	r->t+=dt/2.;
}
void reb_integrator_leapfrog_part2(struct reb_simulation* r){
	const int N = r->N;
	struct reb_particle* restrict const particles = r->particles;
	const double dt = r->dt;
    if (r->ri_whfast.kernel == REB_WHFAST_KERNEL_6ABA363){

        double b1 = 0.2621129352517028;
        reb_leapfrog_apply_C(r, b1, 0.,0); // does NOT recalculate accelerations
     
        double a1 = -0.0682610383918630;
        double a2 =  0.5-a1;
        reb_leapfrog_apply_A(r, a2);
        
        double b2 = 1.-2.*b1;
        double c2 = 0.0164011128160783;
        reb_leapfrog_apply_C(r, b2, c2*2.,1); // recalculates accelerations
        
        reb_leapfrog_apply_A(r, a2);
        
        
        reb_leapfrog_apply_C(r, b1, 0.,1); // recalculates accelerations
        //reb_leapfrog_apply_A(r, a1);
        r->ri_whfast.is_synchronized =0;
    
    }else{
        // Traditional leapfrog
        for (int i=0;i<N;i++){
            particles[i].vx += dt * particles[i].ax;
            particles[i].vy += dt * particles[i].ay;
            particles[i].vz += dt * particles[i].az;
            particles[i].x  += 0.5* dt * particles[i].vx;
            particles[i].y  += 0.5* dt * particles[i].vy;
            particles[i].z  += 0.5* dt * particles[i].vz;
        }
    }
	r->t+=dt/2.;
	r->dt_last_done = r->dt;
}
	
void reb_integrator_leapfrog_synchronize(struct reb_simulation* r){
    if (r->ri_whfast.kernel == REB_WHFAST_KERNEL_6ABA363){
        if (r->ri_whfast.is_synchronized == 0){
            double a1 = -0.0682610383918630;
            reb_leapfrog_apply_A(r, a1);
            double z[6] = { 0.07943288242455420, 0.02974829169467665, -0.7057074964815896, 0.3190423451260838, -0.2869147334299646, 0.};
            z[5] = -(z[0]+z[1]+z[2]+z[3]+z[4]);
            double y[6] = {1.3599424487455264, -0.6505973747535132, -0.033542814598338416, -0.040129915275115030, 0.044579729809902803, 0.};
            y[5] = -(y[0]+y[1]+y[2]+y[3]+y[4]);
            double v[6] = {-0.034841228074994859, 0.031675672097525204, -0.005661054677711889, 0.004262222269023640, 0.005, -0.005};
            for (int i=5;i>=0;i--){
                reb_leapfrog_apply_C(r, -y[i], -v[i]*2.,1); // recalculates accelerations
                reb_leapfrog_apply_A(r, -z[i]);
            }
            r->ri_whfast.is_synchronized =1;
        }
    }
}

void reb_integrator_leapfrog_reset(struct reb_simulation* r){
    r->ri_whfast.is_synchronized =1;
	// Do nothing.
}
