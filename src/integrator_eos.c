/**
 * @file    integrator_eos.c
 * @brief   Embedded Operator Splitting (EOS) method
 * @author  Hanno Rein
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
#include "integrator_eos.h"
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
static const double alpha1 = 1.35120719195965763404768780897;
static double alpha42 = 0.211324865405187117745425609749;
static const double a1_764 = 0.5600879810924619;
static const double b1_764 = 1.5171479707207228;
static const double a2_764 = -0.060087981092461900000;
static const double b2_764 = -2.0342959414414456000;
static const double z1_764 = -0.3346222298730800;
static const double z2_764 = 1.0975679907321640;
static const double z3_764 = -1.0380887460967830;
static const double z4_764 = 0.6234776317921379;
static const double z5_764 = -1.1027532063031910;
static const double z6_764 = -0.0141183222088869;
static const double y1_764 = -1.6218101180868010;
static const double y2_764 = 0.0061709468110142;
static const double y3_764 = 0.8348493592472594;
static const double y4_764 = -0.0511253369989315;
static const double y5_764 = 0.5633782670698199; 
static const double y6_764 = -0.5;
static const double a1_69 = 0.1867; //page 89, [176]
static const double a2_69 = 0.5554970237124784;
static const double a3_69 = 0.1294669489134754;
static const double a4_69 = -0.843265623387734;
static const double a5_69 = 0.9432033015235604;

static const double a1_864 = 0.0711334264982231177779387300061549964174;
static const double a2_864 = 0.241153427956640098736487795326289649618;
static const double a3_864 = 0.521411761772814789212136078067994229991;
static const double a4_864 = -0.333698616227678005726562603400438876027; // ABA(8,6,4)
static const double b1_864 = 0.183083687472197221961703757166430291072;
static const double b2_864 = 0.310782859898574869507522291054262796375;
static const double b3_864 = -0.0265646185119588006972121379164987592663;
static const double b4_864 = 0.0653961422823734184559721793911134363710; // ABA(8,6,4)

static const double a1_8 = 25./194.; 
static const double a2_8 = 0.581514087105251; 
static const double a3_8 = -0.410175371469850; 
static const double a4_8 = 0.1851469357165877; 
static const double a5_8 = -0.4095523434208514; 
static const double a6_8 = 0.1444059410800120; 
static const double a7_8 = 0.2783355003936797; 
static const double a8_8 = 0.3149566839162949; 
static const double a9_8 = -0.6269948254051343979; 
                
static inline void reb_integrator_eos_interaction_shell0(struct reb_simulation* r, double y, double v){
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    struct reb_particle* const particles = r->particles;
    const int testparticle_type   = r->testparticle_type;
    const double G = r->G;
    for (int i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
    
    // Normal force calculation 
    for (int i=1; i<N_active; i++){
        if (reb_sigint) return;
        for (int j=1; j<i; j++){
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = G/(dr*dr*dr);
            const double prefactj = -prefact*particles[j].m;
            particles[i].ax    += prefactj*dx;
            particles[i].ay    += prefactj*dy;
            particles[i].az    += prefactj*dz;
            const double prefacti = prefact*particles[i].m;
            particles[j].ax    += prefacti*dx;
            particles[j].ay    += prefacti*dy;
            particles[j].az    += prefacti*dz;
        }
    }
    for (int i=N_active; i<N; i++){
        if (reb_sigint) return;
        for (int j=1; j<N_active; j++){
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
        for (int i=1; i<N_active; i++){
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
    if (v!=0.){ // is jerk even used?
        for (int j=1; j<N_active; j++){
            if (reb_sigint) return;
            for (int i=1; i<j; i++){
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = v*G /(dr*dr*dr);
                const double prefact2i = prefact2*particles[i].m;
                const double prefact2j = prefact2*particles[j].m;
                particles[j].vx    -= dax*prefact2i;
                particles[j].vy    -= day*prefact2i;
                particles[j].vz    -= daz*prefact2i;
                particles[i].vx    += dax*prefact2j;
                particles[i].vy    += day*prefact2j;
                particles[i].vz    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2/dr *3./dr;
                const double prefact1i = prefact1*particles[i].m;
                const double prefact1j = prefact1*particles[j].m;
                particles[j].vx    += dx*prefact1i;
                particles[j].vy    += dy*prefact1i;
                particles[j].vz    += dz*prefact1i;
                particles[i].vx    -= dx*prefact1j;
                particles[i].vy    -= dy*prefact1j;
                particles[i].vz    -= dz*prefact1j;
            }
        }
        for (int i=N_active; i<N; i++){
            if (reb_sigint) return;
            for (int j=1; j<N_active; j++){
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 
                
                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = v*G /(dr*dr*dr);
                const double prefact2j = prefact2*particles[j].m;
                if (testparticle_type){
                    const double prefact2i = prefact2*particles[i].m;
                    particles[j].vx    -= dax*prefact2i;
                    particles[j].vy    -= day*prefact2i;
                    particles[j].vz    -= daz*prefact2i;
                }
                particles[i].vx    += dax*prefact2j;
                particles[i].vy    += day*prefact2j;
                particles[i].vz    += daz*prefact2j;
                const double prefact1 = alphasum*prefact2/dr*3./dr;
                const double prefact1j = prefact1*particles[j].m;
                if (testparticle_type){
                    const double prefact1i = prefact1*particles[i].m;
                    particles[j].vx    += dx*prefact1i;
                    particles[j].vy    += dy*prefact1i;
                    particles[j].vz    += dz*prefact1i;
                }
                particles[i].vx    -= dx*prefact1j;
                particles[i].vy    -= dy*prefact1j;
                particles[i].vz    -= dz*prefact1j;
            }
        }
    }
    for (int i=0;i<N;i++){
        particles[i].vx += y*particles[i].ax;
        particles[i].vy += y*particles[i].ay;
        particles[i].vz += y*particles[i].az;
    }
}
static inline void reb_integrator_eos_interaction_shell1(struct reb_simulation* r, double y, double v){
    const int N_active = r->N_active==-1?r->N:r->N_active;
    struct reb_particle* const particles = r->particles;
    const double G = r->G;
    
    if (v!=0.){ // is jerk even used?
            // Normal force calculation 
            particles[0].ax = 0;
            particles[0].ay = 0;
            particles[0].az = 0;
            for (int j=1; j<N_active; j++){
                const double dx = particles[0].x - particles[j].x;
                const double dy = particles[0].y - particles[j].y;
                const double dz = particles[0].z - particles[j].z;
                const double dr = sqrt(dx*dx + dy*dy + dz*dz);

                const double prefact = G/(dr*dr*dr);
                const double prefactj = -prefact*particles[j].m;
                particles[0].ax    += prefactj*dx;
                particles[0].ay    += prefactj*dy;
                particles[0].az    += prefactj*dz;
                const double prefacti = prefact*particles[0].m;
                particles[j].ax    = prefacti*dx;
                particles[j].ay    = prefacti*dy;
                particles[j].az    = prefacti*dz;
            }
            // Jerk calculation
            for (int i=1; i<N_active; i++){
                const double dx = particles[0].x - particles[i].x; 
                const double dy = particles[0].y - particles[i].y; 
                const double dz = particles[0].z - particles[i].z; 
                
                const double dax = particles[0].ax - particles[i].ax; 
                const double day = particles[0].ay - particles[i].ay; 
                const double daz = particles[0].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = v*G /(dr*dr*dr);
                const double prefact2i = prefact2*particles[i].m;
                const double prefact2j = prefact2*particles[0].m;
                const double prefact1 = alphasum*prefact2/dr *3./dr;
                const double prefact1i = prefact1*particles[i].m;
                const double prefact1j = prefact1*particles[0].m;
                particles[0].vx    += -dax*prefact2i + dx*prefact1i;
                particles[0].vy    += -day*prefact2i + dy*prefact1i;
                particles[0].vz    += -daz*prefact2i + dz*prefact1i;
                particles[i].vx    += y*particles[i].ax + dax*prefact2j - dx*prefact1j;
                particles[i].vy    += y*particles[i].ay + day*prefact2j - dy*prefact1j;
                particles[i].vz    += y*particles[i].az + daz*prefact2j - dz*prefact1j;
            }
            particles[0].vx += y*particles[0].ax;
            particles[0].vy += y*particles[0].ay;
            particles[0].vz += y*particles[0].az;
    }else{
        // Normal force calculation 
        for (int j=1; j<N_active; j++){
            const double dx = particles[0].x - particles[j].x;
            const double dy = particles[0].y - particles[j].y;
            const double dz = particles[0].z - particles[j].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz);

            const double prefact = y*G/(dr*dr*dr);
            const double prefactj = -prefact*particles[j].m;
            particles[0].vx    += prefactj*dx;
            particles[0].vy    += prefactj*dy;
            particles[0].vz    += prefactj*dz;
            const double prefacti = prefact*particles[0].m;
            particles[j].vx    += prefacti*dx;
            particles[j].vy    += prefacti*dy;
            particles[j].vz    += prefacti*dz;
        }
    }

}
static inline void reb_integrator_eos_preprocessor(struct reb_simulation* const r, double dt, int order, void (*drift_step)(struct reb_simulation* const r, double a), void (*interaction_step)(struct reb_simulation* const r, double y, double v)){
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
        case 764:
            drift_step(r, z1_764*dt);   
            interaction_step(r, y1_764*dt,0);
            drift_step(r, z2_764*dt);   
            interaction_step(r, y2_764*dt,0);
            drift_step(r, z3_764*dt);   
            interaction_step(r, y3_764*dt,0);
            drift_step(r, z4_764*dt);   
            interaction_step(r, y4_764*dt,0);
            drift_step(r, z5_764*dt);   
            interaction_step(r, y5_764*dt,0);
            drift_step(r, z6_764*dt);   
            interaction_step(r, y6_764*dt,0);
            break;
        case 2:
        default:
            break;
    }
}
static inline void reb_integrator_eos_postprocessor(struct reb_simulation* const r, double dt, int order, void (*drift_step)(struct reb_simulation* const r, double a), void (*interaction_step)(struct reb_simulation* const r, double y, double v)){
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
        case 764:
            interaction_step(r, -y6_764*dt,0);
            drift_step(r, -z6_764*dt);   
            interaction_step(r, -y5_764*dt,0);
            drift_step(r, -z5_764*dt);   
            interaction_step(r, -y4_764*dt,0);
            drift_step(r, -z4_764*dt);   
            interaction_step(r, -y3_764*dt,0);
            drift_step(r, -z3_764*dt);   
            interaction_step(r, -y2_764*dt,0);
            drift_step(r, -z2_764*dt);   
            interaction_step(r, -y1_764*dt,0);
            drift_step(r, -z1_764*dt);   
            break;
        case 2:
        default:
            break;
    }
}
static inline void reb_integrator_eos_drift_shell1(struct reb_simulation* const r, double a){
    struct reb_particle* restrict const particles = r->particles;
    unsigned int N = r->N;
    for (int i=0;i<N;i++){  
        particles[i].x += a*particles[i].vx;
        particles[i].y += a*particles[i].vy;
        particles[i].z += a*particles[i].vz;
    } 
}

void reb_integrator_eos_drift_shell0(struct reb_simulation* const r, double a){
    struct reb_simulation_integrator_eos* const rim = &(r->ri_eos);
    const int n = rim->n;
    const double as = a/n;
    reb_integrator_eos_preprocessor(r, as, rim->ordersubsteps, reb_integrator_eos_drift_shell1, reb_integrator_eos_interaction_shell1);
    switch(rim->ordersubsteps){
        case 6:
            for (int i=0;i<n;i++){
                reb_integrator_eos_drift_shell1(r, as*a_6[0]); //TODO combine drift steps
                reb_integrator_eos_interaction_shell1(r, as*b_6[0], as*as*as*c_6[0]*2.); 
                reb_integrator_eos_drift_shell1(r, as*a_6[1]);
                reb_integrator_eos_interaction_shell1(r, as*b_6[1], as*as*as*c_6[1]*2.); 
                reb_integrator_eos_drift_shell1(r, as*a_6[1]);
                reb_integrator_eos_interaction_shell1(r, as*b_6[0], as*as*as*c_6[0]*2.);
                reb_integrator_eos_drift_shell1(r, as*a_6[0]);
            }
            break;
        case 4:
            reb_integrator_eos_drift_shell1(r, as*0.5); 
            for (int i=0;i<n-1;i++){
                reb_integrator_eos_interaction_shell1(r, as, as*as*as/24.*2); 
                reb_integrator_eos_drift_shell1(r, as);
            }
            reb_integrator_eos_interaction_shell1(r, as, as*as*as/24.*2); 
            reb_integrator_eos_drift_shell1(r, as*0.5);
            break;
        case 40:
            reb_integrator_eos_drift_shell1(r, as*alpha1*0.5);
            for (int i=0;i<n-1;i++){
                reb_integrator_eos_interaction_shell1(r, as*alpha1, 0.);
                reb_integrator_eos_drift_shell1(r, as*(1.-alpha1)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*(1.-2.*alpha1), 0.);
                reb_integrator_eos_drift_shell1(r, as*(1.-alpha1)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*alpha1, 0.);
                reb_integrator_eos_drift_shell1(r, as*alpha1);
            }
            reb_integrator_eos_interaction_shell1(r, as*alpha1, 0.);
            reb_integrator_eos_drift_shell1(r, as*(1.-alpha1)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*(1.-2.*alpha1), 0.);
            reb_integrator_eos_drift_shell1(r, as*(1.-alpha1)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*alpha1, 0.);
            reb_integrator_eos_drift_shell1(r, as*alpha1*0.5);
            break;
        case 69:
            reb_integrator_eos_drift_shell1(r, as*a1_69*0.5);
            for (int i=0;i<n-1;i++){
                reb_integrator_eos_interaction_shell1(r, as*a1_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a1_69+a2_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a2_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a2_69+a3_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a3_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a3_69+a4_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a4_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a4_69+a5_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a5_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a4_69+a5_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a4_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a3_69+a4_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a3_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a2_69+a3_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a2_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a1_69+a2_69)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a1_69, 0.);
                reb_integrator_eos_drift_shell1(r, as*a1_69);
            }
            reb_integrator_eos_interaction_shell1(r, as*a1_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a1_69+a2_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a2_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a2_69+a3_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a3_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a3_69+a4_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a4_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a4_69+a5_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a5_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a4_69+a5_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a4_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a3_69+a4_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a3_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a2_69+a3_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a2_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a1_69+a2_69)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a1_69, 0.);
            reb_integrator_eos_drift_shell1(r, as*a1_69*0.5);
            break; 
        case 8: // SS_17^[8] p 91
            reb_integrator_eos_drift_shell1(r, as*a1_8*0.5);
            for (int i=0;i<n-1;i++){
                reb_integrator_eos_interaction_shell1(r, as*a1_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a1_8+a2_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a2_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a2_8+a3_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a3_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a3_8+a4_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a4_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a4_8+a5_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a5_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a5_8+a6_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a6_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a6_8+a7_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a7_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a7_8+a8_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a8_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a8_8+a9_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a9_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a8_8+a9_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a8_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a7_8+a8_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a7_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a6_8+a7_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a6_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a5_8+a6_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a5_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a4_8+a5_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a4_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a3_8+a4_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a3_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a2_8+a3_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a2_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*(a1_8+a2_8)*0.5);
                reb_integrator_eos_interaction_shell1(r, as*a1_8, 0.);
                reb_integrator_eos_drift_shell1(r, as*a1_8);
            }
            reb_integrator_eos_interaction_shell1(r, as*a1_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a1_8+a2_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a2_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a2_8+a3_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a3_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a3_8+a4_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a4_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a4_8+a5_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a5_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a5_8+a6_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a6_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a6_8+a7_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a7_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a7_8+a8_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a8_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a8_8+a9_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a9_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a8_8+a9_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a8_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a7_8+a8_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a7_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a6_8+a7_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a6_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a5_8+a6_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a5_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a4_8+a5_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a4_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a3_8+a4_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a3_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a2_8+a3_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a2_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*(a1_8+a2_8)*0.5);
            reb_integrator_eos_interaction_shell1(r, as*a1_8, 0.);
            reb_integrator_eos_drift_shell1(r, as*a1_8*0.5);
            break;
        case 2:
        default:
            reb_integrator_eos_drift_shell1(r, as*0.5); 
            for (int i=0;i<n-1;i++){
                reb_integrator_eos_interaction_shell1(r, as, 0.);
                reb_integrator_eos_drift_shell1(r, as);
            }
            reb_integrator_eos_interaction_shell1(r, as, 0.);
            reb_integrator_eos_drift_shell1(r, as*0.5);
            break;
    }
    reb_integrator_eos_postprocessor(r, as, rim->ordersubsteps, reb_integrator_eos_drift_shell1, reb_integrator_eos_interaction_shell1);
}

void reb_integrator_eos_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurius does not work with variational equations.");
    }
    
    r->gravity = REB_GRAVITY_NONE;

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurius only works with a direct collision search.");
    }
    
}


void reb_integrator_eos_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_eos* const rim = &(r->ri_eos);
    const double dt = r->dt;

    if (rim->is_synchronized){
        reb_integrator_eos_preprocessor(r, r->dt, rim->order, reb_integrator_eos_drift_shell0, reb_integrator_eos_interaction_shell0);
        switch(rim->order){
            case 6:
                reb_integrator_eos_drift_shell0(r, dt*a_6[0]); 
                break;
            case 4:
                reb_integrator_eos_drift_shell0(r, dt*0.5); 
                break;
            case 40:
                reb_integrator_eos_drift_shell0(r, dt*alpha1*0.5);
                break;
            case 42:
                reb_integrator_eos_drift_shell0(r, dt*alpha42); 
                break;
            case 764:
                reb_integrator_eos_drift_shell0(r, dt*a1_764);   
                break;
            case 864:
                reb_integrator_eos_drift_shell0(r, dt*a1_864);   
                break;
            case 2:
            default:
                reb_integrator_eos_drift_shell0(r, dt*0.5);
                break;
        }
    }else{
        switch(rim->order){
            case 6:
                reb_integrator_eos_drift_shell0(r, 2.*dt*a_6[0]); 
                break;
            case 4:
                reb_integrator_eos_drift_shell0(r, dt); 
                break;
            case 40:
                reb_integrator_eos_drift_shell0(r, dt*alpha1);
                break;
            case 42:
                reb_integrator_eos_drift_shell0(r, 2.*dt*alpha42); 
                break;
            case 764:
                reb_integrator_eos_drift_shell0(r, 2.*dt*a1_764);   
                break;
            case 864:
                reb_integrator_eos_drift_shell0(r, 2.*dt*a1_864);   
                break;
            case 2:
            default:
                reb_integrator_eos_drift_shell0(r, dt);
                break;
        }
    }
    switch(rim->order){
        case 6:
            reb_integrator_eos_interaction_shell0(r, dt*b_6[0], dt*dt*dt*c_6[0]*2.); 
            reb_integrator_eos_drift_shell0(r, dt*a_6[1]);
            reb_integrator_eos_interaction_shell0(r, dt*b_6[1], dt*dt*dt*c_6[1]*2.); 
            reb_integrator_eos_drift_shell0(r, dt*a_6[1]);
            reb_integrator_eos_interaction_shell0(r, dt*b_6[0], dt*dt*dt*c_6[0]*2.);
            break;
        case 4:
            reb_integrator_eos_interaction_shell0(r, dt, dt*dt*dt/24.*2); 
            break;
        case 40:
            reb_integrator_eos_interaction_shell0(r, dt*alpha1, 0.);
            reb_integrator_eos_drift_shell0(r, dt*(1.-alpha1)*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*(1.-2.*alpha1), 0.);
            reb_integrator_eos_drift_shell0(r, dt*(1.-alpha1)*0.5);
            reb_integrator_eos_interaction_shell0(r, dt*alpha1, 0.);
            break;
        case 42: // page 90 
            reb_integrator_eos_interaction_shell0(r, dt*0.5, 0.); 
            reb_integrator_eos_drift_shell0(r, dt*(1.-2.*alpha42));
            reb_integrator_eos_interaction_shell0(r, dt*0.5, 0.); 
            break;
        case 764:
            reb_integrator_eos_interaction_shell0(r, b1_764*dt,0);
            reb_integrator_eos_drift_shell0(r, a2_764*dt);   
            reb_integrator_eos_interaction_shell0(r, b2_764*dt,0);
            reb_integrator_eos_drift_shell0(r, a2_764*dt);   
            reb_integrator_eos_interaction_shell0(r, b1_764*dt,0);
            break;
        case 864:
            reb_integrator_eos_interaction_shell0(r, b1_864*dt,0);
            reb_integrator_eos_drift_shell0(r, a2_864*dt);   
            reb_integrator_eos_interaction_shell0(r, b2_864*dt,0);
            reb_integrator_eos_drift_shell0(r, a3_864*dt);   
            reb_integrator_eos_interaction_shell0(r, b3_864*dt,0);
            reb_integrator_eos_drift_shell0(r, a4_864*dt);   
            reb_integrator_eos_interaction_shell0(r, b4_864*dt,0);
            reb_integrator_eos_drift_shell0(r, a4_864*dt);   
            reb_integrator_eos_interaction_shell0(r, b3_864*dt,0);
            reb_integrator_eos_drift_shell0(r, a3_864*dt);   
            reb_integrator_eos_interaction_shell0(r, b2_864*dt,0);
            reb_integrator_eos_drift_shell0(r, a2_864*dt);   
            reb_integrator_eos_interaction_shell0(r, b1_864*dt,0);
            break;
        case 2:
        default:
            reb_integrator_eos_interaction_shell0(r, dt, 0.);
            break;
    }

    rim->is_synchronized = 0;
    if (rim->safe_mode){
        reb_integrator_eos_synchronize(r);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

void reb_integrator_eos_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_eos* const rim = &(r->ri_eos);
    const double dt = r->dt;
    if (rim->is_synchronized == 0){
        switch(rim->order){
            case 6:
                reb_integrator_eos_drift_shell0(r, dt*a_6[0]); 
                break;
            case 4:
                reb_integrator_eos_drift_shell0(r, dt*0.5); 
                break;
            case 40:
                reb_integrator_eos_drift_shell0(r, 0.5*dt*alpha1);
                break;
            case 42:
                reb_integrator_eos_drift_shell0(r, dt*alpha42); 
                break;
            case 764:
                reb_integrator_eos_drift_shell0(r, dt*a1_764);   
                break;
            case 864:
                reb_integrator_eos_drift_shell0(r, dt*a1_864);   
                break;
            case 2:
            default:
                reb_integrator_eos_drift_shell0(r, dt*0.5);
                break;
        }
        reb_integrator_eos_postprocessor(r, r->dt, rim->order, reb_integrator_eos_drift_shell0, reb_integrator_eos_interaction_shell0);
        rim->is_synchronized = 1;
    }
}

void reb_integrator_eos_reset(struct reb_simulation* r){
    r->ri_eos.n = 2;
    r->ri_eos.Phi0 = REB_EOS_LF;
    r->ri_eos.Phi1 = REB_EOS_LF;
    r->ri_eos.safe_mode = 1;
    r->ri_eos.is_synchronized = 1;
}

