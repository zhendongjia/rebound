/**
 * A string of solid spheres bouncing
 * 
 * This example tests collision detection methods.
 * The example uses a non-square, rectangular box. 10 particles are placed
 * along a line. All except one of the particles are at rest initially.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r){
}

int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_create_simulation();
    // Setup constants
    r->dt             = 1e-0;
    r->integrator        = REB_INTEGRATOR_LEAPFROG;
    r->boundary        = REB_BOUNDARY_OPEN;
    r->collision        = REB_COLLISION_LINETREE;
    r->collision_resolve = reb_collision_resolve_merge;
    r->gravity        = REB_GRAVITY_TREE;
    r->heartbeat = heartbeat;
    r->usleep        = 50000;            // Slow down integration (for visualization only)
    
    reb_configure_box(r,100.,1,1,1);  // boxsize 10., three root boxes in x direction, one in y and z
    r->nghostx = 0; 
    r->nghosty = 0; 
    r->nghostz = 0;

    // Initial conditions
    int N = 1000;
    for(int i=0;i<N;i++){
        struct reb_particle p = {0};
        p.x  = r->boxsize.x*reb_random_uniform(-0.5,0.5); 
        p.y  = r->boxsize.x*reb_random_uniform(-0.5,0.5); 
        p.z  = r->boxsize.x*reb_random_uniform(-0.5,0.5); 
        p.vx  = 0.1*reb_random_uniform(-0.5,0.5); 
        p.vy  = 0.1*reb_random_uniform(-0.5,0.5); 
        p.vz  = 0.1*reb_random_uniform(-0.5,0.5); 
        p.m  = 100./N;
        p.r  = pow(p.m,1./3.);
        reb_add(r, p);
    }

    reb_integrate(r,INFINITY);
}
