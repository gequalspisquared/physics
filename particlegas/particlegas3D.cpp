/* 
A particle based N-Body gas sim
Created by Nick Crane
*/

/* 
Assumes neutral particles, uses Van der Waals to approximate forces
Can be changed to include gravity
*/

// #include <future> // async????
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include "../include/ParticleGas3D.hpp"

#define RUNANIM 0 // used to determine whether to run anim or not

/* const global paramters */
const int N = 1;
const int L = 100;
const int numZone = 100; // number of zones per side of box
const float dt = 0.01f;

int main() 
{
    /* parameters */
    static const double r1 = 1; // critical distance where force is 0
    int ticksPerInfoUpdate = 100;
    int ticksPerFrame = 10;

    /* particle properties */
    double mass = 1, radius = 1;
    double initialv = 1;

    /* derived quantities (do not touch) */
    static const double r2 = r1 * r1;
    static const double r6 = r2 * r2 * r2;
    static const double r7 = r6 * r1;
    static const double r8 = r6 * r2;
    static const double r12 = r6 * r6;
    static const double r13 = r12 * r1;
    static const double r14 = r13 * r1;

    /* variables used by program (do not touch) */
    double momentum = 0, potential = 0;
    float time = 0.0f;
    int tick = 0;

    PARTICLE** particles = new PARTICLE*[N] { 0 };
    // CUBE** cubes = new CUBE*[BoxesPerSide] { 0 };

    InitializeParticlePositionsRandomly(particles, mass, radius, initialv);

    while (time < 10.0f)
    {
        // trying to implement std::async to thread the velocity update (most expensive function)
        // for (int i = 0; i < N - 1; i++)
        // {
        //     p_Futures.push_back(std::async(std::launch::async, UpdateParticleVelocities, particles, i, N, dt, r8, r14));
        // }
        
        /* leapfrog algorithm */
        UpdateParticlePositions(particles, momentum);
        UpdateParticleVelocities(particles, r8, r14);
        UpdateParticlePositions(particles, momentum);
        time += dt; tick++;

        if (tick % ticksPerInfoUpdate == 0)
        {
            // CalculatePotential(particles, r6, r12, potential);
            // UpdateInfo(particles, ticksPerInfoUpdate, momentum, potential);
            // printf("!v = %f, %f, %f\n", particles[0]->v.i, particles[0]->v.j, particles[0]->v.k);
        }

#if RUNANIM==1
        if (tick % ticksPerFrame == 0)
        {
            for(int i = 0; i < N; i++)
                {
                    printf("c3 %f %f %f 1\n", particles[i]->s.i, particles[i]->s.j, particles[i]->s.k);
                }
                printf("F\n");
        }
#endif
    }

    for (int i = 0; i < N; i++)
    {
        delete particles[i];
    }
    delete[] particles;

#if RUNANIM==1
    printf("Q\n");
#endif

    return 0;
}