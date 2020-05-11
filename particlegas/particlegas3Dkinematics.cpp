/*
Created by Nick Crane
*/

/*
Contains all the function definitions necessary to update the:
Velocities, Positions
Of each particle in the simulation
*/

#include "../include/ParticleGas3D.hpp"
#define RANDOM (double)rand()/RAND_MAX // used for initially placing particles



/* uses Van der Waals approximation */
void UpdateParticleVelocities(PARTICLE** pArray, double r8, double r14)
{
    vec3d f, u; // force vector and direction vector
    double mag2, mag8, mag14;

    for (int i = 0; i < N - 1; i++) 
    {
        for (int j = i + 1; j < N; j++)
        {
            u  = pArray[j]->s - pArray[i]->s;
            mag2 = u.mag2();
            mag8 = mag2 * mag2 * mag2 * mag2; mag14 = mag8 * mag8 / mag2;
            f = u * ( -24 * ( ( 2 * r14 / mag14 ) - ( r8 / mag8 ) ) );
            pArray[i]->v += f * ( dt / pArray[i]->m );
            pArray[j]->v -= f * ( dt / pArray[j]->m );
        }
    }
}

/* uses the leapfrog algorithm */
void UpdateParticlePositions(PARTICLE** pArray, double& momentum)
{
    for (int i = 0; i < N; i++)
    {
        pArray[i]->s += pArray[i]->v * (dt / 2);\

        vec3dint zone = GetZone(pArray[i]->s);
        if (zone.i == 0 || zone.i == numZone - 1)
            CheckParticlePosition(pArray[i], momentum);
        if (zone.j == 0 || zone.j == numZone - 1)
            CheckParticlePosition(pArray[i], momentum);
        if (zone.k == 0 || zone.k == numZone - 1)
            CheckParticlePosition(pArray[i], momentum);
    }
}

/* checks to see if particles will leave box */
void CheckParticlePosition(PARTICLE* particle, double& momentum)
{
    double nx, ny, nz;
    for (int i = 0; i < N; i++)
    {
        nx = dabs(particle->s.i + particle->v.i * dt);
        ny = dabs(particle->s.j + particle->v.j * dt);
        nz = dabs(particle->s.k + particle->v.k * dt);

        if (nx > L / 2)
        {
            particle->v.i *= -1;
            momentum += dabs(particle->p().i) * 2;
        }
        if (ny > L / 2)
        {
            particle->v.j *= -1;
            momentum += dabs(particle->p().j) * 2;
        }
        if (nz > L / 2)
        {
            particle->v.k *= -1;
            momentum += dabs(particle->p().k) * 2;
        }
    }
}

/* calculates potential energy of the system */
void CalculatePotential(PARTICLE** pArray, double r6, double r12, double& potential)
{
    vec3d u; // direction vector
    double mag, mag2, mag6, mag12;
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            u  = pArray[j]->s - pArray[i]->s;
            mag = u.mag();
            mag2 = mag * mag; mag6 = mag2 * mag2 * mag2; mag12 = mag6 * mag6;
            potential += 4 * ( ( r12 / mag12 ) - ( r6 / mag6 ) );
        }
    }
    printf("pot = %f\n", potential);
}

/* initializes particles and gives them a random position and velocity */
void InitializeParticlePositionsRandomly(PARTICLE** pArray, double m, double r, double initialv)
{
    for (int i = 0; i < N; i++)
    {
        pArray[i] = new PARTICLE(
            m, r, 
            vec3d(RANDOM * L - L / 2, RANDOM * L - L / 2, RANDOM * L - L / 2),
            vec3d((2 * RANDOM - 1) * initialv, (2 * RANDOM - 1) * initialv, (2 * RANDOM - 1) * initialv),
            vec3d(0, 0, 0) );
    }
}




/* // trying to use std::async to thread the velocity update
void UpdateParticleAccelerations(PARTICLE** pArray, int i, int N, float dt, double r7, double r13)
{
    vec3d f, u;
    double mag, mag2, mag7, mag13;

    for (int j = i + 1; j < N; j++)
    {
        u  = pArray[i]->s - pArray[j]->s;
        mag = u.mag();
        mag2 = mag * mag; mag7 = mag2 * mag2 * mag2 * mag; mag13 = mag7 * mag7 / mag;
        f = u * (24 * ((r13 / mag13) - (r7 / mag7)));
        pArray[i]->a += f * (dt / pArray[i]->m);
        pArray[j]->a -= f * (dt / pArray[j]->m);
    }
}
*/