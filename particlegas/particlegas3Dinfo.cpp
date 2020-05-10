/* 
Created by Nick Crane
*/

/*
Contains all the function definitions necessary to update the:
Information
Of each particle in the simulation
*/

#include "../include/ParticleGas3D.hpp"



void UpdateInfo(PARTICLE** pArray, int ticks, double& momentum, double& potential)
{
    double PV = L * L * L * momentum / ( dt * ticks * L * L * 6 );
    double KE = 0, TEMP = 0;

    for (int i = 0; i < N; i++) // calculate kinetic energy
    {
        double v = pArray[i]->v.mag();
        KE += pArray[i]->ke();
    }
    TEMP = 2 * KE / 3;

    printf("!kT = %f, KE = %f, PV = %f, U = %f, U+KE = %f\n", TEMP, KE, PV, potential, potential + KE);

    momentum = potential = 0;
}