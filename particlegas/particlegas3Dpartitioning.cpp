/* 
Created by Nick Crane
*/

/*
Contains all the function definitions necessary to update the:
Partioning
Of each particle in the simulation
*/

#include "../include/ParticleGas3D.hpp"



// returns zone particle is in
vec3dint GetZone(vec3d pos)
{
    vec3dint zone;
    zone.i = (int)((pos.i + L / 2) * numZone / L);
    zone.j = (int)((pos.j + L / 2) * numZone / L);
    zone.k = (int)((pos.k + L / 2) * numZone / L);
    
    if (zone.i >= numZone)
        zone.i = numZone - 1;
    if (zone.i < 0)
        zone.i = 0;
    if (zone.j >= numZone)
        zone.j = numZone - 1;
    if (zone.j < 0)
        zone.j = 0;
    if (zone.k >= numZone)
        zone.k = numZone - 1;
    if (zone.k < 0)
        zone.k = 0;
    
    return zone;
}




/*
void CreateBoxes(int boxesperside, double l)
{
    for (int i = 0; i < boxesperside; i++)
        for (int j = 0; j < boxesperside; j++)
            for (int k = 0; k < boxesperside; k++)
            {

            }
}


void UpdateParticleBox(CUBE** cubes, PARTICLE** pArray, int n, int x, int y, int z)
{

}
*/