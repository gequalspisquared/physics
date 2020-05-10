#pragma once

#include "PhysicsObjects3D.hpp"

// glabal variables that will be set later
extern const int N; // number of particles
extern const int L; // length of box
extern const int numZone; // total number of zones
extern const float dt;

// particlegas3Dkinematics
void UpdateParticleVelocities(PARTICLE** pArray, double r8, double r14);
void UpdateParticlePositions(PARTICLE** pArray);
void FixParticlePosition(PARTICLE** pArray, double& momentum);
void CalculatePotential(PARTICLE** pArray, double r6, double r12, double& potential);

//particelgas3Dinfo
void UpdateInfo(PARTICLE** pArray, int ticks, double& momentum, double& potential);

//particlegas3Dpartitioning
void CreateBoxes();
void UpdateParticleBox();
void PlaceParticlesInBoxes();
void UpdateBoxParticle();
void ForcesInBoxes();



/*
class NODE 
{
public:
    int data;
    NODE* next;
};



class CUBE 
{
public:
    vec3d center;
    double length;
    int idx;
    int idy;
    int idz;
    int id;
    NODE* phead;

    double volume()
    {
        return length * length * length;
    }

    void updatelist(PARTICLE** pArray, NODE* n)
    {
        NODE* n = phead;

        while (n != NULL)
        {
            if (dabs(pArray[n->data]->s.i - length / 2) > center.i)
            {
                
            }
        }
    }

    CUBE()
    {
        center = vec3d();
        length = 0;
        idx = idy = idz = 0;
        phead->data = NULL;
    }

    CUBE(vec3d pos, double l, int x, int y, int z)
    {
        center = pos;
        length = l;
        idx = x; idy = y; idz = z;
        phead->data = NULL;
    }

    ~CUBE() {}
};
*/


