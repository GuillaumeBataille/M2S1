#ifndef VOXEL_H
#define VOXEL_H
#include "Vec3.h"

class Voxel
{
private:
    unsigned int size;

public:
    Voxel(unsigned int Size);
    ~Voxel();
};

Voxel::Voxel(unsigned int Size)
{
}

Voxel::~Voxel()
{
}

#endif