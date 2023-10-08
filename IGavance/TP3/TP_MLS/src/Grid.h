#ifndef GRID_H
#define GRID_H
#include "Vec3.h"

class Grid
{
private:
    double size;
    unsigned int step;

    std::vector<std::vector<std::vector<Vec3>>> Coord;

public:
    Grid(double Size, unsigned int step);
    ~Grid();
    void draw();
};

Grid::Grid(double Size, unsigned int step)
{
    std::cout << "Created grid " << std::endl;
    this->size = Size;
    this->step = step;

    double step_size = (float)Size / (float)step;

    Coord.resize(step); // X

    // Pour chaque axe x y z : Resize
    for (int y = 0; y < step; y++)
    {
        Coord[y].resize(step); // Y

        for (int z = 0; z < step; z++)
        {
            Coord[y][z].resize(step); // Z
        }
    }

    // Pour chaque position X Y Z : On remplit les données
    for (int x = 0; x < step; x++)
    {
        for (int y = 0; y < step; y++)
        {
            for (int z = 0; z < step; z++)
            {
                Coord[x][y][z] = Vec3(x * step_size, y * step_size, z * step_size);
            }
        }
    }
}

void Grid::draw()
{
    glPointSize(5.0f); // Définir la taille des points, ajustez selon vos besoins

    glBegin(GL_POINTS); // Commencer à dessiner des points

    for (int x = 0; x < this->step; x++)
    {
        for (int y = 0; y < step; y++)
        {
            for (int z = 0; z < step; z++)
            {
                Vec3 vertex = Coord[x][y][z];
                std::cout << "X:" << x << " - Y: " << y << " - Z: " << z << std::endl;

                glVertex3f(vertex[0], vertex[1], vertex[2]); // Dessiner le sommet
            }
        }
    }

    glEnd(); // Terminer de dessiner des points
}

Grid::~Grid()
{
}

#endif