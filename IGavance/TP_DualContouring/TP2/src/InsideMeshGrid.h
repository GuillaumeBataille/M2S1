#ifndef INSIDEMESHGRID_H
#define INSIDEMESHGRID_H

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>

#include "Vec3.h"



class InsideMeshGrid {
public:

	std::vector<Vec3> isoValue_pos;
	std::vector<int> isoValue;

	std::vector<Vec3> generated_vertex;
	std::vector<std::vector<int>> adjacencyList;
	std::vector<bool> validVertex;

	Vec3 bbmin, bbmax;
	float nbStep, stepSize; 
	int X,Y,Z,cX,cY,cZ;

	InsideMeshGrid(){

		bbmin = Vec3(0,0,0);
		bbmax = Vec3(1,1,1);


		stepSize = 1.0f/5.0f;
		
		X = 5;
		Y = 5;
		Z = 5;
		
		cX = X-1;
		cY = Y-1;
		cZ = Z-1;

		this->isoValue_pos = std::vector<Vec3> (X*Y*Z);
		this->isoValue = std::vector<int> (X*Y*Z);

		this->generated_vertex = std::vector<Vec3> (cX*cY*cZ);
		this->adjacencyList = std::vector<std::vector<int>> (cX*cY*cZ);
		this->validVertex = std::vector<bool>(cX*cY*cZ);

		std::fill(isoValue.begin(), isoValue.end(), 0);
		std::fill(adjacencyList.begin(), adjacencyList.end(), std::vector<int>());
		std::fill(validVertex.begin(), validVertex.end(), false);

		for(int x = 0; x < X; ++x){
			for(int y = 0; y < Y; ++y){
				for(int z = 0; z < Z; ++z){
					this->isoValue_pos[z + Z * (y + Y * x)] = Vec3(bbmin[0] + x*stepSize, bbmin[1] + y*stepSize, bbmin[2] + z*stepSize);
				}
			}
		}

		for(int x = 0; x < cX; ++x){
			for(int y = 0; y < cY; ++y){
				for(int z = 0; z < cZ; ++z){
					this->generated_vertex[z + cZ * (y + cY * x)] = Vec3(bbmin[0] + x*stepSize + stepSize/2,
															 	 bbmin[1] + y*stepSize + stepSize/2,
																 bbmin[2] + z*stepSize + stepSize/2);
				}
			}
		}

		isoValue_at(2,2,2) = 1;
	}

    InsideMeshGrid(size_t minPartitionNumber, std::vector<Vec3> & vertex)
	{
        computeBoundingBox(vertex);

		Vec3 boxSize = bbmax - bbmin;
		float minAxisSize = std::min(boxSize[0], std::min(boxSize[1], boxSize[2]));

		stepSize = minAxisSize/minPartitionNumber;
		
		int damp = 2;
		// +damp pour avoir toujours une bordure de point qui ne sont pas dans le mesh
		X = std::ceil(boxSize[0] / stepSize) + damp;
		Y = std::ceil(boxSize[1] / stepSize) + damp;
		Z = std::ceil(boxSize[2] / stepSize) + damp;
		
		cX = X-1;
		cY = Y-1;
		cZ = Z-1;

		this->isoValue_pos = std::vector<Vec3> (X*Y*Z);
		this->isoValue = std::vector<int> (X*Y*Z);

		this->generated_vertex = std::vector<Vec3> (cX*cY*cZ);
		this->adjacencyList = std::vector<std::vector<int>> (cX*cY*cZ);
		this->validVertex = std::vector<bool>(cX*cY*cZ);

		std::fill(isoValue.begin(), isoValue.end(), 0);
		std::fill(adjacencyList.begin(), adjacencyList.end(), std::vector<int>());
		std::fill(validVertex.begin(), validVertex.end(), false);

		Vec3 bbmin_damped = bbmin - Vec3(stepSize * damp/2, stepSize * damp/2, stepSize * damp/2);
		for(int x = 0; x < X; ++x){
			for(int y = 0; y < Y; ++y){
				for(int z = 0; z < Z; ++z){
					this->isoValue_pos[z + Z * (y + Y * x)] = Vec3(bbmin_damped[0] + x*stepSize, bbmin_damped[1] + y*stepSize, bbmin_damped[2] + z*stepSize);
				}
			}
		}

		for(int x = 0; x < cX; ++x){
			for(int y = 0; y < cY; ++y){
				for(int z = 0; z < cZ; ++z){
					this->generated_vertex[z + cZ * (y + cY * x)] = Vec3(bbmin_damped[0] + x*stepSize + stepSize/2,
															 	 bbmin_damped[1] + y*stepSize + stepSize/2,
																 bbmin_damped[2] + z*stepSize + stepSize/2);
				}
			}
		}
    }

    Vec3& operator()(size_t x, size_t y, size_t z) {
        if (x < X && y < Y && z < Z) {
            return isoValue_pos[z + Z * (y + Y * x)];
        }
        throw std::out_of_range("MeshGrid isoValue_pos index out of range.");
    }

    const Vec3& operator()(size_t x, size_t y, size_t z) const {
        if (x < X && y < Y && z < Z) {
            return isoValue_pos[z + Z * (y + Y * x)];
        }
        throw std::out_of_range("MeshGrid isoValue_pos index out of range.");
    }

	int& isoValue_at (size_t x, size_t y, size_t z) {
        if (x < X && y < Y && z < Z) {
            return isoValue[z + Z * (y + Y * x)];
        }
        throw std::out_of_range("MeshGrid isoValue index out of range.");
    }

    const int isoValue_at (size_t x, size_t y, size_t z) const {
        if (x < X && y < Y && z < Z) {
            return isoValue[z + Z * (y + Y * x)];
        }
        throw std::out_of_range("MeshGrid isoValue index out of range.");
    }

	std::string toString() const {
		std::string result = "MeshGrid (" + std::to_string(X) + "x" + std::to_string(Y) + "x" + std::to_string(Z) + "):\n";
		for (size_t x = 0; x < X; ++x) {
			for (size_t y = 0; y < Y; ++y) {
				for (size_t z = 0; z < Z; ++z) {
					const Vec3& vec = (*this)(x, y, z);
					result += "[" + std::to_string(x) + "][" + std::to_string(y) + "][" + std::to_string(z) + "]: (" +
							std::to_string(vec[0]) + ", " + std::to_string(vec[1]) + ", " + std::to_string(vec[2]) + ") isoValue: " + std::to_string(isoValue_at(x,y,z)) + "\n";
				}
			}
		}
		return result;
	}



	void createMesh(std::vector<Vec3> & vertex, std::vector<std::vector<int>> & faces)
	{
		vertex.clear();
		faces.clear();

        // créer les vertex
		for(int x = 0; x < cX; ++x)
		{
			for(int y = 0; y < cY; ++y)
            {
                for(int z = 0; z < cZ; ++z)
                {
                std::vector<int> currentIsoValue = getCellIsoValue(x,y,z);

                if(isDifferentIsoValue(currentIsoValue)){
                    validVertex[z + cZ * (y + cY * x)] = true;
                }
                }
            }
		}

        // créer les faces

		// explore les edges
		// on a pas besoin de check si il y a des vertex entre les edges car l'étapes d'avant les créé obligatoirement
        for(int x = 0; x < cX; ++x)
        {
            for(int y = 0; y < cY; ++y)
            {
                for(int z = 0; z < cZ; ++z)
                {
                    if(y != 0 && z != 0 && isoValue_at(x,y,z) != isoValue_at(x+1,y,z))
					{
						std::vector<int> face1(3), face2(3);

						if(isoValue_at(x,y,z))
						{
							face1[0] = z + cZ * (y + cY * x);
							face1[1] = z + cZ * ((y-1) + cY * x);
							face1[2] = (z-1) + cZ * (y + cY * x);

							face2[0] = z + cZ * ((y-1) + cY * x);
							face2[1] = (z-1) + cZ * ((y-1) + cY * x);
							face2[2] = (z-1) + cZ * (y + cY * x);
						} else
						{
							face1[0] = z + cZ * (y + cY * x);
							face1[2] = z + cZ * ((y-1) + cY * x);
							face1[1] = (z-1) + cZ * (y + cY * x);

							face2[0] = z + cZ * ((y-1) + cY * x);
							face2[2] = (z-1) + cZ * ((y-1) + cY * x);
							face2[1] = (z-1) + cZ * (y + cY * x);
						}
						

						faces.push_back(face1);
						faces.push_back(face2);
					}
					if(x != 0 && z != 0 && isoValue_at(x,y,z) != isoValue_at(x,y+1,z))
					{
						std::vector<int> face1(3), face2(3);

						if(isoValue_at(x,y,z))
						{
							face1[0] = z + cZ * (y + cY * x);
							face1[1] = (z-1) + cZ * (y + cY * x);
							face1[2] = z + cZ * (y + cY * (x-1));

							face2[0] = z + cZ * (y + cY * (x-1));
							face2[1] = (z-1) + cZ * (y + cY * x);
							face2[2] = (z-1) + cZ * (y + cY * (x-1));
						} else 
						{
							face1[0] = z + cZ * (y + cY * x);
							face1[2] = (z-1) + cZ * (y + cY * x);
							face1[1] = z + cZ * (y + cY * (x-1));

							face2[0] = z + cZ * (y + cY * (x-1));
							face2[2] = (z-1) + cZ * (y + cY * x);
							face2[1] = (z-1) + cZ * (y + cY * (x-1));
						}

						faces.push_back(face1);
						faces.push_back(face2);
					}
					if(x != 0 && y != 0 && isoValue_at(x,y,z) != isoValue_at(x,y,z+1))
					{
						std::vector<int> face1(3), face2(3);

						if(isoValue_at(x,y,z))
						{
							face1[0] = z + cZ * (y + cY * x);
							face1[2] = z + cZ * ((y-1) + cY * x);
							face1[1] = z + cZ * ((y-1) + cY * (x-1));

							face2[0] = z + cZ * (y + cY * x);
							face2[2] = z + cZ * ((y-1) + cY * (x-1));
							face2[1] = z + cZ * (y + cY * (x-1));
						}else
						{
							face1[0] = z + cZ * (y + cY * x);
							face1[1] = z + cZ * ((y-1) + cY * x);
							face1[2] = z + cZ * ((y-1) + cY * (x-1));

							face2[0] = z + cZ * (y + cY * x);
							face2[1] = z + cZ * ((y-1) + cY * (x-1));
							face2[2] = z + cZ * (y + cY * (x-1));
						}

						faces.push_back(face1);
						faces.push_back(face2);
					}
                }
            }
        }

		// mise en forme de la structure vertex, faces
		std::map<int, int> vertexNewId;

		int counter = 0;
		for(int i = 0; i < validVertex.size(); ++i)
		{
			if(validVertex[i]){
				vertex.push_back(generated_vertex[i]);
				vertexNewId[i] = counter;
				counter++;
			}
		}

		for(int i = 0; i < faces.size(); ++i)
		{
			faces[i][0] = vertexNewId[faces[i][0]];
			faces[i][1] = vertexNewId[faces[i][1]];
			faces[i][2] = vertexNewId[faces[i][2]];
		}
	}

private:

	// Fonction pour calculer le bbmin et bbmax à partir d'un vecteur de Vec3
	void computeBoundingBox(const std::vector<Vec3>& vertex) {
		if (vertex.empty()) {
			// Si le vecteur est vide, les valeurs de bbmin et bbmax restent à zéro.
			bbmin = Vec3(0, 0, 0);
			bbmax = Vec3(0, 0, 0);
			return;
		}

		// Initialisation des valeurs min et max avec la première valeur du vecteur.
		bbmin = vertex[0];
		bbmax = vertex[0];

		// Parcours du vecteur pour trouver les valeurs min et max.
		for (const Vec3& vec : vertex) {
			// Min
			if (vec[0] < bbmin[0]) bbmin[0] = vec[0];
			if (vec[1] < bbmin[1]) bbmin[1] = vec[1];
			if (vec[2] < bbmin[2]) bbmin[2] = vec[2];

			// Max
			if (vec[0] > bbmax[0]) bbmax[0] = vec[0];
			if (vec[1] > bbmax[1]) bbmax[1] = vec[1];
			if (vec[2] > bbmax[2]) bbmax[2] = vec[2];
		}
	}

	std::vector<int> getCellIsoValue(int x, int y, int z)
	{
		std::vector<int> is(8);
		is[0] = isoValue_at(x,y,z);
		is[1] = isoValue_at(x+1,y,z);
		is[2] = isoValue_at(x,y+1,z);
		is[3] = isoValue_at(x+1,y+1,z);
		is[4] = isoValue_at(x,y,z+1);
		is[5] = isoValue_at(x+1,y,z+1);
		is[6] = isoValue_at(x,y+1,z+1);
		is[7] = isoValue_at(x+1,y+1,z+1);

		return is;
	}

	// check si la cellule à des iso-valeurs différentes
	bool isDifferentIsoValue(std::vector<int> i){
		for(int j = 1; j<8; ++j){
			if(i[j] != i[0]){
				return true;
			}
		}
		return false;
	}

	// check si l'edge à bien 4 vertex autour de lui
	// side : 0 : x+1 | 1 : y+1 | 2 : z+1
	bool isIntersecting4Vertex(int x, int y, int z, int side){

		if(side == 0)
		{
			if(x == 0 || x == X){
				return false;
			}

			if(!validVertex[z + cZ * (y + cY * x)]){
				
			}
		}

		if(side == 1)
		{
			if(y == 0 || y == Y){
				return false;
			}
		}

		if(side == 2)
		{
			if(z == Z){
				return false;
			}
		}
	}


};

#endif