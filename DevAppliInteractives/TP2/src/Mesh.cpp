#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <cmath>

void Mesh::loadOFF(const std::string &filename)
{
    std::ifstream in(filename.c_str());
    if (!in)
        exit(EXIT_FAILURE);
    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    V.resize(sizeV);
    T.resize(sizeT);
    for (unsigned int i = 0; i < sizeV; i++)
        in >> V[i].p;
    int s;
    for (unsigned int i = 0; i < sizeT; i++)
    {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> T[i].v[j];
    }
    in.close();
    recomputeNormals();
}

Vec3 Mesh::projectPointSegment(Vec3 a, Vec3 b, Vec3 c)
{
    double r = Vec3::dot(b - a, b - a);
    if (fabs(r) < 1e-12)
        return a;
    r = Vec3::dot(c - a, b - a) / r;
    if (r < 0)
        return a;
    if (r > 1)
        return b;
    return a + r * (b - a);
}

void Mesh::recomputeNormals()
{
    for (unsigned int i = 0; i < V.size(); i++)
        V[i].n = Vec3(0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < T.size(); i++)
    {
        Vec3 e01 = V[T[i].v[1]].p - V[T[i].v[0]].p;
        Vec3 e02 = V[T[i].v[2]].p - V[T[i].v[0]].p;
        Vec3 n = Vec3::cross(e01, e02);
        n.normalize();
        for (unsigned int j = 0; j < 3; j++)
            V[T[i].v[j]].n += n;
    }
    for (unsigned int i = 0; i < V.size(); i++)
        V[i].n.normalize();
}

void Mesh::compute_skinning_weights(Skeleton &skeleton)
{
    for (size_t i = 0; i < this->V.size(); ++i)
    {
        MeshVertex &vertex = this->V[i];
        // Suppression des vertex précédents
        vertex.w.clear();
        double weights = 0.;
        for (size_t j = 0; j < skeleton.bones.size(); ++j)
        {
            // Calcul du poids
            const Bone &bone = skeleton.bones[j];
            // Les points par lesquels passe l'os
            Vec3 a = skeleton.articulations[bone.joints[0]].p;
            Vec3 b = skeleton.articulations[bone.joints[1]].p;
            Vec3 c = projectPointSegment(a, b, vertex.p); // Projection du vertex courant sur l'os
            double dist = (vertex.p - c).length();        // A quel point le vertex est loin de la projection de l'os
            double w = pow(1 / dist, 5);
            vertex.w.push_back(w);
            weights += w;
        }
        // Normalisation
        for (size_t j = 0; j < vertex.w.size(); ++j)
        {
            vertex.w[j] /= weights;
        }
    }
}

void Mesh::draw(int displayedBone) const
{

    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < T.size(); i++)
        for (unsigned int j = 0; j < 3; j++)
        {
            const MeshVertex &v = V[T[i].v[j]];
            if (displayedBone >= 0 && displayedBone < v.w.size())
            {
                Vec3 rgb = HSVtoRGB(v.w[displayedBone], 0.8, 0.8);
                glColor3f(rgb[0], rgb[1], rgb[2]);
            }
            else
                glColor3f(0.6, 0.6, 0.6);
            glNormal3f(v.n[0], v.n[1], v.n[2]);
            glVertex3f(v.p[0], v.p[1], v.p[2]);
        }

    glEnd();
}

void Mesh::drawTransformedMesh(SkeletonTransformation &transfo) const
{
    std::vector<Vec3> new_positions(this->V.size()); // Liste des nouvelles pos
    std::vector<Vec3> new_normals(this->V.size());   // Liste des normales

    for (size_t i = 0; i < this->V.size(); ++i) // Pour tout les vertex
    {
        Vec3 p = this->V[i].p;
        Vec3 n = this->V[i].n;
        new_positions[i] = Vec3(0, 0, 0);
        new_normals[i] = Vec3(0, 0, 0);
        for (size_t j = 0; j < transfo.bone_transformations.size(); ++j) // Pour chaque transformation qu'on applique
        {
            BoneTransformation &bone = transfo.bone_transformations[j]; // On recupère la transformation de l'os courant
            // wij * (r * p + t)
            // La nouvelle pos = position précedente + poids * (transfo_rota_de_l'os * position de base + transfo_transla_de_l'os)
            new_positions[i] += this->V[i].w[j] *
                                (bone.world_space_rotation * p + bone.world_space_translation);

            Vec3 transformed_normal = (Mat3::inverse(bone.world_space_rotation).getTranspose() * n);
            new_normals[i] += this->V[i].w[j] * transformed_normal;
        }
    }
    // Draw
    glEnable(GL_LIGHTING);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < T.size(); i++)
        for (size_t j = 0; j < 3; j++)
        {
            const MeshVertex &v = this->V[T[i].v[j]];
            Vec3 p = new_positions[T[i].v[j]];
            Vec3 n = new_normals[T[i].v[j]];
            glNormal3f(n[0], n[1], n[2]);
            glVertex3f(p[0], p[1], p[2]);
        }
    glEnd();
}

/*! \brief Convert HSV to RGB color space

  Converts a given set of HSV values `h', `s', `v' into RGB
  coordinates. The output RGB values are in the range [0, 1], and
  the input HSV values are in the ranges h = [0, 360], and s, v =
  [0, 1], respectively.

  \param fH Hue component, used as input, range: [0, 1]
  \param fS Hue component, used as input, range: [0, 1]
  \param fV Hue component, used as input, range: [0, 1]

  \param fR Red component, used as output, range: [0, 1]
  \param fG Green component, used as output, range: [0, 1]
  \param fB Blue component, used as output, range: [0, 1]

*/
Vec3 Mesh::HSVtoRGB(float fH, float fS, float fV) const
{

    fH = (1. - fH) * 0.65 * 360.;

    float fR, fG, fB;
    float fC = fV * fS; // Chroma
    float fHPrime = fmod(fH / 60.0, 6);
    float fX = fC * (1 - fabs(fmod(fHPrime, 2) - 1));
    float fM = fV - fC;

    if (0 <= fHPrime && fHPrime < 1)
    {
        fR = fC;
        fG = fX;
        fB = 0;
    }
    else if (1 <= fHPrime && fHPrime < 2)
    {
        fR = fX;
        fG = fC;
        fB = 0;
    }
    else if (2 <= fHPrime && fHPrime < 3)
    {
        fR = 0;
        fG = fC;
        fB = fX;
    }
    else if (3 <= fHPrime && fHPrime < 4)
    {
        fR = 0;
        fG = fX;
        fB = fC;
    }
    else if (4 <= fHPrime && fHPrime < 5)
    {
        fR = fX;
        fG = 0;
        fB = fC;
    }
    else if (5 <= fHPrime && fHPrime < 6)
    {
        fR = fC;
        fG = 0;
        fB = fX;
    }
    else
    {
        fR = 0;
        fG = 0;
        fB = 0;
    }

    fR += fM;
    fG += fM;
    fB += fM;
    return Vec3(fR, fG, fB);
}
