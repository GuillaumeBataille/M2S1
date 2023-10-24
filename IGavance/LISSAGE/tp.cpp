// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <GL/glut.h>
#include <float.h>
#include "src/Vec3.h"
#include "src/Camera.h"

//NEW
#define LAMBDA 0.330
#define MU -0.331

bool show_tri_quality = false;

float max_curvature;
float min_curvature;

float min_tshape;
float max_tshape;
//ENDNEW

enum DisplayMode
{
    WIRE = 0,
    SOLID = 1,
    LIGHTED_WIRE = 2,
    LIGHTED = 3
};

struct Triangle
{
    inline Triangle()
    {
        v[0] = v[1] = v[2] = 0;
    }
    inline Triangle(const Triangle &t)
    {
        v[0] = t.v[0];
        v[1] = t.v[1];
        v[2] = t.v[2];
    }
    inline Triangle(unsigned int v0, unsigned int v1, unsigned int v2)
    {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
    }
    unsigned int &operator[](unsigned int iv) { return v[iv]; }
    unsigned int operator[](unsigned int iv) const { return v[iv]; }
    bool operator==(const Triangle &other) const
    {
        return (v[0] == other.v[0] && v[1] == other.v[1] && v[2] == other.v[2]);
    }
    inline virtual ~Triangle() {}
    inline Triangle &operator=(const Triangle &t)
    {
        v[0] = t.v[0];
        v[1] = t.v[1];
        v[2] = t.v[2];
        return (*this);
    }
    // membres indices des sommets du triangle:
    unsigned int v[3];
};

struct Mesh
{
    std::vector<Vec3> vertices;                       // array of mesh vertices positions
    std::vector<Vec3> normals;                        // array of vertices normals useful for the display
    std::vector<Triangle> triangles;                  // array of mesh triangles
    std::vector<Vec3> triangle_normals;               // triangle normals to display face normals

    //NEW
    std::vector<Vec3> Laplacien; //Tableau de Laplacien 
    std::vector<Vec3> LaplacienWeighted;
    std::vector<Vec3> Gaussien;

    std::vector<float> vunicurvature;
    std::vector<float> vcurvature;
    std::vector<float> vgausscurvature;

    std::vector <float> tshape;

    std::vector<std::vector<unsigned int>> voisinage; // Liste d'id voisins de chaque vertex
    //ENDNEW

    // Compute face normals for the display
    void computeTrianglesNormals()
    {

        // A faire : implémenter le calcul des normales par face
        // Attention commencer la fonction par triangle_normals.clear();
        // Iterer sur les triangles

        // La normal du triangle i est le resultat du produit vectoriel de deux ses arêtes e_10 et e_20 normalisé (e_10^e_20)
        // L'arete e_10 est représentée par le vecteur partant du sommet 0 (triangles[i][0]) au sommet 1 (triangles[i][1])
        // L'arete e_20 est représentée par le vecteur partant du sommet 0 (triangles[i][0]) au sommet 2 (triangles[i][2])

        // Normaliser et ajouter dans triangle_normales

        triangle_normals.clear();
        for (unsigned int i = 0; i < triangles.size(); i++)
        {
            const Vec3 &e0 = vertices[triangles[i][1]] - vertices[triangles[i][0]];
            const Vec3 &e1 = vertices[triangles[i][2]] - vertices[triangles[i][0]];
            Vec3 n = Vec3::cross(e0, e1);
            n.normalize();
            triangle_normals.push_back(n);
        }
    }

    // Compute vertices normals as the average of its incident faces normals
    void computeVerticesNormals()
    {
        // Utiliser weight_type : 0 uniforme, 1 aire des triangles, 2 angle du triangle

        // A faire : implémenter le calcul des normales par sommet comme la moyenne des normales des triangles incidents
        // Attention commencer la fonction par normals.clear();
        // Initializer le vecteur normals taille vertices.size() avec Vec3(0., 0., 0.)
        // Iterer sur les triangles

        // Pour chaque triangle i
        // Ajouter la normal au triangle à celle de chacun des sommets en utilisant des poids
        // 0 uniforme, 1 aire du triangle, 2 angle du triangle

        // Iterer sur les normales et les normaliser
        normals.clear();
        normals.resize(vertices.size(), Vec3(0., 0., 0.));
        for (unsigned int i = 0; i < triangles.size(); i++)
        {
            for (unsigned int t = 0; t < 3; t++)
                normals[triangles[i][t]] += triangle_normals[i];
        }
        for (unsigned int i = 0; i < vertices.size(); i++)
            normals[i].normalize();
    }

    void computeNormals()
    {
        computeTrianglesNormals();
        computeVerticesNormals();
    }

    void addNoise()
    {
        for (unsigned int i = 0; i < vertices.size(); i++)
        {
            float factor = 0.03;
            const Vec3 &p = vertices[i];
            const Vec3 &n = normals[i];
            vertices[i] = Vec3(p[0] + factor * ((double)(rand()) / (double)(RAND_MAX)) * n[0], p[1] + factor * ((double)(rand()) / (double)(RAND_MAX)) * n[1], p[2] + factor * ((double)(rand()) / (double)(RAND_MAX)) * n[2]);
        }
    }
};

// Transformation made of a rotation and translation
struct Transformation
{
    Mat3 rotation;
    Vec3 translation;
};

bool contain(std::vector<unsigned int> const &i_vector, unsigned int element)
{
    for (unsigned int i = 0; i < i_vector.size(); i++)
    {
        if (i_vector[i] == element)
            return true;
    }
    return false;
}

void collect_one_ring(std::vector<Vec3> const &i_vertices,
                      std::vector<Triangle> const &i_triangles,
                      std::vector<std::vector<unsigned int>> &o_one_ring)
{
    o_one_ring.clear();
    o_one_ring.resize(i_vertices.size()); // one-ring of each vertex, i.e. a list of vertices with which it shares an edge
    // Parcourir les triangles et ajouter les voisins dans le 1-voisinage
    // Attention verifier que l'indice n'est pas deja present
    for (unsigned int i = 0; i < i_triangles.size(); i++)
    {
        // Tous les points opposés dans le triangle sont reliés
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if (j != k)
                {
                    if (!contain(o_one_ring[i_triangles[i][j]], i_triangles[i][k]))
                    {
                        o_one_ring[i_triangles[i][j]].push_back(i_triangles[i][k]);
                    }
                }
            }
        }
    }
}

// Input mesh loaded at the launch of the application
Mesh mesh;
std::vector<float> current_field; //Le champ qu'on veut représenter visuellement

bool display_normals;
bool display_smooth_normals;
bool display_mesh;



DisplayMode displayMode;
int weight_type;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 1600;
static unsigned int SCREENHEIGHT = 900;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX = 0, lastY = 0, lastZoom = 0;
static bool fullScreen = false;

// ------------------------------------
// File I/O
// ------------------------------------
bool saveOFF(const std::string &filename,
             std::vector<Vec3> const &i_vertices,
             std::vector<Vec3> const &i_normals,
             std::vector<Triangle> const &i_triangles,
             std::vector<Vec3> const &i_triangle_normals,
             bool save_normals = false)
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl;

    unsigned int n_vertices = i_vertices.size(), n_triangles = i_triangles.size();
    myfile << n_vertices << " " << n_triangles << " 0" << std::endl;

    for (unsigned int v = 0; v < n_vertices; ++v)
    {
        myfile << i_vertices[v][0] << " " << i_vertices[v][1] << " " << i_vertices[v][2] << " ";
        if (save_normals)
            myfile << i_normals[v][0] << " " << i_normals[v][1] << " " << i_normals[v][2] << std::endl;
        else
            myfile << std::endl;
    }
    for (unsigned int f = 0; f < n_triangles; ++f)
    {
        myfile << 3 << " " << i_triangles[f][0] << " " << i_triangles[f][1] << " " << i_triangles[f][2] << " ";
        if (save_normals)
            myfile << i_triangle_normals[f][0] << " " << i_triangle_normals[f][1] << " " << i_triangle_normals[f][2];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}

void openOFF(std::string const &filename,
             std::vector<Vec3> &o_vertices,
             std::vector<Vec3> &o_normals,
             std::vector<Triangle> &o_triangles,
             std::vector<Vec3> &o_triangle_normals,
             bool load_normals = true)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return;
    }

    std::string magic_s;

    myfile >> magic_s;

    if (magic_s != "OFF")
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        exit(1);
    }

    int n_vertices, n_faces, dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    o_vertices.clear();
    o_normals.clear();

    for (int v = 0; v < n_vertices; ++v)
    {
        float x, y, z;

        myfile >> x >> y >> z;
        o_vertices.push_back(Vec3(x, y, z));

        if (load_normals)
        {
            myfile >> x >> y >> z;
            o_normals.push_back(Vec3(x, y, z));
        }
    }

    o_triangles.clear();
    o_triangle_normals.clear();
    for (int f = 0; f < n_faces; ++f)
    {
        int n_vertices_on_face;
        myfile >> n_vertices_on_face;

        if (n_vertices_on_face == 3)
        {
            unsigned int _v1, _v2, _v3;
            myfile >> _v1 >> _v2 >> _v3;

            o_triangles.push_back(Triangle(_v1, _v2, _v3));

            if (load_normals)
            {
                float x, y, z;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back(Vec3(x, y, z));
            }
        }
        else if (n_vertices_on_face == 4)
        {
            unsigned int _v1, _v2, _v3, _v4;
            myfile >> _v1 >> _v2 >> _v3 >> _v4;

            o_triangles.push_back(Triangle(_v1, _v2, _v3));
            o_triangles.push_back(Triangle(_v1, _v3, _v4));
            if (load_normals)
            {
                float x, y, z;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back(Vec3(x, y, z));
            }
        }
        else
        {
            std::cout << "We handle ONLY *.off files with 3 or 4 vertices per face" << std::endl;
            myfile.close();
            exit(1);
        }
    }
}

// ------------------------------------
// Application initialization
// ------------------------------------
void initLight()
{
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f, -16.0f, -50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv(GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
}

void init()
{
    camera.resize(SCREENWIDTH, SCREENHEIGHT);
    initLight();
    glCullFace(GL_BACK);
    glDisable(GL_CULL_FACE);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    display_normals = false;
    display_mesh = true;
    display_smooth_normals = true;
    displayMode = LIGHTED;
}

// -----------------------TP CODING-----------------------//

//Fonction pour normalize une valeur via un min et max de l'ensemble
float normalizeValue(float &value, float min, float max) {
    return abs(value - min) / (max-min);
}

//Computing du vecteur de laplacien
void computeLaplacien()
{
    mesh.Laplacien.clear();
    mesh.Laplacien.resize(mesh.vertices.size());
    for (int i = 0; i < mesh.vertices.size(); i++) // pour tout les vertex
    {
        Vec3 current_v = mesh.vertices[i];
        Vec3 centroid = Vec3(0, 0, 0);
        float number_of_neighbor = 0.;                     // Nombre de voisin pour normaliser
        for (int j = 0; j < mesh.voisinage[i].size(); j++) // POur tout les voisins du vertex courant
        {
            Vec3 current_neighbour = mesh.vertices[mesh.voisinage[i][j]]; // Le voisin courant du vertex courant
            centroid += current_neighbour;
            number_of_neighbor++;
        }
        mesh.Laplacien[i] = (1.0 / (float)number_of_neighbor) * centroid - current_v; // Definition du vecteur vec_courant vers centroid
    }
}


//Computing du champ Laplacien uniforme normalisé (poids constant 1)
void calc_uniform_mean_curvature()
{
    mesh.vunicurvature.clear();
    mesh.vunicurvature.resize(mesh.Laplacien.size());
    current_field.resize(mesh.Laplacien.size());
    //std::cout << "Uniform mean curvature start " << std::endl;
    min_curvature = FLT_MAX;
    max_curvature = FLT_MIN;
    for (int i = 0; i < mesh.Laplacien.size(); i++)
    {
        mesh.vunicurvature[i] = 0.5 * mesh.Laplacien[i].length();
        max_curvature = std::max(mesh.vunicurvature[i], max_curvature);
        min_curvature = std::min(mesh.vunicurvature[i], min_curvature);
    }

    for (int i = 0; i < mesh.vunicurvature.size(); i++)
    {
        current_field[i] = normalizeValue(mesh.vunicurvature[i],min_curvature,max_curvature);
    }
}


//Calcul de la qualité des triangles en fonction du ratio rayon circonscrit / edge minimale
void calc_triangle_quality()
{
    mesh.tshape.resize(mesh.triangles.size());
    min_tshape = FLT_MAX;
    max_tshape = FLT_MIN;
    for (int i = 0; i < mesh.triangles.size(); i++)
    {
        Triangle current_tri = mesh.triangles[i];
        Vec3 a = mesh.vertices[current_tri.v[0]] - mesh.vertices[current_tri.v[1]];
        Vec3 b = mesh.vertices[current_tri.v[1]] - mesh.vertices[current_tri.v[2]];
        Vec3 c = mesh.vertices[current_tri.v[0]] - mesh.vertices[current_tri.v[2]];

        float a_l = a.length();
        float b_l = b.length();
        float c_l = c.length();

        float edge_lenght_min = std::min(a_l, std::min(b_l, c_l));
        Vec3 ab = Vec3::cross(a, b);
        // std::cout << "Taille de AB : " << ab.length() << std::endl;
        float s = (a_l + b_l + c_l) / 2;
        float area = sqrt(s * (s - a_l) * (s - b_l) * (s - c_l));
        float circumradius = (a_l * b_l * c_l) / (4 * area);

        float r = (a_l * b_l * c_l) / (2 * ab.length());

        if (ab.length() <= 0.0001) // Sil les edges
        {
            mesh.tshape[i] = 0;
        }
        else
            mesh.tshape[i] = r / edge_lenght_min;
        mesh.tshape[i] = std::min(mesh.tshape[i], (float)1);
        max_tshape = std::max(mesh.tshape[i], max_tshape);
        min_tshape = std::min(mesh.tshape[i], min_tshape);

        // std::cout << "Quality min : " << min_tshape << std::endl;
        // std::cout << "Quality max : " << max_tshape << std::endl;
    }

    for (int i = 0; i < mesh.tshape.size(); i++)
    {
        mesh.tshape[i] = normalizeValue(mesh.tshape[i],min_tshape,max_tshape); 
    }

}

//Application du smooth laplacien uniforme
void smoothing_mesh_uniform_mean()
{
    std::vector<Vec3> Temp_vertices;
    Temp_vertices.resize(mesh.vertices.size());

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        Temp_vertices[i] = mesh.vertices[i] + 0.5 * mesh.Laplacien[i];
    }

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        mesh.vertices[i] = Temp_vertices[i];
    }

    mesh.computeNormals();
    computeLaplacien();
    calc_uniform_mean_curvature();
    calc_triangle_quality();
}


//Application du smoothing de taubin uniforme (x2 Laplacien)
void smoothing_taubin(float lambda, float mu)
{
    std::vector<Vec3> Temp_vertices;
    Temp_vertices.resize(mesh.vertices.size());

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        Temp_vertices[i] = mesh.vertices[i] + lambda * mesh.Laplacien[i];
    }

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        mesh.vertices[i] = Temp_vertices[i];
    }
    mesh.computeNormals();
    computeLaplacien();
    calc_uniform_mean_curvature();
    calc_triangle_quality();

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        Temp_vertices[i] = mesh.vertices[i] + mu * mesh.Laplacien[i];
    }

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        mesh.vertices[i] = Temp_vertices[i];
    }
    mesh.computeNormals();
    computeLaplacien();
    calc_uniform_mean_curvature();
    calc_triangle_quality();
}



//Fonction qui recupère les angles opposés d'une arêtes 
std::vector<float> getOppositeAngles(unsigned int vertexIdx, unsigned int neighborIdx) {

    std::vector<float> oppositeAngles;
    size_t numTriangles = mesh.triangles.size();
    for (unsigned int i = 0; i < numTriangles; ++i) // On parcourt tout les triangles
    { // On va chercher a trouver les triangles qui contiennent notre edge courante ( 6 cas possibles car l'edge a 3 cotés possibles et elle peut etre dans un sens ou l'autre)
        
        if (mesh.triangles[i][0] == vertexIdx && mesh.triangles[i][1] == neighborIdx) { // Si l'edge courante est le 0-1 du triangle 0 1 2
            Vec3 oppositeVertex = mesh.vertices[mesh.triangles[i][2]]; // Le sommet opposé est 2 
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][0]] - oppositeVertex; // Edge vers le sommet opposé
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][1]] - oppositeVertex; // Edge vers le sommet opposé
            float cosAngle = Vec3::dot(edge1, edge2); //Cosinus sur l'angle du sommet opposé
            float angle = acos(cosAngle); // L'angle opposé
            oppositeAngles.push_back(angle); // Ajout de l'angle
        }
        else if (mesh.triangles[i][1] == vertexIdx && mesh.triangles[i][0] == neighborIdx) {
            Vec3 oppositeVertex = mesh.vertices[mesh.triangles[i][2]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][1]] - oppositeVertex;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][0]] - oppositeVertex;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            oppositeAngles.push_back(angle);
        }
        else if (mesh.triangles[i][0] == vertexIdx && mesh.triangles[i][2] == neighborIdx) {
            Vec3 oppositeVertex = mesh.vertices[mesh.triangles[i][1]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][0]] - oppositeVertex;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][2]] - oppositeVertex;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            oppositeAngles.push_back(angle);
        }
        else if (mesh.triangles[i][2] == vertexIdx && mesh.triangles[i][0] == neighborIdx) {
            Vec3 oppositeVertex = mesh.vertices[mesh.triangles[i][1]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][2]] - oppositeVertex;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][0]] - oppositeVertex;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            oppositeAngles.push_back(angle);
        }
        else if (mesh.triangles[i][1] == vertexIdx && mesh.triangles[i][2] == neighborIdx) {
            Vec3 oppositeVertex = mesh.vertices[mesh.triangles[i][0]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][1]] - oppositeVertex;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][2]] - oppositeVertex;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            oppositeAngles.push_back(angle);
        }
        else if (mesh.triangles[i][2] == vertexIdx && mesh.triangles[i][1] == neighborIdx) {
            Vec3 oppositeVertex = mesh.vertices[mesh.triangles[i][0]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][2]] - oppositeVertex;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][1]] - oppositeVertex;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            oppositeAngles.push_back(angle);
        }
        if (oppositeAngles.size() == 2) return oppositeAngles; // Dès qu'on a 2 angles opposés, il n'y en a pas plus a chercher
    }
    return oppositeAngles;
}



//Compute du laplacien avec poids sur angles opposés (voir cours)
void computeLaplacienWeighted()
{
    mesh.LaplacienWeighted.clear();
    mesh.LaplacienWeighted.resize(mesh.vertices.size());

    for (int i = 0; i < mesh.vertices.size(); i++) // pour tout les vertex
    {
        Vec3 current_v = mesh.vertices[i];
        Vec3 centroid = Vec3(0, 0, 0);
        float sumofw = 0.;
        for (int j = 0; j < mesh.voisinage[i].size(); j++) // POur tout les voisins du vertex courant
        {
            //Ici application du cours a partir des deux angles opposés
            std::vector<float> angles = getOppositeAngles(i, mesh.voisinage[i][j]);
            float alpha = angles[0];
            float beta = angles[1];

            float cotAlpha = cos(alpha) / sin(alpha);
            float cotBeta = cos(beta) / sin(beta);
            float weight = 0.5 * (cotAlpha + cotBeta);

            Vec3 current_neighbour = mesh.vertices[mesh.voisinage[i][j]]; // Le voisin courant du vertex courant
            centroid += weight * (current_neighbour - current_v); // Centroid pondéré par les poids
            sumofw+= weight;
        }
        mesh.LaplacienWeighted[i] = (1.0 / (float)sumofw) * centroid; //Valeur du laplacien avec les poids
    }

}

//Calcul du camp Laplacien avec poids 
void calc_mean_curvature()
{
    mesh.vcurvature.clear();
    mesh.vcurvature.resize(mesh.LaplacienWeighted.size());
    current_field.resize(mesh.LaplacienWeighted.size());

    min_curvature = FLT_MAX;
    max_curvature = FLT_MIN;
    for (int i = 0; i < mesh.Laplacien.size(); i++)
    {
        mesh.vcurvature[i] = 0.5 * mesh.LaplacienWeighted[i].length();
        max_curvature = std::max(mesh.vcurvature[i], max_curvature);
        min_curvature = std::min(mesh.vcurvature[i], min_curvature);
    }
    for (int i = 0; i < mesh.vcurvature.size(); i++)
    {
        current_field[i] = normalizeValue(mesh.vcurvature[i],min_curvature,max_curvature); // Normalisation
    }
}


//Application du smoothing Laplacien avec les poids
void smoothing_mesh_mean_with_weights()
{
    std::vector<Vec3> Temp_vertices;
    Temp_vertices.resize(mesh.vertices.size());

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        Temp_vertices[i] = mesh.vertices[i] + 0.5 * mesh.LaplacienWeighted[i];
    }

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        mesh.vertices[i] = Temp_vertices[i];
    }

    mesh.computeNormals();
    computeLaplacienWeighted();
    calc_mean_curvature();
    calc_triangle_quality();
}

//Application du smoothing Taubin avec poids (2x Laplacien avec poids)
void smoothing_taubin_weighted(float lambda, float mu)
{
    std::vector<Vec3> Temp_vertices;
    Temp_vertices.resize(mesh.vertices.size());

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        Temp_vertices[i] = mesh.vertices[i] + lambda * mesh.LaplacienWeighted[i];
    }

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        mesh.vertices[i] = Temp_vertices[i];
    }
    mesh.computeNormals();
    computeLaplacienWeighted();
    calc_mean_curvature();
    calc_triangle_quality();

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        Temp_vertices[i] = mesh.vertices[i] + mu * mesh.LaplacienWeighted[i];
    }

    for (int i = 0; i < mesh.vertices.size(); i++)
    {
        mesh.vertices[i] = Temp_vertices[i];
    }
    mesh.computeNormals();
    computeLaplacienWeighted();
    calc_mean_curvature();
    calc_triangle_quality();
}


//Fonction qui retourne la somme des angles des triangles autour d'un vertex en radian
float SumOfAngle(unsigned int vertexIndex) {

    float anglesSum = 0.0;
    for (unsigned int i = 0; i <mesh.triangles.size(); ++i)
    {
        //Même logique que le calcul d'angle opposé, mais pas besoin de vertex opposé
        if(mesh.triangles[i][0] == vertexIndex) { 
            Vec3 current = mesh.vertices[mesh.triangles[i][0]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][2]] - current;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][1]] - current;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            anglesSum += angle;
        }
        else if(mesh.triangles[i][1] == vertexIndex) {
            Vec3 current = mesh.vertices[mesh.triangles[i][1]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][0]] - current;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][2]] - current;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            anglesSum += angle;
        }
        else if(mesh.triangles[i][2] == vertexIndex) {
            Vec3 current = mesh.vertices[mesh.triangles[i][2]];
            Vec3 edge1 = mesh.vertices[mesh.triangles[i][0]] - current;
            Vec3 edge2 = mesh.vertices[mesh.triangles[i][1]] - current;
            float cosAngle = Vec3::dot(edge1, edge2);
            float angle = acos(cosAngle);
            anglesSum += angle;
        }
    }
    return anglesSum;
}

//Calcul du champ Gaussien qui se base sur les angles autour d'un vertex 
void calcGaussCurvature() {
    mesh.vgausscurvature.clear();
    mesh.Gaussien.clear();
    mesh.vgausscurvature.resize(mesh.vertices.size());
    mesh.Gaussien.resize(mesh.vertices.size());
    current_field.resize(mesh.vertices.size());

    min_curvature = FLT_MAX;
    max_curvature = FLT_MIN;

    for (unsigned int i = 0; i < mesh.vertices.size(); ++i) {

        float sumAngles = SumOfAngle(i);
        mesh.vgausscurvature[i] = 2 * M_PI - sumAngles;

        max_curvature = std::max(mesh.vgausscurvature[i], max_curvature);
        min_curvature = std::min(mesh.vgausscurvature[i], min_curvature);
    }       

    for(int i= 0; i < mesh.vertices.size(); i++)
    {
    current_field[i] = normalizeValue(mesh.vgausscurvature[i],min_curvature,max_curvature);
    }

}

//-----------------------TP CODING END-----------------//

// ------------------------------------
// Rendering.
// ------------------------------------

void drawVector(Vec3 const &i_from, Vec3 const &i_to)
{

    glBegin(GL_LINES);
    glVertex3f(i_from[0], i_from[1], i_from[2]);
    glVertex3f(i_to[0], i_to[1], i_to[2]);
    glEnd();
}

void drawAxis(Vec3 const &i_origin, Vec3 const &i_direction)
{

    glLineWidth(4); // for example...
    drawVector(i_origin, i_origin + i_direction);
}

void drawReferenceFrame(Vec3 const &origin, Vec3 const &i, Vec3 const &j, Vec3 const &k)
{

    glDisable(GL_LIGHTING);
    glColor3f(0.8, 0.2, 0.2);
    drawAxis(origin, i);
    glColor3f(0.2, 0.8, 0.2);
    drawAxis(origin, j);
    glColor3f(0.2, 0.2, 0.8);
    drawAxis(origin, k);
    glEnable(GL_LIGHTING);
}

typedef struct
{
    float r; // ∈ [0, 1]
    float g; // ∈ [0, 1]
    float b; // ∈ [0, 1]
} RGB;

RGB scalarToRGB(float scalar_value) // Scalar_value ∈ [0, 1]
{
    RGB rgb;
    float H = scalar_value * 360., S = 1., V = 0.85,
          P, Q, T,
          fract;

    (H == 360.) ? (H = 0.) : (H /= 60.);
    fract = H - floor(H);

    P = V * (1. - S);
    Q = V * (1. - S * fract);
    T = V * (1. - S * (1. - fract));

    if (0. <= H && H < 1.)
        rgb = (RGB){.r = V, .g = T, .b = P};
    else if (1. <= H && H < 2.)
        rgb = (RGB){.r = Q, .g = V, .b = P};
    else if (2. <= H && H < 3.)
        rgb = (RGB){.r = P, .g = V, .b = T};
    else if (3. <= H && H < 4.)
        rgb = (RGB){.r = P, .g = Q, .b = V};
    else if (4. <= H && H < 5.)
        rgb = (RGB){.r = T, .g = P, .b = V};
    else if (5. <= H && H < 6.)
        rgb = (RGB){.r = V, .g = P, .b = Q};
    else
        rgb = (RGB){.r = 0., .g = 0., .b = 0.};

    return rgb;
}

void drawSmoothTriangleMesh(Mesh const &i_mesh, bool draw_field = true)
{
    glBegin(GL_TRIANGLES);
    for (unsigned int tIt = 0; tIt < i_mesh.triangles.size(); ++tIt)
    {   
        RGB color;
        if(show_tri_quality)
            color = scalarToRGB(mesh.tshape[tIt]); // Version pour afficher les trianglequality en couleur
        for (unsigned int i = 0; i < 3; i++)
        {
            const Vec3 &p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; // Vertex position
            const Vec3 &n = i_mesh.normals[i_mesh.triangles[tIt][i]];  // Vertex normal

            if (draw_field && current_field.size() > 0)
            {
                if(!show_tri_quality)
                    color = scalarToRGB(current_field[i_mesh.triangles[tIt][i]]); // Version pour afficher les curvature en couleur
                glColor3f(color.r, color.g, color.b);
            }
            glNormal3f(n[0], n[1], n[2]);
            glVertex3f(p[0], p[1], p[2]);
        }
    }
    glEnd();
}

void drawTriangleMesh(Mesh const &i_mesh, bool draw_field = true)
{
    glBegin(GL_TRIANGLES);
    for (unsigned int tIt = 0; tIt < i_mesh.triangles.size(); ++tIt)
    {
        RGB color;
        if(show_tri_quality)
            color = scalarToRGB(mesh.tshape[tIt]); // Version pour afficher les trianglequality en couleur
        const Vec3 &n = i_mesh.triangle_normals[tIt]; // Triangle normal
        for (unsigned int i = 0; i < 3; i++)
        {
            const Vec3 &p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; // Vertex position

            if (draw_field)
            {
                 if(!show_tri_quality)
                    color = scalarToRGB(current_field[i_mesh.triangles[tIt][i]]); // Version pour afficher les curvature en couleur
            }
            glColor3f(color.r, color.g, color.b);
            glNormal3f(n[0], n[1], n[2]);
            glVertex3f(p[0], p[1], p[2]);
        }
    }
    glEnd();
}

void drawMesh(Mesh const &i_mesh, bool draw_field = true)
{
    if (display_smooth_normals)
        drawSmoothTriangleMesh(i_mesh, draw_field); // Smooth display with vertices normals
    else
        drawTriangleMesh(i_mesh, draw_field); // Display with face normals
}

void drawVectorField(std::vector<Vec3> const &i_positions, std::vector<Vec3> const &i_directions)
{
    glLineWidth(1.);
    for (unsigned int pIt = 0; pIt < i_directions.size(); ++pIt)
    {
        Vec3 to = i_positions[pIt] + 0.02 * i_directions[pIt];
        drawVector(i_positions[pIt], to);
    }
}

void drawNormals(Mesh const &i_mesh)
{

    if (display_smooth_normals)
    {
        drawVectorField(i_mesh.vertices, i_mesh.normals);
    }
    else
    {
        std::vector<Vec3> triangle_baricenters;
        for (const Triangle &triangle : i_mesh.triangles)
        {
            Vec3 triangle_baricenter(0., 0., 0.);
            for (unsigned int i = 0; i < 3; i++)
                triangle_baricenter += i_mesh.vertices[triangle[i]];
            triangle_baricenter /= 3.;
            triangle_baricenters.push_back(triangle_baricenter);
        }

        drawVectorField(triangle_baricenters, i_mesh.triangle_normals);
    }
}

// Draw fonction
void draw()
{

    if (displayMode == LIGHTED || displayMode == LIGHTED_WIRE)
    {

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_LIGHTING);
    }
    else if (displayMode == WIRE)
    {

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable(GL_LIGHTING);
    }
    else if (displayMode == SOLID)
    {
        glDisable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    glColor3f(0.8, 1, 0.8);
    drawMesh(mesh, true);

    if (displayMode == SOLID || displayMode == LIGHTED_WIRE)
    {
        glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glLineWidth(1.0f);
        glPolygonOffset(-2.0, 1.0);

        glColor3f(0., 0., 0.);
        drawMesh(mesh, false);

        glDisable(GL_POLYGON_OFFSET_LINE);
        glEnable(GL_LIGHTING);
    }

    glDisable(GL_LIGHTING);
    if (display_normals)
    {
        glColor3f(1., 0., 0.);
        drawNormals(mesh);
    }

    glEnable(GL_LIGHTING);
}

void changeDisplayMode()
{
    if (displayMode == LIGHTED)
        displayMode = LIGHTED_WIRE;
    else if (displayMode == LIGHTED_WIRE)
        displayMode = SOLID;
    else if (displayMode == SOLID)
        displayMode = WIRE;
    else
        displayMode = LIGHTED;
}

void display()
{
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply();
    draw();
    glFlush();
    glutSwapBuffers();
}

void idle()
{
    glutPostRedisplay();
}

// ------------------------------------
// User inputs
// ------------------------------------
// Keyboard event
void key(unsigned char keyPressed, int x, int y)
{
    switch (keyPressed)
    {
    case 'f':
        if (fullScreen == true)
        {
            glutReshapeWindow(SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        }
        else
        {
            glutFullScreen();
            fullScreen = true;
        }
        break;

    case 'w':
        changeDisplayMode();
        break;

    case 'n': // Press n key to display normals
        display_normals = !display_normals;
        break;

    case '1': // Toggle loaded mesh display
        display_mesh = !display_mesh;
        break;

    case 's': // Switches between face normals and vertices normals
        display_smooth_normals = !display_smooth_normals;
        break;

    case 'i': //Lissage Laplacien avec poids uniforme - Mesh reduit a chaque iteration
        smoothing_mesh_uniform_mean();
        break;

    case 'o': //Lissage de Taubin avec poids uniforme - Mesh équilibré
        smoothing_taubin(LAMBDA, MU);
        break;

    case 'p': //Affichage du la qualité du triangle au lieu de la valeur de     
        //smoothing_taubin(LAMBDA, MU);
        show_tri_quality = !show_tri_quality;
        break;

    case 'k': //Lissage Laplacien avec poids non uniforme - Mesh reduit a chaque iteration
        smoothing_mesh_mean_with_weights();
        break;

    case 'l': //Lissage de Taubin avec poids non uniforme - Mesh équilibré
        smoothing_taubin_weighted(LAMBDA,MU);
        break;

    case 'm'://Courbure gaussienne basée sur l'angle
        calcGaussCurvature();
        break;

    default:
        break;
    }
    idle();
}

// Mouse events
void mouse(int button, int state, int x, int y)
{
    if (state == GLUT_UP)
    {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    }
    else
    {
        if (button == GLUT_LEFT_BUTTON)
        {
            camera.beginRotate(x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        }
        else if (button == GLUT_RIGHT_BUTTON)
        {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        }
        else if (button == GLUT_MIDDLE_BUTTON)
        {
            if (mouseZoomPressed == false)
            {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }

    idle();
}

// Mouse motion, update camera
void motion(int x, int y)
{
    if (mouseRotatePressed == true)
    {
        camera.rotate(x, y);
    }
    else if (mouseMovePressed == true)
    {
        camera.move((x - lastX) / static_cast<float>(SCREENWIDTH), (lastY - y) / static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true)
    {
        camera.zoom(float(y - lastZoom) / SCREENHEIGHT);
        lastZoom = y;
    }
}

void reshape(int w, int h)
{
    camera.resize(w, h);
}

// ------------------------------------
// Start of graphical application
// ------------------------------------
int main(int argc, char **argv)
{
    if (argc > 2)
    {
        exit(EXIT_FAILURE);
    }
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow("TP HAI917I");

    init();
    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutReshapeFunc(reshape);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    key('?', 0, 0);

    // Mesh loaded with precomputed normals
    // openOFF("data/elephant_n.off", mesh.vertices, mesh.normals, mesh.triangles, mesh.triangle_normals);
    // mesh.computeNormals();

    // Version Bruitée
    openOFF("data/elephant_n.off", mesh.vertices, mesh.normals, mesh.triangles,
            mesh.triangle_normals);
    mesh.computeNormals();


    // TP execute start
    mesh.addNoise(); //Ajout de noise
    collect_one_ring(mesh.vertices, mesh.triangles, mesh.voisinage); // Compute tout les voisins d'un vertex
    computeLaplacien(); //Initialisation du Laplacien
    computeLaplacienWeighted();
    //calc_uniform_mean_curvature();
    //calc_mean_curvature();
    //calcGaussCurvature();
    calc_triangle_quality();
    //calc_weights();
    // computeLaplacien();
    //computeLaplacienWithWeights();
    //calc_uniform_mean_curvature();
    //calc_triangle_quality();
    //  TP execute end

    // A faire : normaliser les champs pour avoir une valeur flotante entre 0. et 1. dans current_field
    //***********************************************//

    // current_field.clear();

    glutMainLoop();
    return EXIT_SUCCESS;
}
