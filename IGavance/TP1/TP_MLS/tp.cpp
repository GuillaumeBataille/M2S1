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
#include "src/jmkdtree.h"

std::vector<Vec3> positions; // Vector de positions
std::vector<Vec3> normals;   // Vector de normales

std::vector<Vec3> positions2; // Vector de positions secondaire
std::vector<Vec3> normals2;   // Vector de normales secondaire

std::vector<Vec3> positions3; // Vector de positions secondaire
std::vector<Vec3> normals3;   // Vector de normales secondaire

std::vector<Vec3> positions4; // Vector de positions secondaire
std::vector<Vec3> normals4;   // Vector de normales secondaire
// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

// Parametrage de l'environnement et de la fenetre
static GLint window;
static unsigned int SCREENWIDTH = 640;
static unsigned int SCREENHEIGHT = 480;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX = 0, lastY = 0, lastZoom = 0;
static bool fullScreen = false;
bool showmesh = true;
int mod = 0;
// ------------------------------------------------------------------------------------------------------------
// i/o and some stuff
// ------------------------------------------------------------------------------------------------------------
// Fonction pour charger un nuage de point (Pos + normals) via un path vers un .pn
void loadPN(const std::string &filename, std::vector<Vec3> &o_positions, std::vector<Vec3> &o_normals)
{
    unsigned int surfelSize = 6;
    FILE *in = fopen(filename.c_str(), "rb");
    if (in == NULL)
    {
        std::cout << filename << " is not a valid PN file." << std::endl;
        return;
    }
    size_t READ_BUFFER_SIZE = 1000; // for example...
    float *pn = new float[surfelSize * READ_BUFFER_SIZE];
    o_positions.clear();
    o_normals.clear();
    while (!feof(in))
    {
        unsigned numOfPoints = fread(pn, 4, surfelSize * READ_BUFFER_SIZE, in);
        for (unsigned int i = 0; i < numOfPoints; i += surfelSize)
        {
            o_positions.push_back(Vec3(pn[i], pn[i + 1], pn[i + 2]));
            o_normals.push_back(Vec3(pn[i + 3], pn[i + 4], pn[i + 5]));
        }

        if (numOfPoints < surfelSize * READ_BUFFER_SIZE)
            break;
    }
    fclose(in);
    delete[] pn;
}

// Fonction pour sauvegarder un nuage de point (Pos + normals) via un path vers un .pn
void savePN(const std::string &filename, std::vector<Vec3> const &o_positions, std::vector<Vec3> const &o_normals)
{
    if (o_positions.size() != o_normals.size())
    {
        std::cout << "The pointset you are trying to save does not contain the same number of points and normals." << std::endl;
        return;
    }
    FILE *outfile = fopen(filename.c_str(), "wb");
    if (outfile == NULL)
    {
        std::cout << filename << " is not a valid PN file." << std::endl;
        return;
    }
    for (unsigned int pIt = 0; pIt < o_positions.size(); ++pIt)
    {
        fwrite(&(o_positions[pIt]), sizeof(float), 3, outfile);
        fwrite(&(o_normals[pIt]), sizeof(float), 3, outfile);
    }
    fclose(outfile);
}

// Recentrer et aligner sur l'axe principal
void scaleAndCenter(std::vector<Vec3> &io_positions)
{
    Vec3 bboxMin(FLT_MAX, FLT_MAX, FLT_MAX);
    Vec3 bboxMax(FLT_MIN, FLT_MIN, FLT_MIN);
    // Calcul des bbboxmin et max dans le nuage de points
    for (unsigned int pIt = 0; pIt < io_positions.size(); ++pIt)
    {
        for (unsigned int coord = 0; coord < 3; ++coord)
        {
            bboxMin[coord] = std::min<float>(bboxMin[coord], io_positions[pIt][coord]);
            bboxMax[coord] = std::max<float>(bboxMax[coord], io_positions[pIt][coord]);
        }
    }
    Vec3 bboxCenter = (bboxMin + bboxMax) / 2.f; // Le centre du nuage de point
    float bboxLongestAxis = std::max<float>(bboxMax[0] - bboxMin[0], std::max<float>(bboxMax[1] - bboxMin[1], bboxMax[2] - bboxMin[2]));
    for (unsigned int pIt = 0; pIt < io_positions.size(); ++pIt)
    {
        io_positions[pIt] = (io_positions[pIt] - bboxCenter) / bboxLongestAxis;
    }
}

// Appliquer une matrice de rotation et de translation random sur les normales et les sommets
void applyRandomRigidTransformation(std::vector<Vec3> &io_positions, std::vector<Vec3> &io_normals)
{
    srand(time(NULL));
    Mat3 R = Mat3::RandRotation();
    Vec3 t = Vec3::Rand(1.f);
    for (unsigned int pIt = 0; pIt < io_positions.size(); ++pIt)
    {
        double x = (double)(rand()) / (double)(RAND_MAX);
        double y = (double)(rand()) / (double)(RAND_MAX);
        double z = (double)(rand()) / (double)(RAND_MAX);
        // t = Vec3(0.02*x,0.02*y,0.02*z);
        io_positions[pIt] = R * io_positions[pIt] + t;
        io_normals[pIt] = R * io_normals[pIt];
    }
}

// Subsampling par choix d'un pourcentage d'indices conservés
void subsample(std::vector<Vec3> &i_positions, std::vector<Vec3> &i_normals, float minimumAmount = 0.1f, float maximumAmount = 0.2f)
{
    std::vector<Vec3> newPos, newNormals;
    std::vector<unsigned int> indices(i_positions.size());
    for (unsigned int i = 0; i < indices.size(); ++i)
        indices[i] = i;
    srand(time(NULL));
    std::random_shuffle(indices.begin(), indices.end());
    unsigned int newSize = indices.size() * (minimumAmount + (maximumAmount - minimumAmount) * (float)(rand()) / (float)(RAND_MAX));
    newPos.resize(newSize);
    newNormals.resize(newSize);
    for (unsigned int i = 0; i < newPos.size(); ++i)
    {
        newPos[i] = i_positions[indices[i]];
        newNormals[i] = i_normals[indices[i]];
    }
    i_positions = newPos;
    i_normals = newNormals;
}

// Save la liste des sommets et la liste des triangles dans un file au path donné
bool save(const std::string &filename, std::vector<Vec3> &vertices, std::vector<unsigned int> &triangles)
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl;

    unsigned int n_vertices = vertices.size(), n_triangles = triangles.size() / 3;
    myfile << n_vertices << " " << n_triangles << " 0" << std::endl;

    for (unsigned int v = 0; v < n_vertices; ++v)
    {
        myfile << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << std::endl;
    }
    for (unsigned int f = 0; f < n_triangles; ++f)
    {
        myfile << 3 << " " << triangles[3 * f] << " " << triangles[3 * f + 1] << " " << triangles[3 * f + 2];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}

// Fonction qui projette un point sur un plan avec un poids uniforme (=1) (centroid,normal)
void project(Vec3 point_input, Vec3 &point_output, Vec3 &normal_output, Vec3 centroid, Vec3 normal)
{

    point_output = point_input - ((Vec3::dot((point_input - centroid), normal) / normal.length()) * normal);
    normal_output = normal;
}

// LES NOYAUX (les différents poids qu'on souhaite utiliser):

// Interpolante
double interpole(double r, double d, int s = 2)
{
    return pow(r / d, s);
}
// Gaussienne
double gaussienne(double r, double d)
{
    return (exp(-(d * d) / r * r));
}
// Wendland
double wendland(double r, double d)
{
    return ((1 - pow(d / 2, 4)) * (1 + (4 * (d / r))));
}

// Fonction qui projette un point sur un cercle
void project_to_circle(Vec3 center, double r, Vec3 &point_input, Vec3 &point_output, Vec3 &normal_output)
{
    Vec3 Center_to_Point = (point_input - center) / (point_input - center).length();
    point_output = center + r * Center_to_Point;
    normal_output = Center_to_Point;
}

void compute_weights(std::vector<double> &weight, std::vector<double> &weight_normalized, ANNidxArray id_nearest_neighbors, ANNdistArray square_distances_to_neighbors, int k, double radius = 0.2, int kernel_type = 0)
{
    weight.resize(k);
    weight_normalized.resize(k);
    double sum_w = 0;
    for (int j = 0; j < k; j++) // Pour chaque id_voisin
    {
        double distance = sqrt(square_distances_to_neighbors[j]); // La distance entre le voisin courant et le sommet initial
        weight[j] = gaussienne(radius, distance);
        sum_w += weight[j];
    }

    for (int j = 0; j < k; j++) // Pour chaque id_voisin
    {
        weight_normalized[j] = weight[j] / sum_w;
    }
}
//

// Compute le u4 a partir des positions du voisinage, des normales du maillages et des poids
double compute_u4(std::vector<Vec3> positions, std::vector<Vec3> normals, ANNidxArray id_nearest_neighbors, std::vector<double> weights, std::vector<double> weights_normalized, int k)
{
    double result;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double beta = 0.5;
    for (int j = 0; j < k; j++) // Pour chaque élément du voisinage
    {
        double weight = weights[j];                       // Le poids courant
        double weight_normalized = weights_normalized[j]; // Le poids normalisé courant
        Vec3 pos = positions[id_nearest_neighbors[j]];    // La position du voisin courante
        Vec3 normal = normals[id_nearest_neighbors[j]];   // La normale du voisin courante
        Vec3 pos_T = pos;
        sum1 += weight * (Vec3::dot(pos_T, normal));
        sum2 += Vec3::dot(weight_normalized * pos_T, weight * normal);
        sum3 += weight * (Vec3::dot(pos_T, pos));
        sum4 += Vec3::dot(weight_normalized * pos_T, weight * pos);
    }
    result = beta * ((sum1 - sum2) / (sum3 - sum4));

    return result;
}

Vec3 compute_u1u2u3(std::vector<Vec3> positions, std::vector<Vec3> normals, ANNidxArray id_nearest_neighbors, std::vector<double> weights_normalized, int k, double u4)
{
    Vec3 sum1(0, 0, 0);
    Vec3 sum2(0, 0, 0);
    Vec3 result;

    for (int j = 0; j < k; j++) // Pour chaque élément du voisinage
    {
        double weight_normalized = weights_normalized[j];
        Vec3 pos = positions[id_nearest_neighbors[j]];
        Vec3 normal = normals[id_nearest_neighbors[j]];
        sum1 += weight_normalized * normal;
        sum2 += weight_normalized * pos;
    }
    result = sum1 - 2 * u4 * sum2;
    return result;
}

double compute_u0(std::vector<Vec3> positions, ANNidxArray id_nearest_neighbors, std::vector<double> weights_normalized, int i, int k, Vec3 u1u2u3, double u4)
{
    double result;
    Vec3 sum1(0, 0, 0);
    double sum2 = 0.0;

    for (int j = 0; j < k; j++) // Pour chaque élément du voisinage
    {
        double weight_normalized = weights_normalized[j];
        Vec3 pos = positions[id_nearest_neighbors[j]];
        Vec3 pos_T = pos;
        sum1 += weight_normalized * pos;
        sum2 += weight_normalized * (Vec3::dot(pos_T, pos));
    }
    result = Vec3::dot(-1 * (u1u2u3), sum1) - u4 * sum2;
    return result;
}

void compute_circle(Vec3 &center, double &radius, double u0, Vec3 u1u2u3, double u4)
{
    center = (-1 * u1u2u3) / (2 * u4);
    radius = sqrt(pow(center.length(), 2) - (u0 / u4));
    // normal
}

void APSS(std::vector<Vec3> positions, std::vector<Vec3> normals, std::vector<Vec3> &positions2, std::vector<Vec3> &normals2, BasicANNkdTree const &kdtree, int kernel_type, float radius = 0.2, unsigned int k = 10)
{

    ANNidxArray id_nearest_neighbors = new ANNidx[k];
    ANNdistArray square_distances_to_neighbors = new ANNdist[k];
    double u4, u0;
    Vec3 u1u2u3;
    std::vector<double> weights;
    std::vector<double> weights_normalized;

    for (int i = 0; i < positions2.size(); i++) // Pour chaque point a projeter
    {
        kdtree.knearest(positions2[i], k, id_nearest_neighbors, square_distances_to_neighbors);

        compute_weights(weights, weights_normalized, id_nearest_neighbors, square_distances_to_neighbors, k);
        u4 = compute_u4(positions, normals, id_nearest_neighbors, weights, weights_normalized, k);
        u1u2u3 = compute_u1u2u3(positions, normals, id_nearest_neighbors, weights_normalized, k, u4);
        u0 = compute_u0(positions, id_nearest_neighbors, weights_normalized, i, k, u1u2u3, u4);
        if (sqrt(pow(u4, 2)) < 0.01) // Si u4 est très proche/egal a 0
        {
            // CAS PLAN
            Vec3 centroid(0, 0, 0);
            Vec3 normal(0, 0, 0);

            for (int l = 0; l < k ;l++)
            {   
                int id_neighbor = id_nearest_neighbors[l];
                
                centroid += weights_normalized[l] * positions[id_nearest_neighbors[l]];
                normal += weights_normalized[l] * normals[id_nearest_neighbors[l]];
            }
            project(positions2[i],positions2[i],normals2[i],centroid,normal);
            // TODO
        }
        else
        {
            Vec3 center;
            double radius;
            compute_circle(center, radius, u0, u1u2u3, u4);
            Vec3 position_output;
            Vec3 normal_ouput;
            project_to_circle(center, radius, positions2[i], positions2[i], normals2[i]);
        }
    }
    delete[] id_nearest_neighbors;
    delete[] square_distances_to_neighbors;
}

void HPSS(std::vector<Vec3> positions, std::vector<Vec3> normals, std::vector<Vec3> &positions2, std::vector<Vec3> &normals2, BasicANNkdTree const &kdtree, int kernel_type, float radius, unsigned int nbIterations = 1, unsigned int k = 20)
{

    ANNidxArray id_nearest_neighbors = new ANNidx[k];
    ANNdistArray square_distances_to_neighbors = new ANNdist[k];

    for (int l = 0; l < nbIterations; l++)
    {
        for (int i = 0; i < positions2.size(); i++) // Pour chaque point a projetter
        {
            kdtree.knearest(positions2[i], k, id_nearest_neighbors, square_distances_to_neighbors); // On compute les k voisins
            Vec3 centroid = Vec3(0, 0, 0);
            Vec3 normal = Vec3(0, 0, 0);
            double sum_w = 0.0;
            int weight;
            double max_dist = 0;

            for (int j = 0; j < k; j++) // Pour chaque id_voisin
            {
                int distance = sqrt(square_distances_to_neighbors[j]); // La distance entre le voisin courant et le sommet initial

                if (kernel_type == 0)
                {
                    weight = interpole(radius, distance);
                }

                if (kernel_type == 1)
                {

                    weight = gaussienne(radius, distance);
                }

                if (kernel_type == 2)
                {

                    weight = wendland(radius, distance);
                }
                sum_w += weight;

                centroid += weight * positions[id_nearest_neighbors[j]];
                normal += weight * normals[id_nearest_neighbors[j]];
            }
            // On normalise
            centroid /= sum_w;
            normal /= sum_w;
            project(positions2[i], positions2[i], normals2[i], centroid, normal); // Projeter sur le plan (sommet+normale)
        }
    }

    delete[] id_nearest_neighbors;
    delete[] square_distances_to_neighbors;
}
//

// ------------------------------------------------------------------------------------------------------------
// rendering.
// ------------------------------------------------------------------------------------------------------------

// Initialiser les lights
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

// Initialiser l'environnement du rendu
void init()
{
    camera.resize(SCREENWIDTH, SCREENHEIGHT);
    initLight();
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
}

// Dessiner des triangles a partir des vector de de positions et d'une liste d'indices de position (triangle)
void drawTriangleMesh(std::vector<Vec3> const &i_positions, std::vector<unsigned int> const &i_triangles)
{
    glBegin(GL_TRIANGLES);
    for (unsigned int tIt = 0; tIt < i_triangles.size() / 3; ++tIt)
    {
        Vec3 p0 = i_positions[3 * tIt];
        Vec3 p1 = i_positions[3 * tIt + 1];
        Vec3 p2 = i_positions[3 * tIt + 2];
        Vec3 n = Vec3::cross(p1 - p0, p2 - p0);
        n.normalize();
        glNormal3f(n[0], n[1], n[2]);
        glVertex3f(p0[0], p0[1], p0[2]);
        glVertex3f(p1[0], p1[1], p1[2]);
        glVertex3f(p2[0], p2[1], p2[2]);
    }
    glEnd();
}

// Dessiner a partir de liste de sommets et de normales (nuage de points)
void drawPointSet(std::vector<Vec3> const &i_positions, std::vector<Vec3> const &i_normals)
{
    glBegin(GL_POINTS);
    for (unsigned int pIt = 0; pIt < i_positions.size(); ++pIt)
    {
        glNormal3f(i_normals[pIt][0], i_normals[pIt][1], i_normals[pIt][2]);
        glVertex3f(i_positions[pIt][0], i_positions[pIt][1], i_positions[pIt][2]);
    }
    glEnd();
}

// La fonction draw qui dessiner les nuages de points primaire et secondaire
void draw()
{
    glPointSize(2); // for example...

    glColor3f(0.8, 0.8, 1);
    if (showmesh)
        drawPointSet(positions, normals); // Dessin du mesh initial via nuage de point

    glColor3f(1, 0.5, 0.5);
    if (mod == 0)
        drawPointSet(positions2, normals2);
    if (mod == 1)
        drawPointSet(positions3, normals3);
    if (mod == 2)
        drawPointSet(positions4, normals4);
}

// Affiche dans la frame le draw des nuages de points
void display()
{
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply();
    draw();
    glFlush();
    glutSwapBuffers();
}

// En attente de redraw, de redisplay
void idle()
{
    glutPostRedisplay();
}

// Gestion des inputs keyboard
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
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, polygonMode);
        if (polygonMode[0] != GL_FILL)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        else
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        break;

    case 's':
        showmesh = !showmesh;
        break;

    case 'd':
        mod = (mod + 1) % 3;
        if (mod == 0)
            std::cout << " Mode : Interpoled" << std::endl;
        if (mod == 1)
            std::cout << " Mode : Gaussienne" << std::endl;
        if (mod == 2)
            std::cout << " Mode : WendLand" << std::endl;
        break;

    default:
        break;
    }
    idle();
}

// Gestion de la souris (molette)
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

// Gestion de la souris (rotate when click maintained)
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

// Resize de la camera
void reshape(int w, int h)
{
    camera.resize(w, h);
}

// Main
int main(int argc, char **argv)
{
    if (argc > 2)
    {
        exit(EXIT_FAILURE);
    }
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow("tp point processing");

    init();
    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutReshapeFunc(reshape);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    key('?', 0, 0);

    {
        // Load a first pointset, and build a kd-tree:
        loadPN("pointsets/igea2.pn", positions, normals);
        // applyRandomRigidTransformation(positions,normals);
        BasicANNkdTree kdtree;
        kdtree.build(positions); // Construction du kd tree a partir des positions

        // Create a second pointset that is artificial, and project it on pointset1 using MLS techniques:

        // Le pointset 2 est donc composé d'autant de points que le premier, ils sont generé aléatoirement entre -.6 et 1.2

        positions2.resize(20000);
        normals2.resize(positions2.size());

        for (unsigned int pIt = 0; pIt < positions2.size(); ++pIt)
        {
            positions2[pIt] = Vec3(
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX),
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX),
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX));

            positions2[pIt].normalize();
            positions2[pIt] = 0.6 * positions2[pIt];
        }

        positions3.resize(20000);
        normals3.resize(positions3.size());

        for (unsigned int pIt = 0; pIt < positions3.size(); ++pIt)
        {
            positions3[pIt] = Vec3(
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX),
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX),
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX));

            positions3[pIt].normalize();
            positions3[pIt] = 0.6 * positions3[pIt];
        }

        positions4.resize(20000);
        normals4.resize(positions4.size());
        for (unsigned int pIt = 0; pIt < positions4.size(); ++pIt)
        {
            positions4[pIt] = Vec3(
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX),
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX),
                -0.6 + 1.2 * (double)(rand()) / (double)(RAND_MAX));

            positions4[pIt].normalize();
            positions4[pIt] = 0.6 * positions4[pIt];
        }
        // PROJECT USING MLS (HPSS and APSS):
        // TODO : Il faut faire se projetter les points secondaires sur la surface hypothétiquement definie par les sommet primaire
        // HPSS(positions, normals, positions2, normals2, kdtree, 0, 1, 3);
        // HPSS(positions, normals, positions3, normals3, kdtree, 1, 1, 3);
        // HPSS(positions, normals, positions4, normals4, kdtree, 2, 1, 3);
        APSS(positions, normals, positions2, normals2, kdtree, 0);
    }

    glutMainLoop();
    return EXIT_SUCCESS;
}
/*
Aucune diff entre les kernels?
Sommet random mais?
centroid et normales suffisent pour le calcul de poids

 */