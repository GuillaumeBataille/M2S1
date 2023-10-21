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
#include <random>

#include <algorithm>
#include <GL/glut.h>
#include <float.h>
#include "src/Vec3.h"
#include "src/Camera.h"
#include "src/jmkdtree.h"
#include "src/InsideMeshGrid.h"


// TP2

InsideMeshGrid grid;

std::vector< Vec3 > positions;
std::vector< Vec3 > normals;

std::vector< Vec3 > positions2;
std::vector< Vec3 > normals2;

std::vector< Vec3 > positions_generated;
std::vector< std::vector<int> > faces_generated;
std::vector< unsigned int > faces_generated_formated;
std::vector<Vec3> normals_generated;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 1000;
static unsigned int SCREENHEIGHT = 1000;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static bool fullScreen = false;




// ------------------------------------------------------------------------------------------------------------
// i/o and some stuff
// ------------------------------------------------------------------------------------------------------------
void loadPN (const std::string & filename , std::vector< Vec3 > & o_positions , std::vector< Vec3 > & o_normals ) {
    unsigned int surfelSize = 6;
    FILE * in = fopen (filename.c_str (), "rb");
    if (in == NULL) {
        std::cout << filename << " is not a valid PN file." << std::endl;
        return;
    }
    size_t READ_BUFFER_SIZE = 1000; // for example...
    float * pn = new float[surfelSize*READ_BUFFER_SIZE];
    o_positions.clear ();
    o_normals.clear ();
    while (!feof (in)) {
        unsigned numOfPoints = fread (pn, 4, surfelSize*READ_BUFFER_SIZE, in);
        for (unsigned int i = 0; i < numOfPoints; i += surfelSize) {
            o_positions.push_back (Vec3 (pn[i], pn[i+1], pn[i+2]));
            o_normals.push_back (Vec3 (pn[i+3], pn[i+4], pn[i+5]));
        }

        if (numOfPoints < surfelSize*READ_BUFFER_SIZE) break;
    }
    fclose (in);
    delete [] pn;
}
void savePN (const std::string & filename , std::vector< Vec3 > const & o_positions , std::vector< Vec3 > const & o_normals ) {
    if ( o_positions.size() != o_normals.size() ) {
        std::cout << "The pointset you are trying to save does not contain the same number of points and normals." << std::endl;
        return;
    }
    FILE * outfile = fopen (filename.c_str (), "wb");
    if (outfile == NULL) {
        std::cout << filename << " is not a valid PN file." << std::endl;
        return;
    }
    for(unsigned int pIt = 0 ; pIt < o_positions.size() ; ++pIt) {
        fwrite (&(o_positions[pIt]) , sizeof(float), 3, outfile);
        fwrite (&(o_normals[pIt]) , sizeof(float), 3, outfile);
    }
    fclose (outfile);
}
void scaleAndCenter( std::vector< Vec3 > & io_positions ) {
    Vec3 bboxMin( FLT_MAX , FLT_MAX , FLT_MAX );
    Vec3 bboxMax( FLT_MIN , FLT_MIN , FLT_MIN );
    for(unsigned int pIt = 0 ; pIt < io_positions.size() ; ++pIt) {
        for( unsigned int coord = 0 ; coord < 3 ; ++coord ) {
            bboxMin[coord] = std::min<float>( bboxMin[coord] , io_positions[pIt][coord] );
            bboxMax[coord] = std::max<float>( bboxMax[coord] , io_positions[pIt][coord] );
        }
    }
    Vec3 bboxCenter = (bboxMin + bboxMax) / 2.f;
    float bboxLongestAxis = std::max<float>( bboxMax[0]-bboxMin[0] , std::max<float>( bboxMax[1]-bboxMin[1] , bboxMax[2]-bboxMin[2] ) );
    for(unsigned int pIt = 0 ; pIt < io_positions.size() ; ++pIt) {
        io_positions[pIt] = (io_positions[pIt] - bboxCenter) / bboxLongestAxis;
    }
}

void applyRandomRigidTransformation( std::vector< Vec3 > & io_positions , std::vector< Vec3 > & io_normals ) {
    srand(time(NULL));
    Mat3 R = Mat3::RandRotation();
    Vec3 t = Vec3::Rand(1.f);
    for(unsigned int pIt = 0 ; pIt < io_positions.size() ; ++pIt) {
        io_positions[pIt] = R * io_positions[pIt] + t;
        io_normals[pIt] = R * io_normals[pIt];
    }
}

void subsample( std::vector< Vec3 > & i_positions , std::vector< Vec3 > & i_normals , float minimumAmount = 0.1f , float maximumAmount = 0.2f ) {
    std::vector< Vec3 > newPos , newNormals;
    std::vector< unsigned int > indices(i_positions.size());
    for( unsigned int i = 0 ; i < indices.size() ; ++i ) indices[i] = i;
    srand(time(NULL));
    std::random_shuffle(indices.begin() , indices.end());
    unsigned int newSize = indices.size() * (minimumAmount + (maximumAmount-minimumAmount)*(float)(rand()) / (float)(RAND_MAX));
    newPos.resize( newSize );
    newNormals.resize( newSize );
    for( unsigned int i = 0 ; i < newPos.size() ; ++i ) {
        newPos[i] = i_positions[ indices[i] ];
        newNormals[i] = i_normals[ indices[i] ];
    }
    i_positions = newPos;
    i_normals = newNormals;
}

bool save( const std::string & filename , std::vector< Vec3 > & vertices , std::vector< unsigned int > & triangles ) {
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open()) {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl;

    unsigned int n_vertices = vertices.size() , n_triangles = triangles.size()/3;
    myfile << n_vertices << " " << n_triangles << " 0" << std::endl;

    for( unsigned int v = 0 ; v < n_vertices ; ++v ) {
        myfile << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << std::endl;
    }
    for( unsigned int f = 0 ; f < n_triangles ; ++f ) {
        myfile << 3 << " " << triangles[3*f] << " " << triangles[3*f+1] << " " << triangles[3*f+2];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}



// ------------------------------------------------------------------------------------------------------------
// rendering.
// ------------------------------------------------------------------------------------------------------------

void initLight () {
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f,-16.0f,-50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    glCullFace (GL_BACK);
    glEnable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0, 0, 0, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
}



void drawTriangleMesh( std::vector< Vec3 > const & i_positions , std::vector< unsigned int > const & i_triangles ) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_triangles.size() / 3 ; ++tIt) {
        Vec3 p0 = i_positions[i_triangles[3*tIt]];
        Vec3 p1 = i_positions[i_triangles[3*tIt+1]];
        Vec3 p2 = i_positions[i_triangles[3*tIt+2]];
        Vec3 n = Vec3::cross(p1-p0 , p2-p0);
        n.normalize();
        glNormal3f( n[0] , n[1] , n[2] );
        glVertex3f( p0[0] , p0[1] , p0[2] );
        glVertex3f( p1[0] , p1[1] , p1[2] );
        glVertex3f( p2[0] , p2[1] , p2[2] );
    }
    glEnd();
}

void drawTriangleMesh2( std::vector< Vec3 > const & i_positions , std::vector< Vec3 > const & i_normals, std::vector< unsigned int > const & i_triangles) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_triangles.size() / 3 ; ++tIt) {
        Vec3 p0 = i_positions[i_triangles[3*tIt]];
        Vec3 p1 = i_positions[i_triangles[3*tIt+1]];
        Vec3 p2 = i_positions[i_triangles[3*tIt+2]];
        // Vec3 n = Vec3::cross(p1-p0 , p2-p0);
        Vec3 n0 = i_normals[i_triangles[3*tIt]];
        Vec3 n1 = i_normals[i_triangles[3*tIt+1]];
        Vec3 n2 = i_normals[i_triangles[3*tIt+2]];
        
        Vec3 n = (n0+n1+n2)/3;
        n.normalize();
        glNormal3f( n[0] , n[1] , n[2] );
        glVertex3f( p0[0] , p0[1] , p0[2] );
        glVertex3f( p1[0] , p1[1] , p1[2] );
        glVertex3f( p2[0] , p2[1] , p2[2] );
    }
    glEnd();
}

void drawPointSet( std::vector< Vec3 > const & i_positions , std::vector< Vec3 > const & i_normals ) {
    glBegin(GL_POINTS);
    for(unsigned int pIt = 0 ; pIt < i_positions.size() ; ++pIt) {
        glNormal3f( i_normals[pIt][0] , i_normals[pIt][1] , i_normals[pIt][2] );
        glVertex3f( i_positions[pIt][0] , i_positions[pIt][1] , i_positions[pIt][2] );
    }
    glEnd();
}

	void drawGrid(const InsideMeshGrid& grid) {
    glPointSize(2.0f); // Ajustez la taille des points selon vos besoins
    glBegin(GL_POINTS);

    for (size_t x = 0; x < grid.X; ++x) {
        for (size_t y = 0; y < grid.Y; ++y) {
            for (size_t z = 0; z < grid.Z; ++z) {
                const Vec3& point = grid(x, y, z);
                glVertex3f(point[0], point[1], point[2]);
            }
        }
    }

    glEnd();
}

void draw () {
    glPointSize(0.1); // for example...

    //glColor4f(1,1,1,0.1);
    //drawPointSet(positions , normals);

    //glColor3f(1,1,0);
    //drawPointSet(positions2 , normals2);


    glDisable(GL_LIGHTING);

    // bbmax, min
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    //glColor3f(1,1,0);
    //glVertex3f( grid.bbmin[0] , grid.bbmin[1] , grid.bbmin[2] );
    //glColor3f(0.5,0.5,0);
    //glVertex3f( grid.bbmax[0] , grid.bbmax[1] , grid.bbmax[2] );

    glEnd();

    // isoValue
    // glPointSize(2.0f);
    // glBegin(GL_POINTS);
    
     for(int i = 0; i < grid.isoValue_pos.size(); ++i)
     {
         if(grid.isoValue[i] == 1){
             glColor3f(0,0.9,0.2);
         } else {
             glColor3f(1,0,0.2);
         }
             glVertex3f( grid.isoValue_pos[i][0] , grid.isoValue_pos[i][1] , grid.isoValue_pos[i][2] );
     }
     glEnd();

    // generated_vertex
    // glPointSize(5.0f);
    // glBegin(GL_POINTS);
    // glColor3f(1,0,0);
    
    // for(int i = 0; i < positions_generated.size(); ++i)
    // {
    //     glVertex3f( positions_generated[i][0] , positions_generated[i][1] , positions_generated[i][2] );
    // }
    // glEnd();

    glEnable(GL_LIGHTING);

    // faces_generated
    // glDisable(GL_LIGHTING);
    // glPointSize(1.0f);
    // glBegin(GL_TRIANGLES);

    // for(int i = 0; i < faces_generated.size(); i++)
    // {
    //     glColor3f(1,0,0);
    //     glVertex3f( positions_generated[faces_generated[i][0]][0], positions_generated[faces_generated[i][0]][1], positions_generated[faces_generated[i][0]][2]);
    //     glColor3f(0,1,0);
    //     glVertex3f( positions_generated[faces_generated[i][1]][0], positions_generated[faces_generated[i][1]][1], positions_generated[faces_generated[i][1]][2]);
    //     glColor3f(0,0,1);
    //     glVertex3f( positions_generated[faces_generated[i][2]][0], positions_generated[faces_generated[i][2]][1], positions_generated[faces_generated[i][2]][2]);
    // }
    // glEnd();
    // glEnable(GL_LIGHTING);
    
    // glColor3f(0.870, 0.784, 0.6);
    //drawTriangleMesh2(positions_generated, normals_generated, faces_generated_formated);
    drawTriangleMesh(positions_generated, faces_generated_formated);
    glColor3f(0,0,1);
    //glDisable(GL_LIGHTING);
    //drawPointSet(positions_generated, normals_generated);

    glPointSize(0.1);
    glColor3f(1,0,0);
    //drawGrid(grid);

}

void display () {
    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    draw ();
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    glutPostRedisplay ();
}

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;

    case 'w':
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, polygonMode);
        if(polygonMode[0] != GL_FILL)
            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        else
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        break;

    default:
        break;
    }
    idle ();
}

void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
    idle ();
}

void motion (int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate (x, y);
    }
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH), (lastY-y)/static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}


void reshape(int w, int h) {
    camera.resize (w, h);
}

struct Plan{
    Vec3 point;
    Vec3 normal;
};

void addNoiseAlongNormals(std::vector<Vec3> & pos, std::vector<Vec3> & normals, float amplitude){
    for(int i = 0; i < pos.size(); i++)
    {
        // Créer un générateur de nombres aléatoires
        std::random_device rd;
        std::mt19937 gen(rd());

        // Créer une distribution uniforme dans la plage spécifiée
        std::uniform_real_distribution<float> distribution(-amplitude, amplitude);

        // Générer un nombre aléatoire de type float dans la plage spécifiée
        float randomFloat = distribution(gen);

        pos[i] = pos[i] + ( randomFloat * normals[i]);
    }
}

void project(Vec3 & inputPoint, Plan const & plan, Vec3 & outputPoint, Vec3 & outputNormal){
    outputPoint = inputPoint - Vec3::dot(((inputPoint - plan.point)), plan.normal) * plan.normal;
    outputNormal = plan.normal;
    outputNormal.normalize();
}

// - récupérer les k plus proches voisins
// - calcul du plan de ces voisins
//      - centroide  = somme(pi) / |pi|, moyenne des points
//      - normale = somme(ni) / |ni |, moyenne des normales
//      - Wi = 1, poid
// - projeter inputPoint sur le plan
// - refaire

// kernel_type : 0 : Singulier, 1 : Gaussien, 2 : Wendland
void SPPS(Vec3 inputPoint,
            Vec3 & outputPoint, Vec3 & outputNormal,
            std::vector<Vec3> const & positions, std::vector<Vec3> const & normals,
            BasicANNkdTree const & kdtree,
            int kernel_type, float radius, unsigned int nbIterations = 10, unsigned int knn = 20)
            {
                ANNidxArray id_nearest_neighbors = new ANNidx[ knn ];
                ANNdistArray square_distances_to_neighbors = new ANNdist[ knn ];

                outputPoint = inputPoint;

                for(int iteration = 0; iteration < nbIterations; iteration++)
                {
                    Plan plan;
                    plan.point = Vec3(0,0,0);
                    plan.normal = Vec3(0,0,0);


                    float weightSum = 0;

                    // recherche des k voisins les plus proche du point étudié
                    kdtree.knearest(outputPoint, knn, id_nearest_neighbors, square_distances_to_neighbors);

                    // calcul du plan moyen à partir des voisins
                    for(int i = 0; i < knn; i++){

                        Vec3 pi = positions[id_nearest_neighbors[i]];
                        float distance_x_pi = ( outputPoint - pi).length();
                        //float distance_x_pi = square_distances_to_neighbors[id_nearest_neighbors[i]];
                        float weight = 0;

                        // definition des poids
                        if(kernel_type == 0)
                        {
                            weight = pow(radius/distance_x_pi,1);
                        }
                        else if (kernel_type == 1)
                        {
                            weight = pow(M_E, (pow(-distance_x_pi,2) / pow(radius, 2)) );
                        }
                        else if (kernel_type == 2)
                        {
                            weight = pow( (1 - (distance_x_pi / radius)), 4) * (1 + 4 * (distance_x_pi / radius));
                        }

                        plan.point += weight * pi;
                        plan.normal += weight * normals[id_nearest_neighbors[i]];;

                        weightSum += weight;
                    }

                    // moyennes pondérée
                    plan.point /= weightSum;
                    plan.normal /= weightSum;

                    project(outputPoint, plan, outputPoint, outputNormal);
                }
            }


// - récupérer les k voisins
// - projeter x sur les plans des points voisins
// - calculer le centroide des projections avec des poids
// - ré itérer

void HPPS(Vec3 inputPoint,
          Vec3 & outputPoint, Vec3 & outputNormal,
          std::vector<Vec3> const & positions, std::vector<Vec3> const & normals,
          BasicANNkdTree const & kdtree,
          int kernel_type, float radius, unsigned int nbIterations = 10, unsigned int knn = 20)
{

    ANNidxArray id_nearest_neighbors = new ANNidx[ knn ];
    ANNdistArray square_distances_to_neighbors = new ANNdist[ knn ];

    outputPoint = inputPoint;

    Vec3 point_intermediaire = Vec3(0,0,0);
    Vec3 normale_intermediaire = Vec3(0,0,0);

    for(int iteration = 0; iteration < nbIterations; iteration++)
    {
        float weightSum = 0;

        // recherche des k voisins les plus proche du point étudié
        kdtree.knearest(outputPoint, knn, id_nearest_neighbors, square_distances_to_neighbors);

        // projection de x sur les plans voisins
        for(int i = 0; i < knn; i++){

            Vec3 pi = positions[id_nearest_neighbors[i]];

            Plan planHermiteVoisin;
            planHermiteVoisin.point = Vec3(0,0,0);
            planHermiteVoisin.normal = Vec3(0,0,0);
            planHermiteVoisin.point = pi;
            planHermiteVoisin.normal = normals[id_nearest_neighbors[i]];

            Vec3 outputProjectedPoint, outputProjectedNormals;

            project(pi, planHermiteVoisin, outputProjectedPoint, outputProjectedNormals);

            float distance_x_pip = ( outputPoint - outputProjectedPoint).length();
            float weight = 0;

            // definition des poids
            if(kernel_type == 0)
            {
                weight = pow(radius/distance_x_pip,1);
            }
            else if (kernel_type == 1)
            {
                weight = pow(M_E, (pow(-distance_x_pip,2) / pow(radius, 2)) );
            }
            else if (kernel_type == 2)
            {
                weight = pow( (1 - (distance_x_pip / radius)), 4) * (1 + 4 * (distance_x_pip/ radius));
            }

            point_intermediaire += weight * outputProjectedPoint;
            normale_intermediaire += weight * outputProjectedNormals;

            if(weight != weight || point_intermediaire[0] != point_intermediaire[0]){
                std::cout << "point_intermediaire = Nan" << "\n";
                exit(-1);
            }

            weightSum += weight;
        }

        point_intermediaire /= weightSum;
        normale_intermediaire /= weightSum;

        outputPoint = point_intermediaire;
        outputNormal = normale_intermediaire;
    }

    delete [] id_nearest_neighbors;
    delete [] square_distances_to_neighbors;
}

// Fonction pour calculer le point moyen entre deux points
Vec3 midpoint(const Vec3 &a, const Vec3 &b) {
    return {(a[0] + b[0]) / 2.0f, (a[1] + b[1]) / 2.0f, (a[2] + b[2]) / 2.0f};
}

void subdivideMesh(std::vector<Vec3> &pos, std::vector<std::vector<int>> &faces)
{
    std::vector<Vec3> newPos;
    std::vector< std::vector<int> > newFaces;
    
    for(int f = 0; f < faces.size(); ++f)
    {
        int i_p0 = faces[f][0];
        int i_p1 = faces[f][1];
        int i_p2 = faces[f][2];

        Vec3 p0 = pos[faces[f][0]];
        Vec3 p1 = pos[faces[f][1]];
        Vec3 p2 = pos[faces[f][2]];

        Vec3 n_p0 = p0 + (p1 - p0)/2;
        Vec3 n_p1 = p1 + (p2 - p1)/2;
        Vec3 n_p2 = p2 + (p0 - p2)/2;

        int offset = pos.size();

        pos.push_back(n_p0);
        pos.push_back(n_p1);
        pos.push_back(n_p2);


        std::vector<int> f1 = {offset, offset+1, offset+2};
        std::vector<int> f2 = {i_p0, offset, offset+2};
        std::vector<int> f3 = {offset, i_p1, offset+1};
        std::vector<int> f4 = {offset+1, i_p2, offset+2};


        newFaces.push_back(f1);
        newFaces.push_back(f2);
        newFaces.push_back(f3);
        newFaces.push_back(f4);
    }
    faces = newFaces;
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

    for (int m = 0; m < 2; m++)
    {
    for (int i = 0; i < positions2.size(); i++) // Pour chaque point a projeter
    {
        kdtree.knearest(positions2[i], k, id_nearest_neighbors, square_distances_to_neighbors);

        compute_weights(weights, weights_normalized, id_nearest_neighbors, square_distances_to_neighbors, k);
        u4 = compute_u4(positions, normals, id_nearest_neighbors, weights, weights_normalized, k);

        if (sqrt(pow(u4, 2)) < 0.1) // Si u4 est très proche/egal a 0
        {
            //std::cout<<"0"<<std::endl;
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
            u1u2u3 = compute_u1u2u3(positions, normals, id_nearest_neighbors, weights_normalized, k, u4);
            u0 = compute_u0(positions, id_nearest_neighbors, weights_normalized, i, k, u1u2u3, u4);
            Vec3 center;
            double radius;
            compute_circle(center, radius, u0, u1u2u3, u4);
            project_to_circle(center, radius, positions2[i], positions2[i], normals2[i]);
        }
    }
    }
    delete[] id_nearest_neighbors;
    delete[] square_distances_to_neighbors;
}

int main (int argc, char ** argv) {
    if (argc > 2) {
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ("tp point processing");

    init ();
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);

    
    // Load a first pointset, and build a kd-tree:
    loadPN("pointsets/dino.pn" , positions , normals);

    //addNoiseAlongNormals(positions, normals, 0.02);

    BasicANNkdTree kdtree;
    kdtree.build(positions);

    // Create a second pointset that is artificial, and project it on pointset1 using MLS techniques:
    
    // positions2.resize( 20000 );
    // normals2.resize(positions2.size());
    // for( unsigned int pIt = 0 ; pIt < positions2.size() ; ++pIt ) {
    //     positions2[pIt] = Vec3(
    //                 -0.6 + 1.2 * (double)(rand())/(double)(RAND_MAX),
    //                 -0.6 + 1.2 * (double)(rand())/(double)(RAND_MAX),
    //                 -0.6 + 1.2 * (double)(rand())/(double)(RAND_MAX)
    //                 );
    //     positions2[pIt].normalize();
    //     positions2[pIt] = 0.6 * positions2[pIt];
    // }

    // Vec3 outputPoint;
    // Vec3 outputNormal;

    // for(int i = 0; i < positions2.size(); i++){
    //     HPPS(positions2[i],
    //             outputPoint,
    //             outputNormal,
    //             positions,
    //             normals,
    //             kdtree,
    //             0, 1, 3, 20);

    //     positions2[i] = outputPoint;
    //     normals2[i] = outputNormal;
    // }


    // ------------------ TP2

    // Etapes :
    // - Créer une grille englobant le mesh
    // - diviser la grille en x case
    // - dans chaque case dire si on est dans le mesh ou pas :
    //      - dans chaque cases on a un Vec3 au milieu
    //      - créer le point le plus proche de celui ci sur la surface du mesh avec HPPS
    //      - faire un produit scalaire entre la normale de ce point et le vecteur hpps_point -> point de la case
    //      - si positif ça veut dire que les directions des vecs sont alligné et que le point et hors du mesh.
    //        mettre dans le tableau 3D le signe correspondant 

    Vec3 outputPoint;
    Vec3 outputNormal;

    //grid = InsideMeshGrid();
    int grid_size = 16;
    grid = InsideMeshGrid(grid_size, positions);

    for(int i = 0; i < grid.isoValue.size(); ++i)
    {
        HPPS(grid.isoValue_pos[i],
        outputPoint,
        outputNormal,
        positions,
        normals,
        kdtree,
        0, 1, 3, 20);

        // check si le point de la grille est dedans

        Vec3 hppsPoint_gridPoint = grid.isoValue_pos[i] - outputPoint;
normals_generated.push_back(outputNormal);
        float scal = Vec3::dot(hppsPoint_gridPoint, outputNormal);

        if(scal <= 0){
            grid.isoValue[i] = 1;
        }
    }
    std::cout<<"The grid size :"<< grid_size<<std::endl;
    grid.createMesh(positions_generated,faces_generated);
    std::cout<<"Wait ... APSS x2 takes long time to compute (60-90 sec)..."<<std::endl;
    APSS(positions,normals,positions_generated,normals_generated,kdtree,0);
    std::cout<<"First APSS done, wait for the second..."<<std::endl;
    subdivideMesh(positions_generated, faces_generated);
    APSS(positions,normals,positions_generated,normals_generated,kdtree,0);
    std::cout<<"Done with APSS"<<std::endl;

    // dérouler le vecteur faces
    for(std::vector<int> face : faces_generated)
    {
        faces_generated_formated.push_back(face[0]);
        faces_generated_formated.push_back(face[1]);
        faces_generated_formated.push_back(face[2]);
    }

    glutMainLoop ();
    return EXIT_SUCCESS;
}

