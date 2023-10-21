#ifndef DRAWFRAME_H
#define DRAWFRAME_H

#include <GL/glut.h>

#include "Vec3.h"

static void drawFrame(){
        glBegin(GL_LINES);
        
        glColor3f(1,0,0);
        glVertex3f (0,0,0);
        glVertex3f (1,0,0);
        glColor3f(0,1,0);
        glVertex3f (0,0,0);
        glVertex3f (0,1,0);
        glColor3f(0,0,1);
        glVertex3f (0,0,0);
        glVertex3f (0,0,1);
        glEnd();

        int size = 40;
        for(float i = -size; i < size; i += 0.5){
            glLineWidth(0.1);
            glBegin(GL_LINES);

            glColor3f(0.5,0.5,0.5);
            glVertex3f (i,-size,0);
            glVertex3f (i,size,0);

            glVertex3f (-size,i,0);
            glVertex3f (size,i,0);
            
            glEnd();
        }
}

static void drawPoint(Vec3 p, float r, float g, float b, float size){
        glPointSize(SIZE_MAX);
        glBegin(GL_POINTS);
        glColor3f(r,g,b);
        glVertex3f (p[0],p[1],p[2]);
        glEnd();
}

#endif