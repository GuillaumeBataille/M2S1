#ifndef SKELETON_H
#define SKELETON_H

#include <vector>
#include <queue>
#include <map>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include "Vec3.h"
#include "src/drawFrame.h"

#include <cmath>

#include <GL/glut.h>

// -------------------------------------------
// Basic Skeleton class
// -------------------------------------------

struct Articulation
{
    // membres :
    Vec3 p; // une position

    int fatherBone; // there should be only 1
    std::vector<unsigned int> childBones;

    void setFatherBone(int f)
    {
        if (fatherBone >= 0)
        {
            assert(fatherBone == f);
        }
        fatherBone = f;
    }
    bool isRoot() const
    {
        return fatherBone == -1;
    }

    Articulation() : fatherBone(-1) { childBones.clear(); }
};

struct Bone
{
    // membres :
    unsigned int joints[2];

    int fatherBone; // there should be only 1
    std::vector<unsigned int> childBones;

    void setFatherBone(int f)
    {
        if (fatherBone >= 0)
        {
            assert(fatherBone == f);
        }
        fatherBone = f;
    }
    bool isRoot() const
    {
        return fatherBone == -1;
    }

    static double euclideanWeight(Vec3 vertex, Vec3 projectionOnBone, int n)
    {
        double distance = (projectionOnBone - vertex).length();
        return pow(1.0f / distance, n);
    }

    Bone() : fatherBone(-1) { childBones.clear(); }
};

struct BoneTransformation
{
    // membres :
    Mat3 localRotation;

    Mat3 world_space_rotation;
    Vec3 world_space_translation;

    BoneTransformation() : localRotation(Mat3::Identity()), world_space_rotation(Mat3::Identity()), world_space_translation(0, 0, 0) {}
};

struct SkeletonTransformation
{
    // membres :
    std::vector<BoneTransformation> bone_transformations;
    std::vector<Vec3> articulations_transformed_position;

    void resize(unsigned int n_bones, unsigned int n_articulations)
    {
        bone_transformations.resize(n_bones);
        articulations_transformed_position.resize(n_articulations);
    }
};

struct Skeleton
{
    // membres :
    std::vector<Articulation> articulations;
    std::vector<Bone> bones;
    std::vector<unsigned int> ordered_bone_indices; // process them by order in the hierarchy

    void buildStructure()
    {
        ordered_bone_indices.clear();
        std::vector<unsigned int> rootBones; // why not have several

        for (unsigned int b = 0; b < bones.size(); ++b)
        {
            Articulation &a0 = articulations[bones[b].joints[0]];
            Articulation &a1 = articulations[bones[b].joints[1]];
            a0.childBones.push_back(b);
            a1.setFatherBone(b);
        }

        for (unsigned int aIdx = 0; aIdx < articulations.size(); ++aIdx)
        {
            Articulation &a = articulations[aIdx];
            if (a.isRoot())
            {
                for (unsigned int bIt = 0; bIt < a.childBones.size(); ++bIt)
                {
                    unsigned int b = a.childBones[bIt];
                    rootBones.push_back(b);
                }
            }
            else
            {
                unsigned int bfIdx = a.fatherBone;
                Bone &bf = bones[bfIdx];
                for (unsigned int bIt = 0; bIt < a.childBones.size(); ++bIt)
                {
                    unsigned int bcIdx = a.childBones[bIt];
                    Bone &bc = bones[bcIdx];
                    bc.setFatherBone(bfIdx);
                    bf.childBones.push_back(bcIdx);
                }
            }
        }

        for (unsigned int rIt = 0; rIt < rootBones.size(); ++rIt)
        {
            unsigned int rootboneIdx = rootBones[rIt];
            std::queue<unsigned int> bonesIndices;
            bonesIndices.push(rootboneIdx);
            while (!bonesIndices.empty())
            {
                unsigned int bIdx = bonesIndices.front();
                bonesIndices.pop();
                ordered_bone_indices.push_back(bIdx);
                Bone &b = bones[bIdx];
                for (unsigned int bIt = 0; bIt < b.childBones.size(); ++bIt)
                {
                    unsigned int bcIdx = b.childBones[bIt];
                    bonesIndices.push(bcIdx);
                }
            }
        }

        assert(ordered_bone_indices.size() == bones.size());
    }

    void load(const std::string &filename)
    {
        std::ifstream in(filename.c_str());
        if (!in)
            exit(EXIT_FAILURE);
        std::string tmpString;
        unsigned int sizeA;
        in >> tmpString >> sizeA;
        articulations.resize(sizeA);
        for (unsigned int i = 0; i < sizeA; i++)
            in >> articulations[i].p[0] >> articulations[i].p[1] >> articulations[i].p[2];

        unsigned int sizeB;
        in >> tmpString >> sizeB;
        bones.resize(sizeB);
        for (unsigned int i = 0; i < sizeB; i++)
        {
            for (unsigned int j = 0; j < 2; j++)
                in >> bones[i].joints[j];
        }
        in.close();

        buildStructure();
    }

    // Update les coord de this->articulation avec la transformation s
    void updateCoord(SkeletonTransformation s)
    {
        for (int i = 0; i < s.articulations_transformed_position.size(); i++)
        {
            this->articulations[i].p = s.articulations_transformed_position[i];
        }
    }

    void computeGlobalTransformationParameters(SkeletonTransformation &transfo)
    {
        std::vector<Vec3> &articulations_transformed_position = transfo.articulations_transformed_position;
        articulations_transformed_position.resize(articulations.size());
        for (unsigned int bIt = 0; bIt < ordered_bone_indices.size(); ++bIt)
        {
            unsigned bIdx = ordered_bone_indices[bIt];
            Bone &b = bones[bIdx];

            if (b.isRoot())
            {
                Vec3 a0RestPos = articulations[b.joints[0]].p;
                Vec3 a0TargetPos = a0RestPos;
                articulations_transformed_position[b.joints[0]] = a0TargetPos;
                BoneTransformation &bone_transformation = transfo.bone_transformations[bIdx];
                bone_transformation.world_space_rotation = bone_transformation.localRotation;

                // set the articulation as pivot point :
                bone_transformation.world_space_translation = a0TargetPos - bone_transformation.world_space_rotation * a0RestPos;

                // update the child articulation :
                Vec3 a1RestPos = articulations[b.joints[1]].p;
                Vec3 a1TargetPos = bone_transformation.world_space_rotation * a1RestPos + bone_transformation.world_space_translation;
                articulations_transformed_position[b.joints[1]] = a1TargetPos;
            }
            else
            {
                Vec3 a0RestPos = articulations[b.joints[0]].p;
                Vec3 a0TargetPos = articulations_transformed_position[b.joints[0]];

                BoneTransformation &bone_transformation = transfo.bone_transformations[bIdx];
                BoneTransformation &bFatherTransfo = transfo.bone_transformations[b.fatherBone];
                bone_transformation.world_space_rotation = bFatherTransfo.world_space_rotation * bone_transformation.localRotation;

                // set the articulation as pivot point :
                bone_transformation.world_space_translation = a0TargetPos - bone_transformation.world_space_rotation * a0RestPos;

                // update the child articulation :
                Vec3 a1RestPos = articulations[b.joints[1]].p;
                Vec3 a1TargetPos = bone_transformation.world_space_rotation * a1RestPos + bone_transformation.world_space_translation;
                articulations_transformed_position[b.joints[1]] = a1TargetPos;
            }
        }
    }

    // Fonction pour faire tourner un os ou un groupe d'os
    // Ne fonctionne pas si un bone a plusieurs enfants. Il faudra implenter une pile dans ce cas.
    void rotateBoneAndChilds_UpdateTransfo(int bone_i, SkeletonTransformation &transfo, Mat3 rot)
    {
        int currentBone_i = bone_i;
        int boneRoot_i;
        int boneEnd_i;
        Vec3 rotOrigin = articulations[bone_i].p;

        if (bones[currentBone_i].childBones.size() == 0)
        {
            boneRoot_i = bones[currentBone_i].joints[0];
            boneEnd_i = bones[currentBone_i].joints[1];

            Vec3 boneRoot = articulations[boneRoot_i].p;
            Vec3 boneEnd = articulations[boneEnd_i].p;

            Vec3 boneVec = boneEnd - boneRoot;

            Vec3 newPos = boneRoot + (rot * boneVec);

            // this->articulations[boneEnd_i].p = newPos;
            transfo.articulations_transformed_position[boneEnd_i] = newPos;
            // transfo.bone_transformations[bone_i].world_space_rotation = transfo.bone_transformations[bone_i].world_space_rotation * rot;
            transfo.bone_transformations[bone_i].localRotation = rot;

            // transfo.bone_transformations[bone_i].world_space_translation = newPos;
        }
        else
        {
            while (true)
            {
                boneRoot_i = bones[currentBone_i].joints[0];
                boneEnd_i = bones[currentBone_i].joints[1];

                Vec3 boneRoot = articulations[boneRoot_i].p;
                Vec3 boneEnd = articulations[boneEnd_i].p;

                Vec3 originToPoint = boneEnd - rotOrigin;

                Vec3 newPos = rotOrigin + (rot * originToPoint);

                // this->articulations[boneEnd_i].p = newPos;
                transfo.articulations_transformed_position[boneEnd_i] = newPos;

                // transfo.bone_transformations[bone_i].world_space_rotation = transfo.bone_transformations[bone_i].world_space_rotation * rot;

                // transfo.bone_transformations[bone_i].world_space_translation = newPos;
                transfo.bone_transformations[bone_i].localRotation = rot;

                if (bones[currentBone_i].childBones.size() == 0)
                {
                    break;
                }
                else
                {
                    currentBone_i = bones[currentBone_i].childBones[0];
                }
            }
        }
    }

    void computeProceduralAnim(double t, SkeletonTransformation &transfo)
    {
        transfo.bone_transformations.resize(bones.size());
        for (unsigned int bIt = 0; bIt < ordered_bone_indices.size(); ++bIt)
        {
            unsigned bIdx = ordered_bone_indices[bIt];
            Bone &b = bones[bIdx];
            if (b.isRoot())
            {
                BoneTransformation &bone_transformation = transfo.bone_transformations[bIdx];
                bone_transformation.localRotation = Mat3::Identity();
            }
            else
            {
                BoneTransformation &bone_transformation = transfo.bone_transformations[bIdx];
                Vec3 axis(cos(2 * M_PI * bIt / (double)(bones.size())), sin(2 * M_PI * bIt / (double)(bones.size())), 0.0);
                bone_transformation.localRotation = Mat3::getRotationMatrixFromAxisAndAngle(axis, (0.25 * M_PI) * cos(t));
            }
        }

        // update articulation positions:
        computeGlobalTransformationParameters(transfo);
    }

    void updateIKChain(SkeletonTransformation &transfoIK, unsigned int targetArticulation, Vec3 targetPosition, unsigned int maxIterNumber = 20, double epsilonPrecision = 0.000001)
    {
        // updateCoord(transfoIK);
        bool loopbreak = false;

        int currentArti_i = targetArticulation;
        int currentBone_i = articulations[currentArti_i].fatherBone;
        Vec3 boneEnd = articulations[targetArticulation].p;
        // s'arrêtte à maxIter
        // Cyclic Coordinates Descend
        for (int it = 0; it < maxIterNumber; it++)
        {
            int boneEnd_i = bones[currentBone_i].joints[1];

            while (!loopbreak)
            {
                // ------------- traitement ------------- //

                if ((targetPosition - articulations[targetArticulation].p).length() < epsilonPrecision)
                {
                    std::cout << "Trop colle" << std::endl;
                    return;
                }

                int boneRoot_i = bones[currentBone_i].joints[0];

                Vec3 oldBoneEnd = boneEnd;
                std::cout << "boneEnd : " << boneEnd.toString() << std::endl;
                Vec3 boneRoot = transfoIK.articulations_transformed_position[boneRoot_i];

                Vec3 targetVec = targetPosition - boneRoot;

                Vec3 hingeToTargetArti = boneEnd - boneRoot;

                Mat3 rotMat = Mat3::getRotationMatrixAligning(hingeToTargetArti, targetVec);

                // rotateBoneAndChilds_UpdateTransfo(currentBone_i, transfoIK, rotMat);
                transfoIK.bone_transformations[currentBone_i].localRotation = rotMat;
                computeGlobalTransformationParameters(transfoIK);

                currentArti_i = bones[currentBone_i].joints[0];
                currentBone_i = bones[currentBone_i].fatherBone;

                if (currentBone_i == -1)
                {
                    loopbreak = true;
                }
            }

            currentArti_i = targetArticulation;
            currentBone_i = articulations[currentArti_i].fatherBone;
            loopbreak = false;
        }
    }

    //----------------------------------------------//
    //----------------------------------------------//
    //----------------------------------------------//
    // draw functions :
    //----------------------------------------------//
    //----------------------------------------------//
    //----------------------------------------------//
    void draw(int displayedBone, int displayedArticulation) const
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH);
        glDisable(GL_DEPTH_TEST);
        glLineWidth(3.0);
        glBegin(GL_LINES);
        for (unsigned int i = 0; i < bones.size(); i++)
        {
            glColor3f(1, 0, 0);
            {
                const Articulation &v = articulations[bones[i].joints[0]];
                glVertex3f(v.p[0], v.p[1], v.p[2]);
            }
            glColor3f(1, 1, 1);
            {
                const Articulation &v = articulations[bones[i].joints[1]];
                glVertex3f(v.p[0], v.p[1], v.p[2]);
            }
        }
        glEnd();

        // we highlight the ordered bone number displayedBone
        if (displayedBone >= 0 && displayedBone < ordered_bone_indices.size())
        {
            displayedBone = ordered_bone_indices[displayedBone];
            glLineWidth(8.0);
            glBegin(GL_LINES);
            glColor3f(1, 0, 0);
            {
                const Articulation &v = articulations[bones[displayedBone].joints[0]];
                glVertex3f(v.p[0], v.p[1], v.p[2]);
            }
            glColor3f(1, 0, 0);
            {
                const Articulation &v = articulations[bones[displayedBone].joints[1]];
                glVertex3f(v.p[0], v.p[1], v.p[2]);
            }
            glEnd();
        }

        // draw articulations:
        glPointSize(12.0);
        glBegin(GL_POINTS);
        glColor3f(0.5, 0, 0);
        for (unsigned int i = 0; i < articulations.size(); i++)
        {
            const Articulation &v = articulations[i];
            glVertex3f(v.p[0], v.p[1], v.p[2]);
        }
        glEnd();

        if (displayedArticulation >= 0 && displayedArticulation < articulations.size())
        {
            glPointSize(16.0);
            glBegin(GL_POINTS);
            glColor3f(1, 0, 0);
            {
                const Articulation &v = articulations[displayedArticulation];
                glVertex3f(v.p[0], v.p[1], v.p[2]);
            }
            glEnd();
        }

        glEnable(GL_DEPTH);
        glEnable(GL_DEPTH_TEST);
    }
    void drawTransformedSkeleton(int displayedBone, int displayedArticulation, SkeletonTransformation const &transfo) const
    {
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH);
        glDisable(GL_DEPTH_TEST);
        glLineWidth(3.0);
        glBegin(GL_LINES);
        for (unsigned int i = 0; i < bones.size(); i++)
        {
            glColor3f(1, 0, 0);
            {
                Vec3 p = transfo.articulations_transformed_position[bones[i].joints[0]];
                glVertex3f(p[0], p[1], p[2]);
            }
            glColor3f(1, 1, 1);
            {
                Vec3 p = transfo.articulations_transformed_position[bones[i].joints[1]];
                glVertex3f(p[0], p[1], p[2]);
            }
        }
        glEnd();

        // we highlight the ordered bone number displayedBone
        if (displayedBone >= 0 && displayedBone < ordered_bone_indices.size())
        {
            displayedBone = ordered_bone_indices[displayedBone];
            glLineWidth(8.0);
            glBegin(GL_LINES);
            glColor3f(1, 0, 0);
            {
                Vec3 p = transfo.articulations_transformed_position[bones[displayedBone].joints[0]];
                glVertex3f(p[0], p[1], p[2]);
            }
            glColor3f(1, 0, 0);
            {
                Vec3 p = transfo.articulations_transformed_position[bones[displayedBone].joints[1]];
                glVertex3f(p[0], p[1], p[2]);
            }
            glEnd();
        }

        // draw articulations:
        glPointSize(12.0);
        glBegin(GL_POINTS);
        glColor3f(0.5, 0, 0);
        for (unsigned int i = 0; i < articulations.size(); i++)
        {
            Vec3 p = transfo.articulations_transformed_position[i];
            glVertex3f(p[0], p[1], p[2]);
        }
        glEnd();

        if (displayedArticulation >= 0 && displayedArticulation < articulations.size())
        {
            glPointSize(16.0);
            glBegin(GL_POINTS);
            glColor3f(1, 0, 0);
            {
                Vec3 p = transfo.articulations_transformed_position[displayedArticulation];
                glVertex3f(p[0], p[1], p[2]);
            }
            glEnd();
        }

        glEnable(GL_DEPTH);
        glEnable(GL_DEPTH_TEST);
    }
};

#endif // SKELETON_H
