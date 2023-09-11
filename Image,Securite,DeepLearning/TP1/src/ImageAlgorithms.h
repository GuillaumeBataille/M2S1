//
// Created by bat-portable on 21/01/2023.
//

#ifndef CODAGE_IMAGEALGORITHMS_H
#define CODAGE_IMAGEALGORITHMS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <random>
#include "color.h"


namespace ImageAlgorithms
{
    int findClass(Color &c, std::vector<Color> classes)
    {
        float distanceMin = std::numeric_limits<float>::max();
        int index;

        for(int colorIndex = 0; colorIndex < classes.size(); colorIndex++)
        {
            float d = c.dist(classes[colorIndex]);
            if(d < distanceMin){
                distanceMin = d;
                index = colorIndex;
            }
        }

        return index;
    }

    void findDistantColor(Color &c1, Color &c2, OCTET* ImgIn, int nTaille, int precision)
    {
        int nTaille3 = nTaille*3;

        float distanceMax = 0;
        precision *= 3;

        for (int i = 0; i < nTaille3; i += precision)
        {
            Color temp_c1(ImgIn[i], ImgIn[i+1], ImgIn[i+2]);

            for (int j = 0; j < nTaille3; j += precision)
            {
                Color temp_c2(ImgIn[j], ImgIn[j+1], ImgIn[j+2]);

                float distance = temp_c1.dist(temp_c2);

                if(distance > distanceMax)
                {
                    distanceMax = distance;
                    c1 = temp_c1;
                    c2 = temp_c2;
                }
            }
        }
    }

    std::vector<Color> findRandomColors(OCTET* ImgIn, int number, int width, int height)
    {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distWidth(0,width);
        std::uniform_int_distribution<std::mt19937::result_type> distHeight(0,height);

        std::vector<Color> tab;

        for(int i = 0; i < number; i++)
        {
            int randWidth = distWidth(rng);
            int randHeight = distHeight(rng);

            int index = randWidth*3 + width*randHeight*3;

            Color c(ImgIn[index], ImgIn[index+1], ImgIn[index+2]);
            tab.push_back(c);
        }
        return tab;
    }

    void k_mean(OCTET* ImgIn, std::vector<Color> &classes, int nTaille, int repetition){
        int nTaille3 = nTaille*3;

        for(int k = 0; k < repetition; k++)
        {
            std::vector<Color_float> avgColorPerClasses = std::vector<Color_float>(classes.size());
            std::vector<int> countForAvgClasses = std::vector<int>(classes.size());

            for (int i = 0; i < nTaille3; i += 3)
            {
                Color actualColor(ImgIn[i], ImgIn[i+1], ImgIn[i+2]);
                float distanceMin = std::numeric_limits<float>::max();
                int index;


                for(int colorIndex = 0; colorIndex < classes.size(); colorIndex++)
                {
                    float d = actualColor.dist(classes[colorIndex]);
                    if(d < distanceMin){
                        distanceMin = d;
                        index = colorIndex;
                    }
                }
                avgColorPerClasses[index] = avgColorPerClasses[index] + actualColor;
                countForAvgClasses[index] = countForAvgClasses[index] + 1;
            }

            for(int i = 0; i < avgColorPerClasses.size(); i++)
            {
                classes[i] = (avgColorPerClasses[i] / countForAvgClasses[i]).convertToInt();
            }
        }
    }

    float psnr(OCTET* ImgIn, OCTET* ImgOut, int maxSignal, int arraySize)
    {
        float eqm; // erreur quadratique
        for(int i = 0; i < arraySize; i++)
        {
            eqm += pow(ImgIn[i] - ImgOut[i] , 2);
        }
        eqm /= arraySize;

        float a = 10 * log10(static_cast<float>(pow(maxSignal,2)) / eqm );
        return a;
    }

    std::vector<Color> colorPalette_rgb(OCTET* Img, int nTaille)
    {
        int nTaille3 = nTaille*3;

        std::vector<Color> tab;
        bool isIn;

        for(int i = 0; i < nTaille3; i = i +3)
        {
            Color c(Img[i], Img[i+1],Img[i+2]);

            isIn = false;

            for(Color c2 : tab){
                if(c == c2){
                    isIn = true;
                }
            }

            if(!isIn){
                tab.push_back(c);
            }
        }
        return tab;
    }

    void createColorPalette_rgb(OCTET* &Palette, OCTET* Img, int nTaille, int &nW, int &nH)
    {
        int nTaille3 = nTaille*3;

        std::vector<Color> palette = colorPalette_rgb(Img, nTaille3);

        int tailleImg_x = ceil(sqrt(palette.size()));
        int tailleImg = tailleImg_x*tailleImg_x;
        int tailleImg3 = tailleImg*3;
        nW = tailleImg_x;
        nH = tailleImg_x;

        allocation_tableau(Palette, OCTET, tailleImg3);

        int index = 0;
        for(int i = 0; i < palette.size(); i++)
        {
            Palette[index] = palette[i].r;
            Palette[index+1] = palette[i].g;
            Palette[index+2] = palette[i].b;
            index = index + 3;
        }
    }

    void RGBtoYCbCr(OCTET *ImgIn, OCTET *ImgOut, int taille)
    {
        for (int i = 0; i < taille; i++)
        {

            ImgOut[3 * i] = std::min(255., std::max(0., 0.299 * ImgIn[3 * i] + 0.587 * ImgIn[3 * i + 1] + 0.114 * ImgIn[3 * i + 2]));
            ImgOut[3 * i + 1] = std::min(255., std::max(0., -0.1687 * ImgIn[3 * i] - 0.3313 * ImgIn[3 * i + 1] + 0.5 * ImgIn[3 * i + 2] + 128));
            ImgOut[3 * i + 2] = std::min(255., std::max(0., 0.5 * ImgIn[3 * i] - 0.4187 * ImgIn[3 * i + 1] - 0.0813 * ImgIn[3 * i + 2] + 128));
        }
    }

    void YCbCrtoRGB(OCTET *ImgIn, OCTET *ImgOut, int taille)
    {
        for (int i = 0; i < taille; i++)
        {
            ImgOut[3 * i] = std::min(255., std::max(0., ImgIn[3 * i] + 1.402 * (ImgIn[3 * i + 2] - 128)));
            ImgOut[3 * i + 1] = std::min(255., std::max(0., ImgIn[3 * i] - 0.34414 * (ImgIn[3 * i + 1] - 128) - 0.71414 * (ImgIn[3 * i + 2] - 128)));
            ImgOut[3 * i + 2] = std::min(255., std::max(0., ImgIn[3 * i] + 1.772 * (ImgIn[3 * i + 1] - 128)));
        }
    }

    // erosion : 0 = dilation, 1 = erosion
    void ero_dilat_nonBinary(OCTET *ImgIn, OCTET *ImgOut, int nH, int nW, bool color, bool whiteObject, bool erosion){

        //Check si on est sur tu greyscale ou pas
        int facteur = color ? 3 : 1;
        int boucle = color ? 3 : 1;
        //Boucle sur chaque pixel de l'image
        for (int i = 0; i < nH; i++)
        {
            for (int j = 0; j < nW; j++)
            {
                for (int k = 0; k < boucle; k++)
                {
                    int M;
                    if((whiteObject && erosion) || (!whiteObject && !erosion)){
                        M = 255;
                    } else if((whiteObject && !erosion) || (!whiteObject && erosion)){
                        M = 0;
                    }
                    //boucle sur tout le voisinage de i-1 j -1 a i+1 j+1 en skippant les bords de l'image
                    for (int y = std::max(i - 1, 0); y <= std::min(i + 1, nH - 1); y++)
                    {
                        for (int x = std::max(j - 1, 0); x <= std::min(j + 1, nW - 1); x++)
                        {
                            if((whiteObject && erosion) || (!whiteObject && !erosion)){
                                M = std::min(M, (int)ImgIn[facteur * (y * nW + x) + k]);
                            } else if((whiteObject && !erosion) || (!whiteObject && erosion)){
                                M = std::max(M, (int)ImgIn[facteur * (y * nW + x) + k]);
                            }

                        }
                    }
                    ImgOut[facteur * (i * nW + j) + k] = M;
                }
            }
        }
    }

    // color : 0 = greyScale, 1 = RGB
    // whiteObject : 0 = black object, 1 = white object
    void erosion_nonBinary(OCTET *ImgIn, OCTET *ImgOut, int nH, int nW, bool color, bool whiteObject)
    {
        ero_dilat_nonBinary(ImgIn, ImgOut, nH, nW, color, whiteObject, true);
    }

    // color : 0 = greyScale, 1 = RGB
    // whiteObject : 0 = black object, 1 = white object
    void dilation_nonBinary(OCTET *ImgIn, OCTET *ImgOut, int nH, int nW, bool color, bool whiteObject)
    {
        ero_dilat_nonBinary(ImgIn, ImgOut, nH, nW, color, whiteObject, false);
    }

    // color : 0 = greyScale, 1 = RGB
    // whiteObject : 0 = black object, 1 = white object
    void closing_nonBinary(OCTET *ImgIn, OCTET *ImgOut, int nH, int nW, bool color, bool whiteObject)
    {
        OCTET *ImgTemp;
        allocation_tableau(ImgTemp, OCTET, nW * nH * (color ? 3 : 1));

        dilation_nonBinary(ImgIn, ImgTemp, nH, nW, color, whiteObject);
        erosion_nonBinary(ImgTemp, ImgOut, nH, nW, color, whiteObject);

        free(ImgTemp);
    }

    // color : 0 = greyScale, 1 = RGB
    // whiteObject : 0 = black object, 1 = white object
    void opening_nonBinary(OCTET *ImgIn, OCTET *ImgOut, int nH, int nW, bool color, bool whiteObject)
    {
        OCTET *ImgTemp;
        allocation_tableau(ImgTemp, OCTET, nW * nH * (color ? 3 : 1));

        erosion_nonBinary(ImgIn, ImgTemp, nH, nW, color, whiteObject);
        dilation_nonBinary(ImgTemp, ImgOut, nH, nW, color, whiteObject);

        free(ImgTemp);
    }

    // color : 0 = greyScale, 1 = RGB
    // whiteObject : 0 = black object, 1 = white object
    void boundary_nonBinary(OCTET *ImgIn, OCTET *ImgOut, int nH, int nW, bool color, bool whiteObject)
    {
        OCTET *ImgEro;
        OCTET *ImgDil;
        allocation_tableau(ImgEro, OCTET, nW * nH * (color ? 3 : 1));
        allocation_tableau(ImgDil, OCTET, nW * nH * (color ? 3 : 1));

        erosion_nonBinary(ImgIn, ImgEro, nH, nW, color, whiteObject);
        dilation_nonBinary(ImgIn, ImgDil, nH, nW, color, whiteObject);

        int S = 235;

        for (int i = 0; i < nW * nH * (color ? 3 : 1); ++i)
        {
            int val = 255 - (abs(ImgDil[i] - ImgEro[i])) / 2;
            ImgOut[i] = (val < S ? 0 : 255); // seuil
            //ImgOut[i] = val; // sinon juste la différence
        }

        free(ImgEro);
        free(ImgDil);
    }

    void writeHistoDatFile(OCTET* ImgIn, int arraySize, std::string path, bool isRGB)
    {
        std::ofstream file;

        file.open (path);
        if(!file.is_open()){
            std::cout << "Writer - Open File Error: " << path << std::endl;
            return;
        }

        int max = 255;

        std::vector<int> t_grey, t_r, t_g, t_b;

        if(isRGB)
        {
            t_r = std::vector<int>(max);
            t_g = std::vector<int>(max);
            t_b = std::vector<int>(max);

            for(int i = 0; i < arraySize; i = i+3)
            {
                t_r[ImgIn[i]]++;
                t_g[ImgIn[i+1]]++;
                t_b[ImgIn[i+2]]++;
            }

            for(int i = 0; i < max; i++)
            {
                file << std::to_string(i+1) + " " + std::to_string(t_r[i])
                                                + " " + std::to_string(t_g[i])
                                                + " " + std::to_string(t_b[i]) + "\n";
            }
        } else
        {
            t_grey = std::vector<int>(max);

            for (int i = 0; i < arraySize; i++) {
                t_grey[ImgIn[i]]++;
            }

            for (int i = 0; i < max; i++)
            {
                file << std::to_string(i + 1) + " " + std::to_string(t_grey[i]) + "\n";
            }
        }

        file.close();
    }

    void writeDistibutionDatFile(OCTET* ImgIn, int arraySize, std::string path, bool isRGB)
    {
        std::ofstream file;

        file.open (path);
        if(!file.is_open()){
            std::cout << "Writer - Open File Error: " << path << std::endl;
            return;
        }

        int max = 255;

        std::vector<float> t_grey, t_r, t_g, t_b;

        if(isRGB)
        {
            t_r = std::vector<float>(max);
            t_g = std::vector<float>(max);
            t_b = std::vector<float>(max);

            for(int i = 0; i < arraySize; i = i+3)
            {
                t_r[ImgIn[i]]++;
                t_g[ImgIn[i+1]]++;
                t_b[ImgIn[i+2]]++;
            }

            for(int i = 0; i < max; i++)
            {
                file << std::to_string(i+1) + " " + std::to_string(t_r[i] / arraySize)
                                            + " " + std::to_string(t_g[i] / arraySize)
                                            + " " + std::to_string(t_b[i] / arraySize) + "\n";
            }
        } else
        {
            t_grey = std::vector<float>(max);

            for (int i = 0; i < arraySize; i++) {
                t_grey[ImgIn[i]]++;
            }

            for (int i = 0; i < max; i++)
            {
                file << std::to_string(i + 1) + " " + std::to_string(t_grey[i] / arraySize) + "\n";
            }
        }

        file.close();
    }

    void diffMap(OCTET* ImgIn, OCTET* ImgOut, int nH, int nW)
    {
        std::vector<int> diff(nH*nW);
    
        for(int x = 1; x < nW; x++)
        {
            for(int y = 1; y < nH; y++)
            {
                // A + B / 2
                diff[y * nW + x] = (ImgIn[y * nW + x-1] + ImgIn[(y-1) * nW + x]) / 2;
            }
        }

        for(int i = 0; i < nW*nH; i++)
        {
            ImgOut[i] = diff[i] - ImgIn[i] + 128;
        }
    }

    void reconstructFromDiffMap(OCTET* diffMap, OCTET* ImgOut, int nH, int nW)
    {
        std::vector<int> reconstruct(nH*nW);
    
        for(int x = 1; x < nW; x++)
        {
            for(int y = 1; y < nH; y++)
            {
                // A + B / 2
                reconstruct[y * nW + x] = (diffMap[y * nW + x-1] + diffMap[(y-1) * nW + x]) / 2;
            }
        }

        for(int i = 0; i < nW*nH; i++)
        {
            ImgOut[i] = reconstruct[i] + diffMap[i] + 128;
        }
    }


}

#endif //CODAGE_IMAGEALGORITHMS_H
