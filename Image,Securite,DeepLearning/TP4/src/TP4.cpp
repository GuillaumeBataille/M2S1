

#include "AES.cpp"
#include "ImageAlgorithms.h"
#include "color.h"
#include "image_ppm.h"
#include <iostream>
#include <limits>
#include <random>
#include <stdio.h>
#include <vector>
#include <sstream>

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: ImageIn\n");
        exit(1);
    }

    std::string inputName = argv[1];
    inputName = "../in/" + inputName;

    // unsigned int k = atoi(argv[2]);
    int nH, nW, nTaille;

    OCTET *ImgIn, *ImgOut;

    bool color = !(inputName.substr(inputName.size() - 3, 3)).compare("ppm");
    std::cout << "Image " << (color ? "couleur" : "niveau de gris") << " en input"
              << std::endl;

    (color ? lire_nb_lignes_colonnes_image_ppm(inputName, &nH, &nW)
           : lire_nb_lignes_colonnes_image_pgm(inputName, &nH, &nW));
    nTaille = nH * nW;
    int nTaille3 = nTaille * 3;

    allocation_tableau(ImgIn, OCTET, nTaille3);
    lire_image_ppm(inputName,ImgIn,nTaille);
    //ecrire_image_ppm("../out/test.ppm", ImgIn, nH ,nW);
    // --------------------- Pré Traitement --------------------- //

    std::cout << std::endl;

    int size_filter = 5;
    int pooling_size = 2;
    int width_out = (nW - size_filter + 1);
    int height_out = (nH - size_filter + 1);
    int size_out= width_out * height_out;

    std::vector<OCTET*> ListofImagesPlane = ImageAlgorithms::getImageFromFolder("../in/plane/", 32, 32);
    std::vector<OCTET*> ListofImagesCars = ImageAlgorithms::getImageFromFolder("../in/cars/", 32, 32);

    std::vector<OCTET*> ListofImages = ListofImagesPlane;

    for (int i = 0; i< ListofImagesCars.size();i++)
    {
       ListofImages.push_back(ListofImagesCars[i]);
    }


    ImageAlgorithms::CNN(ListofImages,size_filter, pooling_size,2,nW,nH); // CNN a 2 couches

    OCTET * FlattenedVector = ImageAlgorithms::concatImg(ListofImages,5,5); // Le vecteur applatit

    ecrire_image_ppm("../out/FlattenedVector.ppm", FlattenedVector,1,1000);
    
    std::vector<double> Result = ImageAlgorithms::CNN_Output(FlattenedVector,2,ListofImages.size()*5*5*3);
    std::vector<double> Result_softMaxed = ImageAlgorithms::softmax(Result);
    
    std::cout<<"Classe Plane "<< Result_softMaxed[0]<<std::endl;
    std::cout<<"Classe Cars "<< Result_softMaxed[1]<<std::endl;


    //Test des Convolutions pour Compte rendu
    allocation_tableau(ImgOut, OCTET, 3*size_out);

    ImageAlgorithms::ConvolutionWithFilterRGB(ImgIn, ImgOut, nW, nH, flou_moyen_5x5, size_filter);
    ecrire_image_ppm("../out/flou_moyen5x5.ppm", ImgOut, height_out, height_out);

    ImageAlgorithms::ConvolutionWithFilterRGB(ImgIn, ImgOut, nW, nH, gaussien_5x5, size_filter);
    ecrire_image_ppm("../out/gaussien_5x5.pgm", ImgOut, height_out, height_out);

    ImageAlgorithms::ConvolutionWithFilterRGB(ImgIn, ImgOut, nW, nH, contours_5x5, size_filter);
    ecrire_image_ppm("../out/contours_5x5.pgm", ImgOut, height_out, height_out);

    ImageAlgorithms::ConvolutionWithFilterRGB(ImgIn, ImgOut, nW, nH, embossage_5x5, size_filter);
    ecrire_image_ppm("../out/embossage_5x5.pgm", ImgOut, height_out, height_out);

    ImageAlgorithms::ConvolutionWithFilterRGB(ImgIn, ImgOut, nW, nH, renforcement_bords_5x5, size_filter);
    ecrire_image_ppm("../out/renforcement_bords_5x5.pgm", ImgOut, height_out, height_out);

    ImageAlgorithms::maxPoolingRGB(ImgIn, ImgOut, nW, nH, 2);
    ecrire_image_ppm("../out/PoolingMax_2.pgm", ImgOut, nH / 2, nW / 2);

    ImageAlgorithms::maxPoolingRGB(ImgIn, ImgOut, nW, nH, 4);
    ecrire_image_ppm("../out/PoolingMax_4.pgm", ImgOut, nH / 4, nW / 4);

for(OCTET* img : ListofImages){
        free(img);
    }
    free(ImgIn);
    free(ImgOut);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"
