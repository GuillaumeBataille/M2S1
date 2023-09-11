// test_couleur.cpp : Seuille une image couleur

#include <stdio.h>
#include <iostream>
#include <vector>
#include <limits>
#include <random>
#include "image_ppm.h"
#include "color.h"
#include "ImageAlgorithms.h"

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Usage: ImageIn.pgm Key\n");
        exit(1);
    }

    std::string inputName = argv[1];
    inputName = "../in/"+inputName;
    int nH, nW, nTaille, key;

    key = atoi(argv[2]);
    OCTET *ImgIn, *ImgOut;
    std::cout<<"input name "<< inputName<<std::endl;
    bool color = !(inputName.substr(inputName.size() - 3, 3)).compare("ppm");
    std::cout << "Image " << (color ? "couleur" : "niveau de gris") << " en input" << std::endl;

    (color ? lire_nb_lignes_colonnes_image_ppm(inputName, &nH, &nW) : lire_nb_lignes_colonnes_image_pgm(inputName, &nH, &nW));
    nTaille = nH * nW;
    int nTaille3 = nTaille * 3;
    std::cout<<"Num" <<color <<std::endl;

    allocation_tableau(ImgIn, OCTET, nTaille);

    allocation_tableau(ImgOut, OCTET, nTaille);
    lire_image_ppm(inputName, ImgIn, nH * nW);

    // Traitement
    std::cout<<"Here"<<std::endl;



    ecrire_image_pgm("../out/testtruth.pgm", ImgOut, nH, nW);
    free(ImgIn);
    free(ImgOut);

    return 1;
}