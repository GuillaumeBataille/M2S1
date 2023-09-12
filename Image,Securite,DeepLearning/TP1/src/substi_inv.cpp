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
    std::cout << "argc" << argc << std::endl;
    if (argc != 1)
    {
        printf("Usage:\n");
        exit(1);
    }

    std::string inputName = "../out/Chiffrement_Substitution.pgm";
    int nH, nW, nTaille;

    OCTET *ImgIn, *ImgOut;
    std::cout << "input name " << inputName << std::endl;
    bool color = !(inputName.substr(inputName.size() - 3, 3)).compare("ppm");
    std::cout << "Image " << (color ? "couleur" : "niveau de gris") << " en input" << std::endl;

    (color ? lire_nb_lignes_colonnes_image_ppm(inputName, &nH, &nW) : lire_nb_lignes_colonnes_image_pgm(inputName, &nH, &nW));
    nTaille = nH * nW;
    int nTaille3 = nTaille * 3;

    allocation_tableau(ImgIn, OCTET, nTaille);

    allocation_tableau(ImgOut, OCTET, nTaille);

    lire_image_pgm(inputName, ImgIn, nH * nW);

    //Pré Traitement

    // Traitement

    std::cout << "Substitution inverse par brutforcing" << std::endl;
    ImageAlgorithms::reverse_substitute_brutforce(ImgIn, ImgOut, nTaille);

    ecrire_image_pgm("../out/brutforced.pgm", ImgOut, nH, nW);

    free(ImgIn);
    free(ImgOut);

    return 1;
}
