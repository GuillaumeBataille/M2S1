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
    if (argc != 2)
    {
        printf("Usage: ImageIn.pgm\n");
        exit(1);
    }

    std::string inputName = argv[1];
    inputName = inputName;
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

    std::cout << "Traitement du brutforcing" << std::endl;
    //On teste toute les keys possible
    for (int key = 1; key < 256; key++)
    {
        srand(key); // Init la seed avec la key
        //Pour chaque pixel de l'image
        for (int i = 0; i < nTaille; i++)
        {
            ImgOut[i] = (ImgIn[i] - ImgIn[i - 1] - rand()) % 256;
        }
        double entropy = ImageAlgorithms::entropy(ImgOut, nTaille);
        if (entropy < 7.9)
        {
            std::cout << "Entropy : " << entropy << " associé a la clef :" << key << std::endl;
            ecrire_image_pgm("../out/brutforced.pgm", ImgOut, nH, nW);
        }
    }
    free(ImgIn);
    free(ImgOut);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"