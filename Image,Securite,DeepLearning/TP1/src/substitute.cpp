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
    if (argc != 3)
    {
        printf("Usage: ImageIn.pgm Key\n");
        exit(1);
    }

    std::string inputName = argv[1];
    inputName = "../in/" + inputName;
    int nH, nW, nTaille;
    unsigned int key;
    key = atoi(argv[2]) % 256;
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

    std::cout << std::endl;
    std::cout << "La clef est : " << key << std::endl;
    std::cout << std::endl;

    // Traitement
    std::cout << "Chiffrement par substitution" << std::endl;
    std::cout << std::endl;
    ImageAlgorithms::substitute(ImgIn, ImgOut, nTaille, key);

    //Ecriture
    ImageAlgorithms::writeHistoDatFile(ImgOut, nTaille, "../dat/s_post.dat", false);
    ecrire_image_pgm("../out/Chiffrement_Substitution.pgm", ImgOut, nH, nW);

    //Infos
    std::cout << "Entropy de base: " << ImageAlgorithms::entropy(ImgIn, nTaille) << std::endl;
    std::cout << "PSNR après substitution: " << ImageAlgorithms::psnr(ImgIn, ImgOut, 255, nTaille) << std::endl;
    std::cout << "Entropy après substitution: " << ImageAlgorithms::entropy(ImgOut, nTaille) << std::endl;

    free(ImgIn);
    free(ImgOut);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"