
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
    std::cout << "Chiffrement par permutation" << std::endl;
    std::cout << std::endl;
    ImageAlgorithms::permute(ImgIn, ImgOut, nTaille, key);

    //Histogrammes avant et après chiffrement par permutation
    ImageAlgorithms::writeHistoDatFile(ImgIn, nTaille, "../dat/p_base.dat", false);
    ImageAlgorithms::writeHistoDatFile(ImgOut, nTaille, "../dat/p_post.dat", false);

    //Ecriture
    ecrire_image_pgm("../out/Chiffrement_Permutation.pgm", ImgOut, nH, nW);

    //Infos
    std::cout << "Entropy de base: " << ImageAlgorithms::entropy(ImgIn, nTaille) << std::endl;
    std::cout << "PSNR après permutation: " << ImageAlgorithms::psnr(ImgIn, ImgOut, 256, nTaille) << std::endl;
    std::cout << "Entropy après permutation: " << ImageAlgorithms::entropy(ImgOut, nTaille) << std::endl;
    std::cout << std::endl;

    free(ImgIn);
    free(ImgOut);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"
