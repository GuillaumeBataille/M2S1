

#include "AES.cpp"
#include "ImageAlgorithms.h"
#include "color.h"
#include "image_ppm.h"
#include "Filters.h"
#include <iostream>
#include <limits>
#include <random>
#include <stdio.h>
#include <vector>

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

    allocation_tableau(ImgIn, OCTET, nTaille);

    lire_image_pgm(inputName, ImgIn, nH * nW);

    // --------------------- Pré Traitement --------------------- //

    std::cout << std::endl;

    int size_filter = 5;
    allocation_tableau(ImgOut, OCTET, (nW - size_filter + 1) * (nH - size_filter + 1));

    ImageAlgorithms::ConvolutionWithFilter(ImgIn, ImgOut, nW, nH, flou_moyen_5x5, size_filter);
    ecrire_image_pgm("../out/flou_moyen5x5.pgm", ImgOut, nH - size_filter + 1, nW - size_filter + 1);

    ImageAlgorithms::ConvolutionWithFilter(ImgIn, ImgOut, nW, nH, gaussien_5x5, size_filter);
    ecrire_image_pgm("../out/gaussien_5x5.pgm", ImgOut, nH - size_filter + 1, nW - size_filter + 1);

    ImageAlgorithms::ConvolutionWithFilter(ImgIn, ImgOut, nW, nH, contours_5x5, size_filter);
    ecrire_image_pgm("../out/contours_5x5.pgm", ImgOut, nH - size_filter + 1, nW - size_filter + 1);

    ImageAlgorithms::ConvolutionWithFilter(ImgIn, ImgOut, nW, nH, embossage_5x5, size_filter);
    ecrire_image_pgm("../out/embossage_5x5.pgm", ImgOut, nH - size_filter + 1, nW - size_filter + 1);

    ImageAlgorithms::ConvolutionWithFilter(ImgIn, ImgOut, nW, nH, renforcement_bords_5x5, 5);
    ecrire_image_pgm("../out/renforcement_bords_5x5.pgm", ImgOut, nH - size_filter + 1, nW - size_filter + 1);

    ImageAlgorithms::maxPooling(ImgIn, ImgOut, nW, nH, 2);
    ecrire_image_pgm("../out/PoolingMax_2.pgm", ImgOut, nH / 2, nW / 2);

    ImageAlgorithms::maxPooling(ImgIn, ImgOut, nW, nH, 4);
    ecrire_image_pgm("../out/PoolingMax_4.pgm", ImgOut, nH / 4, nW / 4);

    free(ImgIn);
    free(ImgOut);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"
