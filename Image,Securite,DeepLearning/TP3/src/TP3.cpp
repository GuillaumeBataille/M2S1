

#include "AES.cpp"
#include "ImageAlgorithms.h"
#include "color.h"
#include "image_ppm.h"
#include <iostream>
#include <limits>
#include <random>
#include <stdio.h>
#include <vector>

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Usage: ImageIn key\n");
        exit(1);
    }

    std::string inputName = argv[1];
    inputName = "../in/" + inputName;
    unsigned int k = atoi(argv[2]);
    int nH, nW, nTaille;

    OCTET *ImgIn, *ImgOut, *binPlan;

    bool color = !(inputName.substr(inputName.size() - 3, 3)).compare("ppm");
    std::cout << "Image " << (color ? "couleur" : "niveau de gris") << " en input"
              << std::endl;

    (color ? lire_nb_lignes_colonnes_image_ppm(inputName, &nH, &nW)
           : lire_nb_lignes_colonnes_image_pgm(inputName, &nH, &nW));
    nTaille = nH * nW;
    int nTaille3 = nTaille * 3;

    allocation_tableau(ImgIn, OCTET, nTaille);
    allocation_tableau(ImgOut, OCTET, nTaille);
    allocation_tableau(binPlan, OCTET, nTaille);

    lire_image_pgm(inputName, ImgIn, nH * nW);

    std::cout << "Taille du message : " << nTaille << std::endl;

    // --------------------- Pré Traitement --------------------- //

    std::cout << std::endl;
    unsigned char out[nTaille];

    srand(clock());

    for (int i = 0; i < nTaille; i++)
    {
        out[i] = (rand() % 2 == 0) ? 0 : 255; // Une chance sur deux
    }

    ImageAlgorithms::writeHistoDatFile(ImgIn, nTaille, "../dat/original.dat", false);

    std::string name = "../out/Msg";
    std::string dat_name = "../dat/msg_";
    std::string plan_name = "../out/plan_";

    for (int k = 0; k < 8; k++)
    {
        ImageAlgorithms::getBinaryPlane(ImgIn, binPlan, nTaille, k);
        ecrire_image_pgm(plan_name + std::to_string(k) + ".pgm", binPlan, nH, nW);

        for (int i = 0; i < nTaille; i++)
        {
            ImgOut[i] = ImgIn[i];
            ImageAlgorithms::setBinaryPlane(ImgOut, i, k, out[i] == 255);
        }
        std::cout << "Traitement de l'image " << k << " :"
                  << "PSNR = : " << ImageAlgorithms::psnr(ImgIn, ImgOut, 256, nTaille) << " dB" << std::endl;
        ecrire_image_pgm(name + std::to_string(k) + ".pgm", ImgOut, nH, nW);
        ImageAlgorithms::writeHistoDatFile(ImgOut, nTaille, dat_name + std::to_string(k) + "/" + std::to_string(k) + ".dat", false);
        std::cout << std::endl;
    }
    // Ecriture
    ecrire_image_pgm("../out/BinaryMessage.pgm", out, nH, nW);

    free(ImgIn);
    free(ImgOut);
    free(binPlan);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"
