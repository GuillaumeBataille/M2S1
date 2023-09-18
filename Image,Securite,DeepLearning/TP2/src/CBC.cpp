
#include <stdio.h>
#include <iostream>
#include <vector>
#include <limits>
#include <random>
#include "image_ppm.h"
#include "color.h"
#include "ImageAlgorithms.h"
#include "AES.cpp"

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

    OCTET *ImgIn, *ImgOut;

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
    unsigned char Img[nTaille];
    for (int i = 0; i < nTaille; i++)
    {
        Img[i] = ImgIn[i];
    }
    unsigned char key[] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f}; //key example    AES aes(AESKeyLength::AES_128);
    unsigned int plainLen = nTaille * sizeof(unsigned char);
    bool noise = true;
    unsigned char iv[] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                          0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}; //bytes in plaintext
    unsigned char *out, *innew;
    AES aes(AESKeyLength::AES_128); ////128 - key length, can be 128, 192 or 256
    out = aes.EncryptCBC(Img, plainLen, key, iv);
    srand(out[0]);
    if (noise)
        for (int i = 0; i < nTaille; i++)
        {
            if (rand() % 128 == 0)
                out[i] ^= (1 << rand() % 8);
        }

    for (int i = 0; i < nTaille; i++)
    {
        ImgOut[i] = out[i];
    }

    //Ecriture
    ecrire_image_pgm("../out/CBC_C.pgm", ImgOut, nH, nW);
    std::cout << "Pour chiffrage AES méthode CBC :" << std::endl;
    std::cout << "Entropy de base: " << ImageAlgorithms::entropy(ImgIn, nTaille) << " bits/pixel" << std::endl;
    std::cout << "PSNR après permutation: " << ImageAlgorithms::psnr(ImgIn, ImgOut, 256, nTaille) << " dB" << std::endl;
    std::cout << "Entropy après chiffrement: " << ImageAlgorithms::entropy(ImgOut, nTaille) << " bits/pixel" << std::endl;
    std::cout << std::endl;
    //Histogrammes avant et après chiffrement par permutation
    ImageAlgorithms::writeHistoDatFile(ImgIn, nTaille, "../dat/CBC/original.dat", false);
    ImageAlgorithms::writeHistoDatFile(ImgOut, nTaille, "../dat/CBC/crypt.dat", false);

    innew = aes.DecryptCBC(out, plainLen, key, iv);
    for (int i = 0; i < nTaille; i++)
    {
        ImgOut[i] = innew[i];
    }
    ecrire_image_pgm("../out/CBC_D.pgm", ImgOut, nH, nW);

    //Infos
    std::cout << "Pour déchiffrage AES méthode CBC :" << std::endl;
    std::cout << "Entropy de base: " << ImageAlgorithms::entropy(ImgIn, nTaille) << " bits/pixel" << std::endl;
    std::cout << "PSNR après déchiffrement: " << ImageAlgorithms::psnr(ImgIn, ImgOut, 256, nTaille) << " dB" << std::endl;
    std::cout << "Entropy après déchiffrement: " << ImageAlgorithms::entropy(ImgOut, nTaille) << " bits/pixel" << std::endl;
    ImageAlgorithms::writeHistoDatFile(ImgOut, nTaille, "../dat/CBC/decrypt.dat", false);

    free(ImgIn);
    free(ImgOut);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"
