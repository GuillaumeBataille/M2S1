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
    inputName = inputName;
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
    std::cout << "La clef est : " << key << std::endl;

    OCTET *Map;
    allocation_tableau(Map, OCTET, nTaille);

    // Traitement
    std::cout << "Traitement 1 : Chiffrement par permutation" << std::endl;
    srand(key); // Init la seed avec la key

    //Pour chaque pixel de l'image de base,

    for (int i = 0; i < nTaille; i++)
    {
        int rand_pos = rand() % (nTaille - 1); // Notre position random
        // Si la position random est libre
        if (Map[rand_pos] != 255)
        {
            Map[rand_pos] = 255;         // On colore la map pour dire que c'est occupé
            ImgOut[rand_pos] = ImgIn[i]; // On permute
        }
        else // La pos random est occupée
        {
            int cpt = rand_pos;
            while (Map[cpt] == 255) // Tant qu'on est sur du blanc (donc des pixel occupés)
            {
                cpt = (cpt + 1) % nTaille; // On regarde le pixel suivant
            }
            Map[rand_pos] = 255;         // On colore la map pour dire que c'est occupé
            ImgOut[rand_pos] = ImgIn[i]; // On permute
        }
    }
    //Histogrammes avant et après chiffrement par permutation
    ImageAlgorithms::writeHistoDatFile(ImgIn, nTaille, "../dat/p_base.dat", false);
    ImageAlgorithms::writeHistoDatFile(ImgOut, nTaille, "../dat/p_post.dat", false);

    //Ecriture
    ecrire_image_pgm("../out/Chiffrement_Permutation.pgm", ImgOut, nH, nW);
    ecrire_image_pgm("../out/CartedesDisponibilites.pgm", Map, nH, nW);

    //-----------------//

    std::cout << "Traitement 2 : Chiffrement par substitution" << std::endl;
    //Pour le pixel initial :
    ImgOut[0] = (key + ImgIn[0]) % 256;

    for (int i = 1; i < nTaille; i++)
    {
        ImgOut[i] = (ImgOut[i - 1] + ImgIn[i] + rand()) % 256;
    }
    ImageAlgorithms::writeHistoDatFile(ImgOut, nTaille, "../dat/s_post.dat", false);
    ecrire_image_pgm("../out/Chiffrement_Substitution.pgm", ImgOut, nH, nW);
    free(ImgIn);
    free(ImgOut);

    return 1;
}

// plot "../dat/p_base.dat" lt rgb "blue" with lines title "Histo before"