import cv2
import numpy as np

# Fonction pour effectuer la multiplication C = A.B pour les matrices
def mat_mult(A, B):
    return np.dot(A, B)

# Fonction pour trouver les deux points d'intersection entre une droite et un carré
def intersection(L, Dx, Dy):
    x = [0, Dx - 1]
    y = [0, Dy - 1]
    x_inter = []
    y_inter = []

    if abs(L[0]) > 1e-16:
        b = -L[1] / L[0]
        c = -L[2] / L[0]
        x[0] = b * y[0] + c
        x[1] = b * y[1] + c
    else:
        x[0] = -Dx
        x[1] = -Dx

    if abs(L[1]) > 1e-16:
        a = -L[0] / L[1]
        c = -L[2] / L[1]
        y[0] = a * x[0] + c
        y[1] = a * x[1] + c
    else:
        y[0] = -Dy
        y[1] = -Dy

    for n in range(4):
        if 0 <= x[n] < Dx and 0 <= y[n] < Dy:
            x_inter.append(int(round(x[n])))
            y_inter.append(int(round(y[n])))

    if len(x_inter) == 2:
        return x_inter, y_inter
    else:
        return None

# Charger les deux images TIFF
image1 = cv2.imread("TurtleD.tiff")
image2 = cv2.imread("TurtleG.tiff")

# Créer une fenêtre pour chaque image
cv2.namedWindow("Image 1")
cv2.namedWindow("Image 2")

# Fonction de gestion des clics de souris
def on_mouse_click(event, x, y, flags, param):
    if event == cv2.EVENT_LBUTTONDOWN:
        print(f"Coordonnées dans l'image 1 : ({x}, {y})")

        # Ajoutez ici le code pour effectuer des opérations sur les coordonnées (x, y)

        # Afficher un cercle rouge à l'emplacement du clic dans l'image 2
        cv2.circle(image2, (x, y), 5, (0, 0, 255), -1)
        cv2.imshow("Image 2", image2)

# Associer la fonction de gestion des clics de souris à la fenêtre de l'image 2
cv2.setMouseCallback("Image 2", on_mouse_click)

# Afficher les images
cv2.imshow("Image 1", image1)
cv2.imshow("Image 2", image2)

# Attendre que l'utilisateur appuie sur une touche pour quitter
cv2.waitKey(0)

# Fermer toutes les fenêtres
cv2.destroyAllWindows()