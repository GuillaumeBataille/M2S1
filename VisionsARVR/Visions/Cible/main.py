from PIL import Image
import os
import re
import cv2
import numpy as np
from threading import Timer

# ------------------------------ Image functions ------------------------------------------- #

# [Numpy.Array, (hauteur,largeur)]
def lire_sequence_tif(dossier):
    images = []  # Crée une liste pour stocker les images en tant que tableaux NumPy

    if not os.path.exists(dossier):
        print(f"Le dossier {dossier} n'existe pas.")
        return None

    # Liste tous les fichiers dans le dossier
    fichiers = os.listdir(dossier)

    # Trie les fichiers par ordre numérique en utilisant un tri naturel
    fichiers = sorted(fichiers, key=lambda x: [int(s) if s.isdigit() else s.lower() for s in re.split('(\d+)', x)])

    for fichier in fichiers:
        if fichier.endswith(".tif"):
            chemin_complet = os.path.join(dossier, fichier)
            image = cv2.imread(chemin_complet)
            if image is not None:
                images.append(image)

    return (images, images[0].shape[:2])

def lire_sequence_bmp(dossier):
    images = []  # Crée une liste pour stocker les images en tant que tableaux NumPy

    if not os.path.exists(dossier):
        print(f"Le dossier {dossier} n'existe pas.")
        return None

    # Liste tous les fichiers dans le dossier
    fichiers = os.listdir(dossier)

    # Trie les fichiers par ordre numérique en utilisant un tri naturel
    fichiers = sorted(fichiers, key=lambda x: [int(s) if s.isdigit() else s.lower() for s in re.split('(\d+)', x)])

    for fichier in fichiers:
        if fichier.endswith(".bmp"):
            chemin_complet = os.path.join(dossier, fichier)
            image = cv2.imread(chemin_complet)
            if image is not None:
                images.append(image)

    return (images, images[0].shape[:2])

def displayImgage(image):
    while True:
        cv2.imshow("Image", image)
        key = cv2.waitKey(1) & 0xFF
        if key == 27:  # Appuyez sur la touche "ESC" pour quitter
            break

    cv2.destroyAllWindows()

def displayImgageForTime(image, t):
    
    
    # Fonction pour fermer la fenêtre après un délai
    def closeWindow():
        cv2.destroyAllWindows()

    # Créez un minuteur pour fermer la fenêtre après display_time_ms millisecondes
    timer = Timer(t, closeWindow)
    timer.start()

    while not timer.finished.is_set():
        cv2.imshow("Image", image)



def displayImgages(images):
    for image in images:
        cv2.imshow("Image", image)
        cv2.waitKey(0)  # Attend indéfiniment jusqu'à ce que la touche soit pressée
        cv2.destroyAllWindows()


def selectPointsFromImage(image):

    if image is None:
        print("Impossible de charger l'image.")
        return

    # Copie l'image pour dessiner les points sans modifier l'originale
    image_copy = image.copy()

    points = []

    def clic_souris(event, x, y, flags, param):
        nonlocal points
        if event == cv2.EVENT_LBUTTONDOWN:
            cv2.circle(image_copy, (x, y), 5, (0, 0, 255), -1)  # Dessine un cercle rouge pour le point
            points.append((x, y))
            cv2.imshow("Image", image_copy)

    cv2.namedWindow("Image")
    cv2.setMouseCallback("Image", clic_souris)

    while True:
        cv2.imshow("Image", image_copy)
        key = cv2.waitKey(1) & 0xFF
        if key == 27:  # Appuyez sur la touche "ESC" pour quitter
            break

    cv2.destroyAllWindows()

    return points

def selectRectanglePointsFromImage(image):
    def draw_rectangle(event, x, y, flags, param):
        nonlocal rect_start, rect_end, rect_center, drawing

        if event == cv2.EVENT_LBUTTONDOWN:
            rect_start = (x, y)
            drawing = True
        elif event == cv2.EVENT_MOUSEMOVE and drawing:
            img_copy = image.copy()
            cv2.rectangle(img_copy, rect_start, (x, y), (0, 255, 0), 2)
            cv2.imshow("Image", img_copy)
        elif event == cv2.EVENT_LBUTTONUP:
            rect_end = (x, y)

            if(rect_end[0] < rect_start[0]):
                tmp = rect_start
                rect_start = rect_end
                rect_end = tmp

            rect_center = ( rect_start[0] + (rect_end[0] - rect_start[0]) // 2, rect_start[1] + (rect_end[1] - rect_start[1]) // 2)

            drawing = False
            img_copy = image.copy()
            cv2.rectangle(img_copy, rect_start, rect_end, (0, 255, 0), 2)
            cv2.circle(img_copy, (rect_center[0], rect_center[1]), 4, (0, 0, 255), -1)
            cv2.imshow("Image", img_copy)

    cv2.namedWindow("Image")

    rect_start = None
    rect_end = None
    rect_center = None
    drawing = False

    cv2.imshow("Image", image)
    cv2.setMouseCallback("Image", draw_rectangle)

    print("press Enter to select Area")

    while True:
        key = cv2.waitKey(1) & 0xFF

        if key == ord("c"):
            img_copy = image.copy()
            cv2.imshow("Image", img_copy)
            rect_start = None
            rect_end = None

        if key == 13: # échap
            break

    cv2.destroyAllWindows()

    return rect_start, rect_end, rect_center

# La ya un problème des fois ça fait des tableau de taille x,y-1 sans raison
def extractPattern(img, point1, point2):

    x1, y1 = point1
    x2, y2 = point2

    x1, x2 = min(x1, x2), max(x1, x2)
    y1, y2 = min(y1, y2), max(y1, y2)

    extracted_region = img[int(y1):int(y2), int(x1):int(x2)]

    return extracted_region

def extractRegion(img, point1, point2):

        x1, y1 = point1
        x2, y2 = point2

        x1, x2 = min(x1, x2), max(x1, x2)
        y1, y2 = min(y1, y2), max(y1, y2)

        #extracted_region = img[y1:y2, x1:x2]

        extracted_region = np.zeros((y2-y1, x2-x1), dtype=np.uint8)

        for x in range (x1, x2):
            for y in range(y1, y2):
                x_tab = x - x1
                y_tab = y - y1

                extracted_region[y_tab, x_tab] = img[y,x]

        Ny2, Nx2 = extracted_region.shape[0], extracted_region.shape[1]

        # if Ny2 != 53:
        #     print('b')

        return extracted_region

def drawPointOnImg(img, point):

    image_modifiee = img.copy()

    couleur = (0, 0, 255)
    epaisseur = 2

    x, y = point


    cv2.circle(image_modifiee, (int(x), int(y)), epaisseur, couleur, -1)

    return image_modifiee


def drawRectangleOnImg(img, point1, point2):
    image_modifiee = img.copy()

    couleur = (0, 0, 255)  # Couleur du rectangle (bleu, vert, rouge)
    epaisseur = 2  # Épaisseur des lignes du rectangle

    x1, y1 = point1  # Coin supérieur gauche du rectangle
    x2, y2 = point2  # Coin inférieur droit du rectangle

    cv2.rectangle(image_modifiee, (x1, y1), (x2, y2), couleur, epaisseur)

    return image_modifiee


def translatePoints(points, uv):
    p1 = (points[0][0]+uv[0], points[0][1]+uv[1])
    p2 = (points[1][0]+uv[0], points[1][1]+uv[1])
    p3 = (points[2][0]+uv[0], points[2][1]+uv[1])
    
    return (p1,p2,p3)

def clampBiggestorLittlest(val):
    if(val < 0):
        return np.floor(val)
    else:
        return np.ceil(val)

# ----------------------------------------------------------------------------------------- #

class CorrelationClass:

    base_image = None
    tolerance = None

    windowSize_h = None
    windowSize_l = None

    moy_base = 0
    ecarT_base = 0

    mode = 0

    # mode | 0 : SSD_distance, 1 : SAD_distance, 2 : PearsonCorrelation
    def __init__(self, img, tolerance, mode):
        self.tolerance = tolerance
        self.base_image = img
        self.mode = mode

        self.windowSize_h, self.windowSize_l = img.shape[:2]

        self.moy_base = np.mean(self.base_image)

        Nx, Ny = self.base_image.shape[0], self.base_image.shape[1]
        for x in range(Nx):
            for y in range(Ny):
                self.ecarT_base += (self.base_image[x,y] - self.moy_base) ** 2

        self.ecarT_base /= Nx*Ny

    def CorrelationPearson(self, img):
        Ny, Nx = self.base_image.shape[0], self.base_image.shape[1]

        N = Nx * Ny
        moy1 = np.mean(self.base_image)
        moy2 = np.mean(img)

        var1 = np.sum((self.base_image - moy1) ** 2)
        var2 = np.sum((img - moy2) ** 2)


        if var1 == 0 or var2 == 0:
            return 0.0

        resultat = np.sum((self.base_image - moy1) * (img - moy2))

        correlation = resultat / np.sqrt(var1 * var2)

        return correlation

        # moy1 = self.moy_base
        # moy2 = np.mean(img)

        # ecartT1 = self.ecarT_base
        # ecartT2 = 0

        # for x in range(Nx):
        #     for y in range(Ny):
        #         ecartT2 += (img[x,y] - moy2) ** 2

        
        # ecartT2 /= N

        # correl = 0

        # part1 = 0
        # for x in range(Nx):
        #     for y in range(Ny):
        #         part1 += (self.base_image[x,y] - moy1) * (img[x,y] - moy2) 

        # correl = (part1/np.sqrt(ecartT1) * np.sqrt(ecartT2)) / N

        # return correl

    
    def npCorrelation(self, extractedRegion):
        f_base = self.base_image.flatten()
        f_extracted = extractedRegion.flatten()

        cor = np.corrcoef(f_base, f_extracted)

        retCor = 0

        for c in cor:
            retCor += c
        
        return retCor/len(c)


    def distanceSAD(self, img2):
        Ny, Nx = self.base_image.shape[0], self.base_image.shape[1]
        return np.sum(np.abs(self.base_image - img2)) / Nx * Ny


    def distanceSSD(self, img2):
        Ny, Nx = self.base_image.shape[0], self.base_image.shape[1]
        return np.sum(np.abs(self.base_image - img2)) ** 2 / Nx * Ny
    
    # previewWindowPos = le coin haut gauche de la window
    def findNewPos(self, img, previewWindowPos):
        
        img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

        img_h, img_l = img.shape[:2] 

        h_min = max(0, previewWindowPos[1]-self.tolerance)
        h_max = min(img_h - self.windowSize_h, previewWindowPos[1]+self.tolerance)
        l_min = max(0, previewWindowPos[0]-self.tolerance)
        l_max = min(img_l - self.windowSize_l, previewWindowPos[0]+self.tolerance)

        newPos = previewWindowPos
        correl = float('-inf')
        dist = float('inf')

        for x in range(l_min, l_max):
            for y in range(h_min, h_max):
                
                point1 = (x, y)

                point2 = (self.windowSize_l+x, self.windowSize_h+y)
                extractedRegion = extractRegion(img_gray, point1, point2)

                if self.mode == 0:
                    
                    newDist = self.distanceSSD(extractedRegion)

                    if newDist < dist:

                        newPos = point1
                        dist = newDist
                
                elif self.mode == 1:

                    newDist = self.distanceSAD(extractedRegion)

                    if newDist < dist:
                        newPos = point1
                        dist = newDist

                elif self.mode == 2:

                    newCorrel = self.CorrelationPearson(extractedRegion)

                    if newCorrel > correl and newCorrel >= 0:

                        correl = newCorrel
                        newPos = point1
        
        return newPos





class OpticalFlow:
    
    def __init__(self, img):

        img_gray = img

        self.hauteur, self.largeur = img_gray.shape[:2]
        self.nPix = self.largeur * self.hauteur

        self.A = None
        self.Ainv = None
        self.B = np.zeros((self.nPix, 1))
        self.X = None

        self.h = None
        self.l = None


        self.createMatA(img_gray)


    def gradient(self, img):
        sobelx = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=5)
        sobely = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=5)
        return sobelx, sobely

    def createMatA(self, img):

        self.A = np.zeros((self.nPix, 2))

        Ix, Iy = self.gradient(img)

        for x in range(self.hauteur):
            for y in range(self.largeur):
                self.A[x*2 + y] = [Ix[x,y],Iy[x,y]]

        # calcule la pseudo inverse de A
        self.Ainv = np.linalg.pinv(self.A)

    def setB(self, nextImg):
        img_gray = cv2.cvtColor(nextImg, cv2.COLOR_BGR2GRAY)
        img_gray = img_gray.flatten()

        for i in range(len(img_gray)):
            self.B[i] = img_gray[i]


    def computeX(self):
        self.X = np.dot(self.Ainv, self.B)

        return (self.X[0, 0], self.X[1, 0])


# --------------------------------- Main -------------------------------------

dossier_images = "Ghost4"
sequence_images = lire_sequence_bmp(dossier_images)

hauteur = sequence_images[1][0]
largeur = sequence_images[1][1]


selectedPoints = selectRectanglePointsFromImage(sequence_images[0][0])


# selectedPoints = ((143,51),(184,94),(161,75))
#selectedPoints = ((135,50),(195,100),(161,75))

window_size_l = selectedPoints[1][0] - selectedPoints[0][0]
window_size_h = selectedPoints[1][1] - selectedPoints[0][1]

base_image = cv2.cvtColor(sequence_images[0][0], cv2.COLOR_BGR2GRAY)

base_window = extractRegion(base_image, selectedPoints[0], selectedPoints[1])

Ny, Nx = base_window.shape[0], base_window.shape[1]

# ------------------------------------------ Méthode de la correlation ------------------------------------

corC = CorrelationClass(base_window, 20, 2)

upleftWindowCoord = selectedPoints[0]


count = 0
for img in sequence_images[0]:

    upleftWindowCoord = corC.findNewPos(img, upleftWindowCoord)


    print(count)

    uv = (upleftWindowCoord[0] + window_size_l//2, upleftWindowCoord[1] + window_size_h//2)

    #imgPoint = drawPointOnImg(img, uv)

    point2 = (upleftWindowCoord[0] + window_size_l, upleftWindowCoord[1] + window_size_h)
    imgPoint = drawRectangleOnImg(img, upleftWindowCoord, point2)

    #displayImgage(imgPoint)

    outpath = './out'
    cv2.imwrite(os.path.join(outpath , "p" + str(count) + '_2.jpg'), imgPoint)

    count += 1




# ------------------------------------------ Méthode du flux optique ------------------------------------

# of = OpticalFlow(base_window)

# img0 = sequence_images[0][0]

# count = 0
# for img in sequence_images[0]:
    
#     new_window = extractPattern(img, selectedPoints[0], selectedPoints[1])
#     of.setB(new_window)
#     uv = of.computeX()

#     uvx = clampBiggestorLittlest(uv[0])
#     uvy = clampBiggestorLittlest(uv[1])

#     print("uv: ", uv)
#     uv = (uvx, uvy)

#     print("selectedPoints: ", selectedPoints)
#     print("uv: ", uv)

#     selectedPoints = translatePoints(selectedPoints, uv)

#     print("selectedPoints MODIF: ", selectedPoints)

#     imgPoint = drawPointOnImg(img, selectedPoints[2])

#     outpath = './out/flux'
#     cv2.imwrite(os.path.join(outpath , str(count) + '.jpg'), imgPoint)

#     count += 1
