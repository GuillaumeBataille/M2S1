import numpy as np
from keras.preprocessing import image
from keras.models import load_model
import cv2  # OpenCV for image inversion

# Charger le modèle sauvegardé
loaded_model = load_model('MNIST.keras')

# Charger votre image (assurez-vous qu'elle est prétraitée de la même manière que les données d'entraînement)
img = image.load_img('4_2.png', target_size=(28, 28), color_mode="grayscale")

# Convert the image to a numpy array
img_array = image.img_to_array(img)

# Invert the colors (black becomes white, white becomes black)
inverted_img_array = 255 - img_array

# Normalize the inverted image pixel values to be in the range [0, 1]
inverted_img_array /= 255.0

inverted_img_array = np.expand_dims(inverted_img_array, axis=0)  # Add a batch dimension

# Effectuer la prédiction
predictions = loaded_model.predict(inverted_img_array)

# Obtenir la classe prédite
predicted_class = np.argmax(predictions, axis=1)

print("Classe prédite :", predicted_class)
