import numpy as np
from keras.preprocessing import image
from keras.models import load_model

# Charger le modèle sauvegardé
loaded_model = load_model('MNIST.h5')

# Charger votre image (assurez-vous qu'elle est prétraitée de la même manière que les données d'entraînement)
img = image.load_img('chemin_vers_votre_image.jpg', target_size=(28, 28))
img_array = image.img_to_array(img)
img_array = np.expand_dims(img_array, axis=0)  # Ajouter une dimension batch

# Effectuer la prédiction
predictions = loaded_model.predict(img_array)

# Obtenir la classe prédite
predicted_class = np.argmax(predictions, axis=1)

print("Classe prédite :", predicted_class)