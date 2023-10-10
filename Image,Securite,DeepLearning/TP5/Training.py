#Import des framwork DeepLearning : TensorFlow et KERAS
import tensorflow as tf
from tensorflow import keras

#Data science items
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

#Librairies utiles
import sys,os
from importlib import reload
sys.path.append('..')

#Import LeakyReLu
from keras.layers import LeakyReLU

#Import le save/load model
from keras.models import load_model, save_model

#Récupérer les données CMIST : Des images de chiffres (0-9)
(x_train,y_train),(x_test,y_test) = keras.datasets.mnist.load_data()

#Resize des images pour un format normalisé 28 28
x_train = x_train.reshape(-1,28,28,1)
x_test = x_test.reshape(-1,28,28,1)

#Affichage du train post-shape
print("x_train : ",x_train.shape)
print("y_train : ",y_train.shape)
#Affichage du test post-shape
print("x_test : ",x_test.shape)
print("y_test : ",y_test.shape)

'''
Valeur print :
x_train :  (60000, 28, 28, 1)
y_train :  (60000,)
x_test :  (10000, 28, 28, 1)
y_test :  (10000,)
'''

#Normalisation des pixels (entre 0 et 1) en récupérant le max
print('Avant la normalisation : Min = {}, Max ={}'.format(x_train.min(),x_train.max()))
'''
Valeur print : Avant la normalisation : Min = 0, Max =255
'''
xmax = x_train.max()
x_train = x_train / xmax
x_test = x_test /xmax
print('Après la normalisation : Min = {}, Max ={}'.format(x_train.min(),x_train.max()))
'''
Valeur print : Avant la normalisation : Min = 0, Max =1.0
'''

#Instanciation du modèle
model =keras.models.Sequential()

#--------Les différentes couches qui vont être utilisés dans le modèle

#La tailles des inputs
model.add(keras.layers.Input((28, 28, 1)))

#Définitions des matrices de convolution (2D) et des méthodes de Pooling
model.add(keras.layers.Conv2D(8,(3,3),activation=LeakyReLU(alpha=0.1))) # 8 matrice 3x3 qui s'activent via relu
model.add(keras.layers.MaxPooling2D((2,2))) # Pooling de 2x2
model.add(keras.layers.Conv2D(16, (3, 3), activation=LeakyReLU(alpha=0.1)))# 16 matrice 3x3 qui s'activent via relu
model.add(keras.layers.MaxPooling2D((2, 2)))# Pooling de 2x2

#Flattening des imagettes avec les traitements fullyconnected
model.add(keras.layers.Flatten())
model.add(keras.layers.Dense(100, activation='relu'))
model.add(keras.layers.Dense(10, activation='softmax'))

#Print du resumés du modèle
model.summary()

#Compilation du modèle
model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])

#Taille des batchs et nombres d'époques
batch_size = 256
epochs = 10

#Traitement du modèle
history = model.fit(x_train, y_train, batch_size=batch_size,
                    epochs = epochs,
                    verbose = True,
                    validation_data=(x_test, y_test))

score = model.evaluate(x_test, y_test, verbose=0)
print(f'Test loss   : {score[0]:4.4f}')
print(f'Test accuracy   : {score[1]:4.4f}')

model.metrics_names
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')


class estimator:
    _estimator_type = ''
    classes = []
    def __init__(self, model, classes):
        self.model = model
        self._estimator_type = 'classifier'
        self.classes = classes
    def predict(self, X):
        y_prob = self.model.predict(X)
        y_pred = y_prob.argmax(axis=1)
        return y_pred
    
class_names = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
print("Saving the model")
model.save('MNIST.h5')
classifier = estimator(model, class_names)
ConfusionMatrixDisplay.from_estimator(classifier, x_test, y_test)

plt.show()



'''
Pour 10 epochs
Pour Batch_size 256
Test loss   : 0.0394
Test accuracy   : 0.9877
'''