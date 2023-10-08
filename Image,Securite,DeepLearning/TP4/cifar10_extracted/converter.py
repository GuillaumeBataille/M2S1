import pickle
from PIL import Image
import numpy as np

def unpickle(file):
    with open(file, 'rb') as fo:
        dict = pickle.load(fo, encoding='bytes')
    return dict


imageDict = unpickle("./data_batch_1")


batch_label = [imageDict[b'batch_label']]
labels = []
data = []
filenames = []

for l in imageDict[b'labels']:
    labels.append(l)

for d in imageDict[b'data']:
    data.append(d)

for f in imageDict[b'filenames']:
    filenames.append(d)


for i in range(len(labels)):

    if(labels[i] == 7):

        image_data = data[i].reshape(3,32,32).transpose([1,2,0])

        image_rgb = image_data

        image_pil = Image.fromarray(image_rgb)

        image_pil.save("./batch2/ship/"+str(i)+".ppm")

