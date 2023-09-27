import sys
import matplotlib.pyplot as plt
def plot_histogram_from_file(filename, label=None):
    # Lecture des données
    x = []
    y = []
    with open(filename, 'r') as f:
        for line in f:
            values = line.split()
            x.append(int(values[0]))
            y.append(int(values[1]))

    # Tracer l'histogramme
    plt.bar(x, y, width=1, align='center', label=label, alpha=0.5)


def plot_histogram_from_file_with_i(i, label=None):
    # Générer le chemin du fichier en fonction de la valeur de i
    filename = f'msg_{i}/{i}.dat'
    print(filename)
    # Lecture des données
    x = []
    y = []
    with open(filename, 'r') as f:
        for line in f:
            values = line.split()
            x.append(int(values[0]))
            y.append(int(values[1]))

    # Tracer l'histogramme
    plt.bar(x, y, width=1, align='center', label=label, alpha=0.5)

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <valeur_i>")
        return
    
    i = sys.argv[1]
    plot_histogram_from_file_with_i(i, "Image tatouée "+i)  # Exemple avec un autre fichier
    plot_histogram_from_file("original.dat", "Image originale")
    #plot_histogram_from_file(6, "Autre image tatouée")  # Exemple avec un autre fichier
    
    plt.title("Histogrammes des images tatouées")
    plt.xlabel("Valeur des pixels")
    plt.ylabel("Nombre de pixels")
    plt.legend(loc='upper right')
    plt.show()

if __name__ == '__main__':
    main()
