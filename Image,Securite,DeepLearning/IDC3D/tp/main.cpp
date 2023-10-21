#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h> // Pour read les OFF
#include <cfloat> // Pour le DBL_MAX 


//Fonction pour calculer le barycentre d'un maillage
Eigen::Vector3d getBarycentre(const Eigen::MatrixXd& m) {
    int numRows = m.rows();
    int numCols = m.cols();
    Eigen::Vector3d barycentre(0.0, 0.0, 0.0);

    for (int row = 0; row < numRows; row++) {
        for (int col = 0; col < numCols; col++) {
            barycentre(col) += m(row, col);
        }
    }

    barycentre /= static_cast<double>(numRows);
    return barycentre;
}

//Fonction pour convertir une position carthésienne en une position sphérique
Eigen::Vector3d convertVectorToSpherical(const Eigen::Vector3d& cartesian, const Eigen::Vector3d& bary) {
    double rho = sqrt( pow((cartesian (0)-bary(0)), 2) + pow((cartesian (1)-bary(1)), 2) + pow((cartesian (2)-bary(2)), 2));
    double theta = atan2(cartesian(1) - bary(1), cartesian(0) - bary(0));
    double phi = acos((cartesian(2) - bary(2)) / rho);
    return Eigen::Vector3d(rho, theta, phi);
}

//Fonction pour convertir un maillage carthésien en un maillage sphérique
void convertToSpherical(Eigen::MatrixXd& m, const Eigen::Vector3d& bary) {
    for (int i = 0; i < m.rows(); ++i) {
        Eigen::Vector3d cartesian(m(i, 0), m(i, 1), m(i, 2));
        Eigen::Vector3d spherical = convertVectorToSpherical(cartesian, bary);
        m(i, 0) = spherical(0);
        m(i, 1) = spherical(1);
        m(i, 2) = spherical(2);
    }
}

//Fonction pour convertir une position sphérique en une position carthésienne
Eigen::Vector3d convertSphericalToCartesian(const Eigen::Vector3d& spherical, const Eigen::Vector3d& bary) {
    double x = spherical(0) * cos(spherical(1)) * sin(spherical(2)) + bary(0);
    double y = spherical(0) * sin(spherical(1)) * sin(spherical(2)) + bary(1);
    double z = spherical(0) * cos(spherical(2)) + bary(2);
    return Eigen::Vector3d(x, y, z);
}

//Fonction pour convertir un maillage sphérique en un maillage carthésien
void convertToCartesian(Eigen::MatrixXd& m, const Eigen::Vector3d& bary) {
    for (int i = 0; i < m.rows(); ++i) {
        Eigen::Vector3d spherical(m(i, 0), m(i, 1), m(i, 2));
        Eigen::Vector3d cartesian = convertSphericalToCartesian(spherical, bary);
        m(i, 0) = cartesian(0);
        m(i, 1) = cartesian(1);
        m(i, 2) = cartesian(2);
    }
}

//Fonction qui retourne la valeur radiale min d'un maillage
double getMinRadialFromMatrix(Eigen::MatrixXd m) {
  return m.col(0).minCoeff();
}

//Fonction qui retourne la valeur radiale max d'un maillage
double getMaxRadialFromMatrix(Eigen::MatrixXd m) {
  return m.col(0).maxCoeff();
}



//-----------------Traitement des Bins---------------------//

//Génération des bins
std::vector<std::vector<unsigned int>> GenerateBins(Eigen::MatrixXd m, unsigned int NbrOfBits) {
  std::vector<std::vector<unsigned int>> bins;
  bins.resize(NbrOfBits); // On veut un bin par bits
  double minRadial = getMinRadialFromMatrix(m);
  double maxRadial = getMaxRadialFromMatrix(m);
  double range = (maxRadial-minRadial)/NbrOfBits; //Largeur de la zone sphérique 
  for (int i = 0; i < NbrOfBits; ++i) {// Pour chaque vertex, on va voir si ils tombent ou pas dans la range courante
    for (int j = 0; j < m.rows(); ++j) { 
      double BorneMin = minRadial + range * i; //Borne inf de la zone
      double BorneMax = minRadial + range * (i+1);// Borne sup de la zone
      if(m(j, 0) > BorneMin && m(j, 0) < BorneMax) {
        bins[i].push_back(j); // Si dedans la zone courante i , on push le vertex 
      }
    }
  }
  return bins;
}

//Calcul de la moyenne d'un bin
double MeanBin(Eigen::MatrixXd &m, std::vector<unsigned int> bin) {
  double sum = 0.0;
  size_t k = bin.size();
  for (int i = 0; i < k; ++i){
    sum += m(bin[i], 0);
  }
  return sum / (double)k;
}

//Fonction qui modifie tout les vertex du bin ciblé en vertex^k
void binPowered(Eigen::MatrixXd &m, std::vector<unsigned int> bin, double k) {
  for (int i = 0; i < bin.size(); ++i){
    m(bin[i], 0) = pow( m(bin[i], 0), k );
  }
}


// Calcul de la valeur minimale d'un bin
double getMinRadialFromBin(Eigen::MatrixXd m, std::vector<unsigned int> bin) {
  double min = DBL_MAX;
  for (int i = 0; i < bin.size(); ++i){
    if(m(bin[i], 0) < min) 
    {
      min = m(bin[i], 0);
    }
  } return min;
}

// Calcul de la valeur maximale d'un bin
double getMaxRadialFromBin(Eigen::MatrixXd m, std::vector<unsigned int> bin) {
  double max = DBL_MIN;
  for (int i = 0; i < bin.size(); ++i){
    if(m(bin[i], 0) > max) 
    {
      max = m(bin[i], 0);
    }
  } return max;
}

//Normalisation d'un bin
void normalize(Eigen::MatrixXd &m, std::vector<unsigned int> bin) {
  double minRadial = getMinRadialFromBin(m, bin);
  double maxRadial = getMaxRadialFromBin(m, bin);
  for (int i = 0; i < bin.size(); ++i){
    m(bin[i], 0) = abs(m(bin[i], 0) - minRadial) /  (maxRadial - minRadial);
  }
}

//Dénormalisation d'un bin
void denormalize(Eigen::MatrixXd &m, std::vector<unsigned int> bin, double min, double max) {
  for (int i = 0; i < bin.size(); ++i){
    m(bin[i], 0) = min + (m(bin[i], 0) * (max - min)); // Expansion des valeurs
  }
}

//Fonction d'insertion d'un message dans un mesh
void insertion(Eigen::MatrixXd &m, std::vector<bool> message, unsigned int n, double k, double alpha) {
  Eigen::Vector3d barycentre = getBarycentre(m); 
  convertToSpherical(m, barycentre); 
  std::vector<std::vector<unsigned int>> bins = GenerateBins(m, n); // Générer des bins - Découper le mesh par rayon
  double epsilon = 0.001; // Utilisé pour incrémenter/décrémenter k pour l'ajuster
  
  for (int i = 0; i < n; ++i) { // Pour chaque bit/bin du message
    double minRadial = getMinRadialFromBin(m, bins[i]); 
    double maxRadial = getMaxRadialFromBin(m, bins[i]);
    normalize(m, bins[i]);
    double k_ = k; //k qui être ajusté par epsilon
    double mean = MeanBin(m, bins[i]);
    
    if(message[i]) { // Si le message i est un 1
      
      while(mean < (0.5 + alpha)) { // Tant que mean du bin courant < 0.5+alpha
        k_ -= epsilon; // On décrémente k 
        binPowered(m, bins[i], k_); // On reduit la puissance k pour augmenter les valeurs normalisés
        mean = MeanBin(m, bins[i]); // Augmentation des bins
      }
    } 
    else { // Si le message i est 0
      while(mean > (0.5 - alpha)) { // Tant que mean du bin courant > 0.5+alpha
        k_ += epsilon; // On incrémente k
        binPowered(m, bins[i], k_); // On augmente la puissance k pour diminuer les valeurs normalisés
        mean = MeanBin(m, bins[i]); // Diminution des bins
      }
    }
    // Fin du traitement
    denormalize(m, bins[i], minRadial, maxRadial);
  }
  convertToCartesian(m, barycentre);
}

//Fonction pour extraire un message dans un mesh tatoué
void extraction(Eigen::MatrixXd &m, std::vector<bool> &message, unsigned int n) {
  message.resize(n);
  Eigen::Vector3d barycentre = getBarycentre(m);
  convertToSpherical(m, barycentre);
  std::vector<std::vector<unsigned int>> bins = GenerateBins(m, n);
  for (int i = 0; i < n; ++i) { // Pour tout les bit/bin 
    double minRadial = getMinRadialFromBin(m, bins[i]);
    double maxRadial = getMaxRadialFromBin(m, bins[i]);
    normalize(m, bins[i]);
    double mean = MeanBin(m, bins[i]);
    mean > 0.5 ? message[i] = 1 : message[i] = 0; // Test de quel coté on est de la moyenne dans le bin courant 
    denormalize(m, bins[i], minRadial, maxRadial);
  }
  convertToCartesian(m, barycentre);
}


// Fonction utiles

//Fonction qui retourne une copie d'une matrice en input 
Eigen::MatrixXd copymatrix(Eigen::MatrixXd m) {
  Eigen::MatrixXd copymatrix(m.rows(), m.cols());
  for (int i = 0; i < m.rows(); ++i)
  {
    copymatrix(i, 0) = m(i, 0);
    copymatrix(i, 1) = m(i, 1);
    copymatrix(i, 2) = m(i, 2);
  }
  return copymatrix;
}

    // Génère le message aléatoire de taille n
  void MsgGenerator(int n, std::vector<bool> &Msg)
  {
    for (int i = 0; i < n; i++) {
        bool bit = rand() % 2 == 0;
        Msg.push_back(bit);
    }
  }

//Fonction de calcul du RMSE - Métrique qui mesure la fidélité entre deux maillages
double RMSE(Eigen::MatrixXd &m, Eigen::MatrixXd &t) {
  double eqm = 0.0;
  for(unsigned int i = 0 ; i < m.rows() ; i++){
    for(unsigned int j = 0 ; j < m.cols(); j++){
      eqm += pow(m(i,j) - t(i,j),2);
    }
  }
  eqm/=(m.rows() * m.cols());
  return sqrt(eqm);
}

int main(int argc, char *argv[])
{
 
  // On definit les matrices Vertex et Faces pour récupérer les données d'un .off
  Eigen::MatrixXd Vertices;
  Eigen::MatrixXi Faces;

  igl::readOFF("../tp/bunny.off", Vertices, Faces);

  Eigen::MatrixXd Vertices_tatooed;
  Vertices_tatooed = copymatrix(Vertices);

  unsigned int n = 90;
  std::vector<bool> msg;
  std::vector<bool> msg_extracted;
  MsgGenerator(n,msg);
  double k = 1;
  double alpha = 0.45;

  insertion(Vertices_tatooed, msg, n, k, alpha);
  std::cout << "RMSE : " << RMSE(Vertices, Vertices_tatooed) << std::endl;

  extraction(Vertices_tatooed, msg_extracted, n);

  bool samemsg = true;
  for (int i = 0; i < n; ++i)
  {
    if(msg[i] != msg_extracted[i])
      samemsg = false;
  }

  if(samemsg)

  std::cout << "Le message a bien été recupéré"<<std::endl;
  else
  std::cout << "Le message recupéré est différent du message initial" << std::endl;

  std::cout<<std::endl;
  std::cout<<std::endl;

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(Vertices_tatooed, Faces);
  viewer.data().set_face_based(true);
  viewer.launch();
}

