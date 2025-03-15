### 📄 **README – Simulation de la Propagation du Wi-Fi dans un Appartement**

#### 📌 **Description du projet**

Ce projet simule la **propagation des ondes Wi-Fi** dans un appartement en utilisant **l’équation de Helmholtz** et la **Méthode des Éléments Finis (MEF)**. L’objectif est d’analyser l’effet des **murs**, des **conditions aux limites** et du **placement du routeur** sur la couverture du signal.

#### 🛠 **Technologies utilisées**

- **Python 3**
- **NumPy** – Calcul matriciel et manipulation des données
- **SciPy** – Résolution du système linéaire
- **Matplotlib** – Visualisation des résultats
- **tqdm** – Affichage de la progression du calcul

---

## 🚀 **Installation et exécution**

### **1️⃣ Prérequis**

Avant d’exécuter le code, assurez-vous d’avoir Python et les bibliothèques nécessaires installées :

```bash
pip install numpy scipy matplotlib tqdm
```

### **2️⃣ Exécution du programme**

Lancer la simulation avec :

```bash
python main.py
```

---

## 🏗 **Structure du projet**

📂 **Dossier racine**

- 📜 `main.py` – Fichier principal, exécute la simulation
- 📜 `utils.py` – Fonctions auxiliaires (maillage, conditions aux limites, interpolation...)

---

## 📊 **Explication du code**

1️⃣ **Maillage et discrétisation** (`utils.py`)

- Création d’un maillage triangulaire du domaine
- Détection des **murs** et application de la **fonction de contraste \( n(x) \)**

2️⃣ **Construction des matrices élémentaires**

- **Matrice de rigidité** \( R \)
- **Matrice de masse** \( M \)
- **Matrice de bord** \( B \) pour la condition Robin

3️⃣ **Assemblage et résolution du système \( A U = F \)**

- Résolution avec un **solveur itératif (GMRES)**

4️⃣ **Affichage des résultats** (`main.py`)

- **Carte de chaleur 2D** pour visualiser l’intensité de l’onde
- **Surface 3D** superposée au plan de l’appartement

---

## 📌 **Résultats et analyses**

- 📡 **Propagation de l’onde Wi-Fi** dans l’appartement
- 🏠 **Effet des murs** sur la réflexion et l’atténuation du signal
- 🔍 **Impact du placement du routeur** sur la couverture réseau

---

## 📌 **Améliorations possibles**

✅ Ajout de **matériaux absorbants** pour mieux modéliser les murs  
✅ Simulation en **3D** pour une meilleure précision  
✅ Optimisation du **placement du routeur** avec un algorithme d’apprentissage

---

## 📜 **Auteur**

👨‍💻 **@JudeSeruch** – 🚀
