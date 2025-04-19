# MECA0011-2-Project 
Projet de groupe en langage Python effectué dans le cadre du cours d'Éléments de mécanique des fluides à l'Université de Liège.

## Description du Projet
Ce projet vise à simuler l'écoulement d'un fluide irrotationel autour d'obstacle en utilisant la méthode des différences finies pour résoudre l'équation de Laplace.

## Structure du Projet
```text
MECA0011-2-Project/
├── CL/                     # Conditions Limites
│   ├── 000-README.txt
│   ├── 1-cl.txt           # Conditions limites cas 1
│   ├── 1-dom.txt          # Domaine cas 1
│   ├── 1-num.txt          # Numérotation cas 1
│   ├── 2-contourObj.txt   # Contour de l'obstacle cas 2
│   ├── 2-dom.txt          # Domaine cas 2
│   └── 2-num.txt          # Numérotation cas 2
├── fluid_dynamics/         # Modules de calcul
│   ├── getCoeff.py        # Calcul des coefficients
│   ├── laplacian.py       # Résolution de l'équation de Laplace
│   ├── pressure.py        # Calcul de la pression
│   └── velocity.py        # Calcul des vitesses
├── tools/                  # Outils complémentaires
│   ├── circu.py           # Calcul de la circulation
│   ├── deriv.py           # Calcul des dérivées
│   └── force.py           # Calcul des forces
├── main.py                # Programme principal
├── requirements.txt       # Dépendances Python
└── README.md             # Documentation
```

## Installation
1. Clonez ce dépôt
2. Installez les dépendances :
```bash
pip install -r requirements.txt
```

## Utilisation
Pour exécuter la simulation :
```bash
python main.py
```

## Fonctionnalités
- Résolution de l'équation de Laplace en 2D
- Calcul des champs de vitesse et de pression
- Visualisation des lignes de courant
- Calcul des forces sur l'obstacle
- Analyse de la circulation

## Auteurs
- Benjamin Bock, Baptiste Desmedt & Yazan Saloum

## Licence
M.I.T