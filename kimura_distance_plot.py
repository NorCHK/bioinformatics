import pandas as pd
import matplotlib.pyplot as plt

# Charger le fichier CSV
k_file = "kimura_distance_score.csv"
kimura_distance_score = pd.read_csv(k_file, index_col=0)

# Créer une liste pour stocker les données à tracer
x_labels = []  # Noms des séquences
y_values = []  # Distances (valeurs Kimura)

# Parcourir la matrice et collecter les données pour l'histogramme
for row in kimura_distance_score.index:
    for col in kimura_distance_score.columns:
        if pd.notna(kimura_distance_score.loc[row, col]):  # Vérifier les valeurs non nulles
            x_labels.append(f"{row} vs {col}")
            y_values.append(float(kimura_distance_score.loc[row, col]))

# Créer l'histogramme
plt.figure(figsize=(10, 6))
plt.bar(x_labels, y_values, color='skyblue')

# Ajouter des labels et un titre
plt.xlabel("Paires de séquences")
plt.ylabel("Distances Kimura")
plt.title("Histogramme des distances Kimura entre les séquences")
plt.xticks(rotation=90, fontsize=8)  # Rotation pour améliorer la lisibilité

# Afficher l'histogramme
plt.tight_layout()
plt.show()
