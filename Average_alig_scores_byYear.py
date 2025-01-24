import matplotlib.pyplot as plt
import itertools

# Fonction pour extraire l'année d'un nom de séquence
def extract_year(sequence_name):
    for year in ['2020', '2021', '2022', '2023']:
        if year in sequence_name:
            return year
    return None

# Créer un dictionnaire des scores par combinaison d'années
year_combinations = itertools.combinations_with_replacement(['2020', '2021', '2022', '2023'], 2)
year_combinations = list(year_combinations)  # Liste des combinaisons possibles d'années

# Dictionnaire pour stocker les scores pour chaque combinaison
alignment_scores = {f'{year1}-{year2}': [] for year1, year2 in year_combinations}

# Filtrer les alignements et organiser les scores par combinaison d'années
for align in liste_alignments:
    score = align[2]
    if len(align) >= 6:
        year1 = extract_year(align[5])
        year2 = extract_year(align[6])
        
        if year1 and year2:
            combination = f'{year1}-{year2}'
            if combination in alignment_scores:
                alignment_scores[combination].append(score)

# Calculer la moyenne des scores pour chaque combinaison
average_scores = {}
for combination, scores in alignment_scores.items():
    if scores:  # Si des scores existent pour cette combinaison
        average_scores[combination] = np.mean(scores)
    else:  # Si aucun score, mettre None ou 0
        average_scores[combination] = None

# Vérification des moyennes
print(f"Scores moyens par combinaison d'années: {average_scores}")

# Visualisation des scores moyens par combinaison d'années
plt.figure(figsize=(12, 6))

combinations = list(average_scores.keys())
averages = [score if score is not None else 0 for score in average_scores.values()]

plt.bar(combinations, averages, color='skyblue')

# Configuration du graphique
plt.xlabel('Combinaisons d\'Années')
plt.ylabel('Score moyen d\'Alignement')
plt.title('Scores Moyens d\'Alignement par Combinaisons d\'Années')
plt.xticks(rotation=45)
plt.tight_layout()

plt.show()
