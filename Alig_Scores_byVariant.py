# Fonction pour extraire le variant d'un nom de séquence
def extract_variant(sequence_name):
    return sequence_name.split('_')[0] if '_' in sequence_name else None

# Obtenir toutes les combinaisons possibles de variants (en excluant None et combinaisons identiques)
variants = list(set([extract_variant(align[5]) for align in liste_alignments] +
                    [extract_variant(align[6]) for align in liste_alignments]))
variants = [v for v in variants if v is not None]
variant_combinations = list(itertools.combinations(variants, 2))

# Dictionnaire pour stocker les scores pour chaque combinaison
alignment_scores_variants = {f'{var1}-{var2}': [] for var1, var2 in variant_combinations}
import matplotlib.pyplot as plt
import itertools

# Organiser les alignements et leurs scores par combinaison de variants
for align in liste_alignments:
    score = align[2]
    if len(align) >= 6:
        variant1 = extract_variant(align[5])
        variant2 = extract_variant(align[6])
        if variant1 and variant2 :
            combination = f'{variant1}-{variant2}'
            reverse_combination = f'{variant2}-{variant1}'
            # Ajouter le score à la combinaison correcte
            if combination in alignment_scores_variants:
                alignment_scores_variants[combination].append(score)
            elif reverse_combination in alignment_scores_variants:
                alignment_scores_variants[reverse_combination].append(score)

# Supprimer les combinaisons avec aucun score
alignment_scores_variants = {comb: sco for comb, sco in alignment_scores_variants.items() if sco}

# Visualisation : une combinaison peut avoir plusieurs scores, on les affiche séparément
plt.figure(figsize=(12, 6))
x_labels = []
scores = []

for combination, sco in alignment_scores_variants.items():
    for s in sco:
        x_labels.append(combination)
        scores.append(s)

plt.bar(x_labels, scores, color='skyblue')

# Configuration du graphique
plt.xlabel('Combinaisons de Variants')
plt.ylabel('Score d\'Alignement')
plt.title('Scores d\'Alignement par Combinaisons de Variants')
plt.xticks(rotation=90)
plt.tight_layout()

plt.show()
