import random
import numpy as np
import matplotlib.pyplot as plt

# Fonction pour créer la population initiale
# La population initiale ici est le resultat de l'alignement multiple appliqué aux sequences sars-cov2
def initialize_population(liste_alignments, population_size):
    population = []
    for _ in range(population_size):
        # Sélectionner aléatoirement un alignement de séquences (individu) de la population 
        individual = random.choice(liste_alignments)
        population.append(individual)
    return population

# Fonction pour calculer les scores de fitness pour plusieurs alignements
def calculate_fitness(alignments, substitution_matrix, gap_open_penalty, gap_extension_penalty):
    fitness_scores = []

    for i, alignment in enumerate(alignments):
        seq1_aligned = alignment[3]
        seq2_aligned = alignment[4]
        # Calcul du score d'alignement
        new_score = calculate_alignment_score(
            seq1_aligned, seq2_aligned, blosum62, gap_open_penalty, gap_extension_penalty
        )
        fitness_scores.append((alignment[0], alignment[1], new_score))
    
    # Trier les alignements par score décroissant
    fitness_scores.sort(key=lambda x: x[2], reverse=True)
    return fitness_scores

# Fonction pour effectuer la sélection par roulette wheel
def roulette_wheel_selection(population, fitness_scores, num_selections=2):
    """
    Effectue une sélection par roulette wheel sur les alignements
    
    population (list): La population d'individus (alignements)
    fitness_scores (list): Liste des valeurs de fitness pour chaque alignement
    num_selections (int): Le nombre d'individus à sélectionner (par défaut 2)
    
    Retourne: une liste contenant les individus sélectionnés
    """
    # Extraire uniquement les scores de fitness
    fitness_values = [score[2] for score in fitness_scores]  
    total_fitness = sum(fitness_values)

    # Calculer les probabilités cumulatives
    probabilities = [f / total_fitness for f in fitness_values]
    cumulative_probabilities = [sum(probabilities[:i+1]) for i in range(len(probabilities))]

    selected_individuals = []  # Liste pour stocker les individus sélectionnés
    
    for _ in range(num_selections):  # Effectuer le tirage pour le nombre spécifié
        rand = random.random()  # Tirer un nombre aléatoire entre 0 et 1
        
        # Trouver l'individu correspondant
        for i, cumulative_probability in enumerate(cumulative_probabilities):
            if rand <= cumulative_probability:
                selected_individuals.append(population[i])
                # Retirer la probabilité pour éviter une sélection multiple
                probabilities[i] = 0  # Exclure cet individu
                cumulative_probabilities = [sum(probabilities[:j+1]) for j in range(len(probabilities))]
                break

    return selected_individuals  # Retourne les alignements sélectionnés


# Fonction de crossover uniforme
def uniform_crossover(parent1, parent2):
    """
    Effectue un crossover uniforme entre deux alignements.
    
    Args:
        parent1 (tuple): Premier alignement (seq1, seq2, ...).
        parent2 (tuple): Deuxième alignement (seq1, seq2, ...).

    Returns:
        tuple: Deux nouveaux alignements après crossover.
    """
    seq1_1, seq2_1 = parent1[3], parent1[4]
    seq1_2, seq2_2 = parent2[3], parent2[4]

    length = len(seq1_1)

    # Initialiser les séquences enfants
    child1_seq1, child1_seq2 = [], []
    child2_seq1, child2_seq2 = [], []

    # Parcourir chaque position des séquences
    for i in range(length):
        if random.random() < 0.5:  # Probabilité de 50% pour choisir de chaque parent
            child1_seq1.append(seq1_1[i])
            child1_seq2.append(seq2_1[i])
            child2_seq1.append(seq1_2[i])
            child2_seq2.append(seq2_2[i])
        else:
            child1_seq1.append(seq1_2[i])
            child1_seq2.append(seq2_2[i])
            child2_seq1.append(seq1_1[i])
            child2_seq2.append(seq2_1[i])

    # Combiner les nouvelles séquences pour former deux enfants
    child1 = (parent1[0], parent1[1], None, ''.join(child1_seq1), ''.join(child1_seq2))
    child2 = (parent2[0], parent2[1], None, ''.join(child2_seq1), ''.join(child2_seq2))

    return child1, child2

# Les séquences sont modifiées aléatoirement selon le taux de mutation
# probabilité de mutation par position est 0.01
def mutation(sequence, mutation_rate=0.01):
    """Applique une mutation aléatoire à une séquence alignée."""
    mutated_sequence = list(sequence)  # Convertir la séquence en liste pour modification
    
    for i in range(len(mutated_sequence)):
        # si random number est inferieur a 0.01 alors il y aura une mutation sinon rien
        if random.random() < mutation_rate:
            if mutated_sequence[i] == '-':  # Si c'est un gap, le remplacer par une base aléatoire
                mutated_sequence[i] = random.choice(['A', 'T', 'C', 'G'])
            else:  # Sinon, remplacer la base par un gap ou une autre base
                mutated_sequence[i] = random.choice(['-', 'A', 'T', 'C', 'G'])
    
    return ''.join(mutated_sequence)  # Retourner la séquence mutée


# Paramètres de l'algorithme génétique
population_size = 10       # Taille de la population
num_generations = 100      # Nombre de générations
mutation_rate = 0.01       # Taux de mutation
gap_open_penalty = -2
gap_extension_penalty = -1

# Initialisation de la population
population = initialize_population(liste_alignments, population_size)

# Liste pour stocker l'évolution de la fitness
fitness_progress = []

for generation in range(num_generations):
    print(f"\n=== Génération {generation + 1} ===")
    
    # Étape 1 : Calcul des fitness
    fitness_scores = calculate_fitness(population, blosum62, gap_open_penalty, gap_extension_penalty)
    fitness_scores.sort(key=lambda x: x[2], reverse=True)  # Trier par fitness décroissante
    print(f"Meilleur fitness de la génération : {fitness_scores[0][2]:.2f}")

    # Stocker le meilleur fitness de la génération
    fitness_progress.append(fitness_scores[0][2])

    # Étape 2 : Sélection des parents
    selected_parents = []
    while len(selected_parents) < population_size // 2:  # Nombre de paires
        parents = roulette_wheel_selection(population, fitness_scores, num_selections=2)
        if len(parents) == 2:  # Vérification pour s'assurer que la sélection a produit deux parents
            selected_parents.append(parents)

    # Étape 3 : Croisement pour générer des enfants
    children = []
    for parents in selected_parents:
        parent1, parent2 = parents
        child1, child2 = uniform_crossover(parent1, parent2)
        children.extend([child1, child2])

    # Étape 4 : Mutation des enfants
    mutated_children = []
    for child in children:
        child_seq1_mutated = mutation(child[3], mutation_rate)
        child_seq2_mutated = mutation(child[4], mutation_rate)
        mutated_child = (child[0], child[1], None, child_seq1_mutated, child_seq2_mutated)
        mutated_children.append(mutated_child)

    # Étape 5 : Nouvelle population (fusionner les enfants mutants avec la population existante)
    new_population = population + mutated_children  # Ajouter les enfants mutants à la population

    # Trier la population par fitness décroissante
    new_population = sorted(new_population, key=lambda x: calculate_alignment_score(x[3], x[4], blosum62, gap_open_penalty, gap_extension_penalty), reverse=True)

    # Si la population dépasse la taille, on sélectionne les meilleurs
    new_population = new_population[:population_size]

    # Mettre à jour la population pour la prochaine génération
    population = new_population

# Résultats finaux
print("\n=== Résultats finaux ===")
final_fitness_scores = calculate_fitness(population, blosum62, gap_open_penalty, gap_extension_penalty)
print(f"Meilleur alignement après {num_generations} générations :")
print(f"Alignement entre {final_fitness_scores[0][0]} et {final_fitness_scores[0][1]} - Fitness : {final_fitness_scores[0][2]:.2f}")


# Tracer l'évolution de la fitness
plt.plot(range(1, num_generations + 1), fitness_progress)
plt.xlabel('Générations')
plt.ylabel('Fitness')
plt.title('Évolution de la fitness au fil des générations')
plt.show()

