import math
from math import pi

def kimura_distance(seq1, seq2):
    #assert len(seq1) == len(seq2), "Les séquences doivent être de la même longueur."
    
    # Initialiser les compteurs de transitions et transversions
    transitions = 0
    transversions = 0

    if len(seq1) == len(seq2):
        # Comparer chaque base dans les séquences
        for base1, base2 in zip(seq1, seq2):
            if base1 != base2:
                if (base1 == 'A' and base2 == 'G') or (base1 == 'G' and base2 == 'A') or \
                   (base1 == 'C' and base2 == 'T') or (base1 == 'T' and base2 == 'C'):
                    transitions += 1  # Transition
                else:
                    transversions += 1  # Transversion

    else :
        return float('nan')

    # Calcul des proportions de transitions et transversions
    p = transitions / len(seq1)
    q = transversions / len(seq1)
    
    # Calcul de la distance de Kimura avec protection contre le domaine des logarithmes
    try:
        if p < 0.5 and q < 0.5:  # Vérifier si la distance est calculable
            K = -0.5 * math.log(1 - 2 * p - q) - 0.25 * math.log(1 - 2 * q)
            return K
        else:
            return float('nan')
    except ValueError:  # Si une erreur mathématique se produit (par exemple, log d'un nombre <= 0)
        return float('nan') 

def ktable(alignments):
    # Extraire les noms des séquences à partir des alignments
    seq_names = set()
    for alignment in alignments:
        seq1, seq2 = alignment[5], alignment[6]
        seq_names.add(seq1)
        seq_names.add(seq2)
    
    # Convertir en liste ordonnée pour avoir un ordre déterminé
    seq_names = sorted(seq_names)

    # Créer un DataFrame vide avec les noms des séquences en lignes et en colonnes
    kimura_distance_score = pd.DataFrame(index=seq_names, columns=seq_names)
    

    
    # Remplir le DataFrame avec les scores d'alignement
    for alignment in alignments:
        seq1, seq2, real_name1, real_name2 = alignment[3], alignment[4], alignment[5], alignment[6]

        
        kdistance = kimura_distance(seq1, seq2)               
        # Remplir la matrice kimura_distance_score
        kimura_distance_score.loc[real_name1, real_name2] = kdistance
        #kimura_distance_score.loc[seq2, seq1] = kdistance  # Symétrique (si aligné, le score est le même)

    # Afficher le tableau
    
    return kimura_distance_score



# Créer et afficher le tableau de scores d'alignement
# kimura_distance_score = ktable(liste_alignments)

# k_file = "kimura_distance_score.csv"

# # Enregistrer le DataFrame dans un fichier CSV
# kimura_distance_score.to_csv(k_file, index=True)

# print(" Matrice des distances Kimura \n\n",kimura_distance_score)
