import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
import pandas as pd

# Charger la matrice BLOSUM62
blosum62 = substitution_matrices.load("BLOSUM62")

# Classe Parametre pour definit gap 
# methode check_id verifie si identite ou substitution et retourne le score a partir de la matrice Blosum62
class Parametre:
    def __init__(self, gap):
        self.gap = gap 

    def check_id(self, sq1, sq2):
        if (sq1, sq2) in blosum62:
            return blosum62[(sq1, sq2)]
        elif (sq2, sq1) in blosum62:
            return blosum62[(sq2, sq1)]
        else:
            return self.gap

# Fonction multi_global_alignment qui effectue l'alignement multiple globale 
def multi_global_alignment(x, y, score=Parametre(0)):
    # Definire la matrice des scores dynamic_matrix
    dynamic_matrix = [[0] * (len(x) + 1) for _ in range(len(y) + 1)]

    # Initialisation de la premiere ligne et colonne de dynamic_matrix
    for i in range(len(y) + 1):
        dynamic_matrix[i][0] = score.gap * i
    for j in range(len(x) + 1):
        dynamic_matrix[0][j] = score.gap * j

    # Remplissage de la matrice
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            score_diag = score.check_id(y[i - 1], x[j - 1])
            # on utilise max car on a des sequences de grande taille
            dynamic_matrix[i][j] = max(
                dynamic_matrix[i - 1][j] + score.gap,
                dynamic_matrix[i][j - 1] + score.gap,
                dynamic_matrix[i - 1][j - 1] + score_diag
            )
            
            
    #print("Matrice des scores :")
    #for line in dynamic_matrix:
        #print(line)
    
    i = len(x)
    j= len(y)
    seq1, seq2 = "", ""
    # les derniers nucleotides de chaque sequence sont conservées tels qu'elles sont
    seq1 = ""+ x[i - 1]
    seq2 = ""+ y[i - 1]
    while i > 1 or j > 1:
        current_score = dynamic_matrix[i][j]
        vert = dynamic_matrix[i - 1][j] 
        hori = dynamic_matrix[i][j - 1] 
        diag = dynamic_matrix[i - 1][j - 1] 
            
        max_score = max(vert, hori, diag)
            
        if i > 1 and j > 1 and max_score == diag:
            i -= 1
            j -= 1
            seq1 = x[j - 1] + seq1
            seq2 = y[i - 1] + seq2

        elif i > 1 and max_score == vert:
            i -= 1
            seq1 = "-" + seq1
            seq2 = y[i - 1] + seq2
                
        elif j > 1 and max_score == hori:
            j -= 1
            seq1 = x[j - 1] + seq1
            seq2 = "-" + seq2
                              

    final_score= dynamic_matrix[len(x)][len(y)]

    return seq1, seq2, final_score


def function_msa(labels, sequences, score=Parametre(0)):
    alignments = []
    
    # Dictionnaire pour renommer les labels
    new_labels = {}

    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            seq1, seq2, score_final = multi_global_alignment(sequences[i], sequences[j], score)
            
            # Extraction des informations des labels
            variant1, date1, accession1 = labels[i].split("_")
            variant2, date2, accession2 = labels[j].split("_")
            
            # Créer des clés en concaténant variant, date et accession
            key1 = f"{variant1}_{date1}_{accession1}"
            key2 = f"{variant2}_{date2}_{accession2}"
            
            # Vérifier si la clé1 est déjà dans new_labels
            if key1 not in new_labels:
                new_labels[key1] = f"seq{len(new_labels) + 1}"

            # Vérifier si la clé2 est déjà dans new_labels
            if key2 not in new_labels:
                new_labels[key2] = f"seq{len(new_labels) + 1}"

            # Récupérer les nouveaux noms de séquences
            renamed_variant1 = new_labels[key1]
            renamed_variant2 = new_labels[key2]

            # Ajout à la liste d'alignements avec les labels renommés
            alignments.append((
                renamed_variant1, renamed_variant2, score_final, seq1, seq2,
                f"{variant1}_{accession1}_{date1}", f"{variant2}_{accession2}_{date2}"
            ))

            # Affichage des résultats avec les labels originaux et renommés
            print(f"Alignement : {renamed_variant1} vs {renamed_variant2}")
            print(f"Score de l'alignement : {score_final}")
            print(f"{renamed_variant1} : {variant1}, Date: {date1}, Accession: {accession1}")
            print(f"{renamed_variant2} : {variant2}, Date: {date2}, Accession: {accession2}")
            print()

    return alignments



# Fonction pour rendre toutes les sequences alignées de meme longueur
def funct_equal_length_seq(alignments):
    """
    Fonction pour égaliser la longueur des séquences alignées en ajoutant des gaps ('-') à la fin des séquences plus courtes.
    """
    # Trouver la longueur maximale parmi toutes les séquences alignées (seq1 et seq2)
    max_length = max(len(align[3]) for align in alignments)  # La séquence seq1 (index 3 dans le tuple)
    
    # Après avoir égalisé la longueur des séquences
    for i in range(len(alignments)):
        # Ajustement des séquences alignées à la même longueur
        alignments[i] = (
            alignments[i][0],  # renamed_variant1
            alignments[i][1],  # renamed_variant2
            alignments[i][2],  # score_final
            alignments[i][3].ljust(max_length, '-'),  # seq1 (alignée)
            alignments[i][4].ljust(max_length, '-'),  # seq2 (alignée)
            alignments[i][5],  # variant1
            alignments[i][6]   # variant2
            )
    
        # Afficher la longueur des séquences alignées
        seq1_length = len(alignments[i][3])
        seq2_length = len(alignments[i][4])
    
        print(f"Longueur de {alignments[i][5]} : {seq1_length}")
        print(f"Longueur de {alignments[i][6]} : {seq2_length}")
        print()  # Ligne vide pour séparer les résultats

    return alignments

# Fonction pour calculer le score d'un alignement avec pénalités pour gaps
def calculate_alignment_score(seq1, seq2, substitution_matrix, gap_open_penalty, gap_extension_penalty):
    score = 0
    in_gap = False  # Indique si on est dans un gap continu

    for base1, base2 in zip(seq1, seq2):
        if base1 == '-' or base2 == '-':  # Gap détecté
            if in_gap:
                score += gap_extension_penalty
            else:
                score += gap_open_penalty
                in_gap = True
        else:
            # Fin du gap
            in_gap = False
            # Ajouter le score de la matrice de substitution
            try:
                score += substitution_matrix[base1, base2]
            except KeyError:
                score += substitution_matrix[base2, base1]  # Symétrie

    return score


# Fonction pour lire le fichier FASTA
def read_fasta(filename):
    sequences = []
    labels = []
    for record in SeqIO.parse(filename, "fasta"):
        labels.append(record.id)
        sequences.append(str(record.seq))
    return labels, sequences




# # Chargement des séquences à partir du fichier "data"
# fasta_file = "sars-cov-2_data/data.fasta"
# labels, sequences = read_fasta(fasta_file)

# # Exécution de l'alignement multiple
# alignments = function_msa(labels, sequences)


# # Appliquer cette fonction après avoir obtenu les alignements
# liste_alignments = funct_equal_length_seq(alignments)

# # Définir les pénalités pour gaps
# gap_open_penalty = -2
# gap_extension_penalty = -1

# for i, alignment in enumerate(alignments):
#         seq1_aligned = alignment[3]
#         seq2_aligned = alignment[4]
#         # Calcul du score d'alignement
#         new_score = calculate_alignment_score(
#             seq1_aligned, seq2_aligned, blosum62, gap_open_penalty, gap_extension_penalty
#         )
#         liste_alignments[i] = (alignment[0], alignment[1], new_score, *alignment[3:])
    

# print(liste_alignments)


# seq_names = set()
# for alignment in liste_alignments:
#     seq1, seq2 = alignment[5], alignment[6]
#     seq_names.add(seq1)
#     seq_names.add(seq2)
# # Créer un DataFrame vide avec les noms des séquences en lignes et en colonnes
# alignments_mat = pd.DataFrame(index=seq_names, columns=seq_names)
# score_file = "score_after_alignment.csv"

# # Enregistrer le DataFrame dans un fichier CSV
# alignments_mat.to_csv(score_file, index=True)



