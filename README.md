# SARS-CoV-2  project 

This project includes pairwise sequence alignment (MSA) method, genetic algorithm, and Kimura distance to study SARS-CoV-2 sequence alignments. It explores year-based and variant-based grouping while visualizing alignment scores.
---
### Data Source
The sequence data was retrieved from the NCBI database on SARS-CoV-2 variants. More information can be found at the following link:
[NCBI SARS-CoV-2 Sequence Data](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:2697049)

---

## Collaborators:
This project was developed in collaboration with @LeilaKFZ, who contributed to the analysis and construction of phylogenetic trees for the project.

---

## Tools
- **Languages**: Python
  
- **Biological Tools**: 
  - `BLOSUM62` substitution matrix for scoring alignments.

---

## Genetic Algorithm  

The genetic algorithm is used to optimize sequence alignments over generations. 

### **Methods**
1. **Selection by Roulette Wheel**:
   - **Steps**:
       1_ Calculate the total fitness of the population.
       2_ Calculate the relative fitness for each individual by dividing their fitness by the total fitness.
       3_ Calculate the cumulative probability for each individual by summing up the relative fitnesses.
       4_ Generate a random number between 0 and 1.
       5_ Select the first individual whose cumulative probability is greater than or equal to the random number.

2. **Uniform Crossover**:
   - Combines parents' sequences at each position with a 50% chance of taking the gene from either parent.

3. **Mutation**:
   - Introduces random changes to alignments.
   - **Mutation Rate**: 1% (0.01) probability of mutation per position.


---

## Visualizations

### 1. Fitness Evolution Across Generations
**Description**: Fitness evolution through generations.
![Image](https://github.com/user-attachments/assets/d2f572b5-902e-47bb-ab36-5e4cbc1f233f)

---

### 2. Average Alignment Score by Year
**Description**: Average alignment scores grouped by year.
![Image](https://github.com/user-attachments/assets/a1432f43-ed2c-4b65-93ac-01a88f9cf4d1)

---

### 3. Alignment Score by Variant
**Description**: Alignment scores grouped by variant pairs (e.g., Delta-Omicron).
![Image](https://github.com/user-attachments/assets/7170db60-67c0-4165-a325-389e9afb2efe)

---

### 4. Kimura Distance per Alignment
**Description**: Distribution of Kimura distances for sequence alignments.
![Image](https://github.com/user-attachments/assets/f94bef21-ba9a-4c44-a5c8-2dcf2c9f89f1)

---


