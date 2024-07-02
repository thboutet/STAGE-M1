# STAGE-M1 : Use of gene prediction tools for microsporidian genomes

This github includes all the files,scripts and results generated during my M1 internship at the LMGE.
The aim of this internship was to try to optimize the Microannot tool by comparing different tools for gene prediction.
Here we have 4: Glimmer and prodigal for prokaryotic organisms, and Augustus and Funannotate which are made for eukaryotic organisms.
I will use these tools on 5 microsporidia : _Encephalitozoon cuniculi_, _Nosema ceranae_, _Enterocytozoon bieneusi_ , _Anncaliia algerae_ and _Nematocida parisii_.

# STEP 1 : CREATE A CONDA ENVIRONMENT :

```
curl -O https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh
```  
```
bash ~/Téléchargements/Anaconda3-2024.02-1-Linux-x86_64.sh
```  

```
conda config --add channels defaults  
conda config --add channels bioconda  
conda config --add channels conda-forge


conda create --name micro  
conda activate micro
```

# STEP 2 : USED THE TOOLS (EXAMPLE WITH _E.CUNICULI_)
To use the tools correctly we need a fasta file of the complete genome of microsporidia (found on NCBI)
We also need a data training file generated via Microannot
(All these files are given in this github)

## INSTALL GLIMMER
```
conda install glimmer
```

## TRAIN GLIMMER 
```
build-icm icm_file < data_training_glimmer_sur_e_cuniculi.fa  
```   
=> Creation of the training data file for cuniculi 
## RUN GLIMMER (option = codon start : ATG, gene length > 240 nt)
```
glimmer3 -g 240 --start_codons atg genome_complet/E_cuniculi.fna icm_file glimmer/result_E.cuniculi
```

## CONVERT TO GFF 
Here I use the script "glimmer/Script_gff.py" on the output ('glimmer/result_E.cuniculi.predict') to get a gff file 


## INSTALL PRODIGAL
```
conda install prodigal
```
## TRAIN A NEW SPECIES
```
prodigal -i data_training_glimmer_sur_e_cuniculi.fa -t data_training_glimmer_sur_e_cuniculi.trn -p single 
```
## RUN PRODIGAL 
```
prodigal -i genome_complet/E_cuniculi.fna -t data_training_glimmer_sur_e_cuniculi.trn -f gff > prodigal/result_E.cuniculi
```
> [!IMPORTANT]
> Prodigal does not allow to choose the size of the genes so I will treat the results later with "script_genes.py".


## INSTALL AUGUSTUS 
```
conda install augustus  
```
> [!WARNING]
> Here you may have problems with "scipio.py", so you will also have to download this script and put it in the directory indicated by the error returned by augustus.

## TRAIN A NEW SPECIES
> [!NOTE]
> To train augustus, we must create a new species. For this we need the genome and the proteins predicted for it. For cuniculi, the data training contains proteins from the other microsporidia annotated on Microannot  _Nosema ceranae_, _Enterocytozoon bieneusi_ and _Anncaliia algerae_). It will therefore be worth creating a genome corresponding to all these proteins.
```
cat genome_complet/A_algerae.fna genome_complet/E_bieneusi.fna genome_complet/N_ceranae.fna > genome_complet/all_genome_clear_cuniculi
```
Then we train it :
```
autoAugTrain.pl --species=microsporidie_cuniculi --genome=genome_complet/all_genome_clear_cuniculi --
trainingset=data_training_prot_cuniculi  
```
> [!NOTE]
> The data_training is different, I used the "script.aa.py" in order to create a data training file with amino acid from the nucleotide base file.
## RUN AUGUSTUS 
```
augustus --species=microsporidie_cuniculi --introns=off --stopCodonExcludedFromCDS=False --predictionStart=ATG genome_complet/E_cuniculi.fna > augustus/result_E.cuniculi_augustus.gff
```

## INSTALL FUNANNOTATE 
```
docker pull nextgenusfs/funannotate

wget -O funannotate-docker https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate-docker

chmod +x funannotate-docker
```
```
./funannotate-docker setup -d db/           
```
```
nano ~/.bashrc  
export FUNANNOTATE_DB=/home/path/to/db
source ~/.bashrc
```
> [!WARNING]
> Funannotate will often produce bugs. The only way to train him I found is to reused the gff file produced by augustus. (When I gave him the protein training file he predicted only a hundred genes or don't run).

```
gtf2gff3 --cfg augustus/result_E.cuniculi_augustus.gff augustus/result_E.cuniculi_augustus.gff > augustus/cuniculi_out.gff3       
```
=>Transform the gff from augustus to gff3 so that funannotate can read it 

## RUN FUNANNOTATE
```
./funannotate-docker predict -i genome_complet/E_cuniculi.fna -o funannotate -s "E.cuniculi" --augustus_gff augustus/cuniculi_out.gff3 --max_intronlen 0 --cpus 8
```
# STEP 3: Obtain the files containing the CDS predicted by each tool for each microsporidia

For this step I use the script_gene.py which will allow me to generate 2 files: 1 containing the CDS in nucleotide with a name of this style"microsporidia_CDS_tools" and the CDS in aa "Proteomes_microsporidia_tools"

<details>
<summary> Voir script_gene.py </summary>
```python   

from collections import defaultdict


#Fonction pourtransformer mes séquences nt en aa 
def dna_to_protein(dna_seq):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    protein_seq = ""
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if codon in codon_table:
            aa = codon_table[codon]
            if aa == '_': 	 # Si codon stop j'arrête
                break
            protein_seq += aa
        else:
            protein_seq += 'X'   # Si je trouve pas de correspondance dans la table des codons (à cause des 'N' je mets un 'X'
        
    return protein_seq
    
    
#Dict pour stocker le génome de cuniculi (clé = chromosome, valeur = seq nt)

def dico_genome (fichier) :
	dico_genome= defaultdict(str)
	with open (fichier, 'r') as f1 :
		for lig in f1 :
			lig = lig.rstrip()
			if lig.startswith('>') :
				chromosome = lig[1:]
			else :

				dico_genome[chromosome] += lig.upper()    
	return dico_genome

dict_genome_cuniculi = dico_genome("genome_complet/E_cuniculi.fna")
dict_genome_algerae  = dico_genome("genome_complet/A_algerae.fna")
dict_genome_bieneusi = dico_genome("genome_complet/E_bieneusi.fna")
dict_genome_ceranae = dico_genome("genome_complet/N_ceranae.fna")
dict_genome_ceranae_clear = dico_genome("genome_complet/N_ceranae_clear.fna") #Dico pour funannotate de ceranae car les noms de contigs étaient trop long 
dict_genome_parisii = dico_genome("genome_complet/N_parisii.fna")


#Extraire infos importantes du gff de glimmer (clé = chromosome, valeur = liste de liste contenant start, stop et sens du brin) 
def glimmer (fichier) :
	with open (fichier, 'r') as f1 :
		dict_glimmer = defaultdict(list)
		for lig in f1:
			lig = lig.rstrip() 
			lig = lig.split()
			dict_glimmer[lig[0]] += [[int(lig[3]),int(lig[4]) , lig[6]]]
	return dict_glimmer
	
##ASSOCIER LES FICHIERS GLIMMER A UN DICO
dict_glimmer_cuniculi = glimmer('glimmer/glimmer.gff')
dict_glimmer_algerae = glimmer('glimmer/glimmer_algerae.gff')
dict_glimmer_bieneusi = glimmer('glimmer/glimmer_bieneusi.gff')
dict_glimmer_ceranae = glimmer('glimmer/glimmer_ceranae.gff')
dict_glimmer_parisii = glimmer('glimmer/glimmer_parisii.gff')

#Extraire infos importantes du gff de augustus
def augustus (fichier) :
	dict_augustus = defaultdict(list)
	with open(fichier, 'r') as f1 :
		for lig in f1 :
			lig = lig.rstrip() 
			if not lig.startswith("#") : 				#Ne pas garder les lignes de texte 
				lig=lig.split()
				if lig[2] == 'CDS' and int(lig[4]) - int(lig[3]) >= 240 :				#Ne garder que les gènes prédits qui sont supérieur à 240 nt
					dict_augustus[lig[0]] += [[int(lig[3]), int(lig[4]), lig[6]]]
	return dict_augustus

##ASSOCIER LES FICHIERS AUGUSTUS A UN DICO
dict_augustus = augustus('augustus/result_E.cuniculi_augustus.gff')
dict_augustus_algerae = augustus('augustus/result_A.algerae_augustus.gff')
dict_augustus_bieneusi = augustus('augustus/result_E.bieneusi_augustus.gff') 
dict_augustus_ceranae= augustus('augustus/result_N.ceranae_augustus.gff') 
dict_augustus_parisii= augustus('augustus/result_N.parisii_augustus.gff') 

#Extraire infos importantes du gff de funannotate 
def funannotate (fichier) :
	dict_funannotate = defaultdict(list)
	with open(fichier, 'r') as f1 :
		for lig in f1 :
			lig = lig.rstrip()
			if not lig.startswith("#") : 				#Ne pas garder les lignes de texte 
				lig=lig.split()
				if lig[2] == 'CDS' and int(lig[4]) - int(lig[3]) >= 240 :	#Ne garder que les gènes prédits qui sont supérieur à 240 nt
					dict_funannotate[lig[0]] += [[int(lig[3]), int(lig[4]), lig[6]]]				
	return dict_funannotate

##ASSOCIER LES FICHIERS FUNANNOTATE A UN DICO 
dict_funannotate = funannotate('funannotate/predict_results/E.cuniculi.gff3')
#dict_funannotate_algerae = funannotate('funannotate_algerae/predict_results/A.algerae.gff3')
dict_funannotate_bieneusi = funannotate('funannotate_bieneusi/predict_results/E.bieneusi.gff3')
#dict_funannotate_ceranae = funannotate('funannotate_ceranae/predict_results/N.ceranae.gff3')
dict_funannotate_parisii = funannotate('funannotate_parisii/predict_results/N.parisii.gff3')

#Extraire infos importantes du gff de prodigal 
def prodigal (fichier) :
	dict_prodigal = defaultdict(list)
	with open(fichier, 'r') as f1 :
		for lig in f1:
			lig = lig.rstrip()
			if not lig.startswith("#") : 				#Ne pas garder les lignes de texte
				lig = lig.split()
				if int(lig[4]) - int(lig[3]) >= 240 :					#Ne garder que les gènes prédits qui sont supérieur à 240 nt
					dict_prodigal[lig[0]] += [[int(lig[3]), int(lig[4]), lig[6]]]
	return (dict_prodigal)
	

##ASSOCIER LES FICHIERS FUNANNOTATE A UN DICO 
dict_prodigal 	       = prodigal('prodigal/result_E.cuniculi')
dict_prodigal_algerae  = prodigal('prodigal/result_A.algerae')
dict_prodigal_bieneusi = prodigal('prodigal/result_E.bieneusi')
dict_prodigal_ceranae  = prodigal('prodigal/result_N.ceranae')
dict_prodigal_parisii  = prodigal('prodigal/result_N.parisii')


#Pour obtenir le complément inverse d'une séquence
def reverse_complement (sequence) :
	base = {"A" : "T", "T" : "A", "G" : "C", "C" : "G", "N" : "N"}	#Définition des compléments (NE PAS OUBLIER LE "N" : "N") 
	brin_reverse_comp = "".join(base[nt] for nt in reversed(sequence))
	
	return brin_reverse_comp



#Créer un fichier qui contient les séquences gènes prédites par glimmer 
glimmer_CDS = open("glimmer/E.cuniculi_CDS_glimmer", "w")	
glimmer= open("glimmer/Proteomes_E.cuniculi_glimmer", "w")
glimmer_CDS_algerae = open("glimmer/A.algerae_CDS_glimmer", "w")
glimmer_algerae = open("glimmer/Proteomes_A.algerae_glimmer", "w")
glimmer_CDS_bieneusi = open("glimmer/E.bieneusi_CDS_glimmer", "w")
glimmer_bieneusi = open("glimmer/Proteomes_E.bieneusi_glimmer", "w")
glimmer_CDS_ceranae= open("glimmer/N.ceranae_CDS_glimmer", "w")
glimmer_ceranae = open("glimmer/Proteomes_N.ceranae_glimmer", "w")
glimmer_CDS_parisii= open("glimmer/N.parisii_CDS_glimmer", "w")
glimmer_parisii = open("glimmer/Proteomes_N.parisii_glimmer", "w")

for chromosome, gene_data in dict_glimmer_cuniculi.items() :		#Pour chaque gene prédit
	if chromosome in dict_genome_cuniculi :    			#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_cuniculi[chromosome]	#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_CDS)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_CDS)
				
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer)
				
				
				



#Créer un fichier qui contient les séquences gènes prédites par augustus
augustus_CDS = open("augustus/E.cuniculi_CDS_augustus", "w")
augustus = open("augustus/Proteomes_E.cuniculi_augustus", "w")
augustus_CDS_algerae = open("augustus/A.algerae_CDS_augustus", "w")
augustus_algerae = open("augustus/Proteomes_A.algerae_augustus", "w")
augustus_CDS_bieneusi = open("augustus/E.bieneusi_CDS_augustus", "w")
augustus_bieneusi = open("augustus/Proteomes_E.bieneusi_augustus", "w")
augustus_CDS_ceranae = open("augustus/N.ceranae_CDS_augustus", "w")
augustus_ceranae = open("augustus/Proteomes_N.ceranae_augustus", "w")
augustus_CDS_parisii = open("augustus/N.parisii_CDS_augustus", "w")
augustus_parisii = open("augustus/Proteomes_N.parisii_augustus", "w")

for chromosome, gene_data in dict_augustus.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_cuniculi :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_cuniculi[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 

				gene_seq = reverse_complement(gene_seq)

			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_CDS)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_CDS)
	
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus)		
			



#Créer un fichier qui contient les séquences gènes prédites par funannotate
funannotate_CDS = open("funannotate/E.cuniculi_CDS_funannotate", "w")
funannotate = open("funannotate/Proteomes_E.cuniculi_funannotate", "w")
##funannotate_CDS_algerae = open("funannotate_algerae/A.algerae_CDS_funannotate", "w")
##funannotate_algerae = open("funannotate_algerae/Proteomes_A.algerae_funannotate", "w")
funannotate_CDS_bieneusi = open("funannotate_bieneusi/E.bieneusi_CDS_bieneusi", "w")
funannotate_bieneusi = open("funannotate_bieneusi/Proteomes_E.bieneusi_funannotate", "w")
funannotate_CDS_ceranae = open("funannotate_ceranae/N.ceranae_CDS_bieneusi", "w")
funannotate_ceranae = open("funannotate_ceranae/Proteomes_N.ceranae_funannotate", "w")
funannotate_CDS_parisii = open("funannotate_parisii/N.parisii_CDS_bieneusi", "w")
funannotate_parisii = open("funannotate_parisii/Proteomes_N.parisii_funannotate", "w")


for chromosome, gene_data in dict_funannotate.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_cuniculi :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_cuniculi[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate_CDS)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_CDS)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate)



#Créer un fichier qui contient les séquences gènes prédites par prodigal
prodigal_CDS = open("prodigal/E.cuniculi_CDS_prodigal", "w")
prodigal = open("prodigal/Proteomes_E.cuniculi_prodigal", "w")
prodigal_CDS_algerae = open("prodigal/A.algerae_CDS_prodigal", "w")
prodigal_algerae = open("prodigal/Proteomes_A.algerae_prodigal", "w")
prodigal_CDS_bieneusi = open("prodigal/E.bieneusi_CDS_prodigal", "w")
prodigal_bieneusi = open("prodigal/Proteomes_E.bieneusi_prodigal", "w")
prodigal_CDS_ceranae = open("prodigal/N.ceranae_CDS_prodigal", "w")
prodigal_ceranae = open("prodigal/Proteomes_N.ceranae_prodigal", "w")
prodigal_CDS_parisii = open("prodigal/N.parisii_CDS_prodigal", "w")
prodigal_parisii = open("prodigal/Proteomes_N.parisii_prodigal", "w")


for chromosome, gene_data in dict_prodigal.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_cuniculi :    			#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_cuniculi[chromosome]	#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_CDS)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_CDS)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal)

##=========================================ALGERAE================================================
##GLIMMER
for chromosome, gene_data in dict_glimmer_algerae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_algerae :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_algerae[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_CDS_algerae)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_CDS_algerae)
				
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_algerae)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_algerae)
				
				
##PRODIGAL
for chromosome, gene_data in dict_prodigal_algerae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_algerae :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_algerae[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_CDS_algerae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_CDS_algerae)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)			#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_algerae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_algerae)
				
			
				
##AUGUSTUS				
for chromosome, gene_data in dict_augustus_algerae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_algerae :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_algerae[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 

				gene_seq = reverse_complement(gene_seq)

			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_CDS_algerae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_CDS_algerae)
	
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_algerae)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_algerae)	
				
"""		
##FUNANNOTATE
for chromosome, gene_data in dict_funannotate_algerae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_algerae :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_algerae[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">FUNANNOTATE:{start}-{stop}({brin})", file = funannotate_CDS_algerae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_CDS_algerae)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">FUNANNOTATE:{start}-{stop}({brin})", file = funannotate_algerae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_algerae)
				
				
"""
		
#=========================================BIENEUSI================================================
##GLIMMER
for chromosome, gene_data in dict_glimmer_bieneusi.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_bieneusi :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_bieneusi[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_CDS_bieneusi)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_CDS_bieneusi)
				
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_bieneusi)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_bieneusi)				
			
##PRODIGAL
for chromosome, gene_data in dict_prodigal_bieneusi.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_bieneusi :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_bieneusi[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_CDS_bieneusi)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_CDS_bieneusi)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)			#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_bieneusi)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_bieneusi)
				
				
##AUGUSTUS				
for chromosome, gene_data in dict_augustus_bieneusi.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_bieneusi :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_bieneusi[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 

				gene_seq = reverse_complement(gene_seq)

			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_CDS_bieneusi)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_CDS_bieneusi)
	
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_bieneusi)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_bieneusi)					
				

##FUNANNOTATE
for chromosome, gene_data in dict_funannotate_bieneusi.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_bieneusi :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_bieneusi[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate_CDS_bieneusi)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_CDS_bieneusi)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate_bieneusi)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_bieneusi)			
			

"""		#=========================================CERANAE================================================
##GLIMMER
for chromosome, gene_data in dict_glimmer_ceranae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_ceranae :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_ceranae[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_CDS_ceranae)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_CDS_ceranae)
				
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_ceranae)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_ceranae)		
				
##PRODIGAL
for chromosome, gene_data in dict_prodigal_ceranae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_ceranae :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_ceranae[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_CDS_ceranae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_CDS_ceranae)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)			#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_ceranae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_ceranae)
				
				
##AUGUSTUS				
for chromosome, gene_data in dict_augustus_ceranae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_ceranae :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_ceranae[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 

				gene_seq = reverse_complement(gene_seq)

			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_CDS_ceranae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_CDS_ceranae)
	
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_ceranae)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_ceranae)					
				
			
##FUNANNOTATE
for chromosome, gene_data in dict_funannotate_ceranae.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_ceranae_clear :    				#Utiliser les bons chromosomes 	
		
		chromosome_seq = dict_genome_ceranae_clear[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate_CDS_ceranae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_CDS_ceranae)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate_ceranae)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_ceranae)
"""			
#========================================PARISII=========================================
				
##GLIMMER
for chromosome, gene_data in dict_glimmer_parisii.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_parisii :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_parisii[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_CDS_parisii)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_CDS_parisii)
				
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = glimmer_parisii)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = glimmer_parisii)	
				
				
				
##PRODIGAL
for chromosome, gene_data in dict_prodigal_parisii.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_parisii :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_parisii[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :			#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :		#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :					#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_CDS_parisii)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_CDS_parisii)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)			#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = prodigal_parisii)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = prodigal_parisii)			
				
##AUGUSTUS				
for chromosome, gene_data in dict_augustus_parisii.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_parisii :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_parisii[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides 
				gene_seq = chromosome_seq[start - 1:stop]
			
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 

				gene_seq = reverse_complement(gene_seq)

			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_CDS_parisii)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_CDS_parisii)
	
			#Print en aa
			gene_seq = dna_to_protein(gene_seq)		#Transformation en aa
				
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = augustus_parisii)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = augustus_parisii)
					
##FUNANNOTATE
for chromosome, gene_data in dict_funannotate_parisii.items() :			#Pour chaque gene prédit
	if chromosome in dict_genome_parisii :    				#Utiliser les bons chromosomes 
		chromosome_seq = dict_genome_parisii[chromosome]		#Je récupère la sequence de chaque chromosome
		
		for i, gene_info in enumerate(gene_data, 1) :		#Parcourir les données des gènes
			start, stop, brin = gene_info
			
			if start >= 0 and stop <= len(chromosome_seq) :	#Vérifier que les coordonnées sont valides
				gene_seq = chromosome_seq[start - 1:stop]
					
			if brin == "-" :				#Si mon brin est le "-" alors j'applique la fonction reverse_complement définit plutôt 
				gene_seq = reverse_complement(gene_seq)
				
			#Print en nt
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate_CDS_parisii)		
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_CDS_parisii)
	
				
			
			#Print en aa	
			gene_seq = dna_to_protein(gene_seq)		
			
			
			print(f">{chromosome[:-2]}:{start}-{stop}({brin})", file = funannotate_parisii)
			for i in range(0, len(gene_seq), 60):		
				print(gene_seq[i:i+60], file = funannotate_parisii)	

 	 ```
</details>

# STEP 4 : CLUSTER (EXAMPLE FOR _E.CUNICULI_)
In this step I will make two clusters to see if the genes predicted by my tools correspond to the genes of the initial database

## CD-HIT
```
conda install cd-hit
```
## GLIMMER 
Cluster 1 : database VS predict genes  
```
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 glimmer/Proteomes_E.cuniculi_glimmer -d 0 -o glimmer/cluster_prot100_glimmer -c 1 -A 1 
```
Cluster 2 : unclusterized genes from cluster 1 VS database  
```
cd-hit-2d -i glimmer/cluster_prot100_glimmer -i2 Proteomes_E.cuniculi.txt -d 0 -o glimmer/cluster_supprot100_glimmer -c 0.9 
```
> [!NOTE]
> Doing this two clusters allows to really find all the genes that need to be clustered 

## PRODIGAL
Cluster 1 : database VS predict genes  
```
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 prodigal/Proteomes_E.cuniculi_prodigal -d 0 -o prodigal/cluster_prot100_prodigal -c 1 -A 1 
```
Cluster 2 : unclusterized genes from cluster 1 VS database  
```
cd-hit-2d -i prodigal/cluster_prot100_prodigal -i2 Proteomes_E.cuniculi.txt -d 0 -o glimmer/cluster_supprot100_prodigal -c 0.9 
```
## AUGUSTUS
Cluster 1 : database VS predict genes  
```
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 augustus/Proteomes_E.cuniculi_augustus -d 0 -o augustus/cluster_prot100_augustus -c 1 -A 1 
```
Cluster 2 : unclusterized genes from cluster 1 VS database  
```
cd-hit-2d -i augustus/cluster_prot100_augustus -i2 Proteomes_E.cuniculi.txt -d 0 -o augustus/cluster_supprot100_augustus -c 0.9 
```
## FUNANNOTATE
Cluster 1 : database VS predict genes  
```
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 funannotate/Proteomes_E.cuniculi_funannotate -d 0 -o funannotate/cluster_prot100_funannotate -c 1 -A 1 
```
Cluster 2 : unclusterized genes from cluster 1 VS database  
```
cd-hit-2d -i funannotate/cluster_prot100_funannotate -i2 Proteomes_E.cuniculi.txt -d 0 -o funannotate/cluster_supprot100_funannotate -c 0.9 
```
# STEP 5 : OUTPUT THE RESULTS

To process my results I use the "script_cluster.py"
This will create "result_microsporidia" files which contain :  
- Pie charts on the predictions of the 4 tools ('Pie') :
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/d52dd62d-f00e-4656-b5f8-d34f76cb15ae" width=125% height=125%>
</p>

- A bar charts that compares the set of genes predicted by the set of tools with those of the initial database ('Bar/database_VS_4tools_microsporidia.png') :
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/8343457f-0623-4c43-a3f3-1dd5052ef98e">
</p>
 
- An other bar charts which represents true positives and false positives ('Bar/prediction_microsporidia.png') :
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/a71a0236-d149-462a-9d89-75a48219493e">
</p>

- A histogram of the size of genes predicted and not predicted ('Histogram/microsporidia') :
<p align="center">
<img src= "https://github.com/thboutet/STAGE-M1/assets/174331140/787b71d7-5c25-4432-b762-58469e9bd2b5">
</p>

  
- A Venn Diagramms which represents the clustered genes for all tools  ('Venn/all_tools_microsporidia.png') :
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/34d0d604-6741-4f2d-8623-0aa54d0482ec" width=75% height=75%>
</p>

- The list of genes that are in the database but are never predicted ('unpredicted_genes_microsporidia')  
- The list of genes that are predicted by all tools ('predicted_genes_microsporidia')  
- Csv sheets (one per tool) that will show the different clusters and whether the database genes are predicted correctly, predicted with an error (indication of the type of error), or not predicted. The false positives of the tool will also be indicated. I concatenate these sheets to create a file "MICROSPORIDIA.ods"

#### I complete these results using another script: "script_pandas.py" :
 - An other Venn Diagramm which represents clustered and correct genes for all tools ('Venn/Correct_genes_microsporidia.png') :
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/70cd21e9-cda1-4b99-a682-259f12f06e96" width=75% height=75% >
</p>

 - Two bar charts of errors made in gene prediction,  100% standardized ('Bar/microsporidia100%.png') and number of errors ('Bar/microsporidia.png') :

<p align="center">

<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/180affbf-5a6f-472e-99bb-5e5634f62cab" width=75% height=75% >

<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/6e6b6198-bf6f-4e31-b633-bdcae674dae4" width=75% height=75% >
</p>

