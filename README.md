# STAGE-M1
Use of gene prediction tools for microsporidian genomes

This github includes all the files/scripts and results generated during my M1 internship at the LMGE.
The aim of this internship was to try to optimize the Microannot tool by comparing different tools for gene prediction.
Here we have 4: Glimmer and prodigal for prokaryotic organisms, and Augustus and Funannotate which are made for eukaryotic organisms.
I will use these tools on 5 microsporidia : Encephalitozoon cuniculi, Nosema ceranae, Enterocytozoon bieneusi , Anncaliia algerae and Nematocida parisii.

# STEP 1 : CREATE A CONDA ENVIRONMENT :

curl -O https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh

bash ~/Téléchargements/Anaconda3-2024.02-1-Linux-x86_64.sh

conda config --add channels defaults

conda config --add channels bioconda

conda config --add channels conda-forge

conda create --name micro

conda activate micro

# STEP 2 : USED THE TOOLS (EXAMPLE WITH E.CUNICULI)
To use the tools correctly I need a fasta file of the complete genome of microsporidia (found on NCBI)
I also need a data_training file generated via microannot
(All these files are given in this github)

## INSTALL GLIMMER
conda install glimmer

## TRAIN GLIMMER 
build-icm icm_file < data_training_glimmer_sur_e_cuniculi.fa	   
=> Creation of the training data file for cuniculi 
## RUN GLIMMER (option = codon start : ATG, gene length > 240 nt
glimmer3 -g 240 --start_codons atg genome_complet/E_cuniculi.fna icm_file glimmer/result_E.cuniculi 

## CONVERT TO GFF 
Here I use the script "glimmer/Script_gff.py" to get a gff file 


## INSTALL PRODIGAL
conda install prodigal

## TRAIN A NEW SPECIES
prodigal -i data_training_glimmer_sur_e_cuniculi.fa -t data_training_glimmer_sur_e_cuniculi.trn -p single 
## RUN PRODIGAL 
prodigal -i genome_complet/E_cuniculi.fna -t data_training_glimmer_sur_e_cuniculi.trn -f gff > prodigal/result_E.cuniculi

=> Prodigal does not allow to choose the size of the genes so I will treat the results later.


## INSTALL AUGUSTUS 
conda install augustus
/!\ Here you may have problems with "scipio.py", so you will also have to download this script and put it in the directory indicated by the error returned by augustus

## TRAIN A NEW SPECIES
autoAugTrain.pl --species=microsporidie_cuniculi --genome=genome_complet/all_genome_clear_cuniculi --trainingset=data_training_prot_cuniculi  
=> /!\ The data_training is different, I used the "script.aa.py" in order to create a data training file with amino acid from the nucleotide base file.
## RUN AUGUSTUS 
augustus --species=microsporidie_cuniculi --introns=off --stopCodonExcludedFromCDS=False --predictionStart=ATG genome_complet/E_cuniculi.fna > augustus/result_E.cuniculi_augustus.gff


## INSTALL FUNANNOTATE 
docker pull nextgenusfs/funannotate
wget -O funannotate-docker https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate-docker
chmod +x funannotate-docker

./funannotate-docker setup -d db/            #Repertory for the funannotate database 

nano ~/.bashrc  
export FUNANNOTATE_DB=/home/path/to/db
source ~/.bashrc

/!\ Funannotate will often produce bugs. The only way to train him I found is to reused the gff file produced by augustus. (When I gave him the protein training file he predicted only a hundred genes)

gtf2gff3 --cfg augustus/result_E.cuniculi_augustus.gff augustus/result_E.cuniculi_augustus.gff > augustus/cuniculi_out.gff3        # Transform the gff from augustus to gff3 so that funannotate can read it 

## RUN FUNANNOTATE
./funannotate-docker predict -i genome_complet/E_cuniculi.fna -o funannotate -s "E.cuniculi" --augustus_gff augustus/cuniculi_out.gff3 --max_intronlen 0 --cpus 8

# STEP 3: Obtain the files containing the CDS predicted by each tool for each microsporidia

For this step I use the script_gene.py which will allow me to generate 2 files: 1 containing the CDS in nucleotide with a name of this style"microsporidia_CDS_tools" and the CDS in aa "Proteomes_microsporidia_tools"

# STEP 4 : CLUSTER (EXAMPLE FOR E.CUNICULI)
In this step I will make two clusters to see if the genes predicted by my tools correspond to the genes of the initial database

## CD-HIT
conda install cd-hit

## GLIMMER 
Cluster 1 : database VS predict genes  
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 glimmer/Proteomes_E.cuniculi_glimmer -d 0 -o glimmer/cluster_prot100_glimmer -c 1 -A 1 

Cluster 2 : unclusterized genes from cluster 1 VS database
cd-hit-2d -i glimmer/cluster_prot100_glimmer -i2 Proteomes_E.cuniculi.txt -d 0 -o glimmer/cluster_supprot100_glimmer -c 0.9 

Doing it is two clusters allows me to really find all the genes that need to be clustered 

## PRODIGAL
Cluster 1 : database VS predict genes  
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 prodigal/Proteomes_E.cuniculi_prodigal -d 0 -o prodigal/cluster_prot100_prodigal -c 1 -A 1 

Cluster 2 : unclusterized genes from cluster 1 VS database
cd-hit-2d -i prodigal/cluster_prot100_prodigal -i2 Proteomes_E.cuniculi.txt -d 0 -o glimmer/cluster_supprot100_prodigal -c 0.9 

## AUGUSTUS
Cluster 1 : database VS predict genes  
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 augustus/Proteomes_E.cuniculi_augustus -d 0 -o augustus/cluster_prot100_augustus -c 1 -A 1 

Cluster 2 : unclusterized genes from cluster 1 VS database
cd-hit-2d -i augustus/cluster_prot100_augustus -i2 Proteomes_E.cuniculi.txt -d 0 -o augustus/cluster_supprot100_augustus -c 0.9 

## FUNANNOTATE
Cluster 1 : database VS predict genes  
cd-hit-2d -i Proteomes_E.cuniculi.txt -i2 funannotate/Proteomes_E.cuniculi_funannotate -d 0 -o funannotate/cluster_prot100_funannotate -c 1 -A 1 

Cluster 2 : unclusterized genes from cluster 1 VS database
cd-hit-2d -i funannotate/cluster_prot100_funannotate -i2 Proteomes_E.cuniculi.txt -d 0 -o funannotate/cluster_supprot100_funannotate -c 0.9 

# STEP 5 : OUTPUT THE RESULTS

To process my results I use the "script_cluster.py"
This will create "result_microsporidia" files which contain :
-Pie charts on the predictions of the 4 tools ('Pie')
-A bar charts that compares the set of genes predicted by the set of tools with those of the initial database
('Bar/database_VS_4tools_microsporidia.png')
-An other bar charts which represents true positives and false positives ('Bar/prediction_microsporidia.png')
-A histogram of the size of genes predicted and not predicted ('Histogram/microsporidia')
-2 Venn Diagramms : 1 which represents the clustered genes for all tools ('Venn/all_tools_microsporidia.png') , the other which represents clustered and correct genes for all tools ('Correct_genes_microsporidia.png')
-The list of genes that are in the database but are never predicted ('unpredicted_genes_microsporidia')
-The list of genes that are predicted by all tools ('predicted_genes_microsporidia')


