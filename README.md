# STAGE-M1
Use of gene prediction tools for microsporidian genomes

This github includes all the files/scripts and results generated during my M1 internship at the LMGE.
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
- Pie charts on the predictions of the 4 tools ('Pie')
![Combined_pies](https://github.com/thboutet/STAGE-M1/assets/174331140/d52dd62d-f00e-4656-b5f8-d34f76cb15ae)


- A bar charts that compares the set of genes predicted by the set of tools with those of the initial database ('Bar/database_VS_4tools_microsporidia.png')
<p align="center">
<img srchttps://github.com/thboutet/STAGE-M1/assets/174331140/8343457f-0623-4c43-a3f3-1dd5052ef98e">
</p>
 
- An other bar charts which represents true positives and false positives ('Bar/prediction_microsporidia.png')
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/a71a0236-d149-462a-9d89-75a48219493e">
</p>

- A histogram of the size of genes predicted and not predicted ('Histogram/microsporidia')
file:///home/thomas/Stage/Database/result_cuniculi/Histogram/cuniculi.png

  
- A Venn Diagramms which represents the clustered genes for all tools  ('Venn/all_tools_microsporidia.png')
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/34d0d604-6741-4f2d-8623-0aa54d0482ec" width=75% height=75%>
</p>

- The list of genes that are in the database but are never predicted ('unpredicted_genes_microsporidia')  
- The list of genes that are predicted by all tools ('predicted_genes_microsporidia')  
- Csv sheets (one per tool) that will show the different clusters and whether the database genes are predicted correctly, predicted with an error (indication of the type of error), or not predicted. The false positives of the tool will also be indicated. I concatenate these sheets to create a file "MICROSPORIDIA.ods"

 I complete these results using another script: "script_pandas.py" :
 - An other Venn Diagramm which represents clustered and correct genes for all tools ('Venn/Correct_genes_microsporidia.png')
<p align="center">
<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/70cd21e9-cda1-4b99-a682-259f12f06e96" width=75% height=75% >
</p>

 - Two bar charts of errors made in gene prediction,  100% standardized ('Bar/microsporidia100%.png') and number of errors ('Bar/microsporidia.png')

<p align="center">

<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/180affbf-5a6f-472e-99bb-5e5634f62cab" width=75% height=75% >

<img src="https://github.com/thboutet/STAGE-M1/assets/174331140/6e6b6198-bf6f-4e31-b633-bdcae674dae4" width=75% height=75% >
</p>

