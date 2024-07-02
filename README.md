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
![database_VS_4tools_cuniculi](https://github.com/thboutet/STAGE-M1/assets/174331140/8343457f-0623-4c43-a3f3-1dd5052ef98e)

 
- An other bar charts which represents true positives and false positives ('Bar/prediction_microsporidia.png')
![prediction_cuniculi](https://github.com/thboutet/STAGE-M1/assets/174331140/a71a0236-d149-462a-9d89-75a48219493e)

- A histogram of the size of genes predicted and not predicted ('Histogram/microsporidia')
file:///home/thomas/Stage/Database/result_cuniculi/Histogram/cuniculi.png

  
- A Venn Diagramms which represents the clustered genes for all tools  ('Venn/all_tools_microsporidia.png')

  ![all_tools_cuniculi](https://github.com/thboutet/STAGE-M1/assets/174331140/34d0d604-6741-4f2d-8623-0aa54d0482ec)

- The list of genes that are in the database but are never predicted ('unpredicted_genes_microsporidia')
  [Uploading unpredAL590447:c(124010-124192)
AL590451:75388-75720
AL590449:c(212265-212621)
AL590451:215965-216249
AL590447:c(215534-215809)
AL590445:118662-118880
AL590444:206132-206326
AL590444:209287-209559
AL590446:c(179288-179521)
AL590448:111841-112113
AL590450:53325-53558
AL590451:173388-173762
AL590450:c(28740-28952)
AL590444:join(173787..173861,173886..174011)
AL590447:119537-119776
AL590446:c(7268-7582)
AL590447:c(134766-134984)
AL590450:c(233031-233357)
AL590446:c(143454-143759)
AL590443:c(140713-140883)
AL590451:c(138717-138899)
AL590442:c(91138-91371)
AL391737:111169-111411
AL590448:c(151478-151762)
AL590442:c(101882-102118)
AL590450:53844-54326
AL590445:189237-189464
AL590449:252502-252813
AL590443:95770-99126
AL590446:76983-77234
AL590446:c(26497-26697)
AL590447:c(23135-23443)
AL590446:81124-81423
AL590449:50276-50524
AL590445:16940-17263
AL590443:c(7642-7953)
AL590448:c(43972-44172)
AL590442:168542-168724
AL590448:c(133188-133418)
AL590443:164915-165208
AL590449:join(217935..217937,217983..218159)
AL590443:44793-45059
AL590450:c(11372-11737)
AL590444:join(61296..61466,61490..61666)
AL391737:4117-4593
AL391737:c(11158-11529)
AL590450:c(259745-259954)
AL590448:c(59011-59166)
AL590450:88890-89123
AL590446:130862-131044
AL590451:c(50820-51134)
AL590451:186383-186580
AL590444:c(174018-174362)
AL590449:140750-140905
AL590446:86234-86614
AL590448:c(22412-22636)
AL590446:c(21283-21519)
AL590449:10834-11043
AL590450:124230-124439
AL590450:97808-98251
AL590451:114819-115058
AL590450:c(134479-134709)
AL590447:c(160835-161089)
AL590446:c(182214-182429)
AL590444:c(184432-184650)
AL590446:c(44619-44837)
AL590442:c(67947-68147)
AL590447:join(48136..48246,48271..48462)
AL590450:c(220241-220699)
AL590450:257769-258362
AL391737:71022-71237
AL590447:221421-221828
AL391737:c(44035-44310)
AL590451:56851-57048
AL590445:c(108421-108663)
AL590450:c(261448-261660)
AL590451:c(join(166840..167293,167370..167473))
AL391737:c(83716-83886)
AL590451:80121-80399
AL590444:133980-134153
AL590449:42622-43038
AL590448:204130-204513
AL590445:c(120236-120472)
AL590446:78290-78556
AL590442:86682-86936
AL590447:join(204733..204950,204977..205529)
AL590449:c(243432-244034)
AL590446:c(96059-96640)
AL590447:210847-211260
AL590450:c(182008-182256)
AL590446:c(21585-21782)
AL590444:185748-185942
AL590448:37377-37562
AL590446:c(37433-37804)
AL590449:c(216799-217104)
AL590446:151190-151405
AL590451:c(146280-146477)
AL590446:23744-23968
AL590449:217139-217414
AL590442:c(47845-48162)
AL590451:c(join(52254..52409,52442..52444))
AL590444:90615-90887
AL590448:c(203729-203983)
AL590450:16810-17115
AL590449:c(207056-207220)
AL590446:c(67537-67860)
AL590445:197415-197798
AL590450:72665-72862
AL590446:212714-213025
AL590448:join(140156..140175,140199..140430)
AL590449:233001-233291
AL590448:78332-78475
AL590451:c(237506-238132)
AL590446:c(150884-151105)
AL590448:c(219127-219414)
AL590450:c(49542-49853)
AL590442:c(join(104293..104810,104842..104845))
AL590444:19413-19613
AL590442:36512-36826
AL391737:195149-195454
AL590444:52214-52447
AL590448:68308-68517
AL590448:23526-23717
AL590445:c(27931-28113)
AL590444:c(202852-203085)
icted_genes_cuniculi…]()

- The list of genes that are predicted by all tools ('predicted_genes_microsporidia')
- Csv sheets (one per tool) that will show the different clusters and whether the database genes are predicted correctly, predicted with an error (indication of the type of error), or not predicted. The false positives of the tool will also be indicated. I concatenate these sheets to create a file "MICROSPORIDIA.ods"

 I complete these results using another script: "script_pandas.py" :
 - An other Venn Diagramm which represents clustered and correct genes for all tools ('Venn/Correct_genes_microsporidia.png')
 - Two bar charts of errors made in gene prediction,  100% standardized ('Bar/microsporidia100%.png') and number of errors ('Bar/microsporidia.png')


