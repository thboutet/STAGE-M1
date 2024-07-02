#Fonction pour transformer les fichiers générés par glimmer en gff
def gff(input_file, output_file):
	gene_id = 0
	with open(input_file, 'r') as f1, open(output_file, 'w') as glimmer: 	# J'ouvre le fichier à traiter ('input_file') et en même temps je créer un fichier pour stocker les gff (output_file)
		for lig in f1:
			lig = lig.rstrip()
			if lig.startswith(">"):						# Nom du contig
                		chromosome = lig.strip()[1:]
            		else:
                		lig = lig.split()
                		gene_id += 1						#J'incrème à chaque boucle le nombre pour attribuer à chaque gènes un ID (format gff)
                		start = int(lig[1])
                		end = int(lig[2])

                	if start > end:							#GFF toujours écrit de la position la plus petite à la plus grande
                    		start, end = end, start
                    		brin = "-"						#Indication quand même du brin - si on inverse les valeurs
                	else:
                    		brin = "+"						#Si les valeurs sont déjà dans l'ordre croissant brin +

                	print(f"{chromosome}\tGlimmer\tgene\t{start}\t{end}\t.\t{brin}\t.\tID={gene_id}", file=glimmer)

# Création d'une liste qui contient les différents fichiers input et output par paires 
files = [
    	("result_E.cuniculi.predict", "glimmer.gff"), 
    	("glimmer/result_A.algerae.predict", "glimmer_algerae.gff"),
	("glimmer/result_E.bieneusi.predict", "glimmer_bieneusi.gff"), ("glimmer/result_N.ceranae.predict", "glimmer_ceranae.gff"), 
	("glimmer/result_N.parisii.predict", "glimmer_parisii.gff")]

# J'applique ma fonction à toutes mes "paires" 
for input_file, output_file in files:
	gff(input_file, output_file)
			
					
