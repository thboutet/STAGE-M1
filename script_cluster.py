from collections import defaultdict 
import re
import sys
import matplotlib.pyplot as plt      		#Pour graphique (pie)
import matplotlib.image as mpimg
from PIL import Image 				#Pour traiter des images 
import numpy as np
from matplotlib_venn import venn2		#Pour diagramme de Venn
import matplotlib.patches as patches
from matplotlib.patches import Patch
from venny4py.venny4py import *			#Diagramme de Venn pour comparer les 4 outils 

#Fonction qui récupère les séquences de gènes prédites par les outils d'analyse qui ne corresponde pas à celles de microannot : 

def gene_predict (fichier) :
	dico = defaultdict(str)
	with open (fichier, 'r') as f1 :
		for lig in f1 :
			lig=lig.rstrip()
			if lig.startswith('>') :
				contig = lig[1:]
			else :
				dico[contig] += lig
				
	return dico 
	


print()
print("===========================POUR CUNICULI==============================")
print("Nombre total de gènes prédits par l'outil")
print("Funannotate :", len(gene_predict('funannotate/Proteomes_E.cuniculi_funannotate')))
print("Glimmer :", len(gene_predict('glimmer/Proteomes_E.cuniculi_glimmer')))
print("Augustus :", len(gene_predict('augustus/Proteomes_E.cuniculi_augustus')))
print("Prodigal :", len(gene_predict('prodigal/Proteomes_E.cuniculi_prodigal')))
print()


##Dico des bases de données clean de microannot
dico_cuniculi = gene_predict('Proteomes_E.cuniculi.txt')

##Dico des gènes prédits juste par l'outil
funannotate_cuniculi = gene_predict('funannotate/cluster_prot100_funannotate')
glimmer_cuniculi = gene_predict('glimmer/cluster_prot100_glimmer')
augustus_cuniculi = gene_predict('augustus/cluster_prot100_augustus')
prodigal_cuniculi = gene_predict('prodigal/cluster_prot100_prodigal')

#Fonction qui vérifie que les gènes qui n'ont pas été clusterisé sont bien absent :



#Fonction qui parse les fichiers de cluster de cd-hit, création de dico à partir des fichiers de clustering 
def parse_cdhit (fichier) :
	dico = defaultdict(list)
	with open (fichier, 'r') as f1 :
		for lig in f1 :
			lig = lig.rstrip()
			if lig.startswith('>') :
				cluster = lig[1:]
			else :
				dico[cluster].append(lig.split()[2][1:-3])  #Enlever le chevron et les 3 petits points à la fin 
				dico[cluster].append(lig.split()[1][0:-1])        #Longueur (aa) associé, sans prendre la virgule associé
	return dico

##Création des dicos qui contiennent tous les clusters pour cuniculi			

funannotate_cluster100    = parse_cdhit("funannotate/cluster_prot100_funannotate.clstr")
glimmer_cluster100        = parse_cdhit("glimmer/cluster_prot100_glimmer.clstr")
augustus_cluster100       = parse_cdhit("augustus/cluster_prot100_augustus.clstr")
prodigal_cluster100       = parse_cdhit("prodigal/cluster_prot100_prodigal.clstr")

#Les deuxième cluster contenant les gènes non clusterisés la première fois mais clusterisé la seconde
funannotate_supcluster100 = parse_cdhit("funannotate/cluster_supprot100_funannotate.clstr")
glimmer_supcluster100 	  = parse_cdhit("glimmer/cluster_supprot100_glimmer.clstr")
augustus_supcluster100 	  = parse_cdhit("augustus/cluster_supprot100_augustus.clstr")
prodigal_supcluster100 	  = parse_cdhit("prodigal/cluster_supprot100_prodigal.clstr")



#Le fait d'avoir deux dico de clusters par outils va être handicapant pour la suite, je créer donc une fonction pour me permettre de les fusionner en 1 seul 	
def combined_dico (clusters, supclusters) :
	liste_gene=[]
	compt=0
	for cluster, liste in clusters.items() :		#Ici je récupère les gènes de la database qui n'ont pas été clusterisés par le premier cd-hit  
		if len(liste) == 2 :				#Si la liste est égal à deux le gène est seul dans son cluster (nom + taille) 
			for supcluster, gene_predits in supclusters.items() :	#J'ouvre mon dico du second cluster
				if len(gene_predits) > 2 and liste[0] == gene_predits[2] : #Si j'ai des gènes qui ont était clusterisé (> 2), je regarde s'il s'agit des gènes qui n'avaient pas de cluster à la base. Si l'élément 0 de la liste et égal à l'élément 2 du gène prédits 
					clusters[cluster].append(gene_predits[0])	#J'ajoute dans mon dico de clusters de base le gène prédits par l'outil qui y correspond
					clusters[cluster].append(gene_predits[1])	#J'ajoute également sa taille 
	return clusters
					
#Création des dico complet (contenant l'ensemble des gènes clusterisé)
complete_funannotate_cuniculi = combined_dico(funannotate_cluster100, funannotate_supcluster100)
complete_glimmer_cuniculi = combined_dico(glimmer_cluster100, glimmer_supcluster100)
complete_augustus_cuniculi = combined_dico(augustus_cluster100,augustus_supcluster100)
complete_prodigal_cuniculi = combined_dico(prodigal_cluster100,prodigal_supcluster100)

#Compter le nombre de gène prédit à la fois par l'outil et par microannot 
def compt_gene (dico) :
	compt=0
	for cluster, liste in dico.items() :
		if len(liste) > 2 :			#Je dois avoir plus d'une séquence (1 seq + 1 taille = 2, donc je dois avoir une liste plus longue que 2
			compt += 1	
	return compt		



print("Nombre de clusters > 1 :")
print("Funannotate_100 :", compt_gene(complete_funannotate_cuniculi))
print("Glimmer_100 :", compt_gene(complete_glimmer_cuniculi))
print("Augustus_100 :", compt_gene(complete_augustus_cuniculi))
print("Prodigal_100 :", compt_gene(complete_prodigal_cuniculi))



##Pourcentage de gènes trouvés par rapport à microannot
print()
print("Gènes prédits et clusterisés (%) par rapport à microannot :") 	#Correspond au nombre de clusters > 1 prédit sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_100 :",(compt_gene(complete_funannotate_cuniculi))/len(funannotate_cluster100)*100)	
print("Glimmer_100 :",(compt_gene(complete_glimmer_cuniculi))/len(glimmer_cluster100)*100)	
print("Augustus_100 :",(compt_gene(complete_augustus_cuniculi))/len(augustus_cluster100)*100)
print("Prodigal_100 :",(compt_gene(complete_prodigal_cuniculi))/len(prodigal_cluster100)*100)
#print("All :", compt_gene(all_outil100)/len(all_outil100)*100)

#Fonction pour trouver les gènes identiques qui ont était prédit par Microannot et l'autre outil 
def same_gene (clusters) :
	liste_taille=[]
	for clusters, liste in clusters.items() :
		if len(liste) == 4 :			           #Ici je veux seulement des clusters de 2 et qui ont la même taille 
			for info in liste[1::2] :
			
				liste_taille.append(info)	   #Liste qui contient les tailles des séquences 
			
				compt = 0 
				for i in range(0,len(liste_taille)-1,2):	#Je compare mes tailles deux à deux pour voir si c'est les mêmes 
					if liste_taille[i] == liste_taille[i+1] :
						compt+=1
	return compt
	
print()
print("% gènes prédits à l'identique par l'outil ") #Correspond au nombre de clusters > 1 prédit qui ont la même longueur sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_100 :",same_gene(funannotate_cluster100)/len(funannotate_cluster100)*100)
print("Glimmer_100 :",same_gene(glimmer_cluster100)/len(glimmer_cluster100)*100)	
print("Augustus_100 :",same_gene(augustus_cluster100)/len(augustus_cluster100)*100)
print("Prodigal_100 :",same_gene(prodigal_cluster100)/len(prodigal_cluster100)*100)

	

print()
print("% gènes prédits mais pas identique ") # Nombre de gènes correctement prédits - Nombre de gènes prédits et clusterisés 
print("Funannotate_100 :",(compt_gene(complete_funannotate_cuniculi))/len(funannotate_cluster100)*100 -same_gene(funannotate_cluster100)/len(funannotate_cluster100)*100)	
print("Glimmer_100 :",(compt_gene(complete_glimmer_cuniculi))/len(glimmer_cluster100)*100 - same_gene(glimmer_cluster100)/len(glimmer_cluster100)*100)	
print("Augustus_100 :",(compt_gene(complete_augustus_cuniculi))/len(augustus_cluster100)*100 -same_gene(augustus_cluster100)/len(augustus_cluster100)*100)
print("Prodigal_100 :",(compt_gene(complete_prodigal_cuniculi))/len(prodigal_cluster100)*100 -same_gene(prodigal_cluster100)/len(prodigal_cluster100)*100)


print()
print("% gènes correct non prédits par l'outil ")# 100 - nombres de gènes prédits et clusterisés 
print("Funannotate_100 :",100- (compt_gene(complete_funannotate_cuniculi))/len(funannotate_cluster100)*100)
print("Glimmer_100 :",100- (compt_gene(complete_glimmer_cuniculi))/len(glimmer_cluster100)*100)	
print("Augustus_100 :", 100 -(compt_gene(complete_augustus_cuniculi))/len(augustus_cluster100)*100)
print("Prodigal_100 :", 100 -(compt_gene(complete_prodigal_cuniculi))/len(prodigal_cluster100)*100)
print()


##Print des camemberts de chaque outil pour visualiser leur différents niveaux de prédictions 
##pip install intel-openmp
def pie (unpredicted, misspredicted, predicted, outil, microsporidie) :
	labels = ["Gènes non prédits", "Gènes prédits avec une erreur", "Gènes prédits correctement"] 				#Je définis les noms de mes différentes parties (qui sont les mêmes peut importe l'outil de base
	sizes = [unpredicted, misspredicted, predicted]			#Mettre les valeurs correspondantes
	explode = (0.1,0,0)						#Dissocier le groupe des gènes non prédits et prédits 
	colors=['#FF1500','#fdee00', '#32CD32']
	
	
	plt.figure(figsize=(10,7)) 					#Taille de la figure
	#Print le graphique
	plt.pie(sizes, labels=labels, explode=explode, colors=colors, autopct="%1.1f%%", startangle= 90)
	
	#Titre 
	titre = str("Prédictions de "+ outil + " pour " + microsporidie)
	plt.title(titre, fontweight='bold')
	
	#Enregistrer la figure
	plt.savefig("result_"+microsporidie+"/Pie/"+outil+".png")
	plt.close()

#Faire les graphiques ici pour cuniculi
#pie(100- (compt_gene(complete_funannotate_cuniculi))/len(funannotate_cluster100)*100, (compt_gene(complete_funannotate_cuniculi))/len(funannotate_cluster100)*100 -same_gene(funannotate_cluster100)/len(funannotate_cluster100)*100, same_gene(funannotate_cluster100)/len(funannotate_cluster100)*100, "Funannotate", "cuniculi")	#Faire le graphique de funannotate_cuniculi

#pie(100- (compt_gene(complete_glimmer_cuniculi))/len(glimmer_cluster100)*100, (compt_gene(complete_glimmer_cuniculi))/len(glimmer_cluster100)*100 - same_gene(glimmer_cluster100)/len(glimmer_cluster100)*100, same_gene(glimmer_cluster100)/len(glimmer_cluster100)*100, "Glimmer", "cuniculi")

#pie(100 -(compt_gene(complete_augustus_cuniculi))/len(augustus_cluster100)*100,(compt_gene(complete_augustus_cuniculi))/len(augustus_cluster100)*100 -same_gene(augustus_cluster100)/len(augustus_cluster100)*100, same_gene(augustus_cluster100)/len(augustus_cluster100)*100, "Augustus", "cuniculi")

#pie(100 -(compt_gene(complete_prodigal_cuniculi))/len(prodigal_cluster100)*100, (compt_gene(complete_prodigal_cuniculi))/len(prodigal_cluster100)*100 -same_gene(prodigal_cluster100)/len(prodigal_cluster100)*100, same_gene(prodigal_cluster100)/len(prodigal_cluster100)*100, "Prodigal", "cuniculi")



#Concaténer les pies en un seul :
def combined_pie (microsporidie) :
	files = ['result_'+microsporidie+'/Pie/Funannotate.png','result_'+microsporidie+'/Pie/Glimmer.png', 'result_'+microsporidie+'/Pie/Augustus.png','result_'+microsporidie+'/Pie/Prodigal.png']	#Chargement des différents fichiers 
	
	fig, axs = plt.subplots(2,2, figsize=(48,30))		#Création des sous-graphs en 2x2
	
	axes = axs.flatten()					#Listes des axes
	
	for i, ax in enumerate(axes):				#Boucle pour traiter et afficher toutes les images 
		img = mpimg.imread(files[i])			#Charger l'image
		ax.imshow(img)					#La mettre dans le sous-graph
		ax.axis('off')					#Désactiver les axes
		
	fig.suptitle("  Comparaison des différents outils pour "+ microsporidie, fontsize = 40, fontweight='bold', ha='center', y= 0.95) 				#Titre 
	plt.savefig('result_'+microsporidie+'/Pie/Combined_pies.png', bbox_inches='tight')	#Enregistrer nouveau graph dans le même repertoire 
	
	plt.close()
	

#combined_pie('cuniculi')			#Pour les graphs de cuniculi

	
#Je souhaite faire un diagramme de Venn, pour cela je vais créer la listes des gènes de la base de données qui sont prédits (correctes ou non) pour l'outils, j'en profite pour récupérer les gènes non prédits
def predicted_genes(clusters):
	liste_genes_predicted=[]
	liste_genes_unpredicted=[]
	for cluster, liste in clusters.items() :		
		if len(liste) > 2 :					#Si mon cluster contient un gène prédit par l'outil associé 
			liste_genes_predicted.append(liste[0])		#Je le stock dans cette liste
		else :							#Sinon (gènes non prédits) je le stock dans cette liste 
			liste_genes_unpredicted.append(liste[0])
			
			
	return liste_genes_predicted, liste_genes_unpredicted		#Je retourne les deux

#J'associe mes listes pour cuniculi
predicted_funannotate_cuniculi, unpredicted_funannotate_cuniculi = predicted_genes(complete_funannotate_cuniculi)
predicted_glimmer_cuniculi, unpredicted_glimmer_cuniculi = predicted_genes(complete_glimmer_cuniculi)	
predicted_augustus_cuniculi, unpredicted_augustus_cuniculi = predicted_genes(complete_augustus_cuniculi)
predicted_prodigal_cuniculi, unpredicted_prodigal_cuniculi = predicted_genes(complete_prodigal_cuniculi)

#J'ai besoin pour la suite d'une liste qui contient les gènes predits par l'outil mais qui n'ont pas été clusterisé :
#Pour ça je reprends mon dico contenant mon deuxième cluster et je récupère tous les gènes non clusterisé dedans :
def unclusterized (dico_cluster2) :
	liste_genes = []
	for clusters, liste in dico_cluster2.items():
		if len(liste) == 2 :
			liste_genes.append(liste[0])
			
	return liste_genes
	
unclusterized_funannotate_cuniculi = unclusterized(funannotate_supcluster100)
#print(len(unclusterized_funannotate_cuniculi))
unclusterized_glimmer_cuniculi = unclusterized(glimmer_supcluster100)
#print(len(unclusterized_glimmer_cuniculi))
unclusterized_augustus_cuniculi = unclusterized(augustus_supcluster100)
#print(len(unclusterized_augustus_cuniculi))
unclusterized_prodigal_cuniculi = unclusterized(prodigal_supcluster100)
#print(len(unclusterized_prodigal_cuniculi))



#Je veux print mon diagramme de Venn comparatif des gènes trouvés par mes outils : (besoin de pip install venny4py)
def venn_tools (liste_aug, liste_fun, liste_gli, liste_prod, microsporidie) :
	
	#Je mets toutes mes listes dans un set et venny4py va se charger de comparer et compter les gènes commun entre les différents outils 
	sets = {'Augustus': set(liste_aug),
		'Funannotate' : set(liste_fun),
		'Glimmer': set(liste_gli),
		'Prodigal': set(liste_prod)}
	venny4py(sets=sets)
	
	#Titre
	plt.title("Venn diagram for "+microsporidie, fontweight='bold')
	
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Venn/all_tools_"+microsporidie+".png")
	plt.close()

##Print le diagramme de Venn 
#venn_tools(predicted_augustus_cuniculi, predicted_funannotate_cuniculi, predicted_glimmer_cuniculi, predicted_prodigal_cuniculi, "cuniculi")



#Extraire les gènes qui n'ont jamais été prédit 
def unpredicted_genes(liste_aug, liste_fun, liste_glim, liste_prod) :
	
	#D'abord je les convertis en set :
	set_aug = set(liste_aug)
	set_fun = set(liste_fun)
	set_glim = set(liste_glim)
	set_prod = set(liste_prod)
	
	#Maintenant je peux les comparer pour trouver seulement les éléments qui sont communs à tous ces sets (ce qui correspond aux gènes de la database jamais prédits)
	
	comparaison = set_aug.intersection(set_fun).intersection(set_glim).intersection(set_prod)
	
	return comparaison


#Je stock le nom de ces gènes dans 
unpredicted_genes_cuniculi = open("result_cuniculi/unpredicted_genes_cuniculi", 'w') 

for gene in unpredicted_genes(unpredicted_augustus_cuniculi,unpredicted_funannotate_cuniculi, unpredicted_glimmer_cuniculi, unpredicted_prodigal_cuniculi):
	print(gene, file=unpredicted_genes_cuniculi)
	
#Je stock les gènes non prédits dans cette liste 	
liste_unpredicted_cuniculi = unpredicted_genes(unpredicted_augustus_cuniculi,unpredicted_funannotate_cuniculi, unpredicted_glimmer_cuniculi, unpredicted_prodigal_cuniculi)

#Je peux également me servir de cette fonction pour garder uniquement les gènes prédits par tous les outils:

predicted_genes_cuniculi = open("result_cuniculi/predicted_genes_cuniculi", 'w') 

for gene in unpredicted_genes(predicted_augustus_cuniculi, predicted_funannotate_cuniculi, predicted_glimmer_cuniculi, predicted_prodigal_cuniculi):
	print(gene, file=predicted_genes_cuniculi)
	
	
liste_predicted_cuniculi = unpredicted_genes(predicted_augustus_cuniculi,predicted_funannotate_cuniculi, predicted_glimmer_cuniculi, predicted_prodigal_cuniculi)


#Je veux créer un graphique qui va compter la taille des gènes non prédits pour après le représenter avec un graph
def compt_length (liste) :
	compt = 0
	liste_genes = []
	for gene in liste :
		match = re.search(r'c?\(?(\d+)-(\d+)\)?', gene)	#Je cherche ce motif pour la valeur[0] 
		if match :					#Si je le trouve alors 
			start,end = map(int, match.groups())
			length = end - start
			liste_genes.append(length)
	return(liste_genes)
	
length_unpredicted_cuniculi = compt_length(liste_unpredicted_cuniculi)	#Liste des tailles de gènes non prédits
length_predicted_cuniculi = compt_length(liste_predicted_cuniculi)	#Liste des tailles de gènes prédits par les 4 outils 

#Maintenant je fais mon histogramme 
#Fonction pour le faire :
def hist (unpredicted, predicted, microsporidie) :		#Liste taille unpredicted, liste taille predicted et la microsporidie utilisé

	size_unpredicted = unpredicted
	size_predicted = predicted
	
	fig, ax = plt.subplots()
	
	#Je définis mon échelle pour les abscisses
	bins = np.arange(0, 6001,100)		#Ici je veux aller de 0 à 6000 nt, avec un pas de 100 


	ax.hist(size_unpredicted, bins=bins, alpha=0.5, label = "Unpredicted genes", edgecolor = 'black')
	ax.hist(size_predicted, bins=bins, alpha=0.5, label = "Predicted genes", edgecolor = 'black')
	
	plt.xlabel('Genes size (nt)')
	plt.ylabel('Frequency')
	plt.title('Comparison size of predicted and unpredicted genes for '+microsporidie)
	
	plt.legend(loc = 'upper right')
	
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Histogram/"+microsporidie+".png")
	plt.close()
	
#Faire l'histo
hist(length_unpredicted_cuniculi, length_predicted_cuniculi, "cuniculi")



def bar_charts (dico, cluster_aug, cluster_fun, cluster_glim, cluster_prod, predict_aug, predict_fun, predict_glim, predict_prod, microsporidie) :
	groups = ['Database','Augustus', 'Funannotate', 'Glimmer', 'Prodigal']
	
	#Je fais 2 listes de valeurs, la 1 ère contient le nombre de gènes clusterisés pour chaque outil, la seconde le nombre de gènes prédits par l'outil mais non clusterisés
	values1 =np.array([len(dico),cluster_aug, cluster_fun, cluster_glim, cluster_prod])
	values2 =np.array([0,predict_aug, predict_fun, predict_glim, predict_prod])
	
	
	
	
	fig, ax = plt.subplots()
	#Je vais faire deux graphique, 1 avec les vrais valeurs et l'autre normalisés sur 100 
	#Création des bars 
	ax.bar(groups, values1, color ='green', edgecolor = 'black', linewidth = 1)			#Bar des gènes clusterisés
	ax.bar(groups, values2, bottom = values1, color = '#E73E01', edgecolor = 'black', linewidth = 1)	#Bar des gènes prédits par l'outil mais pas clusterisés 
	
	#Print le nombre de gènes 
	for bar in ax.patches :
		ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() / 2 + bar.get_y(), round(bar.get_height()), ha = 'center', color = 'w', weight = 'bold')
		
	#Je vais print le nombre total d'erreurs trouvés pour chaque outil
	total_values = np.add(values1, values2)			#Nombre de gènes clusterisés + non clusterisés
	#Traiter chaque bars en fonction de son total
	for i, total in enumerate (total_values):
		ax.text(i, total + 0.5, round(total), ha = 'center', weight = 'bold', color = 'black')
	
	#Titre (ici set_title et pad sinon le titre est chevauchant avec la figure)
	ax.set_title("Comparison of tool predictions for "+microsporidie,  pad =20, fontweight = 'bold')
	plt.ylabel("Number of genes")
	plt.legend(["Clusterized genes", "Unclusterized genes"], loc = 'lower right')	#Nom des légendes
	
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Bar/prediction_"+microsporidie+".png")
	plt.close()
	
	
#Faire le bar_charts pour cuniculi 
bar_charts(complete_glimmer_cuniculi,(len(gene_predict('augustus/Proteomes_E.cuniculi_augustus'))-len(unclusterized_augustus_cuniculi)), (len(gene_predict('funannotate/Proteomes_E.cuniculi_funannotate'))- len(unclusterized_funannotate_cuniculi)), (len(gene_predict('glimmer/Proteomes_E.cuniculi_glimmer'))- len(unclusterized_glimmer_cuniculi)), (len(gene_predict('prodigal/Proteomes_E.cuniculi_prodigal'))- len(unclusterized_prodigal_cuniculi)), len(unclusterized_augustus_cuniculi), len(unclusterized_funannotate_cuniculi), len(unclusterized_glimmer_cuniculi), len(unclusterized_prodigal_cuniculi), "cuniculi")
	
#Pour finir je vais faire un petit graphique (bar charts) pour montrer combien de gènes sont prédits (et clusterisés) en combinant les 4 outils vs la database

def bar_charts_genes(nb_genes_totals, nb_genes_predits, microsporidie) :
	groups = ['Database','4 tools']

	
	values1 = [nb_genes_totals, nb_genes_predits]
	colors = ['#3D2B1F', '#476C98']
	
	fig, ax = plt.subplots()
	ax.bar(groups, values1, color= colors)
	for bar in ax.patches :
		ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() / 2 + bar.get_y(), round(bar.get_height()), ha = 'center', color = 'w', weight = 'bold')
	plt.ylabel("Number of genes")
	plt.title("Comparaison of all the predicted genes VS database for "+microsporidie)
	

	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Bar/database_VS_4tools_"+microsporidie+".png")
	plt.close()

#Print ce bar charts
#bar_charts_genes(len(complete_glimmer_cuniculi), (len(complete_glimmer_cuniculi)-len(liste_unpredicted_cuniculi)), "cuniculi")

"""
##Faire des diagrammes de Venn qui répresente le protéome de base, les gènes qui ont été trouvé par au moins un outil à l'intersection, les gènes jamais prédits, et les nouveaux gènes prédits

#Besoin de pip install matplotlib matplotlib-venn
def venn_diagram(A, B, AB, microsporidie) :
	A = len(A)			#unpredicted_genes
	B = len(B)			#new_genes
	AB = int(AB)			#gene trouvés par au moins un outil 
	
	
	venn = venn2(subsets=(1, 1, 1), set_labels=('Microannot database', 'All tools'))	#1,1,1 pour ne pas respecter l'échelle (plus lisible), je nomme mes deux groupes 
	
	#Mettre les valeurs dans les bons groupes
	venn.get_label_by_id('10').set_text(A)
	venn.get_label_by_id('01').set_text(B)
	venn.get_label_by_id('11').set_text(AB)
	
	
	#Titre
	plt.title("Result of all tools on "+microsporidie, fontweight='bold')
	
	#Legende
	plt.xlim(-1, 1)
	plt.ylim(-1, 1)				#Pour que la légende ne chevauche pas le diagramme
	plt.legend(['unpredicted genes', 'new predicted genes', 'predicted genes'], loc='lower right', bbox_to_anchor=(1.05, - 0.05))
	
	
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Venn/"+microsporidie+".png")
	plt.close()

venn_diagram(unpredicted_genes(all_outil100), new_genes(new_all_outil100), compt_gene(all_outil100), "cuniculi")			#POUR CUNICULI 
	
	
##Faire un treemap avec les données des gènes trouvés par au moins un des 4 outils, les gènes trouvés par tous les outils et les gènes correctement trouvés par tous les outils 


def treemap(all_predicted, predicted, same_predicted, microsporidie) :
	values =[int(all_predicted), len(predicted), len(same_predicted)] 	#Associer les valeurs correspondantes
	labels =[f"{values[0]}", f"{values[1]}", f"{values[2]}"]		#Les stocker pour les affichés dans les différents rectangles
	colors = ['#EEAB47','#F6D6A6', '#FCF2E4']				#Couleurs associées
	
	fig, ax = plt.subplots(1, figsize = (12,6))
	
	x, y = 0.1, 0.1
	width, height = 0.9, 0.9		#Données du premier rectangle (Les gènes prédits par au moins 1 outils)
	
	#Boucle pour print mes rectangles
	for i, (value, label, color) in enumerate(zip(values, labels, colors)):
		rect = patches.Rectangle((x,y), width, height, linewidth=1, edgecolor='black', facecolor=color, alpha=0.8)			#Print le rectangle avec les données mises à jours à chaque boucle (3)
		ax.add_patch(rect)
		
		plt.text(x + width * 0.02, y + height * 0.98, label, va='top', ha='left', fontsize=12, color='black', fontweight = 'bold' )			#J'affiche le label correspondant en haut à gauche du rectangle 
		
		#Je diminue la taille pour créer un rectangle plus petit 
		x += width * 0.1	
		y += height * 0.1
		width *= 0.8
		height *= 0.8
	
	
	
	#Je créer mes légends 
	legend_patches = [Patch(facecolor=colors[0], edgecolor='black', label='Genes predicted by at least one tools'), Patch(facecolor=colors[1], edgecolor='black', label='Genes predicted by the 4 tools'), Patch(facecolor=colors[2], edgecolor='black', label='Genes correctly predicted by the 4 tools')]
	
	
	
	#J'applique mes légendes 
	ax.legend(handles = legend_patches, bbox_to_anchor=(1,0.1), ncol = 1) 
	
	ax.axis('off')
	
	#Titre 
	plt.title("Onion diagram of predicted genes for "+microsporidie, fontweight = 'bold')
	
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Onion/"+microsporidie+".png")
	plt.close()
	

 
treemap(compt_gene(all_outil100), predicted_genes(all_outil100), same_predicted_genes(all_outil100), 'cuniculi')
"""


##Je créer une fonction que me permets de parcourir les dicos des clusters pour créer un fichier csv avec en colonne 1 le nom du contig, en colonne 2 le nom du contig trouvé ou "ND" si rien ne correspond et en colonne 3 si les gènes prédits sont les mêmes ou non (erreur start, stop, frameshift...)

def csv (dico) :				#1er dico = clusters, 2ème dico = gène trouvé juste par l'outil (non-clusterisé)
	for c,v in dico.items():			#Pour les clefs et les valeurs de mon dico 
		if len(v) == 4:					#Si v == 2 alors il y a un gène trouvé par l'outil  
			print(v[0],";", v[2],end= "")			#Je print le nom du conting de base puis le nom du contig de l'outil (";" car fichier csv) 
			match = re.search(r'c?\(?(\d+)-(\d+)\)?', v[0])	#Je cherche ce motif pour la valeur[0] 
			if match :					#Si je le trouve alors 
				start,end = map(int, match.groups())	#Start correspond au premier match et end au second map puis int me permettent d'être sûr de convertir chaque élément en integer, 
			
			match1 = re.search(r'c?\(?(\d+)-(\d+)\)?', v[2]) #Pareil pour le deuxième élément 
			if match1 :
				start1, end1 = map(int, match1.groups())
		
			if (start1, end1) == (start, end) :
				print("; OK")
			else :
				if start != start1 and end != end1 :
					print("; error (", start,"/", start1, ")","(", end,"/", end1,")", sep="")
				
				elif start != start1 and v[2][-2] == "+"  :			#Je vérifie à chaque fois le sens du brin pour savoir si l'erreur est au niveau de M ou de STOP
					print("; error - start (", start,"/", start1,")", sep="")
					
				elif start != start1 and v[2][-2] == "-"  :
					print("; error - end (", start,"/", start1,")", sep="")
					
				elif end != end1 and v[2][-2] == "+" :
					print("; error - end (", end,"/", end1, ")", sep="" )
					
				elif end != end1 and v[2][-2] == "-" :
					print("; error - start (", end,"/", end1, ")", sep="" )			
					
					
		if len(v) == 6 :					#Si v == 6 alors 2 gènes ont été trouvés par l'outil 
		
			print(v[0],";", v[2], end="")			#Premier gène
			
			
			match = re.search(r'c?\(?(\d+)-(\d+)\)?', v[0])	#Je cherche ce motif pour la valeur[0] 
			if match :					#Si je le trouve alors 
				start,end = map(int, match.groups())	#Start correspond au premier match et end au second map puis int me permettent d'être sûr de convertir chaque élément en integer, 
			
			match1 = re.search(r'c?\(?(\d+)-(\d+)\)?', v[2]) #Pareil pour le deuxième élément 
			if match1 :
				start1, end1 = map(int, match1.groups())
				
			if (start1, end1) == (start, end) :
				print("; OK")
			else :
				if start != start1 and end != end1 :
					print("; error (", start,"/", start1, ")","(", end,"/", end1,")", sep="")
				
				elif start != start1 and v[2][-2] == "+" :
					print("; error - start (", start,"/", start1,")", sep="")
					
				elif start != start1 and v[2][-2] == "-" :
					print("; error - end (", start,"/", start1,")", sep="")	
					
				elif end != end1 and v[2][-2] == "+" :
					print("; error - end (", end,"/", end1, ")", sep="" )
					
				elif end != end1 and v[2][-2] == "-" :
					print("; error - start (", end,"/", end1, ")", sep="" )
					
					
			print(v[0],";", v[4], end="")			#Deuxième gène 
			

			match2 = re.search(r'c?\(?(\d+)-(\d+)\)?', v[4]) 
			if match2 :
				start2, end2 = map(int, match2.groups())	
			if (start2, end2) == (start, end) :
				print("; OK")
			else :
				if start != start2 and end != end2 :
					print("; error (", start,"/", start2, ")","(", end,"/", end2,")", sep="")
				
				elif start != start2 and v[4][-2] == "+" :
					print("; error - start (", start,"/", start2,")", sep="")
				elif start != start2 and v[4][-2] == "-" :
					print("; error - end (", start,"/", start2,")", sep="")			
					
				elif end != end2 and v[4][-2] == "+" :
					print("; error - end (", end,"/", end2, ")", sep="" )	
					
				elif end != end2 and v[4][-2] == "-" :
					print("; error - start (", end,"/", end2, ")", sep="" )	
					
					
		if len(v) == 8 :					# Trois gène prédits
			print(v[0],";", v[2], end="")			#Premier gène
			
			
			match = re.search(r'c?\(?(\d+)-(\d+)\)?', v[0])	#Je cherche ce motif pour la valeur[0] 
			if match :					#Si je le trouve alors 
				start,end = map(int, match.groups())	#Start correspond au premier match et end au second map puis int me permettent d'être sûr de convertir chaque élément en integer, 
			
			match1 = re.search(r'c?\(?(\d+)-(\d+)\)?', v[2]) #Pareil pour le deuxième élément 
			if match1 :
				start1, end1 = map(int, match1.groups())
				
			if (start1, end1) == (start, end) :
				print("; OK")
			else :
				if start != start1 and end != end1 :
					print("; error (", start,"/", start1, ")","(", end,"/", end1,")", sep="")
				
				elif start != start1 and v[2][-2] == "+" :
					print("; error - start (", start,"/", start1,")", sep="")
					
				elif start != start1 and v[2][-2] == "-" :
					print("; error - end (", start,"/", start1,")", sep="")	
				
				elif end != end1 and v[2][-2] == "+" :
					print("; error - end (", end,"/", end1, ")", sep="" )
					
				elif end != end1 and v[2][-2] == "-" :
					print("; error - start (", end,"/", end1, ")", sep="" )					
					
					
			print(v[0],";", v[4], end="")			#Deuxième gène 
			

			match2 = re.search(r'c?\(?(\d+)-(\d+)\)?', v[4]) 
			if match2 :
				start2, end2 = map(int, match2.groups())	
			if (start2, end2) == (start, end) :
				print("; OK")
			else :
				if start != start2 and end != end2 :
					print("; error (", start,"/", start2, ")","(", end,"/", end2,")", sep="")
				
				elif start != start2 and v[4][-2] == "+":
					print("; error - start (", start,"/", start2,")", sep="")
				elif start != start2 and v[4][-2] == "-":
					print("; error - end (", start,"/", start2,")", sep="")					
					
					
				elif end != end2 and v[4][-2] == "+" :
					print("; error - end (", end,"/", end2, ")", sep="" )
				elif end != end2 and v[4][-2] == "-" :
					print("; error - start (", end,"/", end2, ")", sep="" )


			print(v[0],";", v[6], end="")				#Troisième gène
			
			match3 = re.search(r'c?\(?(\d+)-(\d+)\)?', v[6]) 
			if match3 :
				start3, end3 = map(int, match3.groups())	
			if (start3, end3) == (start, end) :
				print("; OK")
			else :
				if start != start3 and end != end3 :
					print("; error (", start,"/", start3, ")","(", end,"/", end3,")", sep="")
				
				elif start != start3 and v[6][-2] == "+":
					print("; error - start (", start,"/", start3,")", sep="")
					
				elif start != start3 and v[6][-2] == "-":
					print("; error - end (", start,"/", start3,")", sep="")	
						
				elif end != end3 and v[6][-2] == "+":
					print("; error - end (", end,"/", end3, ")", sep="" )
				elif end != end3 and v[6][-2] == "-":
					print("; error - start (", end,"/", end3, ")", sep="" )
			
		if len(v) == 2 :				#Si v == 2 alors le gène n'a pas été prédit par l'outil 
			print(v[0], ";", "ND", "; ND")			#Je print quand même le nom du contig mais il n'y a pas de data associé 
			
		

		
	
##Fonction pour print les fichiers csv 
def stdout (outil, file, dico, liste)	:	
	original_stdout = sys.stdout
	with open (file, 'w') as f :
		sys.stdout = f
		print("PROTEOME ;",outil,"; ERROR",outil)				#Nom de mes colonnes (l'outil est à préciser en temps que str en appelant la fonction 
		csv(dico)
		for gene in liste :					#J'ajoute à mon print les gènes qui ont été trouvé seulement par l'outil 
			print("ND ;",gene,"; ND")
		sys.stdout = original_stdout	
		
##Print les feuilles csv	
"""	
stdout('GLIMMER', 'result_cuniculi/glimmer.csv', complete_glimmer_cuniculi, unclusterized_glimmer_cuniculi)		#Glimmer_cuniculi

stdout('AUGUSTUS', 'result_cuniculi/augustus.csv', complete_augustus_cuniculi, unclusterized_augustus_cuniculi)		#Augustus_cuniculi
stdout('FUNANNOTATE','result_cuniculi/funannotate.csv', complete_funannotate_cuniculi, unclusterized_funannotate_cuniculi)	#Funannotate_cuniculi
stdout('PRODIGAL', 'result_cuniculi/prodigal.csv', complete_prodigal_cuniculi, unclusterized_prodigal_cuniculi)		#Prodigal_cuniculi
"""



print()
print("===========================POUR BIENEUSi===============================")
print("Nombre total de gènes prédits par l'outil")
print("Funannotate :", len(gene_predict('funannotate_bieneusi/Proteomes_E.bieneusi_funannotate')))
print("Glimmer :", len(gene_predict('glimmer/Proteomes_E.bieneusi_glimmer')))
print("Augustus :", len(gene_predict('augustus/Proteomes_E.bieneusi_augustus')))
print("Prodigal :", len(gene_predict('prodigal/Proteomes_E.bieneusi_prodigal')))
print()

##Dico des gènes prédits juste par l'outil
funannotate_bieneusi= gene_predict('funannotate_bieneusi/cluster_bieneusi_funannotate')
glimmer_bieneusi = gene_predict('glimmer/cluster_bieneusi_glimmer')
augustus_bieneusi = gene_predict('augustus/cluster_bieneusi_augustus')
prodigal_bieneusi = gene_predict('prodigal/cluster_bieneusi_prodigal')

##Création des dicos qui contiennent tous les clusters pour bieneusi
funannotate_bieneusi_cluster = parse_cdhit("funannotate_bieneusi/cluster_bieneusi_funannotate.clstr")
glimmer_bieneusi_cluster       = parse_cdhit("glimmer/cluster_bieneusi_glimmer.clstr")
augustus_bieneusi_cluster    = parse_cdhit("augustus/cluster_bieneusi_augustus.clstr")
prodigal_bieneusi_cluster    = parse_cdhit("prodigal/cluster_bieneusi_prodigal.clstr")


#Les deuxième cluster contenant les gènes non clusterisés la première fois mais clusterisé la seconde
funannotate_bieneusi_supcluster   = parse_cdhit("funannotate_bieneusi/cluster_supbieneusi_funannotate.clstr")
glimmer_bieneusi_supcluster 	  = parse_cdhit("glimmer/cluster_supbieneusi_glimmer.clstr")
augustus_bieneusi_supcluster 	  = parse_cdhit("augustus/cluster_supbieneusi_augustus.clstr")
prodigal_bieneusi_supcluster 	  = parse_cdhit("prodigal/cluster_supbieneusi_prodigal.clstr")

#Création des dico complet (contenant l'ensemble des gènes clusterisé)
complete_funannotate_bieneusi = combined_dico(funannotate_bieneusi_cluster, funannotate_bieneusi_supcluster)
complete_glimmer_bieneusi = combined_dico(glimmer_bieneusi_cluster, glimmer_bieneusi_supcluster)
complete_augustus_bieneusi = combined_dico(augustus_bieneusi_cluster,augustus_bieneusi_supcluster)
complete_prodigal_bieneusi = combined_dico(prodigal_bieneusi_cluster,prodigal_bieneusi_supcluster)

#Compter le nombre de clusters 
print("Nombre de clusters > 1 :")
print("Funannotate_bieneusi :", compt_gene(complete_funannotate_bieneusi))
print("Glimmer_bieneusi :", compt_gene(complete_glimmer_bieneusi))
print("Augustus_bieneusi :", compt_gene(complete_augustus_bieneusi))
print("Prodigal_bieneusi :", compt_gene(complete_prodigal_bieneusi))

##Pourcentage de gènes trouvés par rapport à microannot
print()
print("Gènes prédits et clusterisés (%) par rapport à microannot :") 	#Correspond au nombre de clusters > 1 prédit sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_bieneusi :",(compt_gene(complete_funannotate_bieneusi))/len(funannotate_bieneusi_cluster)*100)	
print("Glimmer_bieneusi :",(compt_gene(complete_glimmer_bieneusi))/len(glimmer_bieneusi_cluster)*100)	
print("Augustus_bieneusi :",(compt_gene(complete_augustus_bieneusi))/len(augustus_bieneusi_cluster)*100)
print("Prodigal_bieneusi :",(compt_gene(complete_prodigal_bieneusi))/len(prodigal_bieneusi_cluster)*100)


print()
print("% gènes prédits à l'identique par l'outil ") #Correspond au nombre de clusters > 1 prédit qui ont la même longueur sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_bieneusi :",same_gene(funannotate_bieneusi_cluster)/len(funannotate_bieneusi_cluster)*100)
print("Glimmer_bieneusi :",same_gene(glimmer_bieneusi_cluster)/len(glimmer_bieneusi_cluster)*100)	
print("Augustus_bieneusi :",same_gene(augustus_bieneusi_cluster)/len(augustus_bieneusi_cluster)*100)
print("Prodigal_bieneusi :",same_gene(prodigal_bieneusi_cluster)/len(prodigal_bieneusi_cluster)*100)



print()
print("% gènes prédits mais pas identique ") # Nombre de gènes correctement prédits - Nombre de gènes prédits et clusterisés 
print("Funannotate_bieneusi :",(compt_gene(complete_funannotate_bieneusi))/len(funannotate_bieneusi_cluster)*100 -same_gene(funannotate_bieneusi_cluster)/len(funannotate_bieneusi_cluster)*100)	
print("Glimmer_bieneusi :",(compt_gene(complete_glimmer_bieneusi))/len(glimmer_bieneusi_cluster)*100 - same_gene(glimmer_bieneusi_cluster)/len(glimmer_bieneusi_cluster)*100)	
print("Augustus_bieneusi :",(compt_gene(complete_augustus_bieneusi))/len(augustus_bieneusi_cluster)*100 -same_gene(augustus_bieneusi_cluster)/len(augustus_bieneusi_cluster)*100)
print("Prodigal_bieneusi :",(compt_gene(complete_prodigal_bieneusi))/len(prodigal_bieneusi_cluster)*100 -same_gene(prodigal_bieneusi_cluster)/len(prodigal_bieneusi_cluster)*100)


print()
print("% gènes correct non prédits par l'outil ")# 100 - nombres de gènes prédits et clusterisés 
print("Funannotate_bieneusi :",100- (compt_gene(complete_funannotate_bieneusi))/len(funannotate_bieneusi_cluster)*100)
print("Glimmer_bieneusi :",100- (compt_gene(complete_glimmer_bieneusi))/len(glimmer_bieneusi_cluster)*100)	
print("Augustus_bieneusi :", 100 -(compt_gene(complete_augustus_bieneusi))/len(augustus_bieneusi_cluster)*100)
print("Prodigal_bieneusi :", 100 -(compt_gene(complete_prodigal_bieneusi))/len(prodigal_bieneusi_cluster)*100)
print()

#Faire les pie charts pour bieneusi :

#pie(100- (compt_gene(complete_funannotate_bieneusi))/len(funannotate_bieneusi_cluster)*100, (compt_gene(complete_funannotate_bieneusi))/len(funannotate_bieneusi_cluster)*100 -same_gene(funannotate_bieneusi_cluster)/len(funannotate_bieneusi_cluster)*100, same_gene(funannotate_bieneusi_cluster)/len(funannotate_bieneusi_cluster)*100, "Funannotate", "bieneusi")	#Faire le graphique de funannotate_cuniculi

#pie(100- (compt_gene(complete_glimmer_bieneusi))/len(glimmer_bieneusi_cluster)*100, (compt_gene(complete_glimmer_bieneusi))/len(glimmer_bieneusi_cluster)*100 - same_gene(glimmer_bieneusi_cluster)/len(glimmer_bieneusi_cluster)*100, same_gene(glimmer_bieneusi_cluster)/len(glimmer_bieneusi_cluster)*100, "Glimmer", "bieneusi")

#pie(100 -(compt_gene(complete_augustus_bieneusi))/len(augustus_bieneusi_cluster)*100,(compt_gene(complete_augustus_bieneusi))/len(augustus_bieneusi_cluster)*100 -same_gene(augustus_bieneusi_cluster)/len(augustus_bieneusi_cluster)*100, same_gene(augustus_bieneusi_cluster)/len(augustus_bieneusi_cluster)*100, "Augustus", "bieneusi")

#pie(100 -(compt_gene(complete_prodigal_bieneusi))/len(prodigal_bieneusi_cluster)*100, (compt_gene(complete_prodigal_bieneusi))/len(prodigal_bieneusi_cluster)*100 -same_gene(prodigal_bieneusi_cluster)/len(prodigal_bieneusi_cluster)*100, same_gene(prodigal_bieneusi_cluster)/len(prodigal_bieneusi_cluster)*100, "Prodigal", "bieneusi")

#combined_pie('bieneusi')


#Maintenant je récupère les gènes prédits par tous les outils et ceux jamais prédits
#J'associe mes listes pour bieneusi
predicted_funannotate_bieneusi, unpredicted_funannotate_bieneusi = predicted_genes(complete_funannotate_bieneusi)
predicted_glimmer_bieneusi, unpredicted_glimmer_bieneusi = predicted_genes(complete_glimmer_bieneusi)	
predicted_augustus_bieneusi, unpredicted_augustus_bieneusi = predicted_genes(complete_augustus_bieneusi)
predicted_prodigal_bieneusi, unpredicted_prodigal_bieneusi = predicted_genes(complete_prodigal_bieneusi)


#Faire le diagramme de Venn
#venn_tools(predicted_augustus_bieneusi, predicted_funannotate_bieneusi, predicted_glimmer_bieneusi, predicted_prodigal_bieneusi, "bieneusi")


#Récupérer les gènes prédits par les différents outils qui n'ont jamais était clusterisés
unclusterized_funannotate_bieneusi = unclusterized(funannotate_bieneusi_supcluster)
#print(len(unclusterized_funannotate_bieneusi))
unclusterized_glimmer_bieneusi = unclusterized(glimmer_bieneusi_supcluster)
#print(len(unclusterized_glimmer_bieneusi))
unclusterized_augustus_bieneusi = unclusterized(augustus_bieneusi_supcluster)
#print(len(unclusterized_augustus_bieneusi))
unclusterized_prodigal_bieneusi = unclusterized(prodigal_bieneusi_supcluster)
#print(len(unclusterized_prodigal_bieneusi))


#Trouver les gènes de la base de données qui n'ont jamais été clusterisé
#Je stock le nom de ces gènes dans 
unpredicted_genes_bieneusi = open("result_bieneusi/unpredicted_genes_bieneusi", 'w') 

for gene in unpredicted_genes(unpredicted_augustus_bieneusi,unpredicted_funannotate_bieneusi, unpredicted_glimmer_bieneusi, unpredicted_prodigal_bieneusi):
	print(gene, file=unpredicted_genes_bieneusi)
	
#Je stock les gènes non prédits dans cette liste 	
liste_unpredicted_bieneusi = unpredicted_genes(unpredicted_augustus_bieneusi,unpredicted_funannotate_bieneusi, unpredicted_glimmer_bieneusi, unpredicted_prodigal_bieneusi)

#Je peux également me servir de cette fonction pour garder uniquement les gènes prédits par tous les outils:

predicted_genes_bieneusi = open("result_bieneusi/predicted_genes_bieneusi", 'w') 

for gene in unpredicted_genes(predicted_augustus_bieneusi, predicted_funannotate_bieneusi, predicted_glimmer_bieneusi, predicted_prodigal_bieneusi):
	print(gene, file=predicted_genes_bieneusi)
	
#Je stock les gènes prédits dans cette liste 	
liste_predicted_bieneusi = unpredicted_genes(predicted_augustus_bieneusi, predicted_funannotate_bieneusi, predicted_glimmer_bieneusi, predicted_prodigal_bieneusi)

	
#Pour créer mon histogramme j'ai besoin des listes des tailles	
length_unpredicted_bieneusi = compt_length(liste_unpredicted_bieneusi)	#Liste des tailles de gènes non prédits
length_predicted_bieneusi = compt_length(liste_predicted_bieneusi)	#Liste des tailles de gènes prédits par les 4 outils 

#Faire l'histo
#hist(length_unpredicted_bieneusi, length_predicted_bieneusi, "bieneusi")



#Faire le bar_charts pour bieneusi
bar_charts(complete_glimmer_bieneusi,(len(gene_predict('augustus/Proteomes_E.bieneusi_augustus'))-len(unclusterized_augustus_bieneusi)), (len(gene_predict('funannotate_bieneusi/Proteomes_E.bieneusi_funannotate'))- len(unclusterized_funannotate_bieneusi)), (len(gene_predict('glimmer/Proteomes_E.bieneusi_glimmer'))- len(unclusterized_glimmer_bieneusi)), (len(gene_predict('prodigal/Proteomes_E.bieneusi_prodigal'))- len(unclusterized_prodigal_bieneusi)), len(unclusterized_augustus_bieneusi), len(unclusterized_funannotate_bieneusi), len(unclusterized_glimmer_bieneusi), len(unclusterized_prodigal_bieneusi), "bieneusi")



##Print les feuilles csv
"""
#Faire les feuilles csv
stdout('GLIMMER', 'result_bieneusi/glimmer.csv', complete_glimmer_bieneusi, unclusterized_glimmer_bieneusi)		#Glimmer_bieneusi

stdout('AUGUSTUS', 'result_bieneusi/augustus.csv', complete_augustus_bieneusi, unclusterized_augustus_bieneusi)		#Augustus_bieneusi
stdout('FUNANNOTATE','result_bieneusi/funannotate.csv', complete_funannotate_bieneusi, unclusterized_funannotate_bieneusi)	#Funannotate_bieneusi
stdout('PRODIGAL', 'result_bieneusi/prodigal.csv', complete_prodigal_bieneusi, unclusterized_prodigal_bieneusi)		#Prodigal_bieneusi
"""

#Print le bar_charts du nombre de gènes dans la base de données et le nombre de gènes prédits et clusterisé via les 4 outils
#bar_charts_genes(len(complete_glimmer_bieneusi), (len(complete_glimmer_bieneusi)-len(liste_unpredicted_bieneusi)), "bieneusi")


#======================================================================================================
print()
print("===========================POUR CERANAE==============================")
print("Nombre total de gènes prédits par l'outil")
print("Funannotate :", len(gene_predict('funannotate_ceranae/Proteomes_N.ceranae_funannotate')))
print("Glimmer :", len(gene_predict('glimmer/Proteomes_N.ceranae_glimmer')))
print("Augustus :", len(gene_predict('augustus/Proteomes_N.ceranae_augustus')))
print("Prodigal :", len(gene_predict('prodigal/Proteomes_N.ceranae_prodigal')))
print()

##Dico des bases de données clean de microannot
dico_ceranae = gene_predict('Proteomes_N.ceranae.txt')

##Dico des gènes prédits juste par l'outil
funannotate_ceranae = gene_predict('funannotate_ceranae/cluster_ceranae_funannotate')
glimmer_ceranae = gene_predict('glimmer/cluster_ceranae_glimmer')
augustus_ceranae = gene_predict('augustus/cluster_ceranae_augustus')
prodigal_ceranae = gene_predict('prodigal/cluster_ceranae_prodigal')

##Création des dicos qui contiennent tous les clusters pour bieneusi
funannotate_ceranae_cluster = parse_cdhit("funannotate_ceranae/cluster_ceranae_funannotate.clstr")
glimmer_ceranae_cluster       = parse_cdhit("glimmer/cluster_ceranae_glimmer.clstr")
augustus_ceranae_cluster    = parse_cdhit("augustus/cluster_ceranae_augustus.clstr")
prodigal_ceranae_cluster    = parse_cdhit("prodigal/cluster_ceranae_prodigal.clstr")


#Les deuxième cluster contenant les gènes non clusterisés la première fois mais clusterisé la seconde
funannotate_ceranae_supcluster   = parse_cdhit("funannotate_ceranae/cluster_supceranae_funannotate.clstr")
glimmer_ceranae_supcluster 	  = parse_cdhit("glimmer/cluster_supceranae_glimmer.clstr")
augustus_ceranae_supcluster 	  = parse_cdhit("augustus/cluster_supceranae_augustus.clstr")
prodigal_ceranae_supcluster 	  = parse_cdhit("prodigal/cluster_supceranae_prodigal.clstr")

#Création des dico complet (contenant l'ensemble des gènes clusterisé)
complete_funannotate_ceranae = combined_dico(funannotate_ceranae_cluster, funannotate_ceranae_supcluster)
complete_glimmer_ceranae = combined_dico(glimmer_ceranae_cluster, glimmer_ceranae_supcluster)
complete_augustus_ceranae = combined_dico(augustus_ceranae_cluster,augustus_ceranae_supcluster)
complete_prodigal_ceranae = combined_dico(prodigal_ceranae_cluster,prodigal_ceranae_supcluster)


#Compter le nombre de clusters 
print("Nombre de clusters > 1 :")
print("Funannotate_ceranae :", compt_gene(complete_funannotate_ceranae))
print("Glimmer_ceranae :", compt_gene(complete_glimmer_ceranae))
print("Augustus_ceranae :", compt_gene(complete_augustus_ceranae))
print("Prodigal_ceranae :", compt_gene(complete_prodigal_ceranae))

##Pourcentage de gènes trouvés par rapport à microannot
print()
print("Gènes prédits et clusterisés (%) par rapport à microannot :") 	#Correspond au nombre de clusters > 1 prédit sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_ceranae :",(compt_gene(complete_funannotate_ceranae))/len(funannotate_ceranae_cluster)*100)	
print("Glimmer_ceranae :",(compt_gene(complete_glimmer_ceranae))/len(glimmer_ceranae_cluster)*100)	
print("Augustus_ceranae :",(compt_gene(complete_augustus_ceranae))/len(augustus_ceranae_cluster)*100)
print("Prodigal_ceranae :",(compt_gene(complete_prodigal_ceranae))/len(prodigal_ceranae_cluster)*100)

print()
print("% gènes prédits à l'identique par l'outil ") #Correspond au nombre de clusters > 1 prédit qui ont la même longueur sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_ceranae :",same_gene(funannotate_ceranae_cluster)/len(funannotate_ceranae_cluster)*100)
print("Glimmer_ceranae :",same_gene(glimmer_ceranae_cluster)/len(glimmer_ceranae_cluster)*100)	
print("Augustus_cerane :",same_gene(augustus_ceranae_cluster)/len(augustus_ceranae_cluster)*100)
print("Prodigal_ceranae :",same_gene(prodigal_ceranae_cluster)/len(prodigal_ceranae_cluster)*100)



print()
print("% gènes prédits mais pas identique ") # Nombre de gènes correctement prédits - Nombre de gènes prédits et clusterisés 
print("Funannotate_ceranae :",(compt_gene(complete_funannotate_ceranae))/len(funannotate_ceranae_cluster)*100 -same_gene(funannotate_ceranae_cluster)/len(funannotate_ceranae_cluster)*100)	
print("Glimmer_ceranae :",(compt_gene(complete_glimmer_ceranae))/len(glimmer_ceranae_cluster)*100 - same_gene(glimmer_ceranae_cluster)/len(glimmer_ceranae_cluster)*100)	
print("Augustus_ceranae :",(compt_gene(complete_augustus_ceranae))/len(augustus_ceranae_cluster)*100 -same_gene(augustus_ceranae_cluster)/len(augustus_ceranae_cluster)*100)
print("Prodigal_ceranae :",(compt_gene(complete_prodigal_ceranae))/len(prodigal_ceranae_cluster)*100 -same_gene(prodigal_ceranae_cluster)/len(prodigal_ceranae_cluster)*100)


print()
print("% gènes correct non prédits par l'outil ")# 100 - nombres de gènes prédits et clusterisés 
print("Funannotate_ceranae :",100- (compt_gene(complete_funannotate_ceranae))/len(funannotate_ceranae_cluster)*100)
print("Glimmer_ceranae :",100- (compt_gene(complete_glimmer_ceranae))/len(glimmer_ceranae_cluster)*100)	
print("Augustus_ceranae :", 100 -(compt_gene(complete_augustus_ceranae))/len(augustus_ceranae_cluster)*100)
print("Prodigal_ceranae :", 100 -(compt_gene(complete_prodigal_ceranae))/len(prodigal_ceranae_cluster)*100)
print()


#Faire les pie charts pour ceranae :

pie(100- (compt_gene(complete_funannotate_ceranae))/len(funannotate_ceranae_cluster)*100, (compt_gene(complete_funannotate_ceranae))/len(funannotate_ceranae_cluster)*100 -same_gene(funannotate_ceranae_cluster)/len(funannotate_ceranae_cluster)*100, same_gene(funannotate_ceranae_cluster)/len(funannotate_ceranae_cluster)*100, "Funannotate", "ceranae")	#Faire le graphique de funannotate_cuniculi

pie(100- (compt_gene(complete_glimmer_ceranae))/len(glimmer_ceranae_cluster)*100, (compt_gene(complete_glimmer_ceranae))/len(glimmer_ceranae_cluster)*100 - same_gene(glimmer_ceranae_cluster)/len(glimmer_ceranae_cluster)*100, same_gene(glimmer_ceranae_cluster)/len(glimmer_ceranae_cluster)*100, "Glimmer", "ceranae")

pie(100 -(compt_gene(complete_augustus_ceranae))/len(augustus_ceranae_cluster)*100,(compt_gene(complete_augustus_ceranae))/len(augustus_ceranae_cluster)*100 -same_gene(augustus_ceranae_cluster)/len(augustus_ceranae_cluster)*100, same_gene(augustus_ceranae_cluster)/len(augustus_ceranae_cluster)*100, "Augustus", "ceranae")

pie(100 -(compt_gene(complete_prodigal_ceranae))/len(prodigal_ceranae_cluster)*100, (compt_gene(complete_prodigal_ceranae))/len(prodigal_ceranae_cluster)*100 -same_gene(prodigal_ceranae_cluster)/len(prodigal_ceranae_cluster)*100, same_gene(prodigal_ceranae_cluster)/len(prodigal_ceranae_cluster)*100, "Prodigal", "ceranae")

combined_pie('ceranae')




#Maintenant je récupère les gènes prédits par tous les outils et ceux jamais prédits
#J'associe mes listes pour ceranae
predicted_funannotate_ceranae, unpredicted_funannotate_ceranae = predicted_genes(complete_funannotate_ceranae)
predicted_glimmer_ceranae, unpredicted_glimmer_ceranae = predicted_genes(complete_glimmer_ceranae)	
predicted_augustus_ceranae, unpredicted_augustus_ceranae = predicted_genes(complete_augustus_ceranae)
predicted_prodigal_ceranae, unpredicted_prodigal_ceranae = predicted_genes(complete_prodigal_ceranae)



#Faire le diagramme de Venn
venn_tools(predicted_augustus_ceranae, predicted_funannotate_ceranae, predicted_glimmer_ceranae, predicted_prodigal_ceranae, "ceranae")


#Récupérer les gènes prédits par les différents outils qui n'ont jamais était clusterisés
unclusterized_funannotate_ceranae = unclusterized(funannotate_ceranae_supcluster)
#print(len(unclusterized_funannotate_ceranae))
unclusterized_glimmer_ceranae = unclusterized(glimmer_ceranae_supcluster)
#print(len(unclusterized_glimmer_ceranae))
unclusterized_augustus_ceranae = unclusterized(augustus_ceranae_supcluster)
#print(len(unclusterized_augustus_ceranae))
unclusterized_prodigal_ceranae = unclusterized(prodigal_ceranae_supcluster)
#print(len(unclusterized_prodigal_ceranae))


#Trouver les gènes de la base de données qui n'ont jamais été clusterisé
#Je stock le nom de ces gènes dans 
unpredicted_genes_ceranae = open("result_ceranae/unpredicted_genes_ceranae", 'w') 

for gene in unpredicted_genes(unpredicted_augustus_ceranae,unpredicted_funannotate_ceranae, unpredicted_glimmer_ceranae, unpredicted_prodigal_ceranae):
	print(gene, file=unpredicted_genes_ceranae)
	
#Je stock les gènes non prédits dans cette liste 	
liste_unpredicted_ceranae = unpredicted_genes(unpredicted_augustus_ceranae,unpredicted_funannotate_ceranae, unpredicted_glimmer_ceranae, unpredicted_prodigal_ceranae)

#Je peux également me servir de cette fonction pour garder uniquement les gènes prédits par tous les outils:

predicted_genes_ceranae = open("result_ceranae/predicted_genes_ceranae", 'w') 

for gene in unpredicted_genes(predicted_augustus_ceranae, predicted_funannotate_ceranae, predicted_glimmer_ceranae, predicted_prodigal_ceranae):
	print(gene, file=predicted_genes_ceranae)
	
#Je stock les gènes prédits dans cette liste 	
liste_predicted_ceranae = unpredicted_genes(predicted_augustus_ceranae, predicted_funannotate_ceranae, predicted_glimmer_ceranae, predicted_prodigal_ceranae)

#Pour créer mon histogramme j'ai besoin des listes des tailles	
length_unpredicted_ceranae = compt_length(liste_unpredicted_ceranae)	#Liste des tailles de gènes non prédits
length_predicted_ceranae = compt_length(liste_predicted_ceranae)	#Liste des tailles de gènes prédits par les 4 outils 

##Faire l'histo
#hist(length_unpredicted_ceranae, length_predicted_ceranae, "ceranae")

##Faire le bar_charts pour ceranae


#bar_charts(complete_glimmer_ceranae,(len(gene_predict('augustus/Proteomes_N.ceranae_augustus'))-len(unclusterized_augustus_ceranae)), (len(gene_predict('funannotate_ceranae/Proteomes_N.ceranae_funannotate'))- len(unclusterized_funannotate_ceranae)), (len(gene_predict('glimmer/Proteomes_N.ceranae_glimmer'))- len(unclusterized_glimmer_ceranae)), (len(gene_predict('prodigal/Proteomes_N.ceranae_prodigal'))- len(unclusterized_prodigal_ceranae)), len(unclusterized_augustus_ceranae), len(unclusterized_funannotate_ceranae), len(unclusterized_glimmer_ceranae), len(unclusterized_prodigal_ceranae), "ceranae")

"""
#Faire les feuilles csv
stdout('GLIMMER', 'result_ceranae/glimmer.csv', complete_glimmer_ceranae, unclusterized_glimmer_ceranae)		#Glimmer_ceranae

stdout('AUGUSTUS', 'result_ceranae/augustus.csv', complete_augustus_ceranae, unclusterized_augustus_ceranae)		#Augustus_ceranae
stdout('FUNANNOTATE','result_ceranae/funannotate.csv', complete_funannotate_ceranae, unclusterized_funannotate_ceranae)	#Funannotate_ceranae
stdout('PRODIGAL', 'result_ceranae/prodigal.csv', complete_prodigal_ceranae, unclusterized_prodigal_ceranae)		#Prodigal_ceranae

"""


#Print le bar_charts du nombre de gènes dans la base de données et le nombre de gènes prédits et clusterisé via les 4 outils
#bar_charts_genes(len(complete_glimmer_ceranae), (len(complete_glimmer_ceranae)-len(liste_unpredicted_ceranae)), "ceranae")



#======================================================================================================
print()
print("===========================POUR ALGERAE==============================")
print("Nombre total de gènes prédits par l'outil")
#print("Funannotate :", len(gene_predict('funannotate_algerae/Proteomes_A.algerae_funannotate')))
print("Glimmer :", len(gene_predict('glimmer/Proteomes_A.algerae_glimmer')))
print("Augustus :", len(gene_predict('augustus/Proteomes_A.algerae_augustus')))
print("Prodigal :", len(gene_predict('prodigal/Proteomes_A.algerae_prodigal')))
print()

##Dico des bases de données clean de microannot
dico_algerae = gene_predict('Proteomes_A.algerae.txt')

##Dico des gènes prédits juste par l'outil
#funannotate_algerae = gene_predict('funannotate_algerae/cluster_algerae_funannotate')
glimmer_algerae = gene_predict('glimmer/cluster_algerae_glimmer')
augustus_algerae = gene_predict('augustus/cluster_algerae_augustus')
prodigal_algerae = gene_predict('prodigal/cluster_algerae_prodigal')

##Création des dicos qui contiennent tous les clusters pour bieneusi
#funannotate_algerae_cluster = parse_cdhit("funannotate_algerae/cluster_algerae_funannotate.clstr")
glimmer_algerae_cluster       = parse_cdhit("glimmer/cluster_algerae_glimmer.clstr")
augustus_algerae_cluster    = parse_cdhit("augustus/cluster_algerae_augustus.clstr")
prodigal_algerae_cluster    = parse_cdhit("prodigal/cluster_algerae_prodigal.clstr")


#Les deuxième cluster contenant les gènes non clusterisés la première fois mais clusterisé la seconde
#funannotate_algerae_supcluster   = parse_cdhit("funannotate_algerae/cluster_supalgerae_funannotate.clstr")
glimmer_algerae_supcluster 	  = parse_cdhit("glimmer/cluster_supalgerae_glimmer.clstr")
augustus_algerae_supcluster 	  = parse_cdhit("augustus/cluster_supalgerae_augustus.clstr")
prodigal_algerae_supcluster 	  = parse_cdhit("prodigal/cluster_supalgerae_prodigal.clstr")

#Création des dico complet (contenant l'ensemble des gènes clusterisé)
#complete_funannotate_algerae = combined_dico(funannotate_algerae_cluster, funannotate_algerae_supcluster)
complete_glimmer_algerae = combined_dico(glimmer_algerae_cluster, glimmer_algerae_supcluster)
complete_augustus_algerae = combined_dico(augustus_algerae_cluster,augustus_algerae_supcluster)
complete_prodigal_algerae = combined_dico(prodigal_algerae_cluster,prodigal_algerae_supcluster)


#Compter le nombre de clusters 
print("Nombre de clusters > 1 :")
#print("Funannotate_algerae :", compt_gene(complete_funannotate_algerae))
print("Glimmer_algerae :", compt_gene(complete_glimmer_algerae))
print("Augustus_algerae :", compt_gene(complete_augustus_algerae))
print("Prodigal_algerae :", compt_gene(complete_prodigal_algerae))

##Pourcentage de gènes trouvés par rapport à microannot
print()
print("Gènes prédits et clusterisés (%) par rapport à microannot :") 	#Correspond au nombre de clusters > 1 prédit sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
#print("Funannotate_algerae :",(compt_gene(complete_funannotate_algerae))/len(funannotate_algerae_cluster)*100)	
print("Glimmer_algerae :",(compt_gene(complete_glimmer_algerae))/len(glimmer_algerae_cluster)*100)	
print("Augustus_algerae :",(compt_gene(complete_augustus_algerae))/len(augustus_algerae_cluster)*100)
print("Prodigal_algerae :",(compt_gene(complete_prodigal_algerae))/len(prodigal_algerae_cluster)*100)

print()
print("% gènes prédits à l'identique par l'outil ") #Correspond au nombre de clusters > 1 prédit qui ont la même longueur sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
#print("Funannotate_algerae :",same_gene(funannotate_algerae_cluster)/len(funannotate_algerae_cluster)*100)
print("Glimmer_algerae :",same_gene(glimmer_algerae_cluster)/len(glimmer_algerae_cluster)*100)	
print("Augustus_algerae :",same_gene(augustus_algerae_cluster)/len(augustus_algerae_cluster)*100)
print("Prodigal_algerae :",same_gene(prodigal_algerae_cluster)/len(prodigal_algerae_cluster)*100)



print()
print("% gènes prédits mais pas identique ") # Nombre de gènes correctement prédits - Nombre de gènes prédits et clusterisés 
#print("Funannotate_algerae :",(compt_gene(complete_funannotate_algerae))/len(funannotate_algerae_cluster)*100 -same_gene(funannotate_algerae_cluster)/len(funannotate_algerae_cluster)*100)	
print("Glimmer_algerae :",(compt_gene(complete_glimmer_algerae))/len(glimmer_algerae_cluster)*100 - same_gene(glimmer_algerae_cluster)/len(glimmer_algerae_cluster)*100)	
print("Augustus_algerae :",(compt_gene(complete_augustus_algerae))/len(augustus_algerae_cluster)*100 -same_gene(augustus_algerae_cluster)/len(augustus_algerae_cluster)*100)
print("Prodigal_algerae :",(compt_gene(complete_prodigal_algerae))/len(prodigal_algerae_cluster)*100 -same_gene(prodigal_algerae_cluster)/len(prodigal_algerae_cluster)*100)


print()
print("% gènes correct non prédits par l'outil ")# 100 - nombres de gènes prédits et clusterisés 
#print("Funannotate_algerae :",100- (compt_gene(complete_funannotate_algerae))/len(funannotate_algerae_cluster)*100)
print("Glimmer_algerae :",100- (compt_gene(complete_glimmer_algerae))/len(glimmer_algerae_cluster)*100)	
print("Augustus_algerae :", 100 -(compt_gene(complete_augustus_algerae))/len(augustus_algerae_cluster)*100)
print("Prodigal_algerae :", 100 -(compt_gene(complete_prodigal_algerae))/len(prodigal_algerae_cluster)*100)
print()


#Faire les pie charts pour algerae :

#pie(100- (compt_gene(complete_funannotate_algerae))/len(funannotate_algerae_cluster)*100, (compt_gene(complete_funannotate_algerae))/len(funannotate_algerae_cluster)*100 -same_gene(funannotate_algerae_cluster)/len(funannotate_algerae_cluster)*100, same_gene(funannotate_algerae_cluster)/len(funannotate_algerae_cluster)*100, "Funannotate", "algerae")	#Faire le graphique de funannotate_cuniculi

pie(100- (compt_gene(complete_glimmer_algerae))/len(glimmer_algerae_cluster)*100, (compt_gene(complete_glimmer_algerae))/len(glimmer_algerae_cluster)*100 - same_gene(glimmer_algerae_cluster)/len(glimmer_algerae_cluster)*100, same_gene(glimmer_algerae_cluster)/len(glimmer_algerae_cluster)*100, "Glimmer", "algerae")

pie(100 -(compt_gene(complete_augustus_algerae))/len(augustus_algerae_cluster)*100,(compt_gene(complete_augustus_algerae))/len(augustus_algerae_cluster)*100 -same_gene(augustus_algerae_cluster)/len(augustus_algerae_cluster)*100, same_gene(augustus_algerae_cluster)/len(augustus_algerae_cluster)*100, "Augustus", "algerae")

pie(100 -(compt_gene(complete_prodigal_algerae))/len(prodigal_algerae_cluster)*100, (compt_gene(complete_prodigal_algerae))/len(prodigal_algerae_cluster)*100 -same_gene(prodigal_algerae_cluster)/len(prodigal_algerae_cluster)*100, same_gene(prodigal_algerae_cluster)/len(prodigal_algerae_cluster)*100, "Prodigal", "algerae")


#combined_pie('algerae')
"""
#Maintenant je récupère les gènes prédits par tous les outils et ceux jamais prédits
#J'associe mes listes pour algerae
predicted_funannotate_algerae, unpredicted_funannotate_algerae = predicted_genes(complete_funannotate_algerae)
predicted_glimmer_algerae, unpredicted_glimmer_algerae = predicted_genes(complete_glimmer_algerae)	
predicted_augustus_algerae, unpredicted_augustus_algerae = predicted_genes(complete_augustus_algerae)
predicted_prodigal_algerae, unpredicted_prodigal_algerae = predicted_genes(complete_prodigal_algerae)



#Faire le diagramme de Venn
venn_tools(predicted_augustus_algerae, predicted_funannotate_algerae, predicted_glimmer_algerae, predicted_prodigal_algerae, "algerae")


#Récupérer les gènes prédits par les différents outils qui n'ont jamais était clusterisés
unclusterized_funannotate_algerae = unclusterized(funannotate_algerae_supcluster)
unclusterized_glimmer_algerae = unclusterized(glimmer_algerae_supcluster)
unclusterized_augustus_algerae = unclusterized(augustus_algerae_supcluster))
unclusterized_prodigal_algerae = unclusterized(prodigal_algerae_supcluster)

#Trouver les gènes de la base de données qui n'ont jamais été clusterisé
#Je stock le nom de ces gènes dans 
unpredicted_genes_algerae = open("result_algerae/unpredicted_genes_algerae", 'w') 

for gene in unpredicted_genes(unpredicted_augustus_algerae,unpredicted_funannotate_algerae, unpredicted_glimmer_algerae, unpredicted_prodigal_algerae):
	print(gene, file=unpredicted_genes_algerae)
	
#Je stock les gènes non prédits dans cette liste 	
liste_unpredicted_algerae = unpredicted_genes(unpredicted_augustus_algerae,unpredicted_funannotate_algerae, unpredicted_glimmer_algerae, unpredicted_prodigal_algerae)

#Je peux également me servir de cette fonction pour garder uniquement les gènes prédits par tous les outils:

predicted_genes_algerae = open("result_algerae/predicted_genes_algerae", 'w') 

for gene in unpredicted_genes(predicted_augustus_algerae, predicted_funannotate_algerae, predicted_glimmer_algerae, predicted_prodigal_algerae):
	print(gene, file=predicted_genes_algerae)
	
#Je stock les gènes prédits dans cette liste 	
liste_predicted_algerae = unpredicted_genes(predicted_augustus_algerae, predicted_funannotate_algerae, predicted_glimmer_algerae, predicted_prodigal_algerae)

#Pour créer mon histogramme j'ai besoin des listes des tailles	
length_unpredicted_algerae = compt_length(liste_unpredicted_algerae)	#Liste des tailles de gènes non prédits
length_predicted_algerae = compt_length(liste_predicted_algerae)	#Liste des tailles de gènes prédits par les 4 outils 

#Faire l'histo
hist(length_unpredicted_algerae, length_predicted_algerae, "algerae")

#Faire le bar_charts pour algerae
bar_charts(complete_glimmer_algerae,(len(gene_predict('augustus/Proteomes_A.algerae_augustus'))-len(unclusterized_augustus_algerae)), (len(gene_predict('funannotate_algerae/Proteomes_A.algerae_funannotate'))- len(unclusterized_funannotate_algerae)), (len(gene_predict('glimmer/Proteomes_A.algerae_glimmer'))- len(unclusterized_glimmer_algerae)), (len(gene_predict('prodigal/Proteomes_A.algerae_prodigal'))- len(unclusterized_prodigal_algerae)), len(unclusterized_augustus_algerae), len(unclusterized_funannotate_algerae), len(unclusterized_glimmer_algerae), len(unclusterized_prodigal_algerae), "algerae")


#Faire les feuilles csv
stdout('GLIMMER', 'result_algerae/glimmer.csv', complete_glimmer_algerae, unclusterized_glimmer_algerae)		#Glimmer_algerae

stdout('AUGUSTUS', 'result_algerae/augustus.csv', complete_augustus_algerae, unclusterized_augustus_algerae)		#Augustus_algerae
stdout('FUNANNOTATE','result_algerae/funannotate.csv', complete_funannotate_algerae, unclusterized_funannotate_algerae)	#Funannotate_algerae
stdout('PRODIGAL', 'result_algerae/prodigal.csv', complete_prodigal_algerae, unclusterized_prodigal_algerae)		#Prodigal_algerae


#Print le bar_charts du nombre de gènes dans la base de données et le nombre de gènes prédits et clusterisé via les 4 outils
#bar_charts_genes(len(complete_glimmer_algerae), (len(complete_glimmer_algerae)-len(liste_unpredicted_algerae)), "algerae")
"""

print()
print("===========================POUR PARISII==============================")
print("Nombre total de gènes prédits par l'outil")
print("Funannotate :", len(gene_predict('funannotate_parisii/Proteomes_N.parisii_funannotate')))
print("Glimmer :", len(gene_predict('glimmer/Proteomes_N.parisii_glimmer')))
print("Augustus :", len(gene_predict('augustus/Proteomes_N.parisii_augustus')))
print("Prodigal :", len(gene_predict('prodigal/Proteomes_N.parisii_prodigal')))
print()

##Dico des bases de données clean de microannot
dico_parisii = gene_predict('CDS_nematocida.fasta')

##Dico des gènes prédits juste par l'outil
funannotate_parisii = gene_predict('funannotate_parisii/cluster_parisii_funannotate')
glimmer_parisii = gene_predict('glimmer/cluster_parisii_glimmer')
augustus_parisii = gene_predict('augustus/cluster_parisii_augustus')
prodigal_parisii = gene_predict('prodigal/cluster_parisii_prodigal')

##Création des dicos qui contiennent tous les clusters pour bieneusi
funannotate_parisii_cluster = parse_cdhit("funannotate_parisii/cluster_parisii_funannotate.clstr")
glimmer_parisii_cluster       = parse_cdhit("glimmer/cluster_parisii_glimmer.clstr")
augustus_parisii_cluster    = parse_cdhit("augustus/cluster_parisii_augustus.clstr")
prodigal_parisii_cluster    = parse_cdhit("prodigal/cluster_parisii_prodigal.clstr")



#Les deuxième cluster contenant les gènes non clusterisés la première fois mais clusterisé la seconde
funannotate_parisii_supcluster   = parse_cdhit("funannotate_parisii/cluster_supparisii_funannotate.clstr")
glimmer_parisii_supcluster 	  = parse_cdhit("glimmer/cluster_supparisii_glimmer.clstr")
augustus_parisii_supcluster 	  = parse_cdhit("augustus/cluster_supparisii_augustus.clstr")
prodigal_parisii_supcluster 	  = parse_cdhit("prodigal/cluster_supparisii_prodigal.clstr")

#Création des dico complet (contenant l'ensemble des gènes clusterisé)
complete_funannotate_parisii = combined_dico(funannotate_parisii_cluster, funannotate_parisii_supcluster)
complete_glimmer_parisii = combined_dico(glimmer_parisii_cluster, glimmer_parisii_supcluster)
complete_augustus_parisii = combined_dico(augustus_parisii_cluster,augustus_parisii_supcluster)
complete_prodigal_parisii = combined_dico(prodigal_parisii_cluster,prodigal_parisii_supcluster)


#Compter le nombre de clusters 
print("Nombre de clusters > 1 :")
print("Funannotate_parisii :", compt_gene(complete_funannotate_parisii))
print("Glimmer_parisii :", compt_gene(complete_glimmer_parisii))
print("Augustus_parisii :", compt_gene(complete_augustus_parisii))
print("Prodigal_parisii :", compt_gene(complete_prodigal_parisii))



##Pourcentage de gènes trouvés par rapport à microannot
print()
print("Gènes prédits et clusterisés (%) par rapport à microannot :") 	#Correspond au nombre de clusters > 1 prédit sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_parisii :",(compt_gene(complete_funannotate_parisii))/len(funannotate_parisii_cluster)*100)	
print("Glimmer_parisii :",(compt_gene(complete_glimmer_parisii))/len(glimmer_parisii_cluster)*100)	
print("Augustus_parisii :",(compt_gene(complete_augustus_parisii))/len(augustus_parisii_cluster)*100)
print("Prodigal_parisii :",(compt_gene(complete_prodigal_parisii))/len(prodigal_parisii_cluster)*100)

print()
print("% gènes prédits à l'identique par l'outil ") #Correspond au nombre de clusters > 1 prédit qui ont la même longueur sur le nombre de cluster total (nb de gènes prédits par Microannot = 1978)
print("Funannotate_parisii :",same_gene(funannotate_parisii_cluster)/len(funannotate_parisii_cluster)*100)
print("Glimmer_parisii :",same_gene(glimmer_parisii_cluster)/len(glimmer_parisii_cluster)*100)	
print("Augustus_parisii :",same_gene(augustus_parisii_cluster)/len(augustus_parisii_cluster)*100)
print("Prodigal_parisii :",same_gene(prodigal_parisii_cluster)/len(prodigal_parisii_cluster)*100)


print()
print("% gènes prédits mais pas identique ") # Nombre de gènes correctement prédits - Nombre de gènes prédits et clusterisés 
print("Funannotate_parisii :",(compt_gene(complete_funannotate_parisii))/len(funannotate_parisii_cluster)*100 -same_gene(funannotate_parisii_cluster)/len(funannotate_parisii_cluster)*100)	
print("Glimmer_parisii :",(compt_gene(complete_glimmer_parisii))/len(glimmer_parisii_cluster)*100 - same_gene(glimmer_parisii_cluster)/len(glimmer_parisii_cluster)*100)	
print("Augustus_parisii :",(compt_gene(complete_augustus_parisii))/len(augustus_parisii_cluster)*100 -same_gene(augustus_parisii_cluster)/len(augustus_parisii_cluster)*100)
print("Prodigal_parisii :",(compt_gene(complete_prodigal_parisii))/len(prodigal_parisii_cluster)*100 -same_gene(prodigal_parisii_cluster)/len(prodigal_parisii_cluster)*100)


print()
print("% gènes correct non prédits par l'outil ")# 100 - nombres de gènes prédits et clusterisés 
print("Funannotate_parisii :",100- (compt_gene(complete_funannotate_parisii))/len(funannotate_parisii_cluster)*100)
print("Glimmer_parisii :",100- (compt_gene(complete_glimmer_parisii))/len(glimmer_parisii_cluster)*100)	
print("Augustus_parisii :", 100 -(compt_gene(complete_augustus_parisii))/len(augustus_parisii_cluster)*100)
print("Prodigal_parisii :", 100 -(compt_gene(complete_prodigal_parisii))/len(prodigal_parisii_cluster)*100)
print()



#Faire les pie charts pour parisii :

pie(100- (compt_gene(complete_funannotate_parisii))/len(funannotate_parisii_cluster)*100, (compt_gene(complete_funannotate_parisii))/len(funannotate_parisii_cluster)*100 -same_gene(funannotate_parisii_cluster)/len(funannotate_parisii_cluster)*100, same_gene(funannotate_parisii_cluster)/len(funannotate_parisii_cluster)*100, "Funannotate", "parisii")	#Faire le graphique de funannotate_cuniculi

pie(100- (compt_gene(complete_glimmer_parisii))/len(glimmer_parisii_cluster)*100, (compt_gene(complete_glimmer_parisii))/len(glimmer_parisii_cluster)*100 - same_gene(glimmer_parisii_cluster)/len(glimmer_parisii_cluster)*100, same_gene(glimmer_parisii_cluster)/len(glimmer_parisii_cluster)*100, "Glimmer", "parisii")

pie(100 -(compt_gene(complete_augustus_parisii))/len(augustus_parisii_cluster)*100,(compt_gene(complete_augustus_parisii))/len(augustus_parisii_cluster)*100 -same_gene(augustus_parisii_cluster)/len(augustus_parisii_cluster)*100, same_gene(augustus_parisii_cluster)/len(augustus_parisii_cluster)*100, "Augustus", "parisii")

pie(100 -(compt_gene(complete_prodigal_parisii))/len(prodigal_parisii_cluster)*100, (compt_gene(complete_prodigal_parisii))/len(prodigal_parisii_cluster)*100 -same_gene(prodigal_parisii_cluster)/len(prodigal_parisii_cluster)*100, same_gene(prodigal_parisii_cluster)/len(prodigal_parisii_cluster)*100, "Prodigal", "parisii")

combined_pie('parisii')

#Maintenant je récupère les gènes prédits par tous les outils et ceux jamais prédits
#J'associe mes listes pour parisii
predicted_funannotate_parisii, unpredicted_funannotate_parisii = predicted_genes(complete_funannotate_parisii)
predicted_glimmer_parisii, unpredicted_glimmer_parisii = predicted_genes(complete_glimmer_parisii)	
predicted_augustus_parisii, unpredicted_augustus_parisii = predicted_genes(complete_augustus_parisii)
predicted_prodigal_parisii, unpredicted_prodigal_parisii = predicted_genes(complete_prodigal_parisii)



#Faire le diagramme de Venn
venn_tools(predicted_augustus_parisii, predicted_funannotate_parisii, predicted_glimmer_parisii, predicted_prodigal_parisii, "parisii")

#Récupérer les gènes prédits par les différents outils qui n'ont jamais était clusterisés
unclusterized_funannotate_parisii = unclusterized(funannotate_parisii_supcluster)
unclusterized_glimmer_parisii = unclusterized(glimmer_parisii_supcluster)
unclusterized_augustus_parisii = unclusterized(augustus_parisii_supcluster)
unclusterized_prodigal_parisii = unclusterized(prodigal_parisii_supcluster)

#Trouver les gènes de la base de données qui n'ont jamais été clusterisé
#Je stock le nom de ces gènes dans 
unpredicted_genes_parisii = open("result_parisii/unpredicted_genes_parissi", 'w') 

for gene in unpredicted_genes(unpredicted_augustus_parisii,unpredicted_funannotate_parisii, unpredicted_glimmer_parisii, unpredicted_prodigal_parisii):
	print(gene, file=unpredicted_genes_parisii)
	
#Je stock les gènes non prédits dans cette liste 	
liste_unpredicted_parisii = unpredicted_genes(unpredicted_augustus_parisii,unpredicted_funannotate_parisii, unpredicted_glimmer_parisii, unpredicted_prodigal_parisii)

#Je peux également me servir de cette fonction pour garder uniquement les gènes prédits par tous les outils:

predicted_genes_parisii = open("result_parisii/predicted_genes_parisii", 'w') 

for gene in unpredicted_genes(predicted_augustus_parisii, predicted_funannotate_parisii, predicted_glimmer_parisii, predicted_prodigal_parisii):
	print(gene, file=predicted_genes_parisii)
	
#Je stock les gènes prédits dans cette liste 	
liste_predicted_parisii = unpredicted_genes(predicted_augustus_parisii, predicted_funannotate_parisii, predicted_glimmer_parisii, predicted_prodigal_parisii)


#Pour créer mon histogramme j'ai besoin des listes des tailles	
length_unpredicted_parisii = compt_length(liste_unpredicted_parisii)	#Liste des tailles de gènes non prédits
length_predicted_parisii = compt_length(liste_predicted_parisii)	#Liste des tailles de gènes prédits par les 4 outils 



#Faire l'histo
hist(length_unpredicted_parisii, length_predicted_parisii, "parisii")

#Faire le bar_charts pour parisii
bar_charts(complete_glimmer_parisii,(len(gene_predict('augustus/Proteomes_N.parisii_augustus'))-len(unclusterized_augustus_parisii)), (len(gene_predict('funannotate_parisii/Proteomes_N.parisii_funannotate'))- len(unclusterized_funannotate_parisii)), (len(gene_predict('glimmer/Proteomes_N.parisii_glimmer'))- len(unclusterized_glimmer_parisii)), (len(gene_predict('prodigal/Proteomes_N.parisii_prodigal'))- len(unclusterized_prodigal_parisii)), len(unclusterized_augustus_parisii), len(unclusterized_funannotate_parisii), len(unclusterized_glimmer_parisii), len(unclusterized_prodigal_parisii), "parisii")


#Faire les feuilles csv
stdout('GLIMMER', 'result_parisii/glimmer.csv', complete_glimmer_parisii, unclusterized_glimmer_parisii)		#Glimmer_parisii

stdout('AUGUSTUS', 'result_parisii/augustus.csv', complete_augustus_parisii, unclusterized_augustus_parisii)		#Augustus_parisii
stdout('FUNANNOTATE','result_parisii/funannotate.csv', complete_funannotate_parisii, unclusterized_funannotate_parisii)	#Funannotate_parisii
stdout('PRODIGAL', 'result_parisii/prodigal.csv', complete_prodigal_parisii, unclusterized_prodigal_parisii)		#Prodigal_parisii

#Print le bar_charts du nombre de gènes dans la base de données et le nombre de gènes prédits et clusterisé via les 4 outils
bar_charts_genes(len(complete_glimmer_parisii), (len(complete_glimmer_parisii)-len(liste_unpredicted_parisii)), "parisii")


