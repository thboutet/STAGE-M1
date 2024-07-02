#Packages nécessaire à installer
"""
pip install pandas
pip install matplotlib
pip install numpy
pip install venny4py
"""

import pandas as pd 			#pandas pour traiter les fichiers ods
import matplotlib.pyplot as plt		#Pour bar charts
import numpy as np			#Pour données du bar charts 
from venny4py.venny4py import *		#Diagramme de Venn


##JE VEUX PRINT DES BAR CHARTS DES DIFFERENTES ERREURS 

#Je définis d'abord une fonction pour me permettre de bien compter mes différentes erreur (car sinon elles sont 'unique' et ici je veux regrouper les erreurs starts entre elles, end entre elles...)
def error (df) :
	#Trois types d'erreurs possibles : 
	if ' error - start' in df :
		return ('error - start')
	elif ' error - end' in df :
		return ('error - end')	
	elif ' error ' in df :
		return ('error')	
	else :	
		return None			#None car ici je veux seulement compter et garder mes types d'erreur
	
#Fonction pour récupérer les données de chaques feuilles pour le graphiques bar
def bar (path, microsporidie) :
	#Création des différents df (1 par outil)
	df1 = pd.read_excel(path, sheet_name='augustus')
	df2 = pd.read_excel(path, sheet_name='funannotate')
	df3 = pd.read_excel(path, sheet_name='glimmer')
	df4 = pd.read_excel(path, sheet_name='prodigal')
	
	#J'applique ma fonction à ma df pour créer une nouvelle colonne composé uniquement des erreurs
	df1['normalized_error']= df1.iloc[:,2].apply(error)		
	df2['normalized_error']= df2.iloc[:,2].apply(error)
	df3['normalized_error']= df3.iloc[:,2].apply(error)
	df4['normalized_error']= df4.iloc[:,2].apply(error)
	
	
	#Je m'assure de garder seulement les lignes avec les erreurs, je les comptes en fonction de leur type et je les trie 
	occ_augustus = df1['normalized_error'].dropna().value_counts().sort_index()	
	occ_funannotate = df2['normalized_error'].dropna().value_counts().sort_index()
	occ_glimmer = df3['normalized_error'].dropna().value_counts().sort_index()
	occ_prodigal = df4['normalized_error'].dropna().value_counts().sort_index()

	
	#Faire les figures bars 
	plt.figure(figsize=(20,6))

	#Je souhaite une figure avec des données normalisées sur 100 %
	total_aug = occ_augustus.sum()			#Somme totale de tout les types d'erreurs
	percent_aug = occ_augustus/total_aug *100	#Calcul du pourcentage de chaque type d'erreur
	
	total_fun = occ_funannotate.sum()
	percent_fun = occ_funannotate/total_fun *100

	total_glim = occ_glimmer.sum()
	percent_glim = occ_glimmer/total_glim *100
	
	total_prod = occ_prodigal.sum()
	percent_prod = occ_prodigal/total_prod *100

	##Je fais ma figure normalisé sur 100 %
	x=['Augustus', 'Funannotate', 'Glimmer', 'Prodigal']									#Nom des bars 
	y1=np.array([percent_aug.tolist()[0], percent_fun.tolist()[0], percent_glim.tolist()[0], percent_prod.tolist()[0]])	#Erreur indeterminé
	y2=np.array([percent_aug.tolist()[1],percent_fun.tolist()[1],percent_glim.tolist()[1],percent_prod.tolist()[1]])	#Erreur end
	y3=np.array([percent_aug.tolist()[2],percent_fun.tolist()[2],percent_glim.tolist()[2],percent_prod.tolist()[2]])	#Erreur start
	
	fig, ax = plt.subplots()
	
	ax.bar(x, y1, color='#CC0000', edgecolor = "black", linewidth = 1)			#Pourcentage d'erreur non spécifier ( en rouge, avec des contours noirs)
	ax.bar(x, y2, bottom=y1, color='#007FFF', edgecolor = "black", linewidth = 1)		#Pourcentage d'erreur end (en bleu, avec des contours noirs)
	ax.bar(x, y3, bottom=y1+y2, color='#FCDC12', edgecolor = "black", linewidth = 1)	#Pourcentage d'erreur start (en jaune, avec des contours noirs)
	
	#Print les pourcentages 
	for bar in ax.patches :
		ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() / 2 + bar.get_y(), round(bar.get_height()), ha = 'center', color = 'w', weight = 'bold')
	
	
	plt.ylabel("error (%)")
	plt.legend(["Other", "End error", "Start error"])	#Nom des légendes
	plt.suptitle("Percentage of errors for each tool , "+microsporidie, fontweight='bold')		#Titre (en gras)
	plt.title("(standartized out of 100 %)")		#Sous titre pour spécifier le format des données 

	
		
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Bar/"+microsporidie+"100%.png")
	plt.close()
	
	#=======================================================
	##Je fais une deuxième figure avec les valeurs de bases 
	plt.figure(figsize=(8,6))
	x_=['Augustus', 'Funannotate', 'Glimmer', 'Prodigal']
	y_1=np.array([occ_augustus.tolist()[0], occ_funannotate.tolist()[0], occ_glimmer.tolist()[0], occ_prodigal.tolist()[0]])		#Erreur indeterminé
	y_2=np.array([occ_augustus.tolist()[1],occ_funannotate.tolist()[1],occ_glimmer.tolist()[1],occ_prodigal.tolist()[1]])			#Erreur end
	y_3=np.array([occ_augustus.tolist()[2],occ_funannotate.tolist()[2],occ_glimmer.tolist()[2],occ_prodigal.tolist()[2]])			#Erreur start
	fig, ax = plt.subplots()
	
	ax.bar(x_, y_1, color='#CC0000', edgecolor = "black", linewidth = 1)				#Pourcentage d'erreur non spécifier
	ax.bar(x_, y_2, bottom=y_1, color='#007FFF', edgecolor = "black", linewidth = 1)		#Pourcentage d'erreur end
	ax.bar(x_, y_3, bottom=y_1+y_2, color='#FCDC12', edgecolor = "black", linewidth = 1)		#Pourcentage d'erreur start 
	
	#Print les pourcentages 
	for bar in ax.patches :
		ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() / 2 + bar.get_y(), round(bar.get_height()), ha = 'center', color = 'w', weight = 'bold')
	
	#Je vais print le nombre total d'erreurs trouvés pour chaque outil
	total_values = np.add(y_1, y_2)			#Nombre d'erreur indeterminé + erreur end
	total_values = np.add(total_values, y_3)	#indeterminé + end + start
	
	#Traiter chaque bars en fonction de son total
	for i, total in enumerate (total_values):
		ax.text(i, total + 0.5, round(total), ha = 'center', weight = 'bold', color = 'black')
	
	plt.ylabel("Number of errors")				#Nom axe y
	plt.legend(["Other", "End error", "Start error"])	#Nom des légendes
	plt.title("Number of errors for each tool , "+ microsporidie, fontweight='bold')	#Titre
	
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Bar/"+microsporidie+".png")
	plt.close()
	
##Faire les graphiques bar pour les différentes espèces (enlever le # de la microsporidie souhaité)
#bar('result_cuniculi/CUNICULI.ods', 'cuniculi')
#bar('result_bieneusi/BIENEUSI.ods', 'bieneusi')
#bar('result_ceranae/CERANAE.ods', 'ceranae')
#bar('result_parisii/PARISII.ods', 'parisii')

def venn (path, microsporidie):
	#Création des différents df (1 par outil)
	df1 = pd.read_excel(path, sheet_name='augustus')
	df2 = pd.read_excel(path, sheet_name='funannotate')
	df3 = pd.read_excel(path, sheet_name='glimmer')
	df4 = pd.read_excel(path, sheet_name='prodigal')
	
	
	#Garder seulement les lignes des gènes prédits correctements (Colonne error = OK)	
	df1 = df1[df1[df1.columns[2]] == ' OK']			#Je conserve dans mon df seulement les lignes qui ont "OK" pour la colonnes 3 qui est celle des erreurs					
	df2 = df2[df2[df2.columns[2]] == ' OK']	
	df3 = df3[df3[df3.columns[2]] == ' OK']	
	df4 = df4[df4[df4.columns[2]] == ' OK']
	
	#Maintenant je recupère mes listes des gènes correctement prédits
	liste_aug = df1[df1.columns[0]].tolist()
	liste_fun = df2[df2.columns[0]].tolist()
	liste_gli = df3[df3.columns[0]].tolist()
	liste_prod = df4[df4.columns[0]].tolist()


	#Je attribue les listes des gènes correctement prédits à mes outils et je les transforme en set  
	sets = {'Augustus': set(liste_aug),
		'Funannotate' : set(liste_fun),
		'Glimmer': set(liste_gli),
		'Prodigal': set(liste_prod)}
	
	#Création des diagramme de Venn
	venny4py(sets=sets)		
	
	#Titre
	plt.title("Correctly predicted genes for "+microsporidie, fontweight='bold')
	
	#Enregistrer la figure 
	plt.savefig("result_"+microsporidie+"/Venn/Correct_genes_"+microsporidie+".png")
	plt.close()

##Print les graphiques pour les microsporidie 
#venn('result_cuniculi/CUNICULI.ods', "cuniculi")
#venn('result_bieneusi/BIENEUSI.ods', "bieneusi")
#venn('result_ceranae/CERANAE.ods', "ceranae")
#venn('result_parisii/PARISII.ods', "parisii")




