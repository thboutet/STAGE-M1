#Script pour transformer les génomes nt en aa
from collections import defaultdict

##Fonction que je vais utiliser pour traduire les séquences nt en aa
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
        codon = dna_seq[i:i+3]		#Je créer mes groupes de 3 nt 
        if codon in codon_table:
            aa = codon_table[codon]	#Je cherche si mes 3 nt ont une correspondances en a.a
            if aa == '_':  		# Si codon stop j'arrête
                break
            protein_seq += aa
        else:
            protein_seq += 'X'  	#Si le codon n'est pas dans la table des codons alors je print X
        
    return protein_seq
    
    
    


#Fonction pour créer des dicos de mes génomes 
def dico_aa (file) :
	dico_aa = defaultdict(str)
	with open (file,'r') as f1 :
		for lig in f1 :
			lig = lig.rstrip() 
			if lig.startswith('>') :
				chromosome = lig[1:]	
			else :
	
				dico_aa[chromosome] += dna_to_protein(lig.upper())    #J'applique la fonction pour directement récupérer des séquences protéiques et non en nt 
	return dico_aa

#Dico des génomes en aa : 
dico_aa_cuniculi = dico_aa('data_training_sur_e_cuniculi.fa')
dico_aa_algerae = dico_aa('data_training_A_algerae')
dico_aa_bieneusi = dico_aa('data_training_E_bieneusi')
dico_aa_ceranae = dico_aa('data_training_Nosmea_ceranae')
dico_aa_parisii = dico_aa('data_training_parissi.fr')


##POUR CUNICULI 

prot_cuniculi= open("data_training_prot_cuniculi", 'w')
for id, seq in dico_aa_cuniculi.items():
	print(f">{id}", file = prot_cuniculi)		
	for i in range(0, len(seq), 60):		
		print(seq[i:i+60], file = prot_cuniculi)



##POUR ALGERAE	
prot_algerae= open("data_training_prot_algerae", 'w')

for id, seq in dico_aa_algerae.items():
	print(f">{id}", file = prot_algerae)		
	for i in range(0, len(seq), 60):		
		print(seq[i:i+60], file = prot_algerae)

		

		
#POUR BIENEUSI		
prot_bieneusi= open("data_training_prot_bieneusi", 'w')

for id, seq in dico_aa_bieneusi.items():
	print(f">{id}", file = prot_bieneusi)		
	for i in range(0, len(seq), 60):		
		print(seq[i:i+60], file = prot_bieneusi)		
		
#POUR CERANAE	
prot_ceranae= open("data_training_prot_ceranae", 'w')

for id, seq in dico_aa_ceranae.items():
	print(f">{id}", file = prot_ceranae)		
	for i in range(0, len(seq), 60):		
		print(seq[i:i+60], file = prot_ceranae)



#POUR PARISII	
prot_parisii = open("data_training_prot_parisii", 'w')

for id, seq in dico_aa_parisii.items():
	print(f">{id}", file = prot_parisii)		
	for i in range(0, len(seq), 60):		
		print(seq[i:i+60], file = prot_parisii)



