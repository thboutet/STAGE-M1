"""
#POUR CUNICULI
gene_id = 0 
glimmer = open("glimmer.gff", "w")
with open ('result_E.cuniculi.predict', 'r') as f1 :
	for lig in f1 :
		lig =lig.rstrip()
		
		if lig.startswith(">") : 
			chromosome=lig.strip()[1:]
		else :
			lig=lig.split()	
			gene_id += 1
			start     = int(lig[1])
			end       = int(lig[2])
			
			
			
			if start > end :
				start = lig[2]
				end = lig[1]
				brin = "-"
			else :
				
				brin ="+"
				
				
			
			print(f"{chromosome}\tGlimmer\tgene\t{start}\t{end}\t.\t{brin}\t.\"tID={gene_id}", file = glimmer)
		
			
##POUR ALGERAE		
gene_id_algerae = 0 
glimmer_algerae = open("glimmer_algerae.gff", "w")
with open ('glimmer/result_A.algerae.predict', 'r') as f1 :
	for lig in f1 :
		lig =lig.rstrip()
		
		if lig.startswith(">") : 
			chromosome=lig.strip()[1:]
		else :
			lig=lig.split()	
			gene_id_algerae += 1
			start     = int(lig[1])
			end       = int(lig[2])
			
			
			
			if start > end :
				start = lig[2]
				end = lig[1]
				brin = "-"
			else :
				
				brin ="+"
				
				
			
			print(f"{chromosome}\tGlimmer\tgene\t{start}\t{end}\t.\t{brin}\t.\tID={gene_id_algerae}", file = glimmer_algerae)
		

##POUR BIENEUSI		
gene_id_bieneusi = 0 
glimmer_bieneusi = open("glimmer_bieneusi.gff", "w")
with open ('glimmer/result_E.bieneusi.predict', 'r') as f1 :
	for lig in f1 :
		lig =lig.rstrip()
		
		if lig.startswith(">") : 
			chromosome=lig.strip()[1:]
		else :
			lig=lig.split()	
			gene_id_bieneusi += 1
			start     = int(lig[1])
			end       = int(lig[2])
			
			
			
			if start > end :
				start = lig[2]
				end = lig[1]
				brin = "-"
			else :
				
				brin ="+"
				
				
			
			print(f"{chromosome}\tGlimmer\tgene\t{start}\t{end}\t.\t{brin}\t.\tID={gene_id_bieneusi}", file = glimmer_bieneusi)
						
##POUR CERANAE		
gene_id_ceranae = 0 
glimmer_ceranae = open("glimmer_ceranae.gff", "w")
with open ('glimmer/result_N.ceranae.predict', 'r') as f1 :
	for lig in f1 :
		lig =lig.rstrip()
		
		if lig.startswith(">") : 
			chromosome=lig.strip()[1:]
		else :
			lig=lig.split()	
			gene_id_ceranae += 1
			start     = int(lig[1])
			end       = int(lig[2])
			
			
			
			if start > end :
				start = lig[2]
				end = lig[1]
				brin = "-"
			else :
				
				brin ="+"	
				
			print(f"{chromosome}\tGlimmer\tgene\t{start}\t{end}\t.\t{brin}\t.\tID={gene_id_ceranae}", file = glimmer_ceranae)
"""		
##POUR PARISII		
gene_id_parisii = 0 
glimmer_parisii = open("glimmer_parisii.gff", "w")
with open ('glimmer/result_N.parisii.predict', 'r') as f1 :
	for lig in f1 :
		lig =lig.rstrip()
		
		if lig.startswith(">") : 
			chromosome=lig.strip()[1:]
		else :
			lig=lig.split()	
			gene_id_parisii += 1
			start     = int(lig[1])
			end       = int(lig[2])
			
			
			
			if start > end :
				start = lig[2]
				end = lig[1]
				brin = "-"
			else :
				
				brin ="+"	
				
			print(f"{chromosome}\tGlimmer\tgene\t{start}\t{end}\t.\t{brin}\t.\tID={gene_id_parisii}", file = glimmer_parisii)
					
			
				
					
