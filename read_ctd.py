import csv
import json
from Bio import Entrez

gene2diseases = dict()

Entrez.email = "zhangmx1997@foxmail.com"
cnt = 0
with open('CTD_genes_diseases.csv', 'r') as input_f:
    lines = csv.reader(input_f)
    for line in lines:
        #print(line)
        if len(line) == 9 and 'GeneID' not in line:
            gene = line[1]
            disease = (line[2], line[3])
            if gene not in gene2diseases:
                gene2diseases[gene] = list()
            gene2diseases[gene].append(disease)
            '''
            cnt += 1
            if cnt <= 10:
                handle = Entrez.efetch(db="gene", id=gene, rettype="gb", retmode="text")
                print(handle.read())
            else:
                break
            '''



with open('gene2diseases.json', 'w') as output_f:
    output_f.write(json.dumps(gene2diseases))

print(len(gene2diseases))

'''
for gene, diseases in gene2diseases.items():
    handle = Entrez.efetch(db="gene", id=gene, rettype="gb", retmode="text")
    print(handle.read())
'''
