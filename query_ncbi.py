from Bio import Entrez, SeqIO
import json


Entrez.email = "mxzhang@cse.cuhk.edu.hk"

with open('gene2diseases.json', 'r') as input_f:
    gene2diseases = json.loads(input_f.read())

print(len(gene2diseases))


#gene2seq = dict()
cnt = 0

processed_set = set()
with open('0-gene2seq.txt', 'r') as input_f:
    #output_f.write(json.dumps(gene2seq))
    for line in input_f:
        gene = line.split('\t')[0]
        processed_set.add(gene)

#with open('0-gene2seq.txt', 'w') as output_f:
for gene, diseases in sorted(gene2diseases.items(), key=lambda x:x[0], reverse=False):
    #print(gene)
    cnt += 1
    #if cnt % 100 == 0:
    #    print(cnt)
    #if cnt >= 15200:
    #    break
    if str(gene) in processed_set:
        continue
    try:
        handle = Entrez.efetch(db="nucleotide", id=gene, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        #print(record.seq)
        #gene2seq[gene] = record.seq
        result_str = gene + '\t' + record.seq# + '\n'
        print(result_str)
        #output_f.write(result_str)
    except Exception as e:
        result_str = gene + '\t' + ' ' #' \n'
        #output_f.write(result_str)
        print(result_str)
        continue
'''
handle = Entrez.efetch(db="nucleotide", id="100174880", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
print(type(record), record)
'''

'''
with open('0-gene2seq.json', 'w') as output_f:
    output_f.write(json.dumps(gene2seq))
print(len(gene2seq))
'''
