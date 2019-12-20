import json, os

feature_dir = '/home/mxzhang/Desktop/data_mining/features'
files = os.listdir(feature_dir)
files = [f for f in files if f.endswith('-features.json')]


for f in files:
    f = os.path.join(feature_dir, f)
    with open(f, 'r') as input_f:
        gene2vec = json.loads(input_f.read())
        
    with open('geneIDs.txt', 'w') as output_f:
        for gene, vec in sorted(gene2vec.items(), key=lambda x:x[0], reverse=False):
            line = gene+'\n'
            output_f.write(line)
    break
