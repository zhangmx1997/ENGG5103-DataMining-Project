import json, os, time
from scipy import spatial
import numpy as np

feature_dir = '/home/mxzhang/Desktop/data_mining/features'
files = os.listdir(feature_dir)
files = [f for f in files if f.endswith('-features.json')]



def cosine_similarity(x, y, norm=False):
    if sum(x) == 0 or sum(y) == 0:
        if x == y:
            return 1.0
        else:
            return 0.0
    else:
        res = np.array([[x[i] * y[i], x[i] * x[i], y[i] * y[i]] for i in range(len(x))])
        cos = sum(res[:, 0]) / (np.sqrt(sum(res[:, 1])) * np.sqrt(sum(res[:, 2])))
        if norm:
            return 0.5 * cos + 0.5
        else:
            return cos



def main():
    for f in files:
        f = os.path.join(feature_dir, f)

        print(f)
        start = time.clock()
        with open(f, 'r') as input_f:
            gene2vec = json.loads(input_f.read())
            
        num_gene = len(gene2vec)
        similarity_matrix = np.zeros((num_gene, num_gene))
        #print(num_gene)

        gene_i = -1
        for gene, vec in sorted(gene2vec.items(), key=lambda x:x[0], reverse=False):
            gene_i += 1
            vec1 = vec[0]
            
            gene_j = -1
            for gene_, vec_ in sorted(gene2vec.items(), key=lambda x:x[0], reverse=False):
                gene_j += 1
                if gene == gene_:
                    continue
                vec2 = vec_[0]
                #cosine_similarity = 1 - spatial.distance.cosine(vec1, vec2)
                #cos = cosine_similarity(vec1, vec2, True)
                cos = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
                similarity_matrix[gene_i, gene_j] = cos
        end = time.clock()
        print('finished in %.3f'%((end-start)/60.0))

        feature_name = f.split('/')[-1].split('-features.json')[0]
        output_file = feature_name + '-similarity.txt'
        np.savetxt(output_file, similarity_matrix, fmt='%.5f')
        print(output_file)

if __name__ == '__main__':
    main()
