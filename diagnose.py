import json, os
import numpy as np
import matplotlib.pyplot as plt

def get_gene_to_line_dict(f):
    gene_cnt = 0
    gene2line = dict()
    with open(f, 'r') as input_f:
        for line in input_f:
            gene = line.split('\n')[0]
            gene2line[gene] = gene_cnt 
            gene_cnt += 1
    print('%d gene'%(len(gene2line)))
    return gene2line


def get_disease_to_col_dict(f, gene2line):
    disease_cnt = 0
    disease2col = dict()
    disease_set = set()
    with open(f, 'r') as input_f:
        gene2diseases = json.loads(input_f.read())

    valid_gene2diseases = dict()
    for gene, diseases in gene2diseases.items():
        if gene not in gene2line:
            continue
        if gene not in valid_gene2diseases:
            valid_gene2diseases[gene] = list()
        for d in diseases:
            try:
                disease_set.add(d[0])
                valid_gene2diseases[gene].append(d)
            except Exception as e:
                print(e, d, gene)
                continue
    print('%d diseases'%(len(disease_set)))
    
    for d in disease_set:
        disease2col[d] = disease_cnt
        disease_cnt += 1
    return disease2col, valid_gene2diseases


def construct_gene_disease_matrix(gene2line, disease2col, gene2diseases):
    gene_cnt = len(gene2line)
    disease_cnt = len(disease2col)
    gene2disease_mat = np.zeros((gene_cnt, disease_cnt))

    for gene, diseases in gene2diseases.items():
        line = gene2line[gene]
        for d in diseases:
            d = d[0] # disease name
            col = disease2col[d]
            gene2disease_mat[line, col] = 1
    
    output_file = 'gene2disease_matrix.txt'
    np.savetxt(output_file, gene2disease_mat, fmt='%d')
    print(output_file)
    return gene2disease_mat


# Select top-k sequence for each query, using the similarity in file f
def query(f, k, gene2disease, line2gene, percentage, findTop=False):
    gene2predict = dict()
    # suppose there are M genes and N diseases
    # gene2disease is a M*N 0-1 matrix
    # similarity is a M*M matrix, and the elements on the diagnoal are all 0.00
    gene2gene_sim = np.loadtxt(f, delimiter=' ')
    line, col = gene2gene_sim.shape

    num_disease = len(gene2disease[0])
    for i in range(line):
        weighted_sum = np.zeros((1, num_disease))
        sim_vec = gene2gene_sim[i,:]
        sorted_sim_vec = sorted(enumerate(sim_vec), key=lambda x:x[1], reverse=True)
        top_sim_vec = sorted_sim_vec[:k]
        top_gene_ids = [t[0] for t in top_sim_vec]
        top_gene_sims = [t[1] for t in top_sim_vec]

        for j, sim in enumerate(top_gene_sims):
            #print(np.sum(gene2disease[top_gene_ids[i],:]))
            weighted_sum += sim * gene2disease[top_gene_ids[j],:]
        #print(sim, np.sum(weighted_sum))
        weighted_sum = weighted_sum[0].tolist()
        #print(weighted_sum)
        if findTop:
            predict = sorted(enumerate(weighted_sum), key=lambda x:x[1], reverse=True)[0] # (diseaseID, weighted_sum)
            #predict = (predict[0], predict[1].tolist())
            #print(type(predict), predict)
        else:
            predict = sorted(enumerate(weighted_sum), key=lambda x:x[1], reverse=True)
            top_score = predict[0][1]
            threshold = percentage * top_score# - 0.3
            predict = [p for p in predict if p[1]>=threshold]
            #predict = [(p[0], p[1].tolist()) for p in predict]
        gene2predict[line2gene[i]] = predict
    return gene2predict

    

def main():

    if not os.path.isfile('gene2line.json'):
        input_file = 'geneIDs.txt'
        gene2line = get_gene_to_line_dict(input_file)

        input_file = 'gene2diseases.json'
        disease2col, valid_gene2diseases = get_disease_to_col_dict(input_file, gene2line)
        
        output_file = 'gene2line.json'
        with open(output_file, 'w') as output_f:
            output_f.write(json.dumps(gene2line))
        print(output_file)
            
        output_file = 'disease2col.json'
        with open(output_file, 'w') as output_f:
            output_f.write(json.dumps(disease2col))
        print(output_file)

        gene2disease_mat = construct_gene_disease_matrix(gene2line, disease2col, valid_gene2diseases) # np matrix

        output_file = 'valid_gene2disease.json'
        with open(output_file, 'w') as output_f:
            output_f.write(json.dumps(valid_gene2diseases))
        print(output_file)

    else:
        with open('gene2line.json', 'r') as input_f:
            gene2line = json.loads(input_f.read())
        with open('disease2col.json', 'r') as input_f:
            disease2col = json.loads(input_f.read())
        gene2disease_mat = np.loadtxt('gene2disease_matrix.txt', dtype='i', delimiter=' ')

        line2gene = dict()
        for gene, line in gene2line.items():
            line2gene[line] = gene

        sim_dir = 'similarity_matrix'
        sim_files = os.listdir(sim_dir)
        for sim_file in sim_files:  
            sim_file = os.path.join(sim_dir, sim_file)
            for k in [50]: #, 60, 70, 80, 90]:
            #for k in [100, 300, 500, 700, 900, 1100]:
                #for findTop in [True, False]:
                for findTop in [False]: #[True, False]:
                    precision_list = list()
                    recall_list = list()
                    #x = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
                    #x = [0.8]
                    x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
                    for percentage in x:
                        gene2predict = query(sim_file, k, gene2disease_mat, line2gene, percentage, findTop)
                        '''
                        sim_name = sim_file.split('/')[-1].split('-similarity.txt')[0]
                        output_file = 'threshold-%s-%d-gene2predict.json'%(sim_name, percentage*10)
                        with open(output_file, 'w') as output_f:
                            output_f.write(json.dumps(gene2predict))
                        print(output_file)
                        #print(percentage)
                        '''
                        cnts = list()
                        for g, p in gene2predict.items():
                            cnts.append(len(p))
                        print(percentage, sum(cnts)/len(cnts))
                        
                        ''' 
                        current_precision_list = list()
                        current_recall_list = list()
                        for gene, predicts in gene2predict.items():
                            total_num = len(predicts)
                            expected_num = sum(gene2disease_mat[gene2line[gene],:])
                            correct_num = 0
                            for p in predicts:
                                if gene2disease_mat[gene2line[gene], p[0]] == 1:
                                    correct_num += 1
                            precision = float(correct_num) / float(total_num)
                            recall = float(correct_num) / float(expected_num)
                            current_precision_list.append(precision)
                            current_recall_list.append(recall)
                        precision_list.append(sum(current_precision_list)/len(current_precision_list))
                        recall_list.append(sum(current_recall_list)/len(current_recall_list))
                        '''
                        
                    '''
                    plt.figure() 
                    plt.plot(x, precision_list, marker='x', lw=2, alpha=0.5, label='Accuracy')
                    plt.plot(x, recall_list, marker='o', lw=2, alpha=0.5, label='Recall')

                    print(percentage)
                    print(precision_list)
                    print(recall_list)
                    print('\n\n')
                    plt.legend()
                    plt.xlabel('Threshold')
                    plt.ylabel('Percentage')
                    sim_name = sim_file.split('/')[-1].split('-similarity.txt')[0]
                    fig_name = sim_name + '.png'
                    plt.savefig(fig_name)
                    print(fig_name)
                    #plt.show()
                    '''
                    




if __name__ == '__main__':
    main()
