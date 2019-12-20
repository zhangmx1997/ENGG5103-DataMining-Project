import numpy as np
import json, os
import matplotlib.pyplot as plt

gene2disease_mat = np.loadtxt('gene2disease_matrix.txt', dtype='i', delimiter=' ')
with open('gene2line.json', 'r') as input_f:
    gene2line = json.loads(input_f.read())

prediction_dir = 'prediction'
prediction_dir = 'tmp'

predict_files = os.listdir(prediction_dir)
for f in predict_files:
    f = os.path.join(prediction_dir, f)
    with open(f, 'r') as input_f:
        gene2predict = json.loads(input_f.read())

    precision_list = list()
    recall_list = list()
    #printed_cnt = 0

    for gene, predict in gene2predict.items():
        line = gene2line[gene]

        correct_num = 0
        expected_num = sum(gene2disease_mat[line,:])

        if type(predict[0]) == list:
            total_num = len(predict)
            for p in predict:
                disease_id = p[0]
                similarity_score = p[1]
                if gene2disease_mat[line, disease_id] == 1:
                    correct_num += 1
        else:
            total_num = 1
            disease_id = predict[0]
            similarity_score = predict[1]
            if gene2disease_mat[line, disease_id] == 1:
                correct_num += 1


        precision = float(correct_num) / float(total_num)
        recall = float(correct_num) / float(expected_num)
        precision_list.append(precision)
        recall_list.append(recall)
        
        '''
        if 'False' in f and printed_cnt <= 5:
            print(recall, expected_num, correct_num, total_num)
            printed_cnt += 1 # = True
        '''

    '''
    if '4-DAC-' in f and '-False' in f:
        x = range(len(precision_list))
        fig = plt.figure()
        ax = plt.subplot()
        ax.scatter(x, precision_list, s=0.1, alpha=0.6)
        plt.show()
        print(precision_list)
    '''

    feature_name = f.split('/')[-1].split('-gene2predict.json')[0]
    precision = sum(precision_list) / len(precision_list)
    recall = sum(recall_list) / len(recall_list)

    print('%s\tprecision: %.3f\trecall: %.3f'%(feature_name, precision, recall))
