from repDNA.nac import RevcKmer, Kmer   
from repDNA.ac import DAC, TAC
from repDNA.psenac import PseDNC, PseKNC
import json, time


# The following 2 functions extract Nucleic Acid Composition features
def ReverseComplimentKmer(gene2seq, k):
    X = dict()
    succeed_cnt = 0
    rev_kmer = RevcKmer(k, normalize=True, upto=True)
    for gene, seq in gene2seq.items():
        seq = [seq]
        try:
            pos_vec = rev_kmer.make_revckmer_vec(seq)
            X[gene] = pos_vec
            succeed_cnt += 1
            if succeed_cnt % 100 == 0:
                print(succeed_cnt)
        except Exception as e:
            continue
    print('Reverse Compliment Kmer, succeed for %d gene'%(succeed_cnt))
    with open('%d-revcKmer-features.json'%(k), 'w') as output_f:
        output_f.write(json.dumps(X))
    return X

 
def BasicKmer(gene2seq, k):
    X = dict()
    succeed_cnt = 0
    kmer = Kmer(k)
    for gene, seq in gene2seq.items():
        seq = [seq]
        try:
            vec = kmer.make_kmer_vec(seq)
            X[gene] = vec
            succeed_cnt += 1
            if succeed_cnt % 100 == 0:
                print(succeed_cnt)
        except Exception as e:
            continue
    print('Basic Kmer, succeed for %d gene'%(succeed_cnt))
    with open('%d-basicKmer-features.json'%(k), 'w') as output_f:
        output_f.write(json.dumps(X))
    return X




# The following 2 functions extract Autocorrelation features
def DinucleotideAutoCovariance(gene2seq, k):
    X = dict()
    succeed_cnt = 0
    dac = DAC(k)
    for gene, seq in gene2seq.items():
        seq = [seq]
        try:
            vec = dac.make_dac_vec(seq, all_property=True)
            X[gene] = vec
            succeed_cnt += 1
            if succeed_cnt % 100 == 0:
                print(succeed_cnt)
        except Exception as e:
            continue
    print('DAC, succeed for %d gene'%(succeed_cnt))
    with open('%d-DAC-features.json'%(k), 'w') as output_f:
        output_f.write(json.dumps(X))
    return X


def TrinucleotideAutoCovariance(gene2seq, k):
    X = dict()
    succeed_cnt = 0
    tac = TAC(k)
    for gene, seq in gene2seq.items():
        seq = [seq]
        try:
            vec = tac.make_tac_vec(seq, all_property=True)
            X[gene] = vec
            succeed_cnt += 1
            if succeed_cnt % 100 == 0:
                print(succeed_cnt)
        except Exception as e:
            continue
    print('TAC, succeed for %d gene'%(succeed_cnt))
    with open('%d-TAC-features.json'%(k), 'w') as output_f:
        output_f.write(json.dumps(X))
    return X





# The following 2 functions extract Pseudo Nucleic Acid Composition features
def PseudoDinucleotideComposition(gene2seq):
    X = dict()
    succeed_cnt = 0
    # let's use the default value
    psednc = PseDNC()
    for gene, seq in gene2seq.items():
        seq = [seq]
        try:
            vec = psednc.make_psednc_vec(seq)
            X[gene] = vec
            succeed_cnt += 1
            if succeed_cnt % 100 == 0:
                print(succeed_cnt)
        except Exception as e:
            continue
    print('PseDNC, succeed for %d gene'%(succeed_cnt))
    with open('PseDNC-features.json', 'w') as output_f:
        output_f.write(json.dumps(X))
    return X


def PseudoKTuplerComposition(gene2seq):
    X = dict()
    succeed_cnt = 0
    # let's use the default value
    pseknc = PseKNC()
    for gene, seq in gene2seq.items():
        seq = [seq]
        try:
            vec = pseknc.make_pseknc_vec(seq)
            X[gene] = vec
            succeed_cnt += 1
            if succeed_cnt % 100 == 0:
                print(succeed_cnt)
        except Exception as e:
            continue
    print('PseKNC, succeed for %d gene'%(succeed_cnt))
    with open('PseKNC-features.json', 'w') as output_f:
        output_f.write(json.dumps(X))
    return X



'''
# Given a value k and a sequence list, this function goes through 6 feature-extracting functions
def extract_features(gene2seq, k):
    BasicKmer(gene2seq, k)
'''

def main():
    gene2seq = dict()
    with open('gene2seq.txt', 'r') as input_f:
        for line in input_f:
            split_list = line.split('\n')[0].split('\t')
            gene = split_list[0]
            seq = split_list[1]
            gene2seq[gene] = seq

    print('K=4, RevcKmer start')
    start = time.clock() 
    X1 = ReverseComplimentKmer(gene2seq, 4) 
    end = time.clock()
    print('K=4, RevcKmer finished in %.3f minutes'%((end-start)/60.0))
    
    print('K=5, RevcKmer start')
    start = time.clock() 
    X2 = ReverseComplimentKmer(gene2seq, 5) 
    end = time.clock()
    print('K=5, RevcKmer finished in %.3f minutes'%((end-start)/60.0))


    print('K=4, DAC start')
    start = time.clock() 
    X3 = DinucleotideAutoCovariance(gene2seq, 4)
    end = time.clock()
    print('K=4, DAC finished in %.3f minutes'%((end-start)/60.0))

    print('K=5, DAC start')
    start = time.clock() 
    X4 = DinucleotideAutoCovariance(gene2seq, 5)    
    end = time.clock()
    print('K=5, DAC finished in %.3f minutes'%((end-start)/60.0))

    
    print('PseDNC starts')
    start = time.clock()
    X5 = PseudoDinucleotideComposition(gene2seq)
    end = time.clock()
    print('PseDNC finished in %.3f minutes'%((end-start)/60.0))

    print('PseKNC starts')
    start = time.clock()
    X6 = PseudoKTuplerComposition(gene2seq)
    end = time.clock()
    print('PseKNC finished in %.3f minutes'%((end-start)/60.0))






if __name__ == '__main__':
    main()
