newline_cnt = 0

with open('gene2seq.txt', 'w') as output_f:
    with open('gene2seq_query_logs.txt', 'r') as input_f:
        for line in input_f:
            split_list = line.split('\n')[0].split('\t')
            seq = split_list[1]
            gene = split_list[0]

            validate_list = [elem in ['A', 'G', 'C', 'T'] for elem in list(seq)]
            if seq != ' ' and len(seq) >= 3000 and all(validate_list):
                new_line = gene + '\t' + seq + '\n'
                output_f.write(new_line)
                newline_cnt += 1

print(newline_cnt)
