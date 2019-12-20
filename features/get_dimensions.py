import os, json


feature_dir = '.'
files = os.listdir(feature_dir)
files = [f for f in files if f.endswith('.json')]

for f in files:
    print(f)
    with open(f, 'r') as input_f:
        gene2feature = json.loads(input_f.read())
    for g, f in gene2feature.items():
        print(len(f[0]))
        break
