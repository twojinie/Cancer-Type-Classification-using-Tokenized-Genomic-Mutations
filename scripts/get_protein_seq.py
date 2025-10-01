import requests
import pandas as pd
from tqdm import tqdm 
import pickle

def get_uniprot_id(gene_symbol):
    # Homo sapiens 필터 추가
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_symbol}+AND+organism_id:9606&format=json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if data['results']:
            uniprot_id = data['results'][0]['primaryAccession']
            return uniprot_id
        else:
            return None
    else:
        response.raise_for_status()

def get_protein_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text
        sequence = ''.join(fasta_data.split('\n')[1:])
        return sequence
    else:
        response.raise_for_status()

def get_protein_sequence_by_gene_symbol(gene_symbol):
    uniprot_id = get_uniprot_id(gene_symbol)
    if uniprot_id:
        return get_protein_sequence(uniprot_id)
    else:
        return None
    

train = pd.read_csv("../data/train.csv")
test = pd.read_csv("../data/test.csv")
prot_list = list(train.drop(columns=['SUBCLASS', 'ID']).columns)

print("Protein Sequence Download Start")

prot_dict = {} # {gene_symbol : protein_sequence}
i = 0
for prot in tqdm(prot_list):
    seq = get_protein_sequence_by_gene_symbol(prot)
    prot_dict[prot] = seq
    i += 1
    if i % 500 == 0:
        with open("./protein_sequence_database/gene_seq_dict.pkl", "wb") as f:
            pickle.dump(prot_dict, f)
        print(f'save at {i}th data')
        
with open("./protein_sequence_database/gene_seq_dict.pkl", "wb") as f:
    pickle.dump(prot_dict, f)
    
print("Download Completed")