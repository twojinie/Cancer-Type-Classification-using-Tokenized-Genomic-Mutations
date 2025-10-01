import requests
import pandas as pd
from tqdm import tqdm 
import pickle
import re

def drop_dup_mut_at_same_pos(mutations):
    '''
    같은 포지션에 겹치는 변이들이 여러개 보이는 경우,
    가장 impact가 클 것으로 예상되는 순서로, 하나만 남기기
    * > fs > 다른 성질 작용기 amino acid > 같은 성질 작용기 amino acid > silent 
    '''
    pos_mut_dict = {} # {position : [Original_AA, [mutation list]]}
    # ex. E512K, E512E -> {'512' : [E, [K, E]] 
    
    not_canonical_muts = []
    
    for mutation in mutations:
        if ('>' in mutation) or ('fs' in mutation and '*' in mutation) or ('-' in mutation):
            not_canonical_muts.append(mutation)
        else:
            match = re.match(r'([A-Z\*]*)(\d+)([a-zA-Z\*]*)', mutation)
            if match:
                pos_str = match.groups()[1]
                if pos_str not in pos_mut_dict.keys():
                    pos_mut_dict[pos_str] = [match.groups()[0], [match.groups()[2]]]
                else:
                    if match.groups()[0] == pos_mut_dict[pos_str][0]:
                        pos_mut_dict[pos_str][1].append(match.groups()[2])
                    else: # 같은 site에 다른 original sequence를 가지는 mutation 정보들이 주어진 경우. 
                        # ex. 'FFSGR3341fs', 'FFLV3341fs' -> '-3341fs'
                        pos_mut_dict[pos_str][0] = '-'
                        pos_mut_dict[pos_str][1] = [match.groups()[2]]
    
    muts_trimmed = []
    for pos in pos_mut_dict.keys():
        aa = pos_mut_dict[pos][0]
        muts = pos_mut_dict[pos][1]
        if len(muts) > 1:
            muts = [get_most_impactful_mut(aa, muts)]
        mut_trimmed = aa + pos + muts[0]
        muts_trimmed.append(mut_trimmed)
        
    return muts_trimmed + not_canonical_muts

def get_most_impactful_mut(aa, muts):
    
    aa_property_encodings = {}
    
    def encoding_aa(aa_property_encodings, aas, val):
        for aa in aas:
            aa_property_encodings[aa] = val
        return aa_property_encodings
    
    nonpolar = ['G', 'A', 'V', 'L', 'I', 'F', 'M', 'W']
    proline = ['P'] # nonpolar
    polar = ['S', 'T', 'N', 'Q', 'Y']
    cysteine = ['C']
    acidic = ['D', 'E']
    basic = ['K', 'R', 'H']
    
    aa_kinds = [nonpolar, proline, polar, cysteine, acidic, basic]
    i = 0
    for aas in aa_kinds:
        aa_property_encodings = encoding_aa(aa_property_encodings, aas, i)
        i += 1
    
    if '*' in muts:
        return '*'
    elif 'fs' in muts:
        return 'fs'
    else:
        if aa in muts:
            muts.remove(aa)
        if len(muts) == 1:
            return muts[0]
        else:
            ## 여러 point mutation 중에 선택
         
            # 두 변이의 aa가 같은 group인 경우, 첫번째 변이 선택.
            if aa_property_encodings[muts[0]] == aa_property_encodings[muts[1]]:
                return muts[0]
            else: # 두 변이가 다른 group인 경우, 우선순위가 original의 group에 따라 달라짐.
                if aa_property_encodings[aa] == 0 or aa_property_encodings[aa] == 1: # Pro 혹은 nonpolar
                    if aa_property_encodings[aa] == aa_property_encodings[muts[0]]:
                        return muts[1]
                    elif aa_property_encodings[aa] == aa_property_encodings[muts[1]]:
                        return muts[0]
                    elif aa_property_encodings[muts[0]] in [0, 1]:
                        return muts[1]
                    else:
                        return muts[0]
                
                if aa_property_encodings[aa] == 2 or aa_property_encodings[aa] == 3: # Cys 혹은 polar
                    if aa_property_encodings[aa] == aa_property_encodings[muts[0]]:
                        return muts[1]
                    elif aa_property_encodings[aa] == aa_property_encodings[muts[1]]:
                        return muts[0]
                    elif aa_property_encodings[muts[0]] in [0, 1]:
                        return muts[1]
                    else:
                        return muts[0]
                    
                if aa_property_encodings[aa] == 4: # acidic
                    if aa_property_encodings[muts[0]] == 5: # basic으로 바뀌는 경우가 가장 영향 클 것
                        return muts[0]
                    elif aa_property_encodings[muts[1]] == 5:
                        return muts[1]
                    elif aa_property_encodings[muts[0]] in [0, 1]:
                        return muts[0]
                    elif aa_property_encodings[muts[1]] in [0, 1]:
                        return muts[1]
                    elif aa_property_encodings[muts[0]] in [2, 3]:
                        return muts[0]
                    elif aa_property_encodings[muts[1]] in [2, 3]:
                        return muts[1]
                    else:
                        return muts[0]
                
                if aa_property_encodings[aa] == 5: # basic
                    if aa_property_encodings[muts[0]] == 4: # acidic으로 바뀌는 경우가 가장 영향 클 것 
                        return muts[0]
                    elif aa_property_encodings[muts[1]] == 4:
                        return muts[1]
                    elif aa_property_encodings[muts[0]] in [0, 1]:
                        return muts[0]
                    elif aa_property_encodings[muts[1]] in [0, 1]:
                        return muts[1]
                    elif aa_property_encodings[muts[0]] in [2, 3]:
                        return muts[0]
                    elif aa_property_encodings[muts[1]] in [2, 3]:
                        return muts[1]
                    else:
                        return muts[0]
                    
def get_wt_seq(gene_symbol, prot_seq_dict):
    '''
    gene_symbol : 유전자 이름 (ex. ABCC2)
    prot_seq_dict : {Protein : Amino Acid Sequence}
    '''
    return prot_seq_dict[gene_symbol]

def get_mut_seq(sample_number, gene_symbol, sequence, mutations, pos_err_list, len_err_list, no_seq_err_list):
    '''
    기존 sequence와 mutation 정보가 주어졌을 때, 
    mutation이 일어난 sequence를 return하는 함수 
    '''
    
    # Catch No Sequnece Case 
    # print(mutations)
    try:
        sequence = list(sequence)
    except TypeError:
        # print("No sequence found")
        no_seq_err_list.append([sample_number, gene_symbol, mutations])
        return None
    
    mutations = list(set(mutations.split()))
    
    mutations = drop_dup_mut_at_same_pos(mutations)
    
    fs_mis_del_site = None
    
    for mutation in mutations:
        
        # deletion 
        if 'del' in mutation: 
            match = re.match(r'([a-zA-Z]*)(\d+)del', mutation)
            pos = int(match.groups()[1]) - 1
            
            if fs_mis_del_site is not None:
                if pos <= fs_mis_del_site:
                    fs_mis_del_site = pos
                else:
                    return ''.join(sequence)
            else:
                fs_mis_del_site = pos
            
            if pos < len(sequence):
                if match.groups()[0] != '':
                    if sequence[pos] == match.groups()[0]:
                        sequence = sequence[:pos] + sequence[pos+1:]
                        # return ''.join(sequence)
                    else:
                        pos_err_list.append([sample_number, gene_symbol, mutation])
                        return None
                else:
                    sequence = sequence[:pos] + sequence[pos+1:]
                    # return ''.join(sequence)
            else:
                # print(f"lenght error: Expected >= {pos + 1}, got {len(sequence)}")
                len_err_list.append([sample_number, gene_symbol, mutation])
                return None
            
        # frame shift
        elif 'fs' in mutation: 
            match = re.match(r'([\-]*[a-zA-Z]*[\*]*)(\d+)fs', mutation)
            pos = int(match.groups()[1]) - 1
            
            
            if fs_mis_del_site is not None:
                if pos < fs_mis_del_site:
                    fs_mis_del_site = pos
                else:
                    return ''.join(sequence)
            else:
                fs_mis_del_site = pos
            
            if pos < len(sequence):
                if '-' not in match.groups()[0] and '*' not in match.groups()[0]:
                    # all alphabet case
                    if sequence[pos:pos + len(match.groups()[0])] == list(match.groups()[0]):
                        sequence = sequence[:pos]
                        # return ''.join(sequence)
                    else:
                        pos_err_list.append([sample_number, gene_symbol, mutation])
                        return None
                else:
                    sequence = sequence[:pos]
                    # return ''.join(sequence)
                    
            else:
                len_err_list.append([sample_number, gene_symbol, mutation])
                return None
        
        # *를 그냥 amino acid로 보고 변환한 뒤, *가 포함된 서열을 자르는 방식으로 변환
        # Missense mutation을 따로 코딩하지 않음  
        
        # '>' case 
        elif '>' in mutation: # ex. 1998_1999KR>NW
            chunk, new_seq = mutation.split('>')
            start_pos, chunk = chunk.split('_')
            orig_seq = chunk[-2:]
            pos = int(start_pos) - 1
            if pos < len(sequence):
                if sequence[pos:pos+2] == list(orig_seq):
                    sequence = sequence[:pos] + list(new_seq) + sequence[pos+2:]
                else:
                    pos_err_list.append([sample_number, gene_symbol, mutation])
                    return None
            else:
                len_err_list.append([sample_number, gene_symbol, mutation])
                return None
            
        # Normal mutation case, ex. H336D, K512*
        else:
            match = re.match(r'([A-Z\*-]*)(\d+)([A-Z\*]*)', mutation)
            if match == None:
                print('case not matched')
                return None
            pos = int(match.groups()[1]) - 1
            
            # Remove silent mutation
            if match.groups()[0] == match.groups()[2]:
                continue
                
            if pos < len(sequence):
                if sequence[pos] == match.groups()[0] or match.groups()[0] == '*':
                    sequence = sequence[:pos] + list(match.groups()[2]) + sequence[pos+1:]
                else:
                    pos_err_list.append([sample_number, gene_symbol, mutation])
                    return None
                    
            else:
                len_err_list.append([sample_number, gene_symbol, mutation])
                return None

    # Missense 처리
    if '*' in sequence:
        sequence = sequence[:sequence.index('*')]
    
            
    return ''.join(sequence)

def get_uniprot_id_gene(gene_symbol):
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

def get_num_err_genes(e_gene, mut_list, prot_seq_dict, uniprot_id=None):
    len_err_list = []
    pos_err_list = []
    no_seq_err_list = []
    num_err = 0
    
    if uniprot_id == None:
        for mut in mut_list:
            if mut[1] == e_gene:
                seq = get_wt_seq(mut[1], prot_seq_dict)
                mutated_seq = get_mut_seq(mut[0], mut[1], seq, mut[2], pos_err_list, len_err_list, no_seq_err_list)
                if mutated_seq == None:
                    num_err += 1
    else:
        try:
            seq = get_protein_sequence(uniprot_id)
            for mut in mut_list:
                if mut[1] == e_gene:
                    mutated_seq = get_mut_seq(mut[0], mut[1], seq, mut[2], pos_err_list, len_err_list, no_seq_err_list)
                    if mutated_seq == None:
                        num_err += 1
        
        except:
            return None
    
    return num_err
    