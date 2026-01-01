import os
import subprocess
import numpy as np
import pandas as pd
import warnings

def hmmpfam2(inputFile, outputFile, database):
    """
    Execute the hmmpfam2 command with the specified database, input, and output files.

    Parameters:
    inputFile: str - the path to the input file for hmmpfam2.
    outputFile: str - the path where the hmmpfam2 output will be saved.
    database: str - the path to the hmmpfam2 database to be used.

    Returns:
    bool: True if the hmmpfam2 output file was successfully generated, False otherwise.
    """
    # Construct the command as a list to safely pass to subprocess
    cmd = ["hmmpfam2", database, inputFile]

    try:
        # Open the output file in write mode
        with open(outputFile, 'w') as out:
            # Execute the command and redirect stdout to the output file
            subprocess.run(cmd, stdout=out, check=True)
        
        # Check if the output file was created and is not empty
        if os.path.exists(outputFile) and os.path.getsize(outputFile) > 0:
            print("hmmpfam2 output file generated.")
            return True
        else:
            warnings.warn("hmmpfam2 output file not generated.", UserWarning)
            return False
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing hmmpfam2: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False
    

def remove_multiple_spaces(string, remove_lead_trail=True):
    if remove_lead_trail:
        string = string.strip()
    while '  ' in string:
        string = string.replace('  ', ' ')
    return string

def get_state(state, line, print_state_transitions=False):
    #print(line)
    #print('incoming state ', state)
    newstate = ''
    if state in ['init', 'parsing']:
        if line.startswith('//') or line.startswith('- - - -'):
            newstate = 'recognise'
    elif state == 'recognise':
        if line.startswith('Alignments of top'):
            newstate = 'parsing'

    if print_state_transitions:
        if newstate != '':
            print('Changed state from ',state,' to ', newstate)
        else:
            newstate = state
    
    if newstate == '':
        newstate = state
    return newstate

def parse_hmmsearch_output(lines, hmmfiles):
    align_dict = {}
    align_code_dict = {}
    score_dict = {}
    detail_dict = {}
    line_score_idx = 0
    
    
    line_idx = 0
    line_align_idx = 0
    Id = ''
    state = 'init'
    curr_hmm = ''
    #print(lines)
    
    for line in lines:
        line_idx += 1        
        state = get_state(state, line)
        
        if state == 'recognise':
            if line.startswith('Query sequence'):
                Id = line.split(': ')[1].strip('\n')
                #print('set id to ',Id)
                detail_dict[Id] = {}

        elif state == 'parsing':
            hmmheader = np.asarray([line.startswith(hmmfile) for hmmfile in hmmfiles]).any()
            if hmmheader:
                line_align_idx = line_idx
                curr_hmm = line.split(':')[0]
                #print(curr_hmm)
                split = line.split(' ')
                #print(split)
                score_idx = split.index('score')+1
                from_idx = split.index('from')+1
                to_idx = split.index('to')+1
                detail_dict[Id][curr_hmm] = {'score':float(split[score_idx].strip(',')),'from': int(split[from_idx]),
                                             'to': int(split[to_idx].strip(':')), 'top':'', 'bottom':''}
                #print(detail_dict)
            elif line.startswith(' '):
                #print(line_idx, line_align_idx, line)
                if (line_idx - line_align_idx) % 4 == 1:
                    detail_dict[Id][curr_hmm]['top'] += remove_multiple_spaces(line).strip('*-><')
                elif (line_idx - line_align_idx) % 4 == 3:
                    detail_dict[Id][curr_hmm]['bottom'] += remove_multiple_spaces(line).split()[2]
        
    return detail_dict

def parse_hmmsearch_output_from_file(filename, hmmfile):
    with open(filename, 'r') as file:
        content = file.readlines()
    return parse_hmmsearch_output(content, hmmfile)

def get_best_alignment(mydict):
    #print(mydict)
    ret_dict = {}
    for Id in mydict:
        #print(Id)
        start = False
        score = -1
        best_hmm = ''
        for hmm in mydict[Id]:
            #print(hmm, best_hmm)
            if best_hmm == '' or score < mydict[Id][hmm]['score']:
                best_hmm = hmm
        ret_dict[Id] = mydict[Id][best_hmm].copy()
        ret_dict[Id]['hmm'] = best_hmm
        
    return ret_dict

def get_hmm_alignment(mydict, hmm):
    ret_dict = {}
    for key in mydict:
        try:
            ret_dict[key] = mydict[key][hmm].copy()
        except:
            print('Could not get',hmm,' for Id',key)
    return ret_dict

def removetopindels(indict):
    mydict = indict.copy()
    for Id in mydict:
        top_tmp = ''
        bot_tmp = ''
        idx = mydict[Id]['from']
        idx_list = []
        for a,b in zip(mydict[Id]['top'], mydict[Id]['bottom']):
            if a != '.':
                top_tmp += a
                bot_tmp += b
                if b == '-':
                    idx_list.append(idx-0.5)
                else:
                    idx_list.append(idx)
            if b != '-':
                idx += 1
        if mydict[Id]['top'] != top_tmp:
            print('Id:',Id,' top changed from ',mydict[Id]['top'], 'to', top_tmp)
        if mydict[Id]['bottom'] != bot_tmp:
            print('Id:',Id,' bottom changed from ',mydict[Id]['bottom'], 'to', bot_tmp)
        mydict[Id]['top'] = top_tmp
        mydict[Id]['bottom'] = bot_tmp
        assert(len(mydict[Id]['bottom']) == len(idx_list))
        mydict[Id]['idx_list'] = idx_list.copy()
    return mydict

def extractCharacters(Id, target, source, source_idx_list, pattern, idxs) :
    assert len(source) == len(source_idx_list)
    try:
        start = target.index(pattern)
    except:
        print('Problem at Id ', Id, ' pattern ', pattern, ' target ', target)
    ret = ''
    pos = []
    for idx in idxs:
        ret += source[start+idx]
        pos.append(source_idx_list[start+idx])
    return ret, pos

def extract_sig(Id, top, bottom, idx_list):
    try:
        s1, p1 = extractCharacters(Id, top, bottom, idx_list, "KGVmveHrnvvnlvkwl", [12, 15, 16])
        s2, p2 = extractCharacters(Id, top, bottom, idx_list, "LqfssAysFDaSvweifgaLLnGgt", [3,8,9,10,11,12,13,14,17])
        s3, p3 = extractCharacters(Id, top, bottom, idx_list, "iTvlnltPsl", [4,5])
        s4, p4 = extractCharacters(Id, top, bottom, idx_list, "LrrvlvGGEaL", [4,5,6,7,8])
        s5, p5 = extractCharacters(Id, top, bottom, idx_list, "liNaYGPTEtTVcaTi", [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

        return s1+s2+s3+s4+s5, p1+p2+p3+p4+p5
    except:
        return '', []

def extract_sig_dict(mydict):
    ret_dict = {}
    #ret_dict = mydict.copy()
    for Id in mydict:
        ret_dict[Id] = {}
        ret_dict[Id]['sig'], ret_dict[Id]['pos'] = extract_sig(Id, mydict[Id]['top'], mydict[Id]['bottom'], mydict[Id]['idx_list'])
        #ret_dict[Id]['sig'] = extract_sig(Id, mydict[Id]['top'], mydict[Id]['bottom'])
    return ret_dict

def print_signature_dict_to_file(mydict, filename):
    with open(filename, 'w') as file:
        file.writelines('Id,Extracted Signatures\n')
        for Id in mydict:
            file.writelines(Id + ',' + mydict[Id]['sig'] + '\n')

## hmm data 
database = "/storage/disk1/HRL/Project/SM_in_genome/analysis/tools_1025/34AA//aa-activating.aroundLys.hmm"

# prefix='peptaibol19_A'
prefix='peptaibol39_A'
input = f"{prefix}.fasta"
output = f"{prefix}.hmmout"

run_state = hmmpfam2(input, output,database)

detail_dict = parse_hmmsearch_output_from_file(output, ['aa-activating-core.198-334', 'aroundLys517'])

a_align_dict = get_hmm_alignment(detail_dict, 'aa-activating-core.198-334')
best_align_dict = get_best_alignment(detail_dict)
a_align_dict_no_indel = removetopindels(a_align_dict)

sig_dict = extract_sig_dict(a_align_dict_no_indel)

# rewrite into {XXX for XXX in } form
# key_sig_dict = {key.split("|")[0]: value for key, value in sig_dict.items() if value != {}} # MB

key_sig_dict = {key: value for key, value in sig_dict.items() if value != {}}

out_sig_path = f"{prefix}.sig.csv"
print_signature_dict_to_file(sig_dict, out_sig_path)