#!/usr/bin/env python

# multi-jet


import itertools
import math


def create_all_procs_from_unique_procs(proc,swap_qq):
    all_procs = []
    # Iterate over all permutations of the first two elements
    for perm in itertools.permutations(proc, 2):
        if perm[0] not in proton or perm[1] not in proton : continue
        # Find the remaining elements
        remaining = list(proc)
        remaining.remove(perm[0])
        remaining.remove(perm[1])
        
        # Sort the last three elements since order does not matter
        remaining_sorted = tuple(sorted(remaining))
        remaining_sorted_swap = tuple(sorted(swap_qq.get(i,i) for i in remaining))
        
        # Combine the permuted first two elements with the sorted last three elements
        combined = perm+remaining_sorted
        combined_swap = tuple(swap_qq.get(iperm,iperm) for iperm in combined[0:2]) + remaining_sorted_swap

        if (combined not in all_procs) and (combined_swap not in all_procs):
            # Append to all_procs as a tuple
            all_procs.append(combined)
    return all_procs


def is_valid_permutation(perm):
    # Check the cyclic constraint: negative numbers must be followed by positive ones
    n = len(perm)
    for i in range(n):
        if (perm[i] == 'qbar' or perm[i] == 'qpbar' ):  # If the current element is anti-quark
            next_index = (i + 1) % n  # Cyclic next element
            if (perm[next_index] != 'q' and perm[next_index] != 'qp' ):  # next element must be quark
                return False
        if perm[i] == 'a' :
            next_index = (i + 1) % n  # Cyclic next element
            if (perm[next_index] != 'qbar' and perm[next_index] != 'qpbar'  and perm[next_index] != 'a' ):  # next element must be anti-quark or another singlet
                return False
    return True

def generate_permutations(arr,swap_qq):
    # The first element stays fixed
    first_element = arr[0]
    rest_elements = arr[1:]
    
    valid_permutations = []
    
    # Generate all permutations for the remaining 6 elements
    for perm in itertools.permutations(rest_elements):
        # Form the new list with the fixed first element and permuted rest
        full_perm = (first_element,) + perm
        full_perm_swap = tuple(swap_qq.get(i,i) for i in full_perm)
        
        # Check if the permutation satisfies the cyclic constraint
        if is_valid_permutation(full_perm):
            if full_perm not in valid_permutations:
                if arr[0]=='g' and arr[1]=='g':
                    if full_perm_swap not in valid_permutations:
                        valid_permutations.append(full_perm)
                else:
                    valid_permutations.append(full_perm)
                
    
    return valid_permutations

def particle(part):
    return(canonical_order[swap_ini[part]])
    

def convert_to_input(phase_space_order,perm,swap_ini,pso,unique_procs):
    try:
        shift=perm.index('q')
    except ValueError:
        shift=0

    flavour_scheme={'q':['d','u','s','c','b'],
                    'qbar':['dbar','ubar','sbar','cbar','bbar'],
                    'qp':['d','u','s','c','b'],
                    'qpbar':['dbar','ubar','sbar','cbar','bbar']}
        
    all_proc=[]
    if 'q' in perm or 'qbar' in perm:
        for i,q in enumerate(flavour_scheme['q']):
            new_list1=[item for item in perm]
            new_list1=[flavour_scheme['q'][i] if item == 'q' else item for item in new_list1]
            new_list1=[flavour_scheme['qbar'][i] if item == 'qbar' else item for item in new_list1]
            if 'qp' in new_list1 or 'qpbar' in new_list1:
                for j,qp in enumerate(flavour_scheme['qp']):
                    if flavour_scheme['qp'][j] == flavour_scheme['q'][i] : continue
                    new_list2=[item for item in new_list1]
                    new_list2=[flavour_scheme['qp'][j] if item == 'qp' else item for item in new_list2]
                    new_list2=[flavour_scheme['qpbar'][j] if item == 'qpbar' else item for item in new_list2]
                    all_proc.append(new_list2)
            else:
                all_proc.append(new_list1)
    else:
        all_proc.append(perm)
    
    input=[]
    for perm in all_proc:
        sperm=sorted(perm,key=particle)
        if sperm not in unique_procs:
            unique_procs.append(sperm)
            if sperm[2] != sperm[3] and canonical_order[sperm[2]] < 13 and canonical_order[sperm[3]] < 13:
                unique_procs.append(sperm[0:2]+[sperm[3]]+[sperm[2]]+sperm[4:])
        ini1=swap_ini[perm[0]]
        ini2=swap_ini[perm[pso+1]]
        process=[conversion[ini1]]+[conversion[ini2]]+[conversion[perm[i]] for i in range(1,pso+1)]+[conversion[perm[i]] for i in range(pso+2,len(perm))]
        process+=[' ']+[str(i) for i in phase_space_order[shift:]]+[str(i) for i in phase_space_order[:shift]]# + ['\n']
        input.append(process)

    return input

canonical_order={'dbar':1,
                 'ubar':2,
                 'sbar':3,
                 'cbar':4,
                 'bbar':5,
                 'tbar':6,
                 'd':7,
                 'u':8,
                 's':9,
                 'c':10,
                 'b':11,
                 't':12,
                 'g':13,
                 'a':14}

conversion={'g':'21',
            'd':'1','dbar':'-1',
            'u':'2','ubar':'-2',
            's':'3','sbar':'-3',
            'c':'4','cbar':'-4',
            'b':'5','bbar':'-5',
            't':'6','tbar':'-6',
            'a':'22'}
swap_ini={'g':'g',
          'd':'dbar',
          'dbar':'d',
          'u':'ubar',
          'ubar':'u',
          's':'sbar',
          'sbar':'s',
          'c':'cbar',
          'cbar':'c',
          'b':'bbar',
          'bbar':'b',
          'a':'a'}
swap_qq={'g':'g',
         'q':'qp',
         'qp':'q',
         'qbar':'qpbar',
         'qpbar':'qbar',
         'a':'a'}

nfinal=4

# multi-jet base processes (without gluons):
base_procs=[[],['q','qbar'],['q','qp','qbar','qpbar'],['q','q','qbar','qbar']]
#base_procs=[[],['q','qbar'],['q','qp','qbar','qpbar']]
#base_procs=[['q','qp','qbar','qpbar']]
#base_procs=[['q','qbar']]

# Add a photon
for proc in base_procs:
    proc.append('a')


proton=['g','q','qp','qbar','qpbar']

# extend the base_procs with additional gluons. These are all the unique procs
unique_procs=[]
for proc in base_procs:
    while len(proc) < nfinal+2:
        proc.append('g')
    if (len(proc) == nfinal+2):
        unique_procs.append(proc)

base_proc={}
all_procs=[]
for i,proc in enumerate(unique_procs):
    all_new_procs=create_all_procs_from_unique_procs(proc,swap_qq)
    all_procs=all_procs+all_new_procs
    for iproc in all_new_procs:
        base_proc[iproc]=i

# Get symmetry factor
iden_fac={}
expected_number={}
for proc in all_procs:
    i_fac=1
    i_fac*=max(1,math.factorial(proc[2:].count('g')))
    i_fac*=max(1,math.factorial(proc[2:].count('q')))
    i_fac*=max(1,math.factorial(proc[2:].count('qbar')))
    i_fac*=max(1,math.factorial(proc[2:].count('qp')))
    i_fac*=max(1,math.factorial(proc[2:].count('qpbar')))
    if proc[2:].count('q') == proc[2:].count('qbar') and proc[2:].count('q') == proc[2:].count('qp') and proc[2:].count('q') == proc[2:].count('qpbar') and proc[2:].count('q') == 1:
        i_fac*=2
    
    iden_fac[proc]=i_fac
    if base_proc[proc] == 0:
        n_co=math.factorial(nfinal+1)
    elif base_proc[proc] == 1:
        n_co=math.factorial(nfinal)
    elif base_proc[proc] == 2:
        n_co=math.factorial(nfinal-2)*(nfinal-1)*2
    elif base_proc[proc] == 3:
        n_co=math.factorial(nfinal-2)*(nfinal-1)*2
    if int(n_co/iden_fac[proc])*iden_fac[proc] != n_co:
        print('Error: not an integer')
        quit()
    expected_number[proc]=int(n_co/iden_fac[proc])
    

order=[1]+[i for i in range(3,nfinal+3)]
phase_space_order=[]

for i in range(1,nfinal+2):
    phase_space_order.append(tuple(order[:i]+[2]+order[i:]))

pso_map = [{} for _ in phase_space_order]


to_write=[[] for _ in phase_space_order]


unique_procs=[]
for i,proc in enumerate(all_procs):
    print('looping through processes... ',proc)
    valid_perms=generate_permutations(proc,swap_qq)
    for pso,psorder in enumerate(phase_space_order):
        pso_map[pso][proc]=[]
        for valid_perm in valid_perms:
            if valid_perm[0]==proc[0] and valid_perm[pso+1]==proc[1]:
                pso_map[pso][proc].append(valid_perm)
                fprocs=convert_to_input(psorder,valid_perm,swap_ini,pso,unique_procs)
                for fproc in fprocs:
                    indices1 = sorted([i+1 for i, x in enumerate(fproc) if x == "22"])
                    indices2 = sorted([int(x) for x in fproc[len(fproc)-len(indices1)-1:len(fproc)-1]])
                    if indices1 == indices2 :
                        to_write[pso].append(fproc+[' ']+[str(iden_fac[proc])])


with open('processes.txt','w') as f:
    f.write(str(nfinal+2)+' '+str(len(unique_procs))+'\n')
    for proc in unique_procs:
        f.write(' '.join([conversion[ele] for ele in proc])+'\n')
    f.write('\n')
    f.write('\n')
    f.write(str(len(to_write))+'\n')
    f.write('\n')
    for pso in to_write:
        f.write(str(len(pso))+'\n')
        for ele in pso:
            f.write(' '.join(ele)+'\n')
        f.write('\n')
        f.write('\n')
        f.write('\n')
        f.write('\n')

