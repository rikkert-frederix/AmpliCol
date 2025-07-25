#!/usr/bin/env python

import itertools
import copy
import re
import argparse
from collections import Counter
import multiprocessing
import math

# Global sets (make then 'frozenset' so that they are immutable):
quarks=frozenset({'d','u','s','c','b','t'})
antiquarks=frozenset({'d~','u~','s~','c~','b~','t~'})
singlets=frozenset({'a','z','w+','w-','e+','e-','mu+','mu-','ta+','ta-','ve','ve~','vm','vm~','vt','vt~','h'})
gluons=frozenset({'g'})
#gluons=frozenset({})
flavour_scheme=frozenset({'d','u','s','c','b'}) # all the massless quarks
#flavour_scheme=frozenset({'d','u','s','c'}) # all the massless quarks
all_coloured=quarks | antiquarks | gluons
massless_QCD=flavour_scheme | frozenset([q+'~' for q in flavour_scheme]) | gluons
proton=massless_QCD
jet=massless_QCD
if jet != proton:
    raise ValueError("definition of 'jet' and 'proton' should be the same")
pdgs={'g':'21','d':'1','u':'2','s':'3','c':'4','b':'5','t':'6','d~':'-1','u~':'-2','s~':'-3','c~':'-4','b~':'-5','t~':'-6','a':'22','z':'23','w+':'24','w-':'-24','e+':'-11','e-':'11','mu+':'-13','mu-':'13','ta+':'-15','ta-':'15','ve':'12','ve~':'-12','vm':'14','vm~':'-14','vt':'16','vt~':'-16','h':'25'}
anti_particle={'g':'g','d':'d~','u':'u~','s':'s~','c':'c~','b':'b~','t':'t~','d~':'d','u~':'u','s~':'s','c~':'c','b~':'b','t~':'t','a':'a','z':'z','w+':'w-','w-':'w+','e+':'e-','e-':'e+','mu+':'mu-','mu-':'mu+','ta+':'ta-','ta-':'ta+','ve':'ve~','ve~':'ve','vm':'vm~','vm~':'vm','vt':'vt~','vt~':'vt','h':'h'}
#sort_particles={'g':0,'d':1,'u':2,'s':3,'c':4,'b':5,'t':6,'d~':7,'u~':8,'s~':9,'c~':10,'b~':11,'t~':12,'a':99,'z':99,'w+':99,'w-':99,'e+':99,'e-':99,'mu+':99,'mu-':99,'ta+':99,'ta-':99,'ve':99,'ve~':99,'vm':99,'vm~':99,'vt':99,'vt~':99,'h':99}
sort_particles={'g':13,'d':1,'u':2,'s':3,'c':4,'b':5,'t':6,'d~':7,'u~':8,'s~':9,'c~':10,'b~':11,'t~':12,'a':80,'z':81,'w+':82,'w-':83,'e+':84,'e-':85,'mu+':86,'mu-':87,'ta+':88,'ta-':89,'ve':90,'ve~':91,'vm':92,'vm~':93,'vt':94,'vt~':95,'h':96}



def ProcessProcess(proc):
    """Function to process each 'proc' in parallel"""
    phase_space_orders_local = {}  # Local dictionary to avoid race conditions

    all_possible_color_ord = list(itertools.permutations(range(len(proc))))
    valid_color_ord = [perm for perm in all_possible_color_ord if ValidColorOrd(proc, perm)]
    unique_color_ord = [perm for perm in valid_color_ord if UniqueColorOrd(proc, perm)]

    for perm in unique_color_ord:
        ordered_proc, ordered_perm = OrderProcPerm(proc, perm)
        zero = ordered_perm.index(0)
        perm_mapped = tuple(ordered_perm[zero:] + ordered_perm[:zero])

        if perm_mapped in phase_space_orders_local:
            phase_space_orders_local[perm_mapped].append((ordered_proc, ordered_perm, []))
        else:
            phase_space_orders_local[perm_mapped] = [(ordered_proc, ordered_perm, [])]

    return phase_space_orders_local


def ValidColorOrd(proc,perm):
    # Check if 'perm' is a valid color order for the process 'proc'
    found_quark = found_antiquark = found_singlet = found_gluon = False
    for idx in perm:
        particle = proc[idx]
        if particle in quarks:
            if found_quark or found_gluon:
                return False
            if found_antiquark and idx < perm[0]:  # Two-quark line process, reject one ordering
                return False
            found_quark = True
            found_antiquark = found_singlet = found_gluon = False
        elif particle in antiquarks:
            if found_antiquark or found_singlet or not found_quark:
                return False
            found_antiquark = True
            found_quark = found_singlet = found_gluon = False
        elif particle in gluons:
            if found_antiquark or found_singlet:
                return False
            found_gluon = True
        else:  # Assuming the rest are singlets
            if found_quark or found_gluon or not found_antiquark :
                return False
            found_singlet = True
    if found_gluon: # Only for all-gluon process found_gluon is True here. Use it to remove cyclic permutations
        if perm[0] != 0 :
            return False
    return True

def UniqueColorOrd(proc,perm):
    # Check if 'perm' is a color order in canonical order for the process 'proc'.
    #
    # Start at the first incoming particle. From there, the indentical
    # final state particles should each come in *increasing* order.
    zero=perm.index(0)
    perm_mapped=perm[zero:]+perm[:zero]
    for particle in all_coloured:
        particle_positions=[]
        for i,part in enumerate(proc[2:]):
            if part==particle:
                particle_positions.append(i+2)
        previous_position=0        
        for i in particle_positions:
            if perm_mapped.index(i) < previous_position:
                return False
            else:
                previous_position=perm_mapped.index(i)
    return True

def OrderProcPerm(proc,perm):
    zero=perm.index(0)
    perm_mapped=list(perm[zero:]+perm[:zero])
    proc_mapped=proc[zero:]+proc[:zero]
    elements_to_order=[]
    for i in perm_mapped:
        if proc[i] in massless_QCD and i > 1:
            elements_to_order.append(i)
    indices=[perm_mapped.index(x) for x in elements_to_order]
    sorted_elements=sorted(elements_to_order)
    for i, val in zip(indices,sorted_elements):
        perm_mapped[i]=val
    perm_ordered=perm_mapped[len(perm)-zero:]+perm_mapped[:len(perm)-zero]
    
    proc_ordered=[None]*len(proc)
    for i in range(len(perm_ordered)):
        proc_ordered[perm_ordered[i]]=proc[perm[i]]
        
    # if there are two quark lines there are two options for the order. Pick the right one
    if count_matching_elements(proc,quarks) == 2:
        for q in quarks:
            qi=[(i+1,j) for i,j in enumerate(perm_ordered[1:]) if proc_ordered[j] == q]
            if qi:
                break
        if qi[0][1] < perm_ordered[0]:
            perm_ordered=perm_ordered[qi[0][0]:]+perm_ordered[:qi[0][0]]
    return tuple(proc_ordered),tuple(perm_ordered)


def ParseCollision(input_string):
    input_string=input_string.replace('bar','~')
    parts=input_string.split(">")
    if len(parts) != 2:
        raise ValueError("Invalid collision format. Expected 'p p > ...'.")
    initial_state=parts[0].strip().split()
    for i,p in enumerate(initial_state):
        if p != 'p': 
            initial_state[i]=anti_particle[p]
    final_state=parts[1].strip().split()
    jet_match=re.match(r"(\d+)j",final_state[-1]) if final_state else None
    jet_count=int(jet_match.group(1)) if jet_match else 0
    rest=final_state[:-1] if jet_match else final_state
    return {"initial_state":initial_state,"jet_count":jet_count,"rest":rest}

def count_matching_elements(main_list,check_list):
    main_counts=Counter(main_list)
    count=sum(main_counts[item] for item in check_list if item in main_counts)
    return count

def ValidProc(proc):
    nq=count_matching_elements(proc,quarks)
    naq=count_matching_elements(proc,antiquarks)
    if nq > 2 : return False    # at most two quarks
    if naq > 2 : return False   # at most two anti-quarks
    if nq != naq : return False # same number of quarks and anti-quarks
    # remove flavour changing currents:
#    for q in quarks:
#        if count_matching_elements(proc,[q]) != count_matching_elements(proc,[q+'~']) : return False
    # need at least one quark line if there are colour singlets:
    if nq == 0 and count_matching_elements(proc,singlets) > 0 : return False
    return True

def CompatibleUniqueProc(process,proc):
    mandatory=[]
    proc_local=proc.copy()
    for part in process['initial_state']:
        if part != 'p':
            mandatory.append(part)
    mandatory.extend(process['rest'])
    try:
        for p in mandatory:
            proc_local.remove(p)
        return True
    except:
        return False
    
def CompatibleProc(process,proc):
    proc_local=list(proc[2:])
    for i in [0,1]:
        if process['initial_state'][i] != 'p':
            if proc[i] != process['initial_state'][i]:
                return False
    try:
        for p in process['rest']:
            proc_local.remove(p)
    except:
        return False
    return True
    
def GenerateAllUniqueProcs(process):
    procs=[[]]
    for part in process['initial_state']:
        if part != 'p' and part not in jet:
            raise ValueError("Initial state should be a proton ('p').")
    jp=0
    for part in process['rest']:
        if part not in jet : continue
        jp=jp+1
    for part in range(process['jet_count']+2+jp):
        procs_new=[]
        for proc in procs:
            for p in jet:
                if part > 3 and p not in gluons: continue
                if not procs_new:
                    procs_new=[sorted(proc+[p])]
                else:
                    procs_new.append(sorted(proc+[p]))
        procs=procs_new.copy()
    for part in process['rest']:
        if part in jet : continue
        for proc in procs:
            proc.append(part)
    unique_procs=[]
    for proc in procs:
        if ValidProc(proc) and CompatibleUniqueProc(process,proc):
            unique_procs.append(tuple(proc))
    return set(unique_procs)

def GenerateAllProcs(unique_procs,process):
    procs=set()
    for proc in unique_procs:
        for i,j in itertools.combinations(range(len(proc)),2):
            if proc[i] in jet and proc[j] in jet:
                pair1=[proc[i],proc[j]]
                pair2=[proc[j],proc[i]]
                remaining=[e for k,e in enumerate(proc) if k not in (i,j)]
                proc1=tuple(pair1+remaining)
                if CompatibleProc(process,proc1):
                    procs.add(proc1)
                proc2=tuple(pair2+remaining)
                if CompatibleProc(process,proc2):
                    procs.add(proc2)
    return procs

def CombineResults(results):
    phase_space_orders = {}
    for result in results:
        for key, value in result.items():
            if key in phase_space_orders:
                phase_space_orders[key].extend(value)
            else:
                phase_space_orders[key] = value
    return phase_space_orders
                
def ParseArgument():
    parser=argparse.ArgumentParser(description="Generate the full list of processes, ordered by phase-space order")
    parser.add_argument("process_string",type=str,help="Process to consider (e.g., 'p p > w+ z 4j')")
    args=parser.parse_args()
    return ParseCollision(args.process_string)

def IdenticalParticleSymmetryFactor(proc):
    i_fac=1
    for p in all_coloured:
        i_fac*=max(1,math.factorial(proc[2:].count(p)))
    return i_fac


def MultiChannelPartners(proc,perm,k,l):
    all_possible_perms={perm}
    colour_singlets_in_proc=[perm.index(i) for i,p in enumerate(proc) if p in singlets]
    anti_quarks_in_proc=tuple([perm.index(i) for i,p in enumerate(proc) if p in antiquarks])
    if colour_singlets_in_proc:
        all_singlet_orders=tuple(itertools.permutations(colour_singlets_in_proc))
    else:
        all_singlet_orders=()
    if len(anti_quarks_in_proc) == 1:
        for s in all_singlet_orders:
            order=[]
            for i in range(len(perm)):
                if i == anti_quarks_in_proc[0]:
                    order.extend([perm[p] for p in list(anti_quarks_in_proc+s)])
                elif i not in colour_singlets_in_proc:
                    order.append(perm[i])
            all_possible_perms.add(tuple(order))
    elif len(anti_quarks_in_proc) == 2:
        for j in range(len(all_singlet_orders)+1):
            for s in all_singlet_orders:
                order=[]
                for i in range(len(perm)):
                    if i == anti_quarks_in_proc[0]:
                        order.extend([perm[p] for p in list((anti_quarks_in_proc[0],)+s[:j])])
                    elif i == anti_quarks_in_proc[1]:
                        order.extend([perm[p] for p in list((anti_quarks_in_proc[1],)+s[j:])])
                    elif i not in colour_singlets_in_proc:
                        order.append(perm[i])
                all_possible_perms.add(tuple(order))
    mt=[]
    for o in all_possible_perms:
        found=False
        for i,key in enumerate(all_keys_sorted):
            for (process,order,multichannel) in phase_space_orders[key]:
                if process == proc and order==o:
                    if found:
                        print('FOUND DOUBLE')
                    else:
                        found=True
                        mt.append(i)
#                        print('found',process,proc,order,o,i)
        if not found:
            print('NOT FOUND')
    phase_space_orders[k][l]=(proc,perm,tuple(sorted(mt)))
#    print(proc,':',all_singlet_orders,':',anti_quarks_in_proc,':',all_possible_perms)

def DetermineMultiChannelPartnersAndSymmetryFactor():
    for j,key in enumerate(all_keys_sorted):
        for i,(process,order,multichannel) in enumerate(phase_space_orders[key]):
            MultiChannelPartners(process,order,key,i)
    for key in all_keys_sorted:
        for i,(process,order,multichannel) in enumerate(phase_space_orders[key]):
            phase_space_orders[key][i]=(process,order,multichannel,IdenticalParticleSymmetryFactor(process))

def ConvertProcToString(proc):
    process,order,multi_channel,iden=proc
    crossed=[pdgs[p] if i>1 else pdgs[anti_particle[p]] for i,p in enumerate(process)]
    line=str(len(multi_channel))
    line=line+'   '+' '.join([str(m+1) for m in multi_channel])
    line=line+'   '+' '.join(crossed)
    line=line+'   '+' '.join([str(o+1) for o in order])
    line=line+'   '+str(iden)
    return line

def sort_by_pdg_codes(process):
    nq=count_matching_elements(process,quarks)
    if nq == 2:
        quarks_in_proc=tuple([process[i] for i,p in enumerate(process) if p in quarks])
        same_flavour=quarks_in_proc[0]==quarks_in_proc[1]
    else:
        same_flavour=False
    val=0
    val+=nq*2
    if same_flavour : val=val+1
    return (val,[sort_particles[p] for p in process]) # first sort by 'val', then by (modified) PDG codes.

def sort_by_pdg_codes2(proc):
    process=proc[0]
    return sort_by_pdg_codes(process)

def WriteAllProcsIntoList():
    towrite=[]
    towrite.append(str(len(all_keys_sorted)))
    towrite.append('')
    for i,key in enumerate(all_keys_sorted):
        towrite.append(str(i+1)+'   '+str(len(phase_space_orders[key]))+'   '+str(max(len(proc[2]) for proc in phase_space_orders[key]))+'   '+' '.join([str(k+1) for k in key]))
        process_list=sorted(phase_space_orders[key],key=sort_by_pdg_codes2)
        for proc in process_list:
            process_line=ConvertProcToString(proc)
            towrite.append(process_line)
        towrite.append('')
        towrite.append('')
        towrite.append('')
    return towrite

def Add2qq_dfProcesses(sorted_procs):
    # add the 2qq_df processes with the two incoming particles interchanged:
    i=0
    while i < len(sorted_procs):
        proc = sorted_procs[i]
        if proc[2] in antiquarks and proc[3] in antiquarks and proc[2] != proc[3]:
            swapped_proc=proc[:]
            swapped_proc[2],swapped_proc[3]=swapped_proc[3],swapped_proc[2]
            sorted_procs.insert(i+1,swapped_proc)
            i+=1
        i+=1
    return sorted_procs

def WriteUniqueProcsIntoList(procs):
#    sorted_procs=sorted([sorted(proc,key=lambda x: int(pdgs[x])) for proc in procs],key=sort_by_pdg_codes)
    sorted_procs=sorted([sorted(proc,key=lambda x: sort_particles[x]) for proc in procs],key=sort_by_pdg_codes)
    sorted_procs=Add2qq_dfProcesses(sorted_procs)
    line=[str(len(sorted_procs[0]))+' '+str(len(sorted_procs))]
    for proc in sorted_procs:
        line.append(' '.join(pdgs[p] for p in proc))
    line.append('')
    line.append('')
    line.append('')
    return line
    
if __name__ == "__main__":    
    process=ParseArgument()
    all_unique_procs=GenerateAllUniqueProcs(process)
    all_procs=GenerateAllProcs(all_unique_procs,process)
    
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(ProcessProcess, all_procs)  # Parallelize across procs
    phase_space_orders=CombineResults(results)
    all_keys_sorted=sorted(phase_space_orders.keys())
    DetermineMultiChannelPartnersAndSymmetryFactor() # updates the phase_space_orders dictionary
    towriteunique=WriteUniqueProcsIntoList(all_unique_procs)
    towriteallprocs=WriteAllProcsIntoList() # puts the phase_space_orders dictionary in a writable list
    
    with open('processes.txt','w') as f:
        f.write('\n'.join(towriteunique))
        f.write('\n'.join(towriteallprocs))
        f.write('\n')
