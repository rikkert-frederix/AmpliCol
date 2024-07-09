#!/usr/bin/env python3.8

import copy 

def compare_lists(l1,l2):
    diff=[]
    for i,il in enumerate(l1):
        if (il != l2[i]):
                diff.append(il)
                diff.append(l2[i])
    return diff


list=['3 8 7 1 2 4 5 6','3 7 1 8 2 4 5 6','3 7 1 2 8 4 5 6','3 7 1 2 4 5 8 6',
                    '3 1 8 7 2 4 5 6','3 1 7 2 8 4 5 6','3 1 7 2 4 5 8 6','3 1 2 8 7 4 5 6',
                    '3 1 2 7 4 5 8 6','3 1 2 4 5 8 7 6','3 8 7 1 4 5 2 6','3 7 1 8 4 5 2 6',
                    '3 7 1 4 5 8 2 6','3 7 1 4 5 2 8 6','3 1 8 7 4 5 2 6','3 1 7 4 5 8 2 6',
                    '3 1 7 4 5 2 8 6','3 1 4 5 8 7 2 6','3 1 4 5 7 2 8 6','3 1 4 5 2 8 7 6',
                    '3 8 7 2 1 4 5 6','3 7 2 8 1 4 5 6','3 7 2 1 8 4 5 6','3 7 2 1 4 5 8 6',
                    '3 2 8 7 1 4 5 6','3 2 7 1 8 4 5 6','3 2 7 1 4 5 8 6','3 2 1 8 7 4 5 6',
                    '3 2 1 7 4 5 8 6','3 2 1 4 5 8 7 6'
                   ]

new_list=[]

n='8'

keep_fixed=['3','4','5','6']

symmetric=['7','8','9']

n_new = str(int(n)+1)

for el in list:
    for ip,part in enumerate(el.split()):
        if part == keep_fixed[2]: continue
        if part == keep_fixed[0]: continue
        new_ord=el.split()
        new_ord.insert(ip,n_new)
        new_list.append(new_ord)

out_list = copy.deepcopy(new_list)

for ord in out_list:
    for ord2 in new_list:
        diff=[]
        diff = compare_lists(ord,ord2)
        check=True
        for tr in diff:
            if tr not in symmetric:
                check=False
        if check and diff and ord2 in out_list:
            out_list.remove(ord2)

output=[]

for out in out_list:
    out_string = ' '.join(out)
    output.append(out_string)

print(output)

