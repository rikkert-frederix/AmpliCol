#!/usr/bin/env python

import subprocess
import multiprocessing
import numpy as np
from itertools import chain

def run_program(program, args, output_file):
    """Function to run a program with a list of arguments"""
    with open(output_file, 'w') as f:
        subprocess.call([program] + args, stdout=f)

if __name__ == '__main__':

    flavour={'4':'21 21 21 21',
             '5':'21 21 21 21 21',
             '6':'21 21 21 21 21 21',
             '7':'21 21 21 21 21 21 21',
             '8':'21 21 21 21 21 21 21 21',
             '9':'21 21 21 21 21 21 21 21 21',
            '10':'21 21 21 21 21 21 21 21 21 21',
            '11':'21 21 21 21 21 21 21 21 21 21 21',
            '12':'21 21 21 21 21 21 21 21 21 21 21 21'}
    
    orders={'4':['1 2 3 4','1 3 2 4'],
            '5':['1 2 3 4 5','1 3 2 4 5'],
            '6':['1 2 3 4 5 6','1 3 2 4 5 6','1 3 4 2 5 6'],
            '7':['1 2 3 4 5 6 7','1 3 2 4 5 6 7','1 3 4 2 5 6 7'],
            '8':['1 2 3 4 5 6 7 8','1 3 2 4 5 6 7 8','1 3 4 2 5 6 7 8','1 3 4 5 2 6 7 8'],
            '9':['1 2 3 4 5 6 7 8 9','1 3 2 4 5 6 7 8 9','1 3 4 2 5 6 7 8 9','1 3 4 5 2 6 7 8 9'],
           '10':['1 2 3 4 5 6 7 8 9 10','1 3 2 4 5 6 7 8 9 10','1 3 4 2 5 6 7 8 9 10','1 3 4 5 2 6 7 8 9 10','1 3 4 5 6 2 7 8 9 10'],
           '11':['1 2 3 4 5 6 7 8 9 10 11','1 3 2 4 5 6 7 8 9 10 11','1 3 4 2 5 6 7 8 9 10 11','1 3 4 5 2 6 7 8 9 10 11','1 3 4 5 6 2 7 8 9 10 11'],
           '12':['1 2 3 4 5 6 7 8 9 10 11 12','1 3 2 4 5 6 7 8 9 10 11 12','1 3 4 2 5 6 7 8 9 10 11 12','1 3 4 5 2 6 7 8 9 10 11 12','1 3 4 5 6 2 7 8 9 10 11 12','1 3 4 5 6 7 2 8 9 10 11 12']}

    
    nexternal=['4','5','6','7','8','9']
    imodes=['0','1']
    
    for imode in imodes:
        program='./matrix_integrate_QCD'
        processes=[]
        for n in nexternal:
            for i,color_order in enumerate(orders[n]):
                output_file='log_'+n+'_'+imode+'_'+str(i)+'.txt'
                args=list(chain.from_iterable(['1',imode,n,flavour[n].split(),color_order.split()]))
                print(args,output_file)
                processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))
        
        # Start all processes
        for process in processes:
            process.start()
 
        # Wait for all processes to finish
        for process in processes:
            process.join()
 
#        imodes=['2']
#        program='./matrix_reweight'
#            
#        for imode in imodes:
#            processes=[]
#            for color_order in color_orders:
#                output_file='log_'+n+'_'+imode+'_'+color_order+'_rwgt.txt'
#                args=[n,imode,color_order]
#                processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))
#
#            # Start all processes
#            for process in processes:
#                process.start()
#
#            # Wait for all processes to finish
#            for process in processes:
#                process.join()

