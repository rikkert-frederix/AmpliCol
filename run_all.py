#!/usr/bin/env python

import subprocess
import multiprocessing
import numpy as np

def run_program(program, args, output_file):
    """Function to run a program with a list of arguments"""
    with open(output_file, 'w') as f:
        subprocess.call([program] + args, stdout=f)

if __name__ == '__main__':

    nexternal=['4','5']
#    nexternal=['10','11','12']
#    nexternal=['7']

    for n in nexternal:
        color_orders=[str(i) for i in range(0,int((int(n))/2))]
        imodes=['0','1','2']
        program='./matrix_integrate'
        for imode in imodes:
            processes=[]
            for color_order in color_orders:
                output_file='log_'+n+'_'+imode+'_'+color_order+'.txt'
                args=[n,imode,color_order]
                processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))
        
            # Start all processes
            for process in processes:
                process.start()
 
            # Wait for all processes to finish
            for process in processes:
                process.join()
 
        imodes=['2']
        program='./matrix_reweight'
            
        for imode in imodes:
            processes=[]
            for color_order in color_orders:
                output_file='log_'+n+'_'+imode+'_'+color_order+'_rwgt.txt'
                args=[n,imode,color_order]
                processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))
        
            # Start all processes
            for process in processes:
                process.start()

            # Wait for all processes to finish
            for process in processes:
                process.join()
