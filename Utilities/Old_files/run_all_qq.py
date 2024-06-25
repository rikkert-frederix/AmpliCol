#!/usr/bin/env python

import subprocess
import multiprocessing
import numpy as np

def run_program(program, args, output_file):
    """Function to run a program with a list of arguments"""
    with open(output_file, 'w') as f:
        subprocess.call([program] + args, stdout=f)
        print(program)

if __name__ == '__main__':

    nexternal=[4]
    integ='1'

    for n in nexternal:
        imodes=['0','1','2']
        program='./matrix_integrate_QCD'
        proc_type = [n-2,n-3,n-4]
        for imode in imodes:
            processes=[]
            for ptype in proc_type:
                if ptype==n-2:
                    c_o_i=n+1
                    c_o_j=ptype
                    c_o_k=n+1
                    for type in ['1','2']:
                        output_file='log_'+str(n)+'_'+imode+'_'+type+'_'+str(c_o_i)+\
                            '_'+str(c_o_j)+'_'+str(c_o_k)+'.txt'
                        args=[integ,str(n),imode,type,str(c_o_i),str(c_o_j),str(c_o_k)]
                        processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))
                if ptype==n-3:
                    c_o_i=n+1
                    for i in range(0,n-2):
                        c_o_j=i
                        c_o_k=n-3-c_o_j
                        type='1'
                        output_file='log_'+str(n)+'_'+imode+'_'+type+\
                                    '_'+str(c_o_i)+\
                                    '_'+str(c_o_j)+'_'+str(c_o_k)+'.txt'
                        args=[integ,str(n),imode,type,str(c_o_i),str(c_o_j),str(c_o_k)]
                        processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))
                    c_o_k=n+1
                    for i in range(0,n-2):
                        c_o_j=i
                        c_o_i=n-3-c_o_j
                        type='2'
                        output_file='log_'+str(n)+'_'+imode+'_'\
                                    +type+'_'+str(c_o_i)+\
                                    '_'+str(c_o_j)+'_'+str(c_o_k)+'.txt'
                        args=[integ,str(n),imode,type,str(c_o_i),str(c_o_j),str(c_o_k)]
                        print(args)
                        processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))
                if ptype==n-4:
                    for i in range(0,n-3):
                        c_o_i=i
                        for j in range(0,n-3-c_o_i):
                            c_o_j=j
                            c_o_k=n-4-c_o_i-c_o_j
                            for type in ['1','2']:
                                output_file='log_'+str(n)+'_'+imode+'_'\
                                    +type+'_'+str(c_o_i)+\
                                   '_'+str(c_o_j)+'_'+str(c_o_k)+'.txt'
                                args=[integ,str(n),imode,type,str(c_o_i),str(c_o_j),str(c_o_k)]
                                processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))

            print(processes)
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

