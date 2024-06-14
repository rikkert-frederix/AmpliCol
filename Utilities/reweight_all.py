#!/usr/bin/env python3.8

import concurrent.futures
import time
import os
import subprocess
from itertools import chain
import threading
import math
import shlex
import shutil


def run_fortran_program(executable, argument, output_file):
    """
    Function to run a Fortran program with a given argument.
    """
    try:
        result = subprocess.run([executable]+argument, capture_output=True, text=True)

        # Write the output to the file
        with open(output_file, 'w') as f_out:
            f_out.write(result.stdout)
            f_out.write(result.stderr)

        return (argument, result.stdout, result.stderr)
    except Exception as e:
        print('NO')
        return (argument, str(e), "")
        with open(output_file, 'w') as f_out:
            f_out.write(str(e))

if __name__ == '__main__':

    flavour={'4':'21 21 21 21',
             '5':'21 21 21 21 21',
             '6':'21 21 21 21 21 21',
             '7':'21 21 21 21 21 21 21',
             #'8':'21 21 21 21 21 21 21 21',
             #'9':'21 21 21 21 21 21 21 21 21',
            #'10':'21 21 21 21 21 21 21 21 21 21',
            #'11':'21 21 21 21 21 21 21 21 21 21 21',
            #'2':'21 21 21 21 21 21 21 21 21 21 21 21'
            }
    
    orders={'4':['1 2 3 4','1 3 2 4'],
            '5':['1 2 3 4 5','1 3 2 4 5'],
            '6':['1 2 3 4 5 6','1 3 2 4 5 6','1 3 4 2 5 6'],
            '7':['1 2 3 4 5 6 7','1 3 2 4 5 6 7','1 3 4 2 5 6 7'],
            #'8':['1 2 3 4 5 6 7 8','1 3 2 4 5 6 7 8','1 3 4 2 5 6 7 8','1 3 4 5 2 6 7 8'],
            #'9':['1 2 3 4 5 6 7 8 9','1 3 2 4 5 6 7 8 9','1 3 4 2 5 6 7 8 9','1 3 4 5 2 6 7 8 9'],
           #'10':['1 2 3 4 5 6 7 8 9 10','1 3 2 4 5 6 7 8 9 10','1 3 4 2 5 6 7 8 9 10','1 3 4 5 2 6 7 8 9 10','1 3 4 5 6 2 7 8 9 10'],
           #'11':['1 2 3 4 5 6 7 8 9 10 11','1 3 2 4 5 6 7 8 9 10 11','1 3 4 2 5 6 7 8 9 10 11','1 3 4 5 2 6 7 8 9 10 11','1 3 4 5 6 2 7 8 9 10 11'],
           #'2':['1 2 3 4 5 6 7 8 9 10 11 12','1 3 2 4 5 6 7 8 9 10 11 12','1 3 4 2 5 6 7 8 9 10 11 12','1 3 4 5 2 6 7 8 9 10 11 12','1 3 4 5 6 2 7 8 9 10 11 12','1 3 4 5 6 7 2 8 9 10 11 12']
           }

    nexternal=['4','5','6','7']#,'8','9']
    integrators=['1']#,'2','3','4']
    seeds=['101']#,'102','103','104','105','106','107','108','109','110']

    # Number of workers (adjust to the number of CPU cores or desired level of parallelism)
    max_workers = 8  # Change this to match the number of CPU cores you want to utilize

    executable='../matrix_reweight_QCD'
    # create the argument list
    arguments=[]
    outfiles=[]
    for n in nexternal:
        for i,color_order in enumerate(orders[n]):
            seed=seeds[0]
            integrator=integrators[0]
            add_arg='S'+seed+'I'+integrator
            directory='./Outputs'+add_arg+'/Res_files/'
            if not os.path.exists(directory):
                time.sleep(0.1)
                os.makedirs(directory)
            output_file='./Outputs'+add_arg+'/log_rwgt_'+n+'_'+str(i)+'.txt'
            outfiles.append(output_file)
            args=list(chain.from_iterable([n,flavour[n].split(),color_order.split()]))+[add_arg]
            print(args)
            arguments.append(args)

    # Shared counter for completed jobs
    completed_jobs = 0
    total_jobs = len(arguments)
    counter_lock = threading.Lock()
    print(f" We have {total_jobs} jobs to do")
    
    # Use ProcessPoolExecutor for parallel execution
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks to the executor
        future_to_argument = {executor.submit(run_fortran_program, executable, arg, outfile): arg for arg,outfile in zip(arguments,outfiles)}

        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_argument):
            arg = future_to_argument[future]
            try:
                stdout, stderr = future.result()[1], future.result()[2]
            except Exception as exc:
                print(f"Argument: {arg} generated an exception: {exc}")
                    
            # Update and print the counter
            with counter_lock:
                completed_jobs += 1
                print(f"Completed {completed_jobs}/{total_jobs} jobs.")


# Plot the reweighted event files and collect data
    for n in nexternal:
        for i,color_order in enumerate(orders[n]):
            args=list(chain.from_iterable([n,'2',flavour[n].split()]))
            iseed=seeds[0]
            integrator=integrators[0]
            add_arg='S'+seed+'I'+integrator
            arguments.append(args)
            proc_tag=''
            tag=''
            for el in args:
                proc_tag=proc_tag+'_'+str(el)
            tag=proc_tag
            for el in color_order:
                if (el !=' '):
                    tag=tag+'_'+str(el)
            output_file=' ../Outputs'+add_arg+'/events'+proc_tag
            execut='./plot_events/plot_events.sh'
            #res = subprocess.run(execut+output_file+' '+proc_tag+' '+tag,shell=True)

# Compute secondary unweighting efficiency

    for n in nexternal:
        for i,color_order in enumerate(orders[n]):
            args=list(chain.from_iterable([n,'2',flavour[n].split()]))
            proc_tag=''
            tag=''
            for el in args:
                proc_tag=proc_tag+'_'+str(el)
            tag=proc_tag
            for el in color_order:
                if (el !=' '):
                    tag=tag+'_'+str(el)
            print(tag)
            print(proc_tag)
            iseed=seeds[0]
            integrator=integrators[0]
            execut='./compute_unweighting_efficiency/compute.sh'
            output_file=' ../Outputs'+add_arg+'/events'+tag+'.lhe.rwgt'
            file_src='./Outputs'+add_arg+'/events'+tag+'.lhe.rwgt'
            file_dst='./compute_unweighting_efficiency/events'+tag+'.lhe.rwgt'
            file_name=' events'+tag+'.lhe.rwgt'
            shutil.copyfile(file_src,file_dst)
            res = subprocess.run(execut+file_name,shell=True)
            out_file='./compute_unweighting_efficiency/out_unwgt_'+tag+'.txt'
            out_dst='./plot_events/res_wgt_'+proc_tag+'/out_unwgt_'+tag+'.txt'
            shutil.move(out_file,out_dst)
            os.remove('./compute_unweighting_efficiency/events'+tag+'.lhe.rwgt')









