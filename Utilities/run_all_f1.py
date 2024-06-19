#!/usr/bin/env python3.8

import concurrent.futures
import time
import os
import subprocess
from itertools import chain
import threading
import math


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
        return (argument, str(e), "")
        with open(output_file, 'w') as f_out:
            f_out.write(str(e))

if __name__ == '__main__':

    flavour={'4':'21 21 1 -1',
             '5':'21 21 1 -1 21',
             '6':'21 21 1 -1 21 21',
             '7':'21 21 1 -1 21 21 21',
             '8':'21 21 1 -1 21 21 21 21',
             '9':'21 21 1 -1 21 21 21 21 21',
            #'10':'21 21 1 -1 21 21 21 21 21 21',
            #'11':'21 21 1 -1 21 21 21 21 21 21 21',
            #'12':'21 21 1 -1 21 21 21 21 21 21 21 21'
            }
    
    orders={'4':['3 1 2 4','3 2 1 4'],
            '5':['3 1 2 5 4','3 1 5 2 4','3 2 1 5 4',
                 '3 2 5 1 4','3 5 1 2 4','3 5 2 1 4'],
            '6':['3 1 2 5 6 4','3 1 5 2 6 4','3 1 5 6 2 4','3 5 1 2 6 4',
                 '3 5 6 1 2 4','3 5 1 6 2 4','3 2 1 5 6 4','3 2 5 1 6 4',
                 '3 2 5 6 1 4','3 5 2 1 6 4','3 5 6 2 1 4','3 5 2 6 1 4'],
            '7':['3 1 2 5 6 7 4','3 5 1 2 6 7 4','3 5 6 1 2 7 4','3 5 6 7 1 2 4',
                 '3 1 5 2 6 7 4','3 1 5 6 2 7 4','3 1 5 6 7 2 4','3 5 1 6 2 7 4',
                 '3 5 1 6 7 2 4','3 5 6 1 7 2 4','3 2 1 5 6 7 4','3 5 2 1 6 7 4',
                 '3 5 6 2 1 7 4','3 5 6 7 2 1 4','3 2 5 1 6 7 4','3 2 5 6 1 7 4',
                 '3 2 5 6 7 1 4','3 5 2 6 1 7 4','3 5 2 6 7 1 4','3 5 6 2 7 1 4'],
            '8':['3 1 2 5 6 7 8 4','3 1 5 2 6 7 8 4','3 1 5 6 2 7 8 4','3 1 5 6 7 2 8 4',
                 '3 1 5 6 7 8 2 4','3 5 1 2 6 7 8 4','3 5 1 6 2 7 8 4','3 5 1 6 7 2 8 4',
                 '3 5 1 6 7 8 2 4','3 5 6 1 2 7 8 4','3 5 6 1 7 2 8 4','3 5 6 1 7 8 2 4',
                 '3 5 6 7 1 2 8 4','3 5 6 7 1 8 2 4','3 5 6 7 8 1 2 4','3 2 1 5 6 7 8 4',
                 '3 2 5 1 6 7 8 4','3 2 5 6 1 7 8 4','3 2 5 6 7 1 8 4','3 2 5 6 7 8 1 4',
                 '3 5 2 1 6 7 8 4','3 5 2 6 1 7 8 4','3 5 2 6 7 1 8 4','3 5 2 6 7 8 1 4',
                 '3 5 6 2 1 7 8 4','3 5 6 2 7 1 8 4','3 5 6 2 7 8 1 4','3 5 6 7 2 1 8 4',
                 '3 5 6 7 2 8 1 4','3 5 6 7 8 2 1 4'],
            '9':['3 1 2 5 6 7 8 9 4','3 1 5 2 6 7 8 9 4','3 1 5 6 2 7 8 9 4','3 1 5 6 7 2 8 9 4',
                 '3 1 5 6 7 8 2 9 4','3 1 5 6 7 8 9 2 4','3 5 1 2 6 7 8 9 4','3 5 1 6 2 7 8 9 4',
                 '3 5 1 6 7 2 8 9 4','3 5 1 6 7 8 2 9 4','3 5 1 6 7 8 9 2 4','3 5 6 1 2 7 8 9 4',
                 '3 5 6 1 7 2 8 9 4','3 5 6 1 7 8 2 9 4','3 5 6 1 7 8 9 2 4','3 5 6 7 1 2 8 9 4',
                 '3 5 6 7 1 8 2 9 4','3 5 6 7 1 8 9 2 4','3 5 6 7 8 1 2 9 4','3 5 6 7 8 1 9 2 4',
                 '3 5 6 7 8 9 1 2 4',
                 '3 2 1 5 6 7 8 9 4','3 2 5 1 6 7 8 9 4','3 2 5 6 1 7 8 9 4','3 2 5 6 7 1 8 9 4',
                 '3 2 5 6 7 8 1 9 4','3 2 5 6 7 8 9 1 4','3 5 2 1 6 7 8 9 4','3 5 2 6 1 7 8 9 4',
                 '3 5 2 6 7 1 8 9 4','3 5 1 6 7 8 2 9 4','3 5 1 6 7 8 9 2 4','3 5 6 1 2 7 8 9 4',
                 '3 5 6 2 7 1 8 9 4','3 5 6 2 7 8 1 9 4','3 5 6 2 7 8 9 1 4','3 5 6 7 2 1 8 9 4',
                 '3 5 6 7 2 8 1 9 4','3 5 6 7 2 8 9 1 4','3 5 6 7 8 2 1 9 4','3 5 6 7 8 2 9 1 4',
                 '3 5 6 7 8 9 2 1 4'],
            #'10':[],
            #'11':[],
            #'12':[]
           }

    nexternal=['4','5','6','7','8','9']
    imodes=['0','1','2']
    integrators=['1','2','3','4']
    seeds=['101','102','103','104','105','106','107','108','109','110']

    # Number of workers (adjust to the number of CPU cores or desired level of parallelism)
    max_workers = 8  # Change this to match the number of CPU cores you want to utilize

    for imode in imodes:
        executable='../matrix_integrate_QCD'
        # create the argument list
        arguments=[]
        outfiles=[]
        for n in nexternal:
            for i,color_order in enumerate(orders[n]):
                for integrator in integrators:
                    for seed in seeds:
                        add_arg='S'+seed+'I'+integrator
                        directory='./Outputs'+add_arg+'/Res_files/'
                        if not os.path.exists(directory):
                            time.sleep(0.1)
                            os.makedirs(directory)
                        output_file='./Outputs'+add_arg+'/log_'+n+'_'+imode+'_'+str(i)+'.txt'
                        outfiles.append(output_file)
                        args=list(chain.from_iterable([integrator,imode,n,flavour[n].split(),color_order.split()]))+[add_arg]
                        print(args)
                        arguments.append(args)

        # Shared counter for completed jobs
        completed_jobs = 0
        total_jobs = len(arguments)
        counter_lock = threading.Lock()
        print(f"Doing imode {imode}. We have {total_jobs} jobs to do")
    
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
                    print(f"Completed {completed_jobs}/{total_jobs} jobs for imode {imode}")









