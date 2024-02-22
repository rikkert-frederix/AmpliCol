#!/usr/bin/env python
  
import subprocess
import multiprocessing
import numpy as np
import itertools
import copy

class Integrator:
    def __init__(self, process, integ):
        self.process = process
        self.integ = integ
        
    def convert_proc_line_to_pdg(self):
        n=0
        self.proc=[]
        for l in process.split():
            if l=='d':
                self.proc.append(1)
                n=n+1
            if l=='dx':
                self.proc.append(-1)
                n=n+1
            if l=='u':
                self.proc.append(2)
                n=n+1
            if l=='ux':
                self.proc.append(-2)
                n=n+1
            if l=='s ':
                self.proc=self.proc+'3 '
                n=n+1
            if l=='c':
                self.proc=self.proc+'4 '
                n=n+1
            if l=='b':
                self.proc=self.proc+'5 '
                n=n+1
            if l=='t':
                self.proc=self.proc+'6 '
                n=n+1
            if l=='g':
                self.proc.append(21)
                n=n+1
            if l=='a':
                self.proc.append(22)
                n=n+1
        return n
    
    def set_integrator(self):
        if self.integ=='haag':
            self.int=2
        if self.integ=='gen23':
            self.int=1
        if self.integ=='genpt':
            self.int=3

    def get_all_color_orders(self):
        color_order=[0 for i in range(0,n)]
        for i in [0,1]:
            if abs(self.proc[i])<=6:
                self.proc[i]=-self.proc[i]
        for i,t in enumerate(self.proc):
            if t<0:
                color_order[n-1]=i+1
            if (t>0 and t<7):
                color_order[0]=i+1
        if all([v == 0 for v in color_order]):
            color_order[n-1]=n
        for i,t in enumerate(self.proc):
            if (t != 21 and abs(t)>6):
                color_order[n-2]=i+1

        for i in [0,1]:
            if abs(self.proc[i])<=6:
                self.proc[i]=-self.proc[i]
        to_perm=0
        to_perm_list=[]
        for i,r in enumerate(self.proc):
            if (i+1 not in color_order):
                to_perm=to_perm+1
                to_perm_list.append(i+1)

        self.color_order_list=[]
        for perm in list(itertools.permutations(to_perm_list)):
            color_order_comp=copy.deepcopy(color_order)
            k=0
            for i,t in enumerate(color_order):
                if (t == 0):
                    color_order_comp[i] = perm[k]
                    k=k+1
            self.color_order_list.append(color_order_comp)

        i=0
        print('The color orders needed are:')
        print('(Those in curly brackets are to be permuted)')
        print_co = ''
        for t in color_order:
            if (t != 0):
                print_co = print_co + ' ' + str(t)
            else:
                print_co = print_co + ' {'+ str(to_perm_list[i])+'}'
                i = i+1
        print(print_co)
        print('Explicitly:')
        for i,r in enumerate(self.color_order_list):
            print(str(i+1)+': '+str(r))

    def prompt_on_co(self):
        prompt = 'Pick color order to integrate (type in label)\n'+\
                'Or: type -1 for generating all independent ones.\n'
        self.co_pick = input(prompt)

    def get_ind_co(self):
        self.color_order_list

def run_program(program, args, output_file):
        with open(output_file, 'w') as f:
            subprocess.call([program] + args, stdout=f)


#process = raw_input('Type in process: \n' )
process='u ux > g g'   # TO CHANGE
#integrator = raw_input('Type in integrator to use (haag, genpt or gen23):\n' )
integrator = 'haag'    # TO CHANGE

job = Integrator(process, integrator)
n = job.convert_proc_line_to_pdg()
job.set_integrator()
job.get_all_color_orders()
job.prompt_on_co()

if (job.co_pick != -1):
    color_order = job.color_order_list[job.co_pick-1]
elif (job.co_pick == -1):
    job.get_ind_co()

#nm = input('Pick mode 0,1 or 2 (or -1 for all):\n')
nm = -1  # TO CHANGE
if nm==-1:
    modes = ['0','1','2']
elif nm==0:
    modes = ['0']
elif nm==1:
    modes = ['1']
elif nm==2:
    modes = ['2']

for mode in modes:
    processes=[]
    program='./matrix_integrate_QCD'
    output_file='log_'+str(job.int)+'_'+mode+'_'+str(n)+'_'
    args=[str(job.int),mode,str(n)]
    for r in job.proc:
        args.append(str(r))
        output_file=output_file+str(r)+'_'
    for r in color_order:
        args.append(str(r))
        output_file=output_file+str(r)+'_'
    output_file = output_file+'.txt'
    print('running...')
    print(args)
    processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))

    # Start all processes
    for process in processes:
        process.start()

    # Wait for all processes to finish
    for process in processes:
        process.join()

to_rwgt = raw_input('Would you like to reweight to NLC and FC? (y/n):\n')

if to_rwgt == 'y':
    imodes=['2']
    program='./matrix_reweight_QCD'
            
    for imode in imodes:
        processes=[]
        output_file='log_'+str(n)+'_'
        args=[str(n)]
        for r in job.proc:
            args.append(str(r))
            output_file=output_file+str(r)+'_'
        for r in color_order:
            args.append(str(r))
            output_file=output_file+str(r)+'_'
        output_file = output_file+'.txt'
        print('running rwgt...')
        print(args) 
        processes.append(multiprocessing.Process(target=run_program, args=(program, args, output_file)))

        # Start all processes
        for process in processes:
            process.start()

        # Wait for all processes to finish
        for process in processes:
            process.join()
else:
    print('End of run.')

