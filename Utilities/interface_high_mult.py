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
            if l=='d~':
                self.proc.append(-1)
                n=n+1
            if l=='u':
                self.proc.append(2)
                n=n+1
            if l=='u~':
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
        quarks=0
        for i in range(0,n-1):
            if (abs(self.proc[i])<=6):
                quarks=quarks+1
                if (i <= 1):
                    self.proc[i]=-self.proc[i]
        last=False
        first=False
        quark_pair=[0,0]
        for i,t in enumerate(self.proc):
            if (t<0 and not last):
                color_order[n-1]=i+1
                last = True
            elif (t<0 and last):
                quark_pair[0]=i+1
            if (t>0 and t<=6 and not first):
                color_order[0]=i+1
                first=True
            elif (t>0 and t<=6 and first):
                quark_pair[1]=i+1
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
            if (self.proc[i] == 21):
                to_perm=to_perm+1
                to_perm_list.append(i+1)

        if (quarks > 2):
            to_perm_list.append(n+1)

        self.color_order_list=[]
        if (quarks > 2):
            for perm in list(itertools.permutations(to_perm_list)):
                invalid=False
                color_order_comp=copy.deepcopy(color_order)
                k=0
                for i,t in enumerate(color_order_comp):
                    if (t == 0 and perm[k] != n+1):
                        color_order_comp[i] = perm[k]
                        k=k+1
                    elif ( t==0 and perm[k] == n+1):
                        if (color_order_comp[i+1] == 0):
                            color_order_comp[i] = quark_pair[0]
                            color_order_comp[i+1] = quark_pair[1]
                            k=k+1
                        else:
                            invalid=True
                    if (k > len(perm)-1):
                        break

    
                if (not invalid):
                    self.color_order_list.append(color_order_comp)
                    print(color_order_comp)

            i=0
            print('The color orders needed are:')
            #print('(Those in curly brackets are to be permuted)')
            print_co = ''
            for j,t in enumerate(color_order):
                if (t != 0 and self.proc[j] != 21):
                    print_co = print_co + ' ' + str(t)
                elif ( t != 0 and to_perm_list[i] != n+1):
                    print_co = print_co + ' {'+ str(to_perm_list[i])+'}'
                    i = i+1
            #print(print_co)
            print('Explicitly:')
            for i,r in enumerate(self.color_order_list):
                print(str(i+1)+': '+str(r))

    def prompt_on_co(self):
        prompt = 'Pick color order to integrate (type in label)\n'+\
                'Or: type -1 for generating all independent ones.\n'
        #self.co_pick = input(prompt) 
        self.co_pick = -1

    def get_ind_co(self):
        self.color_order_list

def run_program(program, args, output_file):
        with open(output_file, 'w') as f:
            subprocess.call([program] + args, stdout=f)


#process = raw_input('Type in process: \n' )
process='d~ d > u u~ g g'   # TO CHANGE
#integrator = raw_input('Type in integrator to use (haag, genpt or gen23):\n' )
integrator = 'gen23'    # TO CHANGE

print('Process:')
print(process)
job = Integrator(process, integrator)
n = job.convert_proc_line_to_pdg()
job.set_integrator()
job.get_all_color_orders()
job.prompt_on_co()
color_order=[]

if (job.co_pick != -1):
    color_order.append(job.color_order_list[job.co_pick-1])
elif (job.co_pick == -1):
    job.get_ind_co()
    color_order = job.color_order_list

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

for col_ord in color_order:

 print(col_ord)

 for mode in modes:
    processes=[]
    program='./matrix_integrate_QCD'
    output_file='log_'+str(job.int)+'_'+mode+'_'+str(n)+'_'
    args=[str(job.int),mode,str(n)]
    for r in job.proc:
        args.append(str(r))
        output_file=output_file+str(r)+'_'
    for r in col_ord:
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

print('\nFinished generating all channels for LC events.')
print('####################\n')
to_rwgt = raw_input('Would you like to reweight to NLC and FC? (y/n):\n')

if to_rwgt == 'y':
    imodes=['2']
    make_process = subprocess.Popen("make;", shell=True, stdout=subprocess.PIPE)
    program='./matrix_reweight_QCD'

    for col_ord in color_order:

      for imode in imodes:
        processes=[]
        output_file='log_'+str(n)+'_'
        args=[str(n)]
        for r in job.proc:
            args.append(str(r))
            output_file=output_file+str(r)+'_'
        for r in col_ord:
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
    print('End of reweight.\n')
else:
    print('\nNo reweighting done.\n!!Note: the events are only LC accurate!\nEnd of run.\n')

