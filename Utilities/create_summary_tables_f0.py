#!/usr/bin/env python3.8

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain

def get_string(tag1,tag2,plot):
    string=tag1
#    try:
#        string=string+" & %8.3e"%cross_section[tag]
#    except:
#        string=string+" & --"
#    try:
#        string=string+" & %8.3e"%uncertainty[tag]
#    except:
#        string=string+" & --"
    try:
        if ('s101') in tag1:
             string=string+" & \multirow{10}{*}{\includegraphics[scale=0.2]{save_diag/"+plot+"}}"
        else:
            string=string+" & "
        #string=string+" & %7.2f \%%"%(uncertainty[tag]/cross_section[tag] * 100.)
    except:
        string=string+" & --"
#    try:
#        string=string+" & %2.1f"%(chi2_per_DoF[tag])
#    except:
#        string=string+" & --"
#    try:
#        string=string+" & %ik"%(number_of_events[tag]/1000)
#    except:
#        string=string+" & --"
#    try:
#        string=string+" & %ik"%(number_passing_cuts[tag]/1000)
#    except:
#        string=string+" & --"
#    try:
#        string=string+" & %2.1f \%%"%(number_passing_cuts[tag]/number_of_events[tag] *100.)
#    except:
#        string=string+" & --"
#    try:
#        string=string+" & %i "%(total_time[tag])
#    except:
#        string=string+" & --"
#    try:
#        string=string+" & %2.3f \%%"%(gen_eff[tag] * 100.)
#    except:
#        string=string+" & --"
    try:
        string=string+" & %2.3f \%%"%(gen_eff[tag2]*number_of_events[tag2]/number_passing_cuts[tag2] * 100.)
    except:
        string=string+" & --"
    try:
        string=string+" & %7.4f "%unwgt_eff[tag2]
    except:
        string=string+" & --"
    string=string+" \\\\ \n"
    return string

def compute_averages(list_of_tags):
    ave_tag='average'
    l=len(list_of_tags)
    if l==0:
        try:
            del cross_section[ave_tag]
        except:
            pass
        try:
            del uncertainty[ave_tag]
        except:
            pass
        try:
            del chi2_per_DoF[ave_tag]
        except:
            pass
        try:
            del number_of_events[ave_tag]
        except:
            pass
        try:
            del number_passing_cuts[ave_tag]
        except:
            pass
        try:
            del total_time[ave_tag]
        except:
            pass
        try:
            del gen_eff[ave_tag]
        except:
            pass
        return ave_tag
    try:
        cross_section[ave_tag]=sum(cross_section[tag] for tag in list_of_tags)/l
    except:
        try:
            del cross_section[ave_tag]
        except:
            pass
    try:
        uncertainty[ave_tag]=math.sqrt(sum(uncertainty[tag]**2 for tag in list_of_tags))/l
    except:
        try:
            del uncertainty[ave_tag]
        except:
            pass
#    try:
#        chi2_per_DoF[ave_tag]=sum(chi2_per_DoF[tag] for tag in list_of_tags)/l
#    except:
#        try:
#            del chi2_per_DoF[ave_tag]
#        except:
#            pass
    try:
        number_of_events[ave_tag]=sum(number_of_events[tag] for tag in list_of_tags)/l
    except:
        try:
            del number_of_events[ave_tag]
        except:
            pass
    try:
        number_passing_cuts[ave_tag]=sum(number_passing_cuts[tag] for tag in list_of_tags)/l
    except:
        try:
            del number_passing_cuts[ave_tag]
        except:
            pass
    try:
        total_time[ave_tag]=sum(total_time[tag] for tag in list_of_tags)/l
    except:
        try:
            del total_time[ave_tag]
        except:
            pass
    try:
        gen_eff[ave_tag]=sum(gen_eff[tag] for tag in list_of_tags)/l
    except:
        try:
            del gen_eff[ave_tag]
        except:
            pass
    try:
        unwgt_eff[ave_tag]=sum(unwgt_eff[tag] for tag in list_of_tags)/l
    except:
        try:
            del unwgt_eff[ave_tag]
        except:
            pass
    return ave_tag

if __name__ == '__main__':
    table_header=r"""\noindent The tag has the form \textbf{nX-oX-fX-mX-iX-sXXX}, where X is an integer. The integer after 
\textbf{n} corresponds to the number of external particles,
\textbf{o} the colour order index: minimum number of gluons between the two incoming particles in the colour order,
\textbf{f} flavour index (always 0 for all-gluon),
\textbf{m} corresponds to imode: 0 for grid setup; 1 for upper bounding envelope estimation; 2 for event generation,
\textbf{i} the integrator: 1 for gen23; 2 for haag; 3 for genpt (chili-like); 4 for t-channel
\textbf{s} the random number seed.

\begin{tabularx}{0.95\textwidth}{l r r r}
  \toprule
  tag&rel.unc.&unf.eff(passing)&sec.unwgt.eff\\
  \midrule
"""
    table_footer=r"""\end{tabularx}
\newpage
"""
    latex_header=r"""\documentclass[10pt]{article}
\usepackage{tabularx,booktabs,multirow}
\usepackage[a4paper,left=1cm,right=1cm]{geometry}
\usepackage{graphicx}
\renewcommand\arraystretch{1.3}
\begin{document}
\begin{scriptsize}
"""
    latex_footer=r"""\end{scriptsize}
\end{document}
"""
    
        
    flavours={'4':['21 21 21 21',                        '21 21 1 -1',                        '-1 21 -1 21',                       '-1 1 21 21'                  ],
              '5':['21 21 21 21 21',                     '21 21 1 -1 21',                     '-1 21 -1 21 21',                    '-1 1 21 21 21'               ],
              '6':['21 21 21 21 21 21',                  '21 21 1 -1 21 21',                  '-1 21 -1 21 21 21',                 '-1 1 21 21 21 21'            ],
              '7':['21 21 21 21 21 21 21',               '21 21 1 -1 21 21 21',               '-1 21 -1 21 21 21 21',              '-1 1 21 21 21 21 21'         ],
              '8':['21 21 21 21 21 21 21 21',            '21 21 1 -1 21 21 21 21',            '-1 21 -1 21 21 21 21 21',           '-1 1 21 21 21 21 21 21'      ],
              '9':['21 21 21 21 21 21 21 21 21',         '21 21 1 -1 21 21 21 21 21',         '-1 21 -1 21 21 21 21 21 21',        '-1 1 21 21 21 21 21 21 21'   ],
             '10':['21 21 21 21 21 21 21 21 21 21',      '21 21 1 -1 21 21 21 21 21 21',      '-1 21 -1 21 21 21 21 21 21 21',     '-1 1 21 21 21 21 21 21 21 21']}

    orders=[{},{},{},{}]

    # gg -> (n-2)g
    orders[0]={'4':['1 2 3 4','1 3 2 4'],
               '5':['1 2 3 4 5','1 3 2 4 5'],
               '6':['1 2 3 4 5 6','1 3 2 4 5 6','1 3 4 2 5 6'],
               '7':['1 2 3 4 5 6 7','1 3 2 4 5 6 7','1 3 4 2 5 6 7'],
               '8':['1 2 3 4 5 6 7 8','1 3 2 4 5 6 7 8','1 3 4 2 5 6 7 8','1 3 4 5 2 6 7 8'],
               '9':['1 2 3 4 5 6 7 8 9','1 3 2 4 5 6 7 8 9','1 3 4 2 5 6 7 8 9','1 3 4 5 2 6 7 8 9'],
              '10':['1 2 3 4 5 6 7 8 9 10','1 3 2 4 5 6 7 8 9 10','1 3 4 2 5 6 7 8 9 10','1 3 4 5 2 6 7 8 9 10','1 3 4 5 6 2 7 8 9 10']
              }

    # gg -> qqbar + (n-4)g
    orders[1]={'4':['3 1 2 4'],
               '5':['3 1 2 5 4','3 1 5 2 4'],
               '6':['3 1 2 5 6 4','3 1 5 2 6 4','3 1 5 6 2 4',
                    '3 5 1 2 6 4'],
               '7':['3 1 2 5 6 7 4','3 5 1 2 6 7 4','3 1 5 2 6 7 4','3 1 5 6 2 7 4','3 1 5 6 7 2 4','3 5 1 6 2 7 4'],
               '8':['3 1 2 5 6 7 8 4','3 1 5 2 6 7 8 4','3 1 5 6 2 7 8 4','3 1 5 6 7 2 8 4','3 1 5 6 7 8 2 4',
                    '3 5 1 2 6 7 8 4','3 5 1 6 2 7 8 4','3 5 1 6 7 2 8 4',
                    '3 5 6 1 2 7 8 4'],
               '9':['3 1 2 5 6 7 8 9 4','3 1 5 2 6 7 8 9 4','3 1 5 6 2 7 8 9 4','3 1 5 6 7 2 8 9 4','3 1 5 6 7 8 2 9 4','3 1 5 6 7 8 9 2 4',
                    '3 5 1 2 6 7 8 9 4','3 5 1 6 2 7 8 9 4','3 5 1 6 7 2 8 9 4','3 5 1 6 7 8 2 9 4',
                    '3 5 6 1 2 7 8 9 4','3 5 6 1 7 2 8 9 4'],
              '10':['3 1 2 5 6 7 8 9 10 4','3 1 5 2 6 7 8 9 10 4','3 1 5 6 2 7 8 9 10 4','3 1 5 6 7 2 8 9 10 4','3 1 5 6 7 8 2 9 10 4','3 1 5 6 7 8 9 2 10 4','3 1 5 6 7 8 9 10 2 4',
                    '3 5 1 2 6 7 8 9 10 4','3 5 1 6 2 7 8 9 10 4','3 5 1 6 7 2 8 9 10 4','3 5 1 6 7 8 2 9 10 4','3 5 1 6 7 8 9 2 10 4',
                    '3 5 6 1 2 7 8 9 10 4','3 5 6 1 7 2 8 9 10 4','3 5 6 1 7 8 2 9 10 4',
                    '3 5 6 7 1 2 8 9 10 4']
              }

    # qbar g -> qbar + (n-3)g
    orders[2]={'4':['1 2 4 3','1 4 2 3'],
               '5':['1 2 4 5 3','1 4 2 5 3','1 4 5 2 3'],
               '6':['1 2 4 5 6 3','1 4 2 5 6 3','1 4 5 2 6 3','1 4 5 6 2 3'],
               '7':['1 2 4 5 6 7 3','1 4 2 5 6 7 3','1 4 5 2 6 7 3',
                    '1 4 5 6 2 7 3','1 4 5 6 7 2 3'],
               '8':['1 2 4 5 6 7 8 3','1 4 2 5 6 7 8 3','1 4 5 2 6 7 8 3',
                    '1 4 5 6 2 7 8 3','1 4 5 6 7 2 8 3','1 4 5 6 7 8 2 3'],
               '9':['1 2 4 5 6 7 8 9 3','1 4 2 5 6 7 8 9 3','1 4 5 2 6 7 8 9 3',
                    '1 4 5 6 2 7 8 9 3','1 4 5 6 7 2 8 9 3','1 4 5 6 7 8 2 9 3','1 4 5 6 7 8 9 2 3'],
              '10':['1 2 4 5 6 7 8 9 10 3','1 4 2 5 6 7 8 9 10 3','1 4 5 2 6 7 8 9 10 3','1 4 5 6 2 7 8 9 10 3','1 4 5 6 7 2 8 9 10 3','1 4 5 6 7 8 2 9 10 3','1 4 5 6 7 8 9 2 10 3','1 4 5 6 7 8 9 10 2 3']
              }

    # qbar q -> (n-2)g
    orders[3]={'4':['1 3 4 2'],
               '5':['1 3 4 5 2'],
               '6':['1 3 4 5 6 2'],
               '7':['1 3 4 5 6 7 2'],
               '8':['1 3 4 5 6 7 8 2'],
               '9':['1 3 4 5 6 7 8 9 2'],
              '10':['1 3 4 5 6 7 8 9 10 2']
              }

            
    nexternal=['4','5','6','7','8']#,'9']#,'10','11','12']
    imodes=['0','1','2']
    integrators=['1','2','3','4']
    seeds=['101','102','103','104','105','106','107','108','109','110']
        
    
    
    # collect all the results.
    cross_section={}
    uncertainty={}
    chi2_per_DoF={}
    number_of_events={}
    number_passing_cuts={}
    total_time={}
    gen_eff={}
    tags=[]
    unwgt_eff={}
    for n in nexternal:
        for iflav,flavour in enumerate(flavours[n]):
            for io,order in enumerate(orders[iflav][n]):
                for imode in imodes:
                    for integrator in integrators:
                        for seed in seeds:
                            log_file='./OutputsS'+seed+'I'+integrator+'/log_'+n+'_'+imode+'_'+str(io)+'.txt'
                            args=list(chain.from_iterable([n,'2',flavour.split()]))
                            proc_tag=''
                            for el in args:
                                proc_tag=proc_tag+'_'+str(el)
                            tag=proc_tag
                            for el in order:
                                if (el !=' '):
                                    tag=tag+'_'+str(el)
                            log_file2='./res_wgt_'+proc_tag+'/out_unwgt_'+tag+'.txt'
                            try:
                                with open(log_file) as file:
                                    tag='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+imode+'-i'+integrator+'-s'+seed
                                    for line in file:
                                        if 'Final result:'in line:
                                            cross_section[tag]=float(line.split()[2])
                                            uncertainty[tag]=float(line.split()[4])
                                        if 'per D.o.F.:' in line:
                                            chi2_per_DoF[tag]=float(line.split()[3])
                                        if 'Number of events:' in line:
                                            number_of_events[tag]=int(line.split()[3])
                                        if 'Number passing cuts:' in line:
                                            number_passing_cuts[tag]=int(line.split()[3])
                                        if 'Total time:' in line:
                                            total_time[tag]=float(line.split()[2])
                                            tags.append(tag) 
                                        if 'Generation efficiencies:' in line:
                                            gen_eff[tag]=float(line.split()[2])
                                with open(log_file2) as file:
                                    tag='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+imode+'-i'+integrator+'-s'+seed
                                    for line in file:
                                        if 'efficiency'in line:
                                            unwgt_eff[tag]=float(line.split()[3])
                            except:
                                pass

    for n in nexternal:
        for iflav,flavour in enumerate(flavours[n]):
            for io,order in enumerate(orders[iflav][n]):
                for imode in ['1']:
                    for i,integrator in enumerate(integrators):
                        rel_unc_data = np.array([])
                        for seed in seeds:
                            tag='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+imode+'-i'+integrator+'-s'+seed
                            rel_unc_data = np.append(rel_unc_data,[uncertainty[tag]/cross_section[tag]])
                        fig,ax = plt.subplots()
                        plt.boxplot(rel_unc_data,whis=100)
                        fig.patch.set_visible(False)
                        ax.axis('off')
                        plt.savefig('save_diag/box_unc_'+str(n)+'_'+str(io)+'_'+str(i)+'.pdf',transparent=True,format="pdf", bbox_inches="tight")

    with open('tables.tex','w') as f:
        f.write(latex_header)
        for n in nexternal:
            for iflav,flavour in enumerate(flavours[n]):
                for io,order in enumerate(orders[iflav][n]):
                    f.write(table_header)
                    for i,integrator in enumerate(integrators):
                        tags_to_average=[]
                        for seed in seeds:
                            tag1='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+'1'+'-i'+integrator+'-s'+seed
                            tag2='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+'2'+'-i'+integrator+'-s'+seed
                            plot='box_unc_'+str(n)+'_'+str(io)+'_'+str(i)+'.pdf'
                            string=get_string(tag1,tag2,plot)
                            f.write(string)
                            tags_to_average.append(tag2)
                        f.write('\\midrule \n')
                        tag=compute_averages(tags_to_average)
                        string=get_string(tag,tag,plot)
                        f.write(string)
                        f.write('\\bottomrule \n')
                    f.write(table_footer)
        f.write(latex_footer)

