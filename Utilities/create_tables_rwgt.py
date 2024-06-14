#!/usr/bin/env python

import sys
from itertools import chain

def get_string(tag):
    string=tag
    try:
        string=string+" & %7.4f "%unwgt_eff[tag]
    except:
        string=string+" & --"
    string=string+" \\\\ \n"
    return string


if __name__ == '__main__':
    table_header=r"""\noindent The tag has the for nX-oX-fX-mX-iX-sXXX, where X is an integer. The integer after 
\begin{description}
\item[n] corresponds to the number of external particles,
\item[o] the colour order index: minimum number of gluons between the two incoming particles in the colour order,
\item[f] flavour index (always 0 for all-gluon),
\item[m] corresponds to imode: 0 for grid setup; 1 for upper bounding envelope estimation; 2 for event generation,
\item[i] the integrator: 1 for gen23; 2 for haag; 3 for genpt (chili-like); 4 for t-channel
\item[s] the random number seed.
\end{description}
\begin{tabularx}{0.95\textwidth}{l r}
  \toprule
  tag&unwgt.eff.\\
  \midrule
"""
    table_footer=r"""  \bottomrule
\end{tabularx}
\newpage
"""
    latex_header=r"""\documentclass[10pt]{article}
\usepackage{tabularx,booktabs,multirow}
\usepackage[a4paper,left=1cm,right=1cm]{geometry}
\renewcommand\arraystretch{1.3}
\begin{document}
\begin{scriptsize}
"""
    latex_footer=r"""\end{scriptsize}
\end{document}
"""
    
    
    flavours={'4':['21 21 21 21'],
              '5':['21 21 21 21 21'],
              '6':['21 21 21 21 21 21'],
              '7':['21 21 21 21 21 21 21'],
              '8':['21 21 21 21 21 21 21 21'],
              '9':['21 21 21 21 21 21 21 21 21'],
             '10':['21 21 21 21 21 21 21 21 21 21'],
             '11':['21 21 21 21 21 21 21 21 21 21 21'],
              '12':['21 21 21 21 21 21 21 21 21 21 21 21']}
    
    orders={'4':['1 2 3 4','1 3 2 4'],
            '5':['1 2 3 4 5','1 3 2 4 5'],
            '6':['1 2 3 4 5 6','1 3 2 4 5 6','1 3 4 2 5 6'],
            '7':['1 2 3 4 5 6 7','1 3 2 4 5 6 7','1 3 4 2 5 6 7'],
            '8':['1 2 3 4 5 6 7 8','1 3 2 4 5 6 7 8','1 3 4 2 5 6 7 8','1 3 4 5 2 6 7 8'],
            '9':['1 2 3 4 5 6 7 8 9','1 3 2 4 5 6 7 8 9','1 3 4 2 5 6 7 8 9','1 3 4 5 2 6 7 8 9'],
           '10':['1 2 3 4 5 6 7 8 9 10','1 3 2 4 5 6 7 8 9 10','1 3 4 2 5 6 7 8 9 10','1 3 4 5 2 6 7 8 9 10','1 3 4 5 6 2 7 8 9 10'],
           '11':['1 2 3 4 5 6 7 8 9 10 11','1 3 2 4 5 6 7 8 9 10 11','1 3 4 2 5 6 7 8 9 10 11','1 3 4 5 2 6 7 8 9 10 11','1 3 4 5 6 2 7 8 9 10 11'],
           '12':['1 2 3 4 5 6 7 8 9 10 11 12','1 3 2 4 5 6 7 8 9 10 11 12','1 3 4 2 5 6 7 8 9 10 11 12','1 3 4 5 2 6 7 8 9 10 11 12','1 3 4 5 6 2 7 8 9 10 11 12','1 3 4 5 6 7 2 8 9 10 11 12']}

    nexternal=['4','5','6','7']#,'8']#,'9','10','11','12']
    imodes=['1']
    integrators=['1']#,'2','3','4']
    seeds=['101']#,'102','103','104','105','106','107','108','109','110']
        
    
    # collect all the results.
    cross_section={}
    uncertainty={}
    chi2_per_DoF={}
    number_of_events={}
    number_passing_cuts={}
    total_time={}
    gen_eff={}
    unwgt_eff={}
    tags=[]
    for n in nexternal:
        for io,order in enumerate(orders[n]):
            for iflav,flavour in enumerate(flavours[n]):
                imode='2'
                seed=seeds[0]
                integrator=integrators[0]
                args=list(chain.from_iterable([n,'2',flavour.split()]))
                proc_tag=''
                for el in args:
                    proc_tag=proc_tag+'_'+str(el)
                tag=proc_tag
                for el in order:
                    if (el !=' '):
                        tag=tag+'_'+str(el)
                log_file='./plot_events/res_wgt_'+proc_tag+'/out_unwgt_'+tag+'.txt'
                tag='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+imode+'-i'+integrator+'-s'+seed
                tags.append(tag)
                try:
                    with open(log_file) as file:
                        for line in file:
                            if 'efficiency'in line:
                                unwgt_eff[tag]=float(line.split()[3])
                except:
                    continue

    with open('tables_rwgt.tex','w') as f:
        f.write(latex_header)
        for n in nexternal:
            f.write(table_header)
            for io,order in enumerate(orders[n]):
                for iflav,flavour in enumerate(flavours[n]):
                    for integrator in integrators:
                        for seed in seeds:
                            tag='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+imode+'-i'+integrator+'-s'+seed
                            if tag in tags:
                                string=get_string(tag)
                                f.write(string)
                        f.write('\\midrule \n')
            f.write(table_footer)

        f.write(latex_footer)
