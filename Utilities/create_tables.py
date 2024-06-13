#!/usr/bin/env python

import sys
    
def get_string(tag):
    string=tag
    try:
        string=string+" & %8.3e"%cross_section[tag]
    except:
        string=string+" & --"
    try:
        string=string+" & %8.3e"%uncertainty[tag]
    except:
        string=string+" & --"
    try:
        string=string+" & %7.4f \%%"%(uncertainty[tag]/cross_section[tag] * 100.)
    except:
        string=string+" & --"
    try:
        string=string+" & %2.1f"%(chi2_per_DoF[tag])
    except:
        string=string+" & --"
    try:
        string=string+" & %ik"%(number_of_events[tag]/1000)
    except:
        string=string+" & --"
    try:
        string=string+" & %ik"%(number_passing_cuts[tag]/1000)
    except:
        string=string+" & --"
    try:
        string=string+" & %2.1f \%%"%(number_passing_cuts[tag]/number_of_events[tag] *100.)
    except:
        string=string+" & --"
    try:
        string=string+" & %i "%(total_time[tag])
    except:
        string=string+" & --"
    try:
        string=string+" & %2.3f \%%"%(gen_eff[tag] * 100.)
    except:
        string=string+" & --"
    try:
        string=string+" & %2.3f \%%"%(gen_eff[tag]*number_of_events[tag]/number_passing_cuts[tag] * 100.)
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
\begin{tabularx}{0.95\textwidth}{l r r r r r r r r r r}
  \toprule
  tag&Xsec&unc&rel.unc.&$\chi^2$/&events&events&fraction&time&unw.eff.&unw.eff\\
    &(in pb)&&&D.o.F.&(total)&(pass)&passing&(in s)&(total)&(pass)\\
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

    nexternal=['4','5','6','7','8']#,'9','10','11','12']
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
    for n in nexternal:
        for io,order in enumerate(orders[n]):
            for iflav,flavour in enumerate(flavours[n]):
                for imode in imodes:
                    for integrator in integrators:
                        for seed in seeds:
                            log_file='./OutputsS'+seed+'I'+integrator+'/log_'+n+'_'+imode+'_'+str(io)+'.txt'
                            tag='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+imode+'-i'+integrator+'-s'+seed
                            tags.append(tag)
                            try:
                                with open(log_file) as file:
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
                                        if 'Generation efficiencies:' in line:
                                            gen_eff[tag]=float(line.split()[2])
                            except:
                                continue

    with open('tables.tex','w') as f:
        f.write(latex_header)

        for n in nexternal:
            for io,order in enumerate(orders[n]):
                for iflav,flavour in enumerate(flavours[n]):
                    for imode in imodes:
                        f.write(table_header)
                        for integrator in integrators:
                            for seed in seeds:
                                tag='n'+n+'-o'+str(io)+'-f'+str(iflav)+'-m'+imode+'-i'+integrator+'-s'+seed
                                if tag in tags:
                                    string=get_string(tag)
                                    f.write(string)
                            f.write('\\midrule \n')
                        f.write(table_footer)
                
        f.write(latex_footer)
