#!/usr/bin/env python


# Some useful bash one-liners:
# from main directory:
# for dir in OutputsS*/ ; do cd $dir ; for file in *.lhe.rwgt ; do ../Utilities/compute_unweighting_efficiency/compute_unw_eff $file ; done ; cd .. ; done
# for dir in OutputsS*/ ; do cd $dir ; for file in *.lhe.rwgt ; do ../Utilities/plot_events/plot_events $file ; done ; cd .. ; done

# from with HwU_plots:
# for file in ../OutputsS101I1/*.HwU ; do ../Utilities/plot_events/internal/histograms.py ../OutputsS1*I1/${file:17} --average --out=${file:17:-4} ; done


import sys
import math
import os.path
    
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
        string=string+" & %7.2f \%%"%(uncertainty[tag]/cross_section[tag] * 100.)
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

def convert_proc_to_string(flavour):
    string=''
    for i,iflav in enumerate(flavour.split()):
        if i==0:
            string+=pdg_to_str(iflav)+' '
        elif i==1:
            string+=pdg_to_str(iflav)+' -> '
        else:
            if iflav == '21' :
                string+=' '+str(len(flavour.split())-i)+pdg_to_str(iflav)
                return string
            else:
                string+=pdg_to_str(iflav)
    return string

def pdg_to_str(ipdg):
    string=''
    if ipdg == '21' : string='g'
    if ipdg == '1' : string='d'
    if ipdg == '-1' : string='~d{.55-}'
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
    return ave_tag

if __name__ == '__main__':
    table_header=r"""\noindent The tag has the form \textbf{nX-oX-fX-mX-iX-sXXX}, where X is an integer. The integer after 
\textbf{n} corresponds to the number of external particles;
\textbf{f} flavour index: 0 for $gg \to (n-2)g$, 1 for $gg\to d\bar{d}+(n-4)g$, 2 for $\bar{d} g\to \bar{d}+(n-3)g$ , 3 for $\bar{d}d \to (n-2)g$ );
\textbf{o} the colour order index: minimum number of gluons between the two incoming particles in the colour order;
\textbf{m} corresponds to imode: 0 for grid setup, 1 for upper bounding envelope estimation, 2 for event generation;
\textbf{i} the integrator: 1 for gen23, 2 for haag, 3 for genpt (chili-like), 4 for t-channel;
\textbf{s} the random number seed.

\begin{tabularx}{0.95\textwidth}{l r r r r r r r r r r}
  \toprule
  tag&Xsec&unc&rel.unc.&$\chi^2$/&events&events&fraction&time&unw.eff.&unw.eff\\
    &(in pb)&&&D.o.F.&(total)&(pass)&passing&(in s)&(total)&(pass)\\
  \midrule
"""
    table_footer=r"""\end{tabularx}
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

    gnuplot_header=r"""reset
set lmargin 5
set rmargin 0
set bmargin 0.5
set tmargin 0.5
set terminal postscript portrait enhanced color "Helvetica" 10
set output "plots.ps"
set linetype 1 lc rgb '#a9e5bb' lw 3
set linetype 2 lc rgb '#f7b32b' lw 3
set linetype 3 lc rgb '#f72c25' lw 3
set linetype 4 lc rgb '#533756' lw 3
set linetype 5 lc rgb '#a4d3cf' lw 3
set linetype 6 lc rgb '#89b487' lw 3
set linetype 7 lc rgb '#7c7056' lw 3
unset xtics
"""
    gnuplot_footer=r"""unset multiplot
!ps2pdf "plots.ps" &> /dev/null
"""
    gnuplot_middle="""unset multiplot
set multiplot layout 3,3
set xrange [-2:51]
set size 0.4,0.5
set origin 0,0.5
set label 'Cross section [pb]' at graph 0.1, graph 0.97
plot 'event_numbers.dat' i %s u 0:($2):3:(int($0/13)+1) w yerrorbar pointtype 7 pointsize 0 linecolor variable notitle
unset label
set logscale y
set size 0.4,0.25
set origin 0,0.25
set label 'Unweighting eff.' at graph 0.1, graph 0.95
plot 'event_numbers.dat' i %s u 0:($4):(int($0/13)+1) w points pointtype 7 linecolor variable notitle
unset label
set size 0.25,0.25
set origin 0.0,0.0
set label '2nd. Unw. eff. (LC->full)' at graph 0.1, graph 0.95
set xrange [-2:12]
plot 'event_numbers.dat' i %s u 0:($5):(int($0/13)+5) w points pointtype 7 linecolor variable notitle
unset label
set size 0.25,0.25
set origin 0.3,0.0
set label '2nd. Unw. eff. (LC->NLC)' at graph 0.1, graph 0.95
set xrange [-2:12]
plot 'event_numbers.dat' i %s u 0:($6):(int($0/13)+6) w points pointtype 7 linecolor variable notitle
unset label
set size 0.25,0.25
set origin 0.6,0.0
set label '2nd. Unw. eff. (NLC->full)' at graph 0.1, graph 0.95
set xrange [-2:12]
plot 'event_numbers.dat' i %s u 0:($7):(int($0/13)+7) w points pointtype 7 linecolor variable notitle
unset label
unset logscale y
"""
    gnuplot_middle2="""set xrange[0:3]
set size 0.5,0.5
set origin 0.5,0.4
set label 'weight distribution' at graph 0.1, graph 0.95
set xtics
set logscale y
set yrange [5e-7:1.1]
set format y '10^{%%T}'
plot '%s' i 0 u ($1+$2)/2:3 w histep linetype 5 t 'LC->full',\
     '%s' i 1 u ($1+$2)/2:3 w histep linetype 6 t 'LC->NLC',\
     '%s' i 2 u ($1+$2)/2:3 w histep linetype 7 t 'NLC->full'
unset xtics
unset logscale y
unset yrange
unset format y
unset label
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

    nexternal=['4','5','6','7','8','9','10']
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
    sec_unw_eff_LF={}
    sec_unw_eff_LN={}
    sec_unw_eff_NF={}
    tags=[]
    for n in nexternal:
        for iflav,flavour in enumerate(flavours[n]):
            for io,order in enumerate(orders[iflav][n]):
                for imode in imodes:
                    for integrator in integrators:
                        for seed in seeds:
                            log_file='./OutputsS'+seed+'I'+integrator+'/log_'+n+'_'+imode+'_'+str(io)+'_'+str(iflav)+'.txt'
                            try:
                                with open(log_file) as file:
                                    tag='n'+n+'-f'+str(iflav)+'-o'+str(io)+'-m'+imode+'-i'+integrator+'-s'+seed
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
                                    if imode=='2':
                                        try:
                                            log_file='./OutputsS'+seed+'I'+integrator+'/out_unwgt__'+n+'_2_'+'_'.join(flavour.split())+'_'+'_'.join(order.split())+'.txt'
                                            with open(log_file) as file:
                                                for line in file:
                                                    if ' LC->full' in line:
                                                        sec_unw_eff_LF[tag]=float(line.split()[3])
                                                    if ' LC->NLC' in line:
                                                        sec_unw_eff_LN[tag]=float(line.split()[3])
                                                    if ' NLC->full' in line:
                                                        sec_unw_eff_NF[tag]=float(line.split()[3])
                                        except:
                                            sec_unw_eff_LF[tag]=1.
                                            sec_unw_eff_LN[tag]=1.
                                            sec_unw_eff_NF[tag]=1.
                            except:
                                pass

    with open('tables.tex','w') as f:
        f.write(latex_header)

        for n in nexternal:
            for iflav,flavour in enumerate(flavours[n]):
                for io,order in enumerate(orders[iflav][n]):
                    for imode in imodes:
                        f.write(table_header)
                        for i,integrator in enumerate(integrators):
                            tags_to_average=[]
                            for seed in seeds:
                                tag='n'+n+'-f'+str(iflav)+'-o'+str(io)+'-m'+imode+'-i'+integrator+'-s'+seed
                                if tag in tags:
                                    string=get_string(tag)
                                    f.write(string)
                                    tags_to_average.append(tag)
                            f.write('\\midrule \n')
                            tag=compute_averages(tags_to_average)
                            string=get_string(tag)
                            f.write(string)
                            f.write('\\bottomrule \n')
                        f.write(table_footer)
        f.write(latex_footer)


    with open('event_numbers.gnuplot','w') as f1:
        with open('event_numbers.dat','w') as f2:
            f1.write(gnuplot_header)
            icount=0
            for n in nexternal:
                for iflav,flavour in enumerate(flavours[n]):
                    for io,order in enumerate(orders[iflav][n]):
                        to_write=[]
                        for i,integrator in enumerate(integrators):
                            for seed in seeds:
                                line=''
                                for imode in ['1','2']:
                                    tag='n'+n+'-f'+str(iflav)+'-o'+str(io)+'-m'+imode+'-i'+integrator+'-s'+seed
                                    if tag in tags:
                                        if imode=='1':
                                            line+=tag+'  '+str(cross_section[tag])+'  '+str(uncertainty[tag])
                                        elif imode=='2':
                                            line+='  '+str(gen_eff[tag]*number_of_events[tag]/number_passing_cuts[tag])+'  '+str(sec_unw_eff_LF[tag])+'  '+str(sec_unw_eff_LN[tag])+'  '+str(sec_unw_eff_NF[tag])
                                if line: to_write.append(line)
                            to_write.append('""')
                            to_write.append('""')
                            to_write.append('""')
                        to_write.append('')
                        to_write.append('')
                        to_write.append('')
                        if (len(to_write)==55):
                            f2.write('\n'.join(to_write))
                            string='set label "'+convert_proc_to_string(flavour)+', '+''.join(order)+'" at graph 0,graph 1.05\n'
                            f1.write(string)
                            f1.write(gnuplot_middle%(icount,icount,icount,icount,icount))
                            icount=icount+1

                            HwU_file='HwU_plots/events_'+n+'_2_'+'_'.join(flavour.split())+'_'+'_'.join(order.split())+'.HwU'
                            if os.path.isfile(HwU_file):
                                f1.write(gnuplot_middle2%(HwU_file,HwU_file,HwU_file))

                                
            f1.write(gnuplot_footer)            

