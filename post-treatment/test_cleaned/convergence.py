import numpy as numpy
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sp
import array
import matplotlib.ticker as ticker

from matplotlib import gridspec

plt.rc('text', usetex=True)
font = {'family' : 'sans-serif', 'sans-serif' : 'Helvetica',
          'weight' : 'bold',
          'size'   : 15.}
plt.rc('font', **font )
 
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')

mpl.rcParams['grid.linestyle'] = ':'

gs = gridspec.GridSpec(1, 1,
                       #width_ratios=[1,1]
                       )


fig = plt.figure(figsize=(4.8,2.8),dpi=500)

p1 = fig.add_subplot(gs[0])



f1 = open('test_clean.txt_RHS_0.txt', 'r')
f2 = open('test_clean.txt_RHS_1.txt', 'r')

val1=[]
val2=[]


for line in f1:
    val1.append(float(line.split()[0]))
for line1 in f2:
    val2.append(float(line1.split()[0]))

p1.set_yscale('log')
p1.yaxis.grid()
p1.xaxis.grid()
X = numpy.linspace(0, 3800, 750, endpoint=True)
Y = numpy.linspace(1e-10, 1e-10, 500, endpoint=True)
p1.plot(val1, label='1st Rhs',lw='1.',marker='', markevery=20,markersize='3',mew=1.3)
p1.plot(val2, label='2nd Rhs', linestyle='-.')

#plt.plot(val4, label='GMRES without preconditioning', lw='2', color='b')

#p1.plot(X,Y, '--', lw='.6', color='r')
y1 = [1e-8,1e-4,1,1e4,1e8,1e12,1e16,1e20]
ylabels1 = ['1e-8','1e-4','1','1e4','1e8','1e12','1e16','1e20']
p1.set_xlabel('Iteration step number of Classic Block GMRES with 2 Rhs', size='10')
p1.set_ylabel('Relative Residual', size='10')
formatter = ticker.ScalarFormatter()
formatter.set_scientific(True)
#formatter.set_powerlimits((0,0))
p1.yaxis.set_major_formatter(formatter)
p1.set_xticks(numpy.arange(0,1080,120))
#p1.set_yticks(y1,ylabels1)
p1.set_yticks(y1)
p1.set_yticklabels(ylabels1, rotation=0)
#xlabels = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20']
#p1.set_xticklabels(xlabels)
#p1.set_title(r'$\textbf{matBlock}$', size='10')
p1.tick_params(axis='y', labelsize=7.8,which='both')
p1.tick_params(axis='x', labelsize=7.8)
p1.set_xlim(0,1080)
p1.set_ylim(1e-6,1e8)

p1.legend(loc='upper right',
        #bbox_to_anchor=(1.1, -1.595), 
        prop={'size':7},ncol=2,frameon=True)

###########################

#plt.tight_layout(pad=.8,h_pad=.3)
plt.tight_layout(w_pad=.8, h_pad=0.2)

plt.savefig("test_convergence.eps",dpi=500)
#plt.show()