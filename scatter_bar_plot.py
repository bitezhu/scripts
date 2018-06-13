import pysam
import sys
import seqio
import re
import math
import copy
import pdb
from matplotlib import pyplot

if len(sys.argv)!=4:
	sys.stderr.write("[INFO]we assume that your bamfile has already been sorted and indexed\n")
	sys.stderr.write("usage:py xxx.py fasta sorted.bamfile interval\n")
	exit(1)
fasta = sys.argv[1]
bamfile = sys.argv[2]
interval = int(sys.argv[3])

gc_dict={}
fa_dict={}
pt = re.compile("[G|C]",re.I)
g = []
for seq in seqio.fasta_read(fasta):
	lentmp = len(seq.seq)
	fa_dict[seq.id]=lentmp
	n = int(math.ceil(lentmp*1.0/interval))
	GClist = []
	for i in xrange(n):
		subseq = seq.seq[i*interval:(i+1)*interval]
		GC = len(pt.findall(str(subseq)))
		gc = GC*1.0/interval
		GClist.append(gc)
	gc_dict[seq.id] = GClist
	
	

dict = {}
#pysam.sort(sys.argv[2], "sorted")
f = file("tmp.txt","w")
samfile = pysam.AlignmentFile(bamfile, "rb" )
for chrID in fa_dict.keys():
	n = int(math.ceil(fa_dict[chrID]*1.0/interval))
	C = []
	for i in xrange(n):
		count = samfile.count(chrID,i*interval,(i+1)*interval,until_eof=True)
		C.append(count)
	dict[chrID] = C


samfile.close()


f.write("chrID\tGCper\tcount\n")
for each in sorted(dict.keys()):
	for idx in range(len(dict[each])):
		f.write(str(each) + "\t" + str(gc_dict[each][idx]) + "\t" + str(dict[each][idx]) + "\n")


f.close()()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

f = file("tmp.txt","r")
X = []
Y = []
for line in f:
	if line.startswith("chr"):continue
	X.append((line.split("\t"))[1])
	Y.append((line.split("\t"))[2])

f.close()

x = np.asarray(X,dtype=np.float)
y = np.asarray(Y,dtype=np.float)
#x = np.log10(x+1)
y = np.log10(y+1)

fig, axScatter = plt.subplots(figsize=(8,6),dpi=300)
# the scatter plot:
axScatter.scatter(x,y,c='m',marker = u'+',alpha=0.01 )
axScatter.set_xlabel("GC percent")
axScatter.set_ylabel("segment count(log)")
axScatter.set_aspect('auto')
#x_density = x/float(np.sum(x))
#y_density = y/float(np.sum(y))
# create new axes on the right and on the top of the current axes
# The first argument of the new_vertical(new_horizontal) method is
# the height (width) of the axes to be created in inches.
divider = make_axes_locatable(axScatter)
axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=axScatter )
#axHistx.set_xticks(minor=False)
axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey=axScatter)
#axHisty.set_xticks(minor=False)
axHistx.spines['right'].set_visible(False)
axHistx.spines['top'].set_visible(False)
axHisty.spines['right'].set_visible(False)
axHisty.spines['top'].set_visible(False)
axHistx.yaxis.set_ticks_position('left')
axHistx.xaxis.set_ticks_position('bottom')
axHisty.yaxis.set_ticks_position('left')
axHisty.xaxis.set_ticks_position('bottom')

# make some labels invisible
plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
         visible=True)

# now determine nice limits by hand:
xbins = np.arange(0,1, 0.01)
ybins = np.arange(0,10, 1)

axHistx.hist(x, bins=xbins,color="g",log=True)
#axHistx.semilogy()
axHisty.hist(y, bins=ybins, color="g",log=True,orientation='horizontal')

# the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
# thus there is no need to manually adjust the xlim and ylim of these
# axis.
#axHistx.axis["bottom"].major_ticklabels.set_visible(False)

for tl in axHistx.get_xticklabels():
    tl.set_visible(False)
#axHistx.set_yticks([0, 30000, 60000])

#axHisty.axis["left"].major_ticklabels.set_visible(False)
for tl in axHisty.get_yticklabels():
    tl.set_visible(False)
#axHisty.set_xticks([0, 150000, 300000])

#plt.draw()
#plt.show()

plt.savefig("scatter.jpg")

	


