import sys
import gzip
import numpy as np
fh = gzip.GzipFile(sys.argv[1],)
idx = sorted(np.random.randint(1000000,size=int(sys.argv[2])))
i=0
while 1:
    seqID = fh.readline();sequence = fh.readline();	mark = fh.readline();qual = fh.readline()
    i+=1
    for randomnum in idx:
        if i==randomnum:
	    print ">%s"%seqID,
	    print sequence,
    if i > 1000000:break
fh.close()
	
