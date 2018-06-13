import sys
import pysam

'''
use pysam module to fetch fragment sequnece from a fasta file
'''

gtf=file(sys.argv[1])

faidx = pysam.FastaFile(sys.argv[2])

chridx,sidx,eidx=map(int,sys.argv[3].split(','))

for line in gtf:
	if line.startswith("#"):
		continue
	arr=line.rstrip().split("\t")
	chr=arr[chridx]
	s=int(arr[sidx])
	e=int(arr[eidx])
	seqID=arr[0]
	try:
		print ">%s"%seqID
		print faidx.fetch(chr,s-1,e)
	except:
		pass
gtf.close()


