import sys
import pysam

bamfile = pysam.AlignmentFile(sys.argv[1],'rb')
foutbamfile = pysam.Samfile( sys.argv[2],'wb',template=bamfile)
bedfile = file(sys.argv[3],'r')
chr,start,end = map(int,sys.argv[4].split(":"))

regions = []

for line in bedfile:
	if line.startswith("#"):
		continue
	arr=line.rstrip().split("\t")
	chrstr=arr[chr]
	rstart=int(arr[start])-600
	rend=int(arr[end])+600
	regions.append([chrstr,rstart,rend])
bedfile.close()

for region in regions:
	chr,start,end = region	
	for segment in bamfile.fetch(chr,start,end):
		if segment.is_supplementary or segment.is_secondary:
			continue
		if segment.is_unmapped:
			continue
		foutbamfile.write(segment)
bamfile.close()
foutbamfile.close()
