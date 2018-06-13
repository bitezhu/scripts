import sys
import pysam


if len(sys.argv) != 3:
	sys.stderr.write("usage:\npython XX.py test.bam test.annotate.meaningful.xls\n")
	exit(0)

bamfile=pysam.AlignmentFile(sys.argv[1])
fh=file(sys.argv[2])

def formatDict(tmpdict,totalDepth):
	outstr=[]
	for basestr in ['A','C','G','T']:
		baseformat=''
		basedepth,[revDepth,forDepth] = tmpdict[basestr]
		depthRatio = float(basedepth)/totalDepth*100
		formatstr=';'.join(map(str,[basedepth,'%.0f%%'%depthRatio,'%d+,%d-'%(forDepth,revDepth)]))
		outstr.append(formatstr)
	return '\t'.join(outstr)

def consensusGenerate(context):
	ret=[]
	for element in context:
		if len(element)==0:continue
		if len(element) == 6:
			ret.append(element)
		else:
			for tmpelement in ret:
				if element in tmpelement:
					break
				else:
					ret.append(element)
	return ret

for line in fh:
	if line.startswith("#"):
		continue
	arr=line.rstrip().split("\t")
	if arr[8] != "SNP":
		continue
	chr,start,end,samplename,ref,alt,owner = [arr[4],arr[5],arr[6],arr[0],arr[9],arr[10],'smone']
	
	#chr,start,end,samplename,ref,alt,owner,other=line.rstrip().split("\t",7)
	start,end=map(int,[start,end])
	pileupcolumns=bamfile.pileup(chr,start-1,start,truncate=True)
	for pileupcolumn in pileupcolumns:
		print samplename+'\t'+pileupcolumn.reference_name + '\t' + str(pileupcolumn.pos) + '\t'+ ref+'\t'+alt +'\t'+ str(pileupcolumn.n) + "\t",
		
		baseQarr=[]
		baseChararr=[]
		mappingOrient=[]	
		upstreamContext=[]
		downstreamContext=[]
		borderDistances = Readtotallen=0
		readMapq=[]	
		for pileupread in pileupcolumn.pileups:
			readalignment = pileupread.alignment
			#print line
			#if not pileupread.query_position:
			if pileupread.is_del or pileupread.is_refskip:
				#print readalignment
				continue
			baseChararr.append(readalignment.query_sequence[pileupread.query_position])
			baseQarr.append(readalignment.query_qualities[pileupread.query_position])
			readMapq.append(readalignment.mapping_quality)
			qlen=len(readalignment.query_sequence)
			borderDistance = qlen-pileupread.query_position if pileupread.query_position > 0.5*qlen else pileupread.query_position 
			#print str(qlen) + '\t'+str(pileupread.query_position)+ '\t' + str(borderDistance)
			borderDistances+=borderDistance
			Readtotallen+=qlen
			if readalignment.is_reverse:
				mappingOrient.append(0)
			else:
				mappingOrient.append(1)
			'''
			try:
				upstreamContext.append(readalignment.query_sequence[pileupread.query_position-10:pileupread.query_position])
				downstreamContext.append(readalignment.query_sequence[pileupread.query_position+1:pileupread.query_position+10])
			except:
				continue
			'''
			#print readalignment.get_tag('NM')
			if qlen==readalignment.query_alignment_length and readalignment.get_tag('NM')==0:
				#print readalignment
				#print str(qlen) + '\t'+str(readalignment.query_alignment_length)
				upstreamContext.append(readalignment.query_sequence[pileupread.query_position-6:pileupread.query_position])
				downstreamContext.append(readalignment.query_sequence[pileupread.query_position+1:pileupread.query_position+7])

		tmpdict={'A':[0,[0,0]],'C':[0,[0,0]],'G':[0,[0,0]],'T':[0,[0,0]]}
		for i,tmpstr in  enumerate(baseChararr):
			strandflag=mappingOrient[i]
			if tmpstr=="A":
				tmpdict[tmpstr][0]+=1
				tmpdict[tmpstr][1][strandflag]+=1
			elif tmpstr=='T':
				tmpdict[tmpstr][0]+=1
				tmpdict[tmpstr][1][strandflag]+=1
			elif tmpstr=='G':
				tmpdict[tmpstr][0]+=1
				tmpdict[tmpstr][1][strandflag]+=1
			elif tmpstr=='C':
				tmpdict[tmpstr][0]+=1
				tmpdict[tmpstr][1][strandflag]+=1
			else:
				sys.stderr.write("[Warning]There is an ambiguous Base here!! Plz check %s in %dth read\n"%(tmpstr,i))
				continue

		#print ''.join(map(str,baseQarr)) + "\t" ''.join(baseChararr)

		print formatDict(tmpdict,pileupcolumn.n) + '\t',
		print '%.2f'%(float(borderDistances)/Readtotallen)+'\t',
		upstreamContext=consensusGenerate(upstreamContext)
		downstreamContext=consensusGenerate(downstreamContext)
		print ','.join(list(set(upstreamContext))) + "\t"+ ','.join(list(set(downstreamContext))) + "\t",
		print ';'.join(map(str,baseQarr))+'\t'+''.join(baseChararr)+"\t",
		print ";".join(map(str,readMapq)),
		print "\t"+owner
bamfile.close()
fh.close()


