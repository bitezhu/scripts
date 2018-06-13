import sys
import pysam
import bamio
import numpy as np
import time
import os

from scipy.stats import fisher_exact
from sklearn import svm
from sklearn import preprocessing 
from sklearn.grid_search import GridSearchCV

vcf='/biocluster/data/biobk/user_test/zhusihui/mydata/dbsnp_138.b37.vcf'
vcfidx=bamio.Tabix(vcf,chrpos=0,start=1,end=1)

def penaltyScore(altbase,queryStr):
	gapflag = 0
	polymerScore = 0
	for i in range(len(queryStr)):
		if queryStr[i] == altbase:
			print i
			if gapflag == 1:
				polymerScore -= 2
			elif gapflag == 1 and i == 0:
				polymerScore -= 6
			else:
				polymerScore -= 1
			gapflag = 1
		else:
			polymerScore += 2
	return polymerScore


class SNPfactory(object):
	def __init__(self,sample,chr,pos,refbase,altbase,totaldepth,basefreqDict,tailDistance,upbases,downbases,basequality,basestr,mapQ,):
		self.sample       = sample
		self.chr          = chr
		self.pos          = int(pos)
		self.ref          = refbase
		self.alt          = altbase
		self.totaldepth   = int(totaldepth)
		self.basefreqDict = basefreqDict
		self.tailDistance = float(tailDistance)
		self.upbases      = upbases
		self.downbases    = downbases
		self.basequality  = basequality
		self.basestr      = basestr
		self.mapQ         = mapQ
		self.SBscore      = 0
		self.transition   = 1
		self.polymerScore = 0

	def TiTv(self,):
		# to determine a mutation is a transversion or transition; transition mutations are generated at higher frequency than transversions
		refbase = self.ref
		altbase = self.alt
		tmpstr = refbase+altbase
		transitionFlag = 1
		if tmpstr in ['AG','GA','CT','TC']:
			transitionFlag = 1
		else:  #['AT','TA','AC','CA','CG','GC','GT','TG']
			transitionFlag = 0 
		self.transition = transitionFlag
		return 0

	
	def calSBscore(self):
		majorFreq = self.basefreqDict[self.ref]
		minorFreq = self.basefreqDict[self.alt]
		self.majorFreq = majorFreq
		self.minorFreq = minorFreq
		# caculate strand bias score based GATK-SB score ,referred in Yan Guo,BMC Genomics 2012
		majordepth,majorfrac,other=majorFreq.split(";")
		majorFor,majorRev = other.split(",")
		self.majordepth = int(majordepth)
		self.majorfrac = float(majorfrac.strip("%"))/100.0

		minordepth,minorfrac,other=minorFreq.split(";")
		minorFor,minorRev = other.split(",")
		self.minordepth = int(minordepth)
		self.minorfrac = float(minorfrac.strip("%"))/100.0

		majorFor=majorFor.strip("+")
		majorRev=majorRev.strip("-")
		minorFor=minorFor.strip("+")
		minorRev=minorRev.strip("-")
		majorFor,majorRev,minorFor,minorRev = map(int,[majorFor,majorRev,minorFor,minorRev])
		commonDom=float(minorFor+minorRev)/(majorFor+majorRev+minorFor+minorRev)
		'''
		S1=minorFor*majorRev/(float(majorFor+minorFor)*(majorRev+minorRev))
		S2=majorFor*minorRev/(float(majorFor+minorFor)*(majorRev+minorRev))
		score=np.max([S1/commonDom,S2/commonDom])
		'''
		S1 = np.abs(minorFor/float(majorFor+minorFor) - minorRev/float(majorRev+minorRev))
		score = S1/commonDom
		score = 1 - fisher_exact(np.asarray([[majorFor,minorFor],[majorRev,minorRev]]))[1]
		self.SBscore = score
		return 0

	def homopolymerPenalty(self):
		queryStr = ''.join(reversed(list(self.upbases)))
		uppolymerScore = 0
		downpolymerScore = 0
		altbase = self.alt
		uppolymerScore = penaltyScore(altbase,queryStr)
		queryStr = self.downbases
		downpolymerScore = penaltyScore(altbase,queryStr)
		self.polymerScore = min(uppolymerScore,downpolymerScore)
		return 0

	def known_novel_anno(self,vcfidx):
		chr = self.chr
		pos = self.pos
		items = vcfidx.fetch(chr,pos,pos+1)
		knownflag = 0
		#print str(pos) + "\t",
		try:
			items.next()
			knownflag = 1
		except StopIteration,e:
			pass
		'''
		items = vcfidx.fetch(chr,pos,pos+1)
		for item in items:
			print  item
		'''
		self.knownflag = knownflag
		return 0

	def sequenceErr(self):
		altbasescore = 0
		altreadscore = 0
		altbaseNum = 0.0
		self.basequality = map(int,self.basequality.split(";"))
		self.mapQ = map(int,self.mapQ.split(";"))
		for i,basestr in enumerate(self.basestr):
			if basestr == self.alt:
				altbasescore += self.basequality[i]
				altreadscore += self.mapQ[i]
				altbaseNum += 1
		self.altbaseAveScore = altbasescore/altbaseNum
		self.altreadAveScore = altreadscore/altbaseNum
		return 0


def generateData(fh,traindata=0):
	backgroudDict = {}
	f=file("/biocluster/data/biobk/user_test/zhusihui/SNP_classification/database_file",'r')
	for line in f:
		if line.startswith("#"):
			continue
		arr=line.rstrip("\n").split("\t")
		chr,start,end,sample,ref,alt,owner,gene,panel,flagval=arr
		backgroudDict["%s-%s-%s-%s-%s"%(arr[0],arr[1],arr[3],arr[4],arr[5])] = arr[-1]
	f.close()
	trainDatasetsX = []
	trainDatasetsY = []
	keystrs = []
	for line in fh:
		if line.startswith("#"):
			continue
		sample,chr,pos,refbase,altbase,totaldepth,Afreq,Cfreq,Gfreq,Tfreq,tailDistance,upbases,downbases,basequality,basestr,mapQ,owner = line.rstrip("\n").split("\t")
		basefreqDict = {'A':Afreq,'C':Cfreq,'G':Gfreq,'T':Tfreq}
		SNPobj = SNPfactory(sample,chr,pos,refbase,altbase,totaldepth,basefreqDict,tailDistance,upbases,downbases,basequality,basestr,mapQ)
		SNPobj.TiTv()
		SNPobj.calSBscore()
		SNPobj.homopolymerPenalty()
		SNPobj.known_novel_anno(vcfidx)
		SNPobj.sequenceErr()
		keystr = "\t".join([SNPobj.sample,SNPobj.chr,str(SNPobj.pos+1),SNPobj.ref,SNPobj.alt])

		outstr = [SNPobj.sample,SNPobj.chr,SNPobj.pos,SNPobj.ref,SNPobj.alt,SNPobj.transition ,SNPobj.majorfrac,SNPobj.minorfrac,SNPobj.polymerScore ,SNPobj.knownflag,SNPobj.altbaseAveScore,SNPobj.altreadAveScore]
		#if "%s-%s-%s-%s-%s"%(SNPobj.chr,str(SNPobj.pos+1),SNPobj.sample,SNPobj.ref,SNPobj.alt) not in backgroudDict:
		#	continue
		if not traindata:
			trainDatasetsX.append(map(float,[SNPobj.transition,SNPobj.majorfrac,SNPobj.minorfrac,SNPobj.polymerScore,SNPobj.knownflag,SNPobj.altbaseAveScore,SNPobj.altreadAveScore,SNPobj.SBscore]))
			trainDatasetsY.append(0)
			keystrs.append(keystr)
		else:
			try :
				if backgroudDict["%s-%s-%s-%s-%s"%(SNPobj.chr,str(SNPobj.pos+1),SNPobj.sample,SNPobj.ref,SNPobj.alt)]=="True":
					trainDatasetsY.append(1)
				else:
					trainDatasetsY.append(0)
				trainDatasetsX.append(map(float,[SNPobj.transition,SNPobj.majorfrac,SNPobj.minorfrac,SNPobj.polymerScore,SNPobj.knownflag,SNPobj.altbaseAveScore,SNPobj.altreadAveScore,SNPobj.SBscore]))
			except KeyError,e:
				pass
	return trainDatasetsX,trainDatasetsY,keystrs

def fmtresult(rawTabfile,keystrs,predictArr):
	fh = file(rawTabfile,'r')
	fout = file(rawTabfile.rsplit(".",1)[0]+"_classed_pos.xls",'w')
	fouttmp = file(rawTabfile.rsplit(".",1)[0]+"_classed_neg.xls",'w')
	for line in fh:
		if line.startswith("#"):
			fout.write(line)
			continue
		arr=line.rstrip("\n").split("\t")
		if arr[8] != "SNP":
			fout.write(line)
			continue
		#chr,start,samplename,ref,alt = arr[4],arr[5],arr[0],arr[9],arr[10]
		subarr = [arr[0],arr[4],arr[5],arr[9],arr[10]]
		if predictArr[keystrs.index("\t".join(subarr))]:
			fout.write(line)
		else:
			fouttmp.write(line)
	fout.close()
	fouttmp.close()
	return 0

def fmtstr(x):
	return "%.3f"%x

def scaleFeature(mat,colidx,vmax,vmin,traindata=0):
	# scale to [-1,1]
	if traindata:
		vmin = np.min(mat[:,colidx])
		vmax = np.max(mat[:,colidx])
	mat[:,colidx] = -1.0 + (1.0 - -1.0) * (mat[:,colidx] - vmin)/(vmax-vmin)
	return mat,vmax,vmin

def __main():
	start_time = time.time()
	ret = 0
	if len(sys.argv) != 4:
		sys.stderr.write("usage:\npython SNPfilter.py sampleCollection.xls test_predict_SNP_data test.annotate.meaningful.xls\n")
		sys.exit(0)
	fh=file(sys.argv[1],'r')
	trainDatasetsX,trainDatasetsY,keystrs = generateData(fh,traindata=1)
	fh.close()
	fh=file(sys.argv[2],'r')
	predictDatasetsX,predictDatasetsY,keystrs = generateData(fh)
	fh.close()
	### classification 
	trainDatasetsX = np.asarray(trainDatasetsX)
	predictDatasetsX = np.asarray(predictDatasetsX)
	
	clf = svm.SVC(kernel='linear',C=1000.,)#class_weight={0: 10, 1: 10, 2: 10, 3: 30, 4: 10, 5: 30, 6: 30, 7: 30})
	
	if len(predictDatasetsX) == 0:
		sys.stderr.write("empty array found here,please check your input data first!!!\nThis program can only be applied for SNP F/T determenation for now, so it'll ignore InDel record.\n ")
		ret = 1
	else:
		### scale data,only for basequal and mappingQual
		# ref http://suanfazu.com/t/is-scaling-of-feature-values-in-libsvm-necessary/2030/2
		'''	
		trainDatasetsX,vmax,vmin = scaleFeature(trainDatasetsX,-3,0,0,traindata=1)
		predictDatasetsX,vmax,vmin = scaleFeature(predictDatasetsX,-3,vmax,vmin,traindata=0)
		trainDatasetsX,vmax,vmin = scaleFeature(trainDatasetsX,-2,0,0,traindata=1)
		predictDatasetsX,vmax,vmin = scaleFeature(predictDatasetsX,-2,vmax,vmin,traindata=0)
		'''
		min_max_scaler = preprocessing.MinMaxScaler()
		trainDatasetsX = min_max_scaler.fit_transform(trainDatasetsX)
		predictDatasetsX = min_max_scaler.fit_transform(predictDatasetsX)
		clf.fit(trainDatasetsX, trainDatasetsY)
		predictArr = clf.predict(predictDatasetsX)
		f=file('%s_log.xls'%(os.path.basename(sys.argv[3]).split(".",1)[0]),'w')
		f.write("#TiTv\tmajorF\tminorF\tpolymer\tknownFlag\tbaseQ\tmappingQ\tfishScore\n")
		for i in range(len(predictArr)):
			f.write("\t".join(map(fmtstr,predictDatasetsX[i,:].tolist()))+'\t'+str(predictArr[i])+"\n")
		f.close()
		fmtresult(sys.argv[3],keystrs,predictArr)
	end_time = time.time()
	sys.stderr.write("Task cost %ds\n"%(end_time-start_time)) 
	return ret

if __name__ == "__main__":
	ret= __main()
	if ret != 0:
		sys.stderr.write("Task Failed !!!\n")
	else:
		sys.stderr.write("Task Done ! \n ") 
