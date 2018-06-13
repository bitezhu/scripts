import sys
import glob
import time
import os
import pwd


panellist=['panel18','panel68','panle88','panel203','panel509']
resultDict={}
for panel in panellist:
	targetfiles=glob.glob('/leostore/analysis/clinical/%s/LAA*/'%(panel))
	for targetfile in targetfiles:
		panel=targetfile.rstrip("\/").split("/")[4]
		#print os.path.split(targetfile)
		Mon=time.ctime(os.path.getctime(targetfile)).split(" ")
		Year=Mon[-1]
		if Year !='2017':
			continue
		Mon=Mon[1]
		if Mon not in ['Jun','Jul',]:#,'Aug']: # 'May'
			continue
		#print targetfile
		
		try:
			unchangedFile=glob.glob(targetfile+"/02_var/vep/*.annotate.meaningful.xls")[0]
			resultFile=glob.glob(targetfile+"/03_report/qc_report/*.annotate.meaningful.xls")[0]
			tumorBam=glob.glob(targetfile+"/01_aln/*T1_recal.bam")[0]
			docx=glob.glob(targetfile+"/03_report/qc_report/*docx")[0]
		except:
			continue
		if not os.path.exists(unchangedFile):
			continue
		tmpdict={}
		tmpfh=file("SNP_pos.txt",'w')
		f1=file(unchangedFile,'r')
		f2=file(resultFile,'r')
		owner=pwd.getpwuid(os.stat(docx).st_uid).pw_name
		for line in f1:
			if line.startswith("#"):
				continue
			arr=line.rstrip("\n").split("\t")
			if arr[8] != 'SNP':
				continue
			#keyval="\t".join([arr[4],arr[5],arr[6], arr[0],arr[9],arr[10],owner,arr[1],panel])
			keyval = line.rstrip()
			resultDict[keyval]='False'
			tmpdict[keyval]='False'
		f1.close()
		for line in f2:
			if line.startswith("#"):
				continue
			arr=line.rstrip("\n").split("\t")
			if arr[8] != 'SNP':
				continue
			#keyval="\t".join([arr[4],arr[5],arr[6], arr[0],arr[9],arr[10],owner,arr[1],panel])
			keyval=line.rstrip()
			if keyval in resultDict.keys():
				resultDict[keyval]='True'
				tmpdict[keyval]='True'
		f2.close()
		print unchangedFile + "\t" + owner
		for k,v in tmpdict.items():
			tmpfh.write(k+'\t'+v+'\n')
		tmpfh.close()
		os.system("python  /home/zhusihui/scripts/sam2pileup.py %s SNP_pos.txt  >> sampleCollection.xls "%tumorBam)
'''
for k,v in resultDict.items():
	print k + "\t" + v
	
'''
