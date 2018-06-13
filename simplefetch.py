import sys
import bamio


vcf='/biocluster/data/biobk/user_test/zhusihui/mydata/dbsnp_138.b37.vcf'
vcfidx=bamio.Tabix(vcf,chrpos=0,start=1,end=1)

for iten in vcfidx.fetch('6',160500619,160500622):
	print iten
