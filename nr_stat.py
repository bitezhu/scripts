#-*- coding: UTF-8 -*-
#!/usr/bin/python
import sys
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
def plot_pie(labels,fracs,figname):
	###设置图布长宽比为8:6
	fig = plt.figure(figsize=(8,6),dpi=300)
	###设置图形占整个画图的区域是111
	ax5 = fig.add_subplot(111)
	###从cm中选取颜色，否则默认颜色很丑
	colors = cm.Accent(np.linspace(0, 1, len(labels)))
	###圆里面的文本格式，%3.1f%%表示小数有三位，整数有一位的浮点数
	autopct='%1.1f%%'
	###给每个labels加一个数字标识，以免有些labels太长，图例和图形上都不好展示
	sign=range(1, len(labels)+1)
	###画图 startangle表示从第一块从哪个角度开始逆时针旋转，startangle=0表示第一块从90度开始逆时针旋转
	patches, texts = ax5.pie(fracs,labels=sign,labeldistance=1.1,shadow=False,colors=colors,startangle=0)
	#for t in texts2:
	#	t.set_size = 300
	tmplabels = []
	total = sum(fracs)
	a=0
	###得到图例legend展示所需的labels即tmplabels
	for i in xrange(len(labels)):
		a+=1
		lable = labels[i]
		size = float(fracs[i])/total*100
		tmplabels.append(str(a)+": "+lable+"  %1.1f%%"%size)
	###图例的绘制
	###bbox_to_anchor表示图例位置；fontsize表示图例大小；frameon表示图例最外框线
	legend = ax5.legend(patches,tmplabels,bbox_to_anchor=(1.6, 1.1),borderaxespad=0,fontsize=7,frameon=True,shadow=True,fancybox=True,)
	frame = legend.get_frame()
	frame.set_facecolor('0.90')
	#frame = legend.get_frame() 
	#frame.set_alpha(1) 
	#frame.set_facecolor('black')
	#frame.set_edgecolor('red')

	###设置图形位置
	pie = ax5.get_position()
	###pie.width表示图形占宽从左边开始的0.7，这是为了更好的展示出图例，否则有时图例会与图形重合
	ax5.set_position([pie.x0, pie.y0, pie.width*0.7, pie.height])
	###设置图形中各个部分之前的线条宽度及颜色，白色会好看些，默认为黑色
	for w in patches:
		w.set_linewidth(0.2)
		w.set_edgecolor('white')
	###为了保持圆形为正规圆形，如果不设置，当figsize不成比例，图形大小也会改变
	plt.axis('equal')
	plt.savefig(figname+".png",format='png',dpi=300)
	plt.savefig(figname+".svg",format='svg',dpi=300)
	plt.close(0)
	return(0)
identity_a = [0,0,0,0,0,0]
evalue_a   = [0,0,0,0,0,0]
identity_l = ["0~40%","40%~70%","70%~80%","80%~88%","88%~95%%","95%~100%"]
evalue_l   = ["0~1e-100","1e-100~1e-80","1e-80~1e-60","1e-60~1e-40","1e-40~1e-20","1e-20~"]
hs = {}
f = file(sys.argv[1],"r")
for line in f:
	if line.startswith("#"):continue
	geneid,nranno,identity,o1,o2,o3,o4,o5,o6,o7,evalue,genename= line.strip("\n").split("\t")
	if ":" in genename:
		specie=" ".join(genename.split(": ")[1].split(" ")[0:2])
	if ":" not in genename:
		specie=" ".join(genename.split(" ")[0:2])
	try:
		identity = float(identity)
	except:
		continue
	evalue = float(evalue)
	if identity < 40:identity_a[0] += 1
	elif identity < 70:identity_a[1] += 1
	elif identity < 80:identity_a[2] += 1
	elif identity < 88:identity_a[3] += 1
	elif identity < 95:identity_a[4] += 1
	elif identity <=101:identity_a[5] += 1
	else:sys.stderr.write("[WARN] identity is '%f'\n"%identity)
	if evalue < 1e-100:evalue_a[0] += 1
	elif evalue < 1e-80:evalue_a[1] += 1
	elif evalue < 1e-60:evalue_a[2] += 1
	elif evalue < 1e-40:evalue_a[3] += 1
	elif evalue < 1e-20:evalue_a[4] += 1
	elif evalue < 1:evalue_a[5] += 1
	else:sys.stderr.write("[WARN] evalue is '%f'\n"%evalue)
	if specie in hs:
		hs[specie] += 1
	else:
		hs[specie] = 1
f.close()

f_NR_dat = file("NR_identity_pie.dat","w")
f_NR_dat.write("\t".join(identity_l)+"\n")
f_NR_dat.write("\t".join(map(str,identity_a))+"\n")
f_NR_dat.close()
plot_pie(identity_l,identity_a,"NR_identity")

f_NR_dat = file("NR_evalue_pie.dat","w")
f_NR_dat.write("\t".join(evalue_l)+"\n")
f_NR_dat.write("\t".join(map(str,evalue_a))+"\n")
f_NR_dat.close()
plot_pie(evalue_l,evalue_a,"NR_evalue")

tmpdata = []
tmplabel = []
for specie in hs:
	if specie.strip() == "-":continue
	tmpdata.append(hs[specie])
	tmplabel.append(specie)
f_NR_dat = file("NR_species_pie.dat","w")
f_NR_dat.write("\t".join(tmplabel)+"\n")
f_NR_dat.write("\t".join(map(str,tmpdata))+"\n")
f_NR_dat.close()
plot_pie(tmplabel,tmpdata,"NR_species")
