import pyzdde.zdde as pyz
import glob
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.optimize
from scipy.optimize import leastsq
from timeit import default_timer as timer
import re

fdirs = [r"D:\Temp\A+TMAG_FINE\TMAGS"]
ddir = r"D:\Temp\A+TMAG_FINE"
customfile = 'DENC_F4_WLFINALTEST_0p1wv.txt'
AOIs = [0,15,30,45]
FNUMs = [1,2,3,4]
WLs = [0.633]

getfiles = 1
cuton = 0
cutoff = 0

zdir = r"D:\AOI+Fnum_Test"


def pmSPO(geospo_list,RMSspo_list,RMSspoX_list,RMSspoY_list):
	# Matrix print assist, SPO
	for i in range(len(geospo_list)):
		print(geospo_list[i],end='')
		print(' ',end='')
		print(RMSspo_list[i],end='')
		print(' ',end='')
		print(RMSspoX_list[i],end='')
		print(' ',end='')
		print(RMSspoY_list[i])
	return
def pm(M):
	# Matrix print assist
	tmpstr=''
	for i in M:
		print(i)
	return
	
def zStartLink():
	try:
		ln
	except(NameError):
		ln = pyz.createLink() # DDE link object
	return ln


if __name__ == "__main__":
	nstatus = 0
	filelist = glob.glob(zdir + '\\*.zmx')
	if getfiles == 1:
		ln = zStartLink()

	for nf,file in enumerate(filelist):
		if getfiles == 1:
			ln.zLoadFile(file)
			ln.zSetConfig(2)
		for nfdir,fdir in enumerate(fdirs):
			datfiles = glob.glob(fdir + '\\**\\*.dat',recursive=True)
			for nAOI,AOI in enumerate(AOIs):
				if len(datfiles)==0:
					raise('No Files Found')
				
				wfms = []
				if cutoff == 0:
					cutoff = len(datfiles)
				
				##norm = colors.Normalize(vmin=0, vmax=len(datfiles), clip=True)
				norm = colors.Normalize(vmin=0, vmax=cutoff-cuton, clip=True)
				mapper = cm.ScalarMappable(norm=norm, cmap=cm.jet)
				
				fig, ax = plt.subplots(nrows=1,ncols=2)
				meanT_list = []
				meanS_list = []
				minT_list = []
				maxT_list = []
				minS_list = []
				maxS_list = []
				error_list = []
				RMSspo_list = []
				geospo_list = []
				enc_list = []
				RMSspoX_list = []
				RMSspoY_list = []
				
				print(file.split('\\')[-1])
				if getfiles==1:
					ln.zSetSurfaceParameter(2,3,AOI)
					ln.zPushLens()
				for nd,datfile in enumerate(datfiles[cuton:cutoff]):
					#nstatus = (nf+1)*(nd+1)+nd
					#print("----- "+str(nstatus)+"/"+str(len(filelist)*len(datfiles))+" -----")
					print("----- "+str(nstatus+1)+"/"+str(len(filelist)*len(datfiles[cuton:cutoff])*len(fdirs)*len(AOIs)*len(FNUMs))+" -----",end='')
					if nd == 0:
						t = timer()
					else:
						dt = (timer()-t)/(nstatus) * (len(datfiles[cuton:cutoff])*len(filelist)*len(fdirs)*len(AOIs)*len(FNUMs)-nstatus)
						print(' ETA: '+str(int(dt/60))+'m'+str(int(dt%60))+'s -----')
					nstatus+=1
					#if nstatus < 3504:
					#	print('SKIPPING...')
					#	continue
					print(datfile.split('\\')[-1])
					if getfiles == 1:
						#ln.zSetWave(1,0.6328,1)
						print("Appending .dat file...", end = '')
						ln.zImportExtraData(3,datfile)
						#print("Getting SPO...", end = '')
						#ln.zGetTextFile(datfile[:-4]+'_SPO_'+str(AOI)+'deg.txt','Spt',settingsFile=r"D:\Documents\School\Master's Report\SPO_Centroid.CFG")
						print("Wfm...", end = '')
						ln.zGetTextFile(datfile[:-4]+'_Wfm_'+str(AOI)+'deg.txt','Wfm',settingsFile=r"D:\Documents\School\Master's Report\Wfm.CFG")
						
						print("Getting DENC...", end = '')
						fname = datfile.split('\\')[-1]
						fn = re.sub('[^0-9,_,.]','',fname)
						fn = [float(i) for i in fn.strip('__').split('_')]
						sf = [1/100, 1/18, 1/100, 1/10000, 1/1000, 1/1000, 1/1000, 1/1000]
						fn = [a * b for a,b in zip(fn,sf)]
						denc = []
						dencx = []
						dency = []
						densq = []
						for fnum in FNUMs:
							ln.zSetSurfaceParameter(1,1,fnum*18)
							for WL in WLs:
								ln.zSetWave(1,WL,1)
								ln.zPushLens()
								ln.zOptimize(-1)
								denc.append(ln.zOperandValue('DENC',5,1,1,0.9,1,0,0,0))
								dencx.append(ln.zOperandValue('DENC',5,1,1,0.9,2,0,0,0))
								dency.append(ln.zOperandValue('DENC',5,1,1,0.9,3,0,0,0))
								densq.append(ln.zOperandValue('DENC',5,1,1,0.9,4,0,0,0))
								print(ln.zGetSurfaceParameter(1,1))
					with open(datfile[:-4]+'_Wfm_'+str(AOI)+'deg.txt','r',encoding='utf-16') as tmp:
						tmp2 = tmp.read().split('\n')
						wfm = float(tmp2[9].split(' ')[-2])
						wfms.append(wfm)
						fn[3] = wfm
					if customfile == '':
						customfile == 'DENC.txt'
					with open(ddir + '\\' + customfile,'a') as fil:
						##fil.write('\t'.join([str(i) for i in fn]) + '\t'+str(denc) + '\t'+str(dencx) + '\t'+str(dency) + '\t'+str(densq) + '\t'+str(AOI) + '\n')
						##fil.write('\t'.join([str(i) for i in denc]) + '\t' + '\t'.join([str(i) for i in dencx]) + '\t' + '\t'.join([str(i) for i in dency])+ '\t' + '\t'.join([str(i) for i in densq])+'\n')
						fil.write('\t'.join([str(i) for i in fn]) + '\t'+'\t'.join([str(i) for i in denc]) + '\t' + '\t'.join([str(i) for i in dencx]) + '\t' + '\t'.join([str(i) for i in dency])+ '\t' + '\t'.join([str(i) for i in densq]) + '\t'+str(AOI) + '\n')
						"""
						ln.zOptimize(-1)
						print("Getting RMS SPO...", end = '')
						RMSspo_list.append(ln.zGetOperand(1,10))
						print("Getting Geo SPO...", end = '')
						RMSspo_list.append(ln.zGetOperand(2,10))
						print("Getting Encircled Energy...", end = '')
						enc_list.append(ln.zGetOperand(3,10))
			`			"""
						print("Done!")
					"""
					with open(datfile[:-4]+'_SPO_'+str(AOI)+'deg.txt','r',encoding='utf-16') as tmp:
						tmp2 = tmp.read().split('\n')
						geospo_list.append(float(tmp2[-3].split(' ')[-2]))
						RMSspoY_list.append(float(tmp2[-4].split(' ')[-2]))
						RMSspoX_list.append(float(tmp2[-5].split(' ')[-2]))
						RMSspo_list.append(float(tmp2[-6].split(' ')[-2]))
					"""
			#flist_sorted = [float(i.split('DecX')[-1].split('_')[0]) for i in datfiles]