# -*- coding: utf-8 -*-
import os
import csv
import numpy as np
from scipy.interpolate import griddata
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from timeit import default_timer as timer

"""
Note:
	"A" = Angular. New convention: "Tangential"
	"T" = Thermal. New convention: "Depth"
"""

#------- INPUTS --------\
toolfreqs = [10]		# Tooling Mark Frequencies [Cyc/ap]
#toolfreqs = [i/5 for i in range(1,int(500*1/3))]	# Tooling Mark Frequencies [Cyc/ap]
tooltype = 'sin'		# Tooling Mark Waveform ('sin', 'trap', 'square')
dia = 18				# Diameter of part [in]
offsy = 0				# Y Offset [in]
offsx = 0				# X Offset [in]
res = 0.05				# Sampling resolution [in]
rmsSurfSpecWVs = [.05]	# RMS Surface Specs [λ @ 632.8nm]
dutyCycles = [1]		# Waveform Duty Cycles
#dutyCycles = [(i+1)/100 for i in range(100)]		# Duty Cycle of Waveform
DX = [0]				# Decenter from on-axis [in]
#DX = [i*5/20 for i in range(200)]				# Decenter from on-axis [in]
afreqs = [0]			# Angular ("Tan") TM Seed Frequencies [cy/360°]
amags = [0]				# Angular ("Tan") TM overlap magnitude into next track
thermalMags = [0]		# Thermal ("Depth") Magnitude
thermalFreqs = [5]		# Thermal ("Depth") Seed Frequencies
#-----------------------/
output = 1		 		# Output grid files?
showfig = 0				# Draw 2D plots [1=yes, 0=no]
showfig3D = 0			# Draw 3D figures [1=yes, 0=no]
n = 1					# Overscale? (visualization ONLY) (?)
verbose = 0				# Output progress during scaling?
#-----------------------/
fdir = r"D:\Temp\A+TMAG_FINE"
#-----------------------/

wfms = [] # Used for 2nd pass refinement. Import .dat into Zemax -> Calculate RMS WFE map → Copy here and re-run.

#----- General Functions ---------------------------------------\
#																#
#																#
def calcRMS(Z):
	# Calculate Root Mean Square of matrix
	RMS = np.sqrt(1/len(Z)*np.mean(Z**2))
	return RMS

def calcRMSSlope(Z,res):
	return calcRMS(gradMag(Z,res))

def calcPSS(Z,res):
	return np.max(gradMag(Z,res))

def gradMag(Z,res):
	# Calculate gradient magnitude
	gZ = np.gradient(Z,res)
	return np.sqrt(gZ[0]**2+gZ[1]**2)

def plotGrad(X,Y,Z,res):
	# Plot gradient magnitude of matrix
	return plot3D(X,Y,gradMag(Z,res))

def plotProfile(X,Z,nd):
	# Plot central cross-sectional profile of matrix
	#fig = plt.figure()
	resx=len(X[0])
	[r, g, b, a] = mapper.to_rgba(nd)
	return plt.plot(X[int(resx/2)][:],Z[int(resx/2)][:],color=[r,g,b])

def plot3D(X,Y,Z):
	# Plot 3D Surface Plot of Meshgrid
	plt.close('all')
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	ax.set_zlim(-1.01, 1.01)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)
	ax.view_init(azim=0, elev=90)
	return plt.show()

def manualGrad(Z,resX,resY):
	# Manual Implementation of gradient (for error checking)
	dZx=np.zeros([len(Z)-2,len(Z)-2])
	dZy=np.zeros([len(Z)-2,len(Z)-2])
	dZz=np.zeros([len(Z)-2,len(Z)-2])
	dZ=[dZx,dZy,dZz]
	for i in range(1,len(Z)-1):
		for j in range(1,len(Z[0])-1):
			tmp = np.mean([Z[i+1,j]-Z[i,j],Z[i,j]-Z[i-1,j]])
			dZ[0][i-1,j-1] = 0 if (tmp == 0) else tmp*resx
			tmp = np.mean([Z[i,j+1]-Z[i,j],Z[i,j]-Z[i,j-1]])
			dZ[1][i-1,j-1] = 0 if (tmp == 0) else tmp*resy
			dZ[2][i-1,j-1] = np.sqrt(dZ[0][i-1,j-1]**2+dZ[1][i-1,j-1]**2)
	return dZ
#																#
#----- Shape Generation ----------------------------------------\
#																#
def sin_signal(X,Y,R, toolfreq, dia, A, dutyCycle, thermalMag, thermalFreq):
	sinWave = A*(np.cos(np.pi*toolfreq/(dia/2)*R)-1+dutyCycle*2)
	#sinWave +=A/5*(np.cos(np.pi*toolfreq/(dia/2)*thermalVar)-1+dutyCycle*2)
	# Add in thermal variation which does not depend with constant spatial frequency
	sinWave +=A*thermalMag*(np.cos(np.pi*thermalFreq/(dia/2)*X)-1+dutyCycle*2) 
	sinWave +=A*thermalMag*(np.cos(np.pi*thermalFreq/(dia/2)*Y)-1+dutyCycle*2)
	return (sinWave > 0) * sinWave

def square_signal(R, toolfreq, dia, A, dutyCycle):
	squareWave = sin_signal(R, toolfreq, dia, A, dutyCycle)
	squareWave = (squareWave > 0) * A
	return squareWave.astype(float)

def trapzoid_signal(t, width=2., slope=1., amp=1., offs=0):
	# Generate Trapezoidal Waveform, First Attempt
	a = slope*width*signal.sawtooth(2*np.pi*t/width, width=0.5)/4.
	a[a>amp/2.] = amp/2.
	a[a<-amp/2.] = -amp/2.
	return a + amp/2. + offs

#																#
#----- Processing for Zemax ------------------------------------\
#																#
def grid(x, y, z, resX=200, resY=200):
	# Convert 3 column data to matplotlib grid
	grid_x, grid_y = np.mgrid[min(x):max(x):1j * resX, min(y): max(y):1j * resY]
	grd = np.vstack((x,y))
	Z = griddata(grd.T, np.array(z), (grid_x, grid_y), method='cubic', )

	return grid_x, grid_y, Z

def surf_slopes(Zgrid):
	# Output first and second derivatives for surface slopes
	Fx, Fy = np.gradient(Zgrid)
	Fxx, Fxy = np.gradient(Fx)
	Fyx, Fyy = np.gradient(Fy)
	return Fx, Fy, Fxy, Fxx, Fyy

#																#
#																#
#---------------------------------------------------------------/



if __name__ == "__main__":
	for rmsSurfSpecWV in rmsSurfSpecWVs:
		norm = colors.Normalize(vmin=0, vmax=len(dutyCycles)*len(toolfreqs), clip=True)
		mapper = cm.ScalarMappable(norm=norm, cmap=cm.jet)
		
		rmsSurfSpecM = rmsSurfSpecWV*632.8*10**(-9)/0.0254*n
		rmsSurfSpecM *= 0.4382120946538124*2 # Scale Factor to correct for final RMS in Zemax...TBR. Use 2nd pass correction via wfms variable.

		#----- Containers -------------------------------------------\
		rmsSurfs = []
		rmsSlopes = []
		PSSs = []
		nstatus = 1
		nlen = len(toolfreqs)*len(dutyCycles)*len(DX)*len(afreqs)*len(amags)*len(thermalFreqs)*len(thermalMags)
		for ndx,decx in enumerate(DX):
			for nduty,dutyCycle in enumerate(dutyCycles):
				for nfreq,toolfreq in enumerate(toolfreqs):
					for nafreq,atoolfreq in enumerate(afreqs):
						for namag,atoolmag in enumerate(amags):
							for ntfreq,thermalFreq in enumerate(thermalFreqs):
								for ntmag,thermalMag in enumerate(thermalMags):
									#----- Initial Setup -------------------------------------------\
									fileout = r"ToolingMarks_{}_DC{}_DecX{}_Freq{}_RMS{}_Afreq{}_Amag{}_Tfreq{}_Tmag{}.dat".format(tooltype,str(round(100*dutyCycle)).zfill(5),str(decx).zfill(5),str(round(100*toolfreq)).zfill(5),str(1000*rmsSurfSpecWV),str(round(1000*atoolfreq)).zfill(5),str(round(1000*atoolmag)).zfill(5),str(round(1000*thermalFreq)).zfill(5),str(round(1000*thermalMag)).zfill(5))
									x = np.arange(-dia/2*1.1,dia/2*1.1,res)
									y = np.arange(-dia/2*1.1,dia/2*1.1,res)
									X, Y = np.meshgrid(x, y)
																	
									R = np.sqrt((X-decx)**2 + Y**2)*(1+atoolmag/np.sqrt((X-decx)**2 + Y**2)*np.cos(atoolfreq/np.pi**2*(X-decx))*np.cos(atoolfreq/np.pi**2*Y))
									
									A = 0.1*632.8*10**(-9)/0.0254*n

									#----- Shape Generation ----------------------------------------\
									print("")
									print("----- "+str(nstatus)+"/"+str(nlen)+" -----",end='')
									if nstatus == 1:
										t = timer()
									else:
										dt = (timer()-t)/(nstatus-1) * (nlen-nstatus-1)
										print(' ETA: '+str(int(dt/60))+'m'+str(int(dt%60))+'s -----')
									nstatus+=1
									print("Tooling Mark Type: "+str(tooltype))
									print("Tooling Mark Frequncy: "+str(toolfreq))
									print("Duty Cycle: "+str(dutyCycle))
									print("Decenter: "+str(decx))
									print("")
									if showfig:
										print("Plot overscale factor: "+str(n))
										print("")
									if tooltype in ['tri','triangle','trap','trapezoid']:
										if toolfreq == 10:
											Z = trapzoid_signal(R, width=2, slope=.5, amp=0.5, offs=0) # 10 cyc/ap
										elif toolfreq == 20:
											Z = trapzoid_signal(R, width=.95, slope=1000, amp=0.5, offs=0) # 20 cyc/ap
										else:
											throw('Invalid freq.')
										Z*=A
									elif tooltype in ['sin','sine']:
										Z = sin_signal(X,Y,R, toolfreq, dia, 1, dutyCycle,thermalMag,thermalFreq)
									elif tooltype in ['square','sq']:
										Z = square_signal(R, toolfreq, dia, 1, dutyCycle)
									else:
										throw('Invalid tooling mark type.')
									
									# Normalize
									Z-=np.min(Z)
									Z/=np.max(Z)


									#----- Scale to Specs ------------------------------------------\
									Z *= A
									rmsSurf = calcRMS(Z)
									print("RMS Spec = " + str(rmsSurfSpecWV)+"wv @ 632.8nm")
									print("Scaling Waveform RMS...",end="")
									
									# Check at end to force at least one iteration of surface figure optimization
									while np.abs(rmsSurf-rmsSurfSpecM) > rmsSurfSpecM*0.001:
										if rmsSurf < rmsSurfSpecM:
											if tooltype in ['tri','triangle','trap','trapezoid']:
												Z*=1.01
											elif tooltype in ['sin','sine','sq','square']:
												Z *= 1.01
										elif rmsSurf > rmsSurfSpecM:
											if tooltype in ['tri','triangle','trap','trapezoid']:
												Z*=0.99
											elif tooltype in ['sin','sine','sq','square']:
												Z*=0.99
										rmsSurf = calcRMS(Z)
										
										if verbose:
											print("Surf: "+str(rmsSurf))
										#plotProfile(X,Z)
										rmsSurf = calcRMS(Z)

									print("Done!")
									print("RMS Surf Difference = "+str((rmsSurf-rmsSurfSpecM)/rmsSurfSpecM*100)+"%")
									print("dutyCycle = "+str(dutyCycle))
									print("")
									
									if len(wfms)>0: # Correct RMS WFE discrepency to Zemax (2nd Run)
										Z/=wfms[nstatus-1]
										Z*=rmsSurfSpecWV
										
									rmsSurfs.append(rmsSurf/A/0.4382120946538124/2/10)
									rmsSlopes.append(calcRMSSlope(Z,res)/A/0.4382120946538124/2/10)
									PSSs.append(calcPSS(Z,res)/A/0.4382120946538124/2/10)
									
									#----- Processing for Zemax ------------------------------------\
									Xval = X.ravel()
									Yval = Y.ravel()
									Zval1 = Z.ravel()
									if showfig:
										plotProfile(X,Z,nstatus)
									if showfig3D:
										plot3D(X,Y,Z)
									fx, fy, fxy, fxx, fyy = surf_slopes(Z)

									# Create Grid Array
									print("Generating grid sag file...", end="")
									output_data = []
									nx = len(X)
									ny = len(Y)
									delx = res
									dely = res
									unitflag = 2 # 0 for mm, 1 for cm, 2 for in, 3 for m
									xdec = offsx
									ydec = offsy

									# Format for Zemax Grid Sag
									output_data.append([nx, ny, delx, dely, unitflag, xdec, ydec])
									for idr, row in enumerate(Z):
										for idc, col in enumerate(row):
											out_set = [Z[idr][idc], fx[idr][idc], fy[idr][idc], fxy[idr][idc]]
											output_data.append(out_set)

									# Write Grid Sag to Zemax .dat file
									if output:
										with open(os.path.join(fdir, fileout), 'w') as fil:
											outwrite = csv.writer(fil, delimiter=" ")
											outwrite.writerows(output_data)

									print("Done!")
									print("")
									if output:
										print("Fileout: ")
										print(fdir+'\\'+fileout)
	if showfig == 1:
		plt.show()
		