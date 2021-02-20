import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)



##### Computes post-strain twist angle function as a function of psi0, lambda, and zeta

def Psi(psi0,lambda_value,zeta):
    x = ((zeta+1)*(lambda_value**3 -1) + (zeta-1)*(lambda_value**3 +1) * np.cos(2*psi0))/( 2*lambda_value**(3/2)*(zeta-1)*np.sin(2*psi0) )
    return (1/2) * np.arctan(1/x)

def ExactPsi(psi0,lambda_array,zeta):
	#Returns correct branch of the solution
	psi = Psi(psi0,lambda_array,zeta)
	for i in range(len(lambda_array)):
		if psi[i]<0:
			psi[i]=psi[i]+ (np.pi/2)
	return psi



##### Computes D-band Strain as a function of psi0, lambda, and zeta

def DBandStrain(psi0,lambda_strain, zeta):







def PlotFig1():

	#Parameters
	zeta = 1.3
	N=1000
	lambda_array = np.linspace(1,1.2,num=N)
	epsilonF_array = lambda_array-1

	# Figure Setup
	fig=plt.figure()
	gs=gridspec.GridSpec(1,3,width_ratios=[1,1,1],height_ratios=[1])
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])
	ax3=plt.subplot(gs[2])
	ax1.minorticks_on()
	ax2.minorticks_on()
	ax3.minorticks_on()




