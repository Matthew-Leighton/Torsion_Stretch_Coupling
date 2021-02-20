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

def Psi_Array(psi0_array,lambda_value,zeta):
	#Returns correct branch of the solution
	psi = Psi(psi0_array,lambda_value,zeta)
	for i in range(len(psi0_array)):
		if psi[i]<0:
			psi[i]=psi[i]+ (np.pi/2)
	return psi



##### Computes D-band Strain as a function of psi0, lambda, and zeta

def DBandStrain(psi0_array,r_array,lambda_array, zeta):
	n = len(lambda_array)
	epsilonD_array = np.zeros(n)

	for j in range(n):
		psi_array = Psi_Array(psi0_array,lambda_array[j],zeta)

		mean_cos4_psi = 0
		mean_cos2_psi = 0
		mean_cos4_psi0 = 0
		mean_cos2_psi0 = 0

		for i in range(len(r_array)-1):
			mean_cos2_psi += r_array[i+1]*(np.cos(psi_array[i+1])**2)* (r_array[i+1] - r_array[i])
			mean_cos4_psi += r_array[i+1]*(np.cos(psi_array[i+1])**4)* (r_array[i+1] - r_array[i])
			mean_cos2_psi0 += r_array[i+1]*(np.cos(psi0_array[i+1])**2)* (r_array[i+1] - r_array[i])
			mean_cos4_psi0 += r_array[i+1]*(np.cos(psi0_array[i+1])**4)* (r_array[i+1] - r_array[i])

		mean_cos2_psi *= 2/(r_array[-1]**2)
		mean_cos4_psi *= 2/(r_array[-1]**2)
		mean_cos2_psi0 *= 2/(r_array[-1]**2)
		mean_cos4_psi0 *= 2/(r_array[-1]**2)

		x = mean_cos4_psi * mean_cos2_psi0 / (mean_cos2_psi * mean_cos4_psi0)
		epsilonD_value = np.sqrt(x)-1

		epsilonD_array[j] = epsilonD_value


	return epsilonD_array





def PlotFig1():

	#Parameters
	N=1000
	n=100
	lambda_array = np.linspace(1,1.1,num=n)
	epsilonF_array = lambda_array-1

	#Unstrained Twist Functions:
	r_array = np.linspace(0,1,num=N)
	psi0_const_array = np.ones(N)*0.1
	psi0_lin_array = 0.3*r_array
	psi0_tanh_array = 0.3*np.tanh(2*r_array)

	# For each twist function we'll show three zeta values:
	zeta_array = np.array([1.1,1.3,1.5])
	ls_list = ['-','--','-.']
	lw_list = [2,2.5,3]
	color_list = ['black','blue','xkcd:orange']


	# Figure Setup
	fig=plt.figure()
	gs=gridspec.GridSpec(1,3,width_ratios=[1,1,1],height_ratios=[1])
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])
	ax3=plt.subplot(gs[2])
	ax1.minorticks_on()
	ax2.minorticks_on()
	ax3.minorticks_on()


	for i in range(len(zeta_array)):
		zeta = zeta_array[i]
		ax1.plot(100*(lambda_array-1),100*DBandStrain(psi0_const_array,r_array,lambda_array,zeta),label='$\zeta = $'+str(zeta),lw=lw_list[i],ls=ls_list[i],color = color_list[i])
		ax2.plot(100*(lambda_array-1),100*DBandStrain(psi0_lin_array,r_array,lambda_array,zeta),lw=lw_list[i],ls=ls_list[i],color = color_list[i])
		ax3.plot(100*(lambda_array-1),100*DBandStrain(psi0_tanh_array,r_array,lambda_array,zeta),lw=lw_list[i],ls=ls_list[i],color = color_list[i])



	ax1.set_xlabel('Fibril Strain (\%)',fontsize=16)
	ax1.set_ylabel('D-Band Strain (\%)',fontsize=16)
	ax1.set_xlim(0,10)
	ax1.set_ylim(0,)
	ax1.set_title('$\psi_0(r) = 0.1$',fontsize=16)
	ax1.legend(loc='best',fontsize=14)

	ax2.set_xlabel('Fibril Strain (\%)',fontsize=16)
	#ax2.set_ylabel('D-Band Strain (\%)',fontsize=16)
	ax2.set_xlim(0,10)
	ax2.set_ylim(0,)
	ax2.set_title('$\psi_0(r) = 0.3r/R$',fontsize=16)

	ax3.set_xlabel('Fibril Strain (\%)',fontsize=16)
	#ax3.set_ylabel('D-Band Strain (\%)',fontsize=16)
	ax3.set_xlim(0,10)
	ax3.set_ylim(0,)
	ax3.set_title(r'$\psi_0(r) = 0.3\tanh(2r/R)$',fontsize=16)


	ax1.tick_params(axis='x', labelsize=14)
	ax1.tick_params(axis='y', labelsize=14)
	ax2.tick_params(axis='x', labelsize=14)
	ax2.tick_params(axis='y', labelsize=14)
	ax3.tick_params(axis='x', labelsize=14)
	ax3.tick_params(axis='y', labelsize=14)

	plt.tight_layout(pad=0.5)
	plt.show()


PlotFig1()
