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
 
def DBandStrain(mean_psi0_squared,lambda_array,zeta):
	a = ((zeta-1)/(zeta*lambda_array**(3/2) - lambda_array**(-3/2)))
	return 1/2 * (1 - a**2)*mean_psi0_squared

def psi_D(epsilon_D,psi0,psi0squared):

	return psi0 * np.sqrt(1 - 2*epsilon_D/psi0squared)


def DBandStrain_ZetaInfty(lambda_array,mean_psi0_squared):
	return (1/2)*(1 - 1/lambda_array**3)*mean_psi0_squared

def DBandStrain_ZetaOne(lambda_array,mean_psi0_squared):
	N = len(lambda_array)
	Dbandstrain = np.zeros(N)

	for i in range(N):
		if i>0:
			Dbandstrain[i] = mean_psi0_squared/2

	return Dbandstrain





def PlotFig1():

	#Parameters
	N=1000
	n=1000
	lambda_array = np.linspace(1,1.1,num=n)
	epsilonF_array = lambda_array-1

	# A) Stuff
	bellepsilondata=np.array([0,0.1,0.24,0.3])
	belltwistdata = np.array([16,14,12,11])*np.pi/180
	psierror = [np.pi/180,np.pi/180,np.pi/180,np.pi/180]
	epsilonerror = [0,0.065,0.035,0.055]

	bellfdata = np.array([1,14/16,12/16,11/16])
	bellferror = np.array([np.sqrt(2)/16,np.sqrt((14/16)**2 + 1)/16,np.sqrt((12/16)**2 + 1)/16,np.sqrt((11/16)**2 + 1)/16])

	#B) data
	bell_eD_data = np.array([0,0.1,0.24,0.3,0.588]) #percent
	bell_eD_error = np.array([0,0.065,0.035,0.055,0.084]) #percent
	bell_eTissue_data = 100*np.array([0,0.014,0.028,0.05,0.08]) # percent
	bell_eTissue_error = np.array([0,0.3,0.3,0.3,0.3]) #percent

	epsilon_D_list = np.linspace(0,0.01,num=10000)
	psi0=16*np.pi/180
	mean_psi0_squaredlist = [0.005,0.01,0.02,0.04]#,0.12,0.16]


	# For each twist function we'll show three zeta values:
	zeta_array = np.array([1.1,1.3,1.5])
	ls_list = ['-','--','-.',':']
	lw_list = [2,2.5,3,3.5]
	color_list = ['black','blue','xkcd:orange','xkcd:red']


	# Figure Setup
	fig=plt.figure()
	gs=gridspec.GridSpec(1,2,width_ratios=[1,1],height_ratios=[1])
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])
	ax1.minorticks_on()
	ax2.minorticks_on()

	# Fig 1A lines and data

	for i in range(len(mean_psi0_squaredlist)):
		ax1.plot(100*epsilon_D_list,psi_D(epsilon_D_list,psi0,mean_psi0_squaredlist[i])/psi0,label=r'$\langle\psi_0^2\rangle =$'+str(mean_psi0_squaredlist[i]),ls=ls_list[i],color=color_list[i],lw=lw_list[i])
	#ax1.errorbar(bellepsilondata,belltwistdata,yerr = psierror,xerr = epsilonerror,color='tab:green',lw=3,zorder=1)
	ax1.errorbar(bellepsilondata,bellfdata,yerr = bellferror,xerr = epsilonerror,color='tab:green',lw=3,zorder=1,ls='none')


	# Best fit <psi_0^2>
	mean_psi0_squared=0.01


	# Fig 1B lines and data

	ax2.plot(100*epsilonF_array,100*DBandStrain_ZetaOne(lambda_array,mean_psi0_squared),label=r'$\zeta\to1^+$',lw=1.5,ls = '-.',color='xkcd:green')

	for i in range(len(zeta_array)):
		ax2.plot(100*epsilonF_array,100*DBandStrain(mean_psi0_squared,lambda_array,zeta_array[i]),label='$\zeta=$'+str(zeta_array[i]),ls=ls_list[i],color=color_list[i],lw=lw_list[i])

	ax2.plot(100*epsilonF_array,100*DBandStrain_ZetaInfty(lambda_array,mean_psi0_squared),label=r'$\zeta\to\infty$',lw=1.5,ls = ':',color='xkcd:red')

	ax2.errorbar(bell_eTissue_data,bell_eD_data,yerr = bell_eD_error,xerr = bell_eTissue_error,color='tab:green',lw=3,zorder=1,ls='none')


	# Figure Niceness Stuff

	ax1.set_ylabel(r'$\langle \psi\rangle/\langle\psi_0\rangle$',fontsize=16)
	ax1.set_xlabel('D-Band Strain (\%)',fontsize=16)
	ax1.set_xlim(0,1)
	ax1.set_ylim(0,1.3)
	ax1.set_title('A)',loc='left',fontsize=16)
	ax1.legend(loc='best',fontsize=14,frameon=False)

	ax2.set_xlabel('Fibril Strain (\%)',fontsize=16)
	#ax2.set_ylabel(r'$\varepsilon_D/\langle\psi_0^2\rangle$',fontsize=16)
	ax2.set_ylabel('D-Band Strain (\%)',fontsize=16)
	#ax2.set_ylabel('D-Band Strain (\%)',fontsize=16)
	ax2.set_xlim(0,10)
	ax2.set_ylim(0,0.75)
	ax2.legend(loc='best',fontsize=14,frameon=False)
	ax2.set_title('B)',loc='left',fontsize=16)


	ax1.tick_params(axis='x', labelsize=14)
	ax1.tick_params(axis='y', labelsize=14)
	ax2.tick_params(axis='x', labelsize=14)
	ax2.tick_params(axis='y', labelsize=14)

	plt.tight_layout(pad=0.5)
	plt.show()


PlotFig1()
