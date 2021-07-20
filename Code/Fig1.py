import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

plt.style.use('Leighton_Style_Colors2')


##### Functions for computing twist and dband strain in the small-angle limit:


# Computes D-band Strain as a function of psi0, lambda, and zeta in the small psi0 limit
def DBandStrain(mean_psi0_squared,lambda_array,zeta):
	a = ((zeta-1)/(zeta*lambda_array**(3/2) - lambda_array**(-3/2)))
	return 1/2 * (1 - a**2)*mean_psi0_squared

# Calculates <psi> as a function of Dband strain in the small psi_0 limit
def psi_D(epsilon_D,psi0,psi0squared):

	return psi0 * np.sqrt(1 - 2*epsilon_D/psi0squared)

# Dband strain in the zeta\to\infty limit
def DBandStrain_ZetaInfty(lambda_array,mean_psi0_squared):
	return (1/2)*(1 - 1/lambda_array**3)*mean_psi0_squared

# Dband strain in the zeta\to1 limit
def DBandStrain_ZetaOne(lambda_array,mean_psi0_squared):
	N = len(lambda_array)
	Dbandstrain = np.zeros(N)

	for i in range(N):
		if i>0:
			Dbandstrain[i] = mean_psi0_squared/2

	return Dbandstrain

##############################################################

##### Parameters, setup, and data

# Number of gridpoints for model curves:
N=10000

# Linestyle stuff:
ls_list = ['-','--','-.',':']
lw_list = [2,2.5,3,3.5]
color_list = ['black','blue','xkcd:orange','xkcd:red']

# Array for fibril strain:
lambda_array = np.linspace(1,1.1,num=N)
epsilonF_array = lambda_array-1

# Array for Dband strain:
epsilon_D_list = np.linspace(0,0.01,num=N)

# psi0 value in Bell:
psi0=16*np.pi/180

# List of <psi_0^2> values to show:
mean_psi0_squaredlist = [0.005,0.01,0.02,0.04]

# List of zeta values to show:
zeta_array = np.array([1.1,1.2,1.3])


##### Digitized data from Bell et al, 2018:

# Dband strain measurements:
bell_eD_data = np.array([0,0.1,0.24,0.3]) #percent
bell_eD_error = np.array([0,0.065,0.035,0.055]) #percent

# Tissue strain measurements:
bell_eTissue_data = 100*np.array([0,0.014,0.028,0.05]) # percent
bell_eTissue_error = np.array([0,0.3,0.3,0.3]) #percent

# Twist measurements <psi>/<psi_0>:
bell_twist_data= np.array([1,14/16,12/16,11/16]) 
bell_twist_error = np.array([0,np.sqrt((14/16)**2 + 1)/16,np.sqrt((12/16)**2 + 1)/16,np.sqrt((11/16)**2 + 1)/16])

##############################################################

##### This function plots Figure 1:

def PlotFig1():

	##### Figure Setup

	width = 8.6
	height = 9


	fig=plt.figure(figsize=(2*width/2.54,height/2.54))
	gs=gridspec.GridSpec(1,2,width_ratios=[1,1],height_ratios=[1])#,left=0.083,bottom=0.112,right=0.985,top=0.926,wspace=0.315)
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])
	ax1.minorticks_on()
	ax2.minorticks_on()


	##### Fig 1A (twist vs Dband strain) lines and data

	# Model curves for different <psi_0^2> values:
	for i in range(len(mean_psi0_squaredlist)):
		ax1.plot(100*epsilon_D_list,psi_D(epsilon_D_list,psi0,mean_psi0_squaredlist[i])/psi0,label=r'$\langle\psi_0^2\rangle =$'+str(mean_psi0_squaredlist[i]),ls=ls_list[i],color=color_list[i],lw=lw_list[i])

	# Bell Data (Twist vs Dband strain_)
	ax1.errorbar(bell_eD_data,bell_twist_data,yerr = bell_twist_error,xerr = bell_eD_error,color='tab:green',lw=3,zorder=1,ls='none',fmt='o')


	# approximate fit for <psi_0^2> (by eye)
	mean_psi0_squared=0.01


	##### Fig 1B (Dband strain (normalized by <psi_0^2>) vs Fibril Strain) lines and data

	# Zeta to 1 limit:
	ax2.plot(100*epsilonF_array,DBandStrain_ZetaOne(lambda_array,mean_psi0_squared)/mean_psi0_squared,lw=1.5,ls = '-.',color='xkcd:purple')#,label=r'$\zeta\to1^+$')
	
	# Model curves for different zeta values
	for i in range(len(zeta_array)):
		ax2.plot(100*epsilonF_array,DBandStrain(mean_psi0_squared,lambda_array,zeta_array[i])/mean_psi0_squared,label='$\zeta=$'+str(zeta_array[i]),ls=ls_list[i],color=color_list[i],lw=lw_list[i])

	# Zeta to infty limit:
	ax2.plot(100*epsilonF_array,DBandStrain_ZetaInfty(lambda_array,mean_psi0_squared)/mean_psi0_squared,lw=1.5,ls = ':',color='xkcd:red')#,label=r'$\zeta\to\infty$')

	# Bell Data:
	ax2.errorbar(bell_eTissue_data,bell_eD_data/100 /mean_psi0_squared,yerr = bell_eD_error/100 /mean_psi0_squared,xerr = bell_eTissue_error,color='tab:green',lw=3,zorder=10,ls='none',fmt='o')

	#ax2.errorbar(bell_eTissue_data*0.5,bell_eD_data/100 /mean_psi0_squared,yerr = bell_eD_error/100 /mean_psi0_squared,xerr = bell_eTissue_error,color='tab:green',lw=1.5,zorder=1,ls='none',fmt='o')
	ax2.scatter(bell_eTissue_data*0.5,bell_eD_data/100 /mean_psi0_squared,color='lawngreen',marker='*',zorder=9,s=140)



	##### Figure Labels/Legends/Formatting Stuff

	ax1.set_ylabel(r'$\langle \psi\rangle/\langle\psi_0\rangle$')
	ax1.set_xlabel('D-Band Strain (\%)')
	ax1.set_xlim(0,1)
	ax1.set_ylim(0,1.3)
	ax1.set_title('A)',loc='left')
	ax1.legend(loc='best',frameon=False)

	ax2.set_xlabel('Fibril Strain (\%)')
	ax2.set_ylabel(r'$\varepsilon_D/\langle\psi_0^2\rangle$')
	ax2.set_xlim(0,10)
	ax2.set_ylim(0,0.51)
	ax2.legend(loc='best',frameon=False)
	ax2.set_title('B)',loc='left')

	ax2.text(6,0.475,r'$\zeta\to1^+$',color='xkcd:purple')
	ax2.text(6,0.09,r'$\zeta\to\infty$',color='xkcd:red',rotation=15)

	ax1.tick_params(axis='x')
	ax1.tick_params(axis='y')
	ax2.tick_params(axis='x')
	ax2.tick_params(axis='y')

	plt.tight_layout(pad=0)
	plt.show()





PlotFig1()




