##### This file contains example functions to compute numerical solutions to key equations from the paper
##### "D-band strain underestimates collagen fibril strain" by Matthew P. Leighton, Andrew D. Rutenberg, and Laurent Kreplak.


# Required package
import numpy as np



##### These functions compute post-strain twist angle function as a function of psi0, lambda, and zeta:

def Psi(psi0,lambda_value,zeta):
    x = ((zeta+1)*(lambda_value**3 -1) + (zeta-1)*(lambda_value**3 +1) * np.cos(2*psi0))/( 2*lambda_value**(3/2)*(zeta-1)*np.sin(2*psi0) )
    return (1/2) * np.arctan(1/x)

# Computes an array of psi values given an array of psi_0 values and a single value of each of lambda and zeta
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


##### Small-angle limit:

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


	


