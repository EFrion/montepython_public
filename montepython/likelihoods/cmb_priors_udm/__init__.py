import os
import numpy as np
import warnings
import montepython.io_mp as io_mp
from montepython.likelihood_class import Likelihood
import scipy.constants as cons


class cmb_priors_udm(Likelihood):

	# initialization routine

	def __init__(self, path, data, command_line):

		Likelihood.__init__(self, path, data, command_line)

		# read covariance matrix
		self.cov_data = np.loadtxt(os.path.join(self.data_directory, self.cov_file))

	# initialisation of the class is done within the parent Likelihood_prior. For
	# this case, it does not differ, actually, from the __init__ method in
	# Likelihood class.
	def loglkl(self, cosmo, data):

		omega_b = cosmo.omega_b()
		n_s = cosmo.n_s()
		h = cosmo.h()
		z_star = cosmo.get_current_derived_parameters(['z_star'])['z_star']
		rs_star = cosmo.get_current_derived_parameters(['rs_star'])['rs_star']

		l_a = (1.+z_star)*np.pi*cosmo.angular_distance(z_star)/rs_star
		Rshift =  100*h*(1.+z_star)*cosmo.angular_distance(z_star)*np.sqrt(cosmo.Omega0_m())/(cons.c/1000)

		theo_cmb = np.array([Rshift,l_a,omega_b,n_s])
		data_cmb = np.array([self.Rshift,self.l_a,self.omega_b,self.n_s])
		diff = theo_cmb - data_cmb
		
		# compute chi squared
		inv_cov_data = np.linalg.inv(self.cov_data)
		chi2 = np.dot(np.dot(diff,inv_cov_data),diff)

		# return ln(L)
		loglkl = -0.5 * chi2

		data.derived_lkl = {'Omega_DE':0.0}

		#log10beta_udm = data.mcmc_parameters['log10beta_udm']['current'] * data.mcmc_parameters['log10beta_udm']['scale']
		#beta_udm = 1.0
		#beta_udm = 10**(log10beta_udm)
		beta_udm = data.mcmc_parameters['beta_udm']['current'] * data.mcmc_parameters['beta_udm']['scale']
		a_t = data.mcmc_parameters['a_t']['current'] * data.mcmc_parameters['a_t']['scale']
		#log10a_t = data.mcmc_parameters['log10a_t']['current'] * data.mcmc_parameters['log10a_t']['scale']
		#a_t = 10**(log10a_t)
		Omega0_fld_fake = data.mcmc_parameters['Omega0_fld']['current'] * data.mcmc_parameters['Omega0_fld']['scale']
		Omega0_fld = cosmo.get_current_derived_parameters(['Omega0_fld'])['Omega0_fld']

		Omega_DE = Omega0_fld / (0.5 + 1.5 * np.log(np.cosh(beta_udm * (1.-pow(a_t,3.))/3.))/beta_udm)

		data.derived_lkl.update({'Omega_DE':Omega_DE})

		#print(Omega0_fld_fake,Omega0_fld)
		return loglkl
