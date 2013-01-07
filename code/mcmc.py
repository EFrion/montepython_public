import os,sys
import math
import random as rd
import numpy as np
import io
import data

# COMPUTE LKL
def compute_lkl(_cosmo,data):
  # Prepare the cosmological module with the new set of parameters
  _cosmo.set(data.cosmo_arguments)

  # Compute the model, keeping track of the errors
  failure=False
  # In classy.pyx, we made use of two type of python errors, to handle two
  # different situations.
  # - AttributeError is returned if a parameter was not properly set during the
  # initialisation (for instance, you entered Ommega_cdm instead of Omega_cdm).
  # Then, the code exits, to prevent running with imaginary parameters. This
  # behaviour is also used in case you want to kill the process.
  # - NameError is returned if Class fails to compute the output given the
  # parameter values. It will display the Class error, the code will register a
  # new failure, and start again with a new point.
  try:
    _cosmo._compute(["lensing"])
  except NameError :
    return True,0 
  except (AttributeError,KeyboardInterrupt):
    exit()

  # For each desired likelihood, compute its value against the theoretical model
  loglike=0
  flag_wrote_fiducial = 0

  for likelihood in data.lkl.itervalues():
    value = likelihood.loglkl(_cosmo,data)
    loglike += value
    # In case the fiducial file was written, store this information
    if value==1:
      flag_wrote_fiducial += 1

  # Compute the derived parameters if relevant
  if data.get_mcmc_parameters(['derived']) != []:
    try:
      _cosmo.get_current_derived_parameters(data)
    except NameError:
      print('Terminating now')
      exit()
  for elem in data.get_mcmc_parameters(['derived']):
    data.mcmc_parameters[elem]['current'] /= data.mcmc_parameters[elem]['initial'][4]


  # Clean the cosmological strucutre
  _cosmo._struct_cleanup(set(["lensing","nonlinear","spectra","primordial","transfer","perturb","thermodynamics","background","bessel"]))

  if flag_wrote_fiducial > 0:
    if flag_wrote_fiducial == len(data.lkl):
      print('--> Fiducial file(s) was(were) created, please start a new chain')
      exit()
    else:
      print('--> Some previously non-existing fiducial files were created, but potentially not all of them')
      print('--> Please check now manually on the headers of the corresponding that all')
      print('--> parameters are coherent for your tested models')
      exit()

  return failure,loglike


# Function used only when the restart flag was set. It will simply pick up the
# last accepted values as a starting point.
def read_args_from_chain(data,chain):
  # Chain is defined as a special File (that inherits from File, and has one
  # new method, tail). The class is defined in code/io.py
  Chain = io.File(chain,'r')
  parameter_names = data.get_mcmc_parameters(['varying'])

  # Warning, that feature was not tested since the adding of derived paramaters

  # BE CAREFUL: Here it works because of the particular presentation of the
  # chain, and the use of tabbings (not spaces). Please keep this in mind if
  # having difficulties
  i = 1
  for elem in parameter_names:
    data.mcmc_parameters[elem]['last_accepted'] = float(Chain.tail(1)[0].split('\t')[i])
    i+=1

# Will deduce the starting covariance matrix, either from the prior, or from an
# existing matrix. Reordering of the names and scaling take place here, in a
# serie of potentially hard to read methods. For the sake of clarity, and to
# avoid confusions, the code will, by default, print out a succession of 4
# covariance matrices at the beginning of the run, if starting from an existing
# one. This way, you can control that the paramters are set properly.
def get_cov(data,command_line):
  np.set_printoptions(precision=2,linewidth=150)
  parameter_names = data.get_mcmc_parameters(['varying'])
  i=0

  # if the user wants to use a covmat file
  if command_line.cov is not None:
    cov=open('{0}'.format(command_line.cov),'r')
    for line in cov:
      if line.find('#')!=-1:
	covnames = line.strip('#').replace(' ','').replace('\n','').split(',')
	M=np.zeros((len(covnames),len(covnames)),'float64')
	rot=np.zeros((len(covnames),len(covnames)))
      else:
	line=line.split()
	for j in range(len(line)):
	  M[i][j]=np.array(line[j],'float64')
	i+=1

    # First print out
    print('\nInput covariance matrix:')
    print(covnames)
    print(M)
    # Deal with the all problematic cases. 
    # First, adjust the scales between stored parameters and the ones used in mcmc
    scales = []
    for elem in covnames:
      if elem in parameter_names:
	scales.append(data.mcmc_parameters[elem]['initial'][4])
      else:
	scales.append(1)
    scales = np.diag(scales)
    invscales = np.linalg.inv(scales)
    M = np.dot(invscales.T,np.dot(M,invscales))
    
    # Second print out, after having applied the scale factors
    print('\nFirst treatment (scaling)')
    print(covnames)
    print(M)

    # Then, rotate M for the parameters to be well ordered, even if some names
    # are missing or some are in extra.
    temp_names = []
    for elem in parameter_names:
      if elem in covnames:
	temp_names.append(elem)

    # Trick if parameter_names contains less things than covnames:
    temp_names_2 = []
    h = 0
    not_in = [elem for elem in covnames if elem not in temp_names]
    for k in range(len(covnames)):
      if covnames[k] not in not_in:
        temp_names_2.append(temp_names[h])
	h+=1
      else:
	temp_names_2.append('')

    for k in range(len(covnames)):
      for h in range(len(covnames)):
	try:
	  if covnames[k]==temp_names_2[h]:
	    rot[h][k] = 1.
	  else:
	    rot[h][k] = 0.
	except IndexError:
	  rot[h][k] = 0.
    M=np.dot(rot,np.dot(M,np.transpose(rot)))

    # Third print out
    print('\nSecond treatment (partial reordering and cleaning)')
    print(temp_names_2)
    print(M)
    
    M_temp    = np.ones((len(parameter_names),len(parameter_names)),'float64')
    indices_1 = np.zeros(len(parameter_names))
    indices_2 = np.zeros(len(covnames))
    #Remove names that are in param_names but not in covnames
    for k in range(len(parameter_names)):
      if parameter_names[k] in covnames:
	indices_1[k]=1
    for zeros in np.where(indices_1 == 0)[0]:
      M_temp[zeros,:] = 0
      M_temp[:,zeros] = 0
    # Remove names that are in covnames but not in param_names
    for h in range(len(covnames)):
      if covnames[h] in parameter_names:
	indices_2[h]=1
    # There, put zeros in the final matrix (the sigmas will be used there)
    for zeros in np.where(indices_2 == 0)[0]:
      M[zeros,:] = 0
      M[:,zeros] = 0
    # super trick, now everything is in place
    M_temp[M_temp == 1]=M[M!=0]
    M = np.copy(M_temp)
    # on all other lines, just use sigma^2
    for zeros in np.where(indices_1 == 0)[0]:
      M[zeros,zeros] = np.array(data.mcmc_parameters[parameter_names[zeros]]['initial'][3],'float64')**2

  # else, take sigmas^2.
  else:
    M = np.identity(len(parameter_names),'float64')
    for elem in parameter_names:
      M[i][i]=np.array(data.mcmc_parameters[elem]['initial'][3],'float64')**2
      i+=1

  # Final print out, the actually used covariance matrix
  sys.stdout.write('\nDeduced starting covariance matrix:\n')
  print(parameter_names)
  print(M)

  #inverse, and diagonalization
  eigv,eigV=np.linalg.eig(np.linalg.inv(M))
  return eigv,eigV


# Routine to obtain a new position in the parameter space from the eigen values
# of the inverse covariance matrix.
def get_new_pos(data,eigv,U,k):
  
  parameter_names = data.get_mcmc_parameters(['varying'])
  vector_new=np.zeros(len(parameter_names),'float64')
  sigmas=np.zeros(len(parameter_names),'float64')

  # Write the vector of last accepted points, 
  vector = np.zeros(len(parameter_names),'float64')
  try:
    for elem in parameter_names:
      vector[parameter_names.index(elem)] = data.mcmc_parameters[elem]['last_accepted']
  except KeyError: # If it does not exist yet (initialization routine), take the mean value
    for elem in parameter_names:
      vector[parameter_names.index(elem)] = data.mcmc_parameters[elem]['initial'][0]

  rd.seed()

  # Choice here between sequential and global change of direction
  if data.jumping == 'global':
    for i in range(len(vector)):
      sigmas[i]=(math.sqrt(1/eigv[i]/len(vector)))*rd.gauss(0,1)*data.jumping_factor
  elif data.jumping == 'sequential':
    i = k%len(vector)
    sigmas[i] = (math.sqrt(1/eigv[i]))*rd.gauss(0,1)*data.jumping_factor
  else:
    print('\n\n  Jumping method unknown (accepted : global (default), sequential)')

  # Fill in the new vector
  vector_new = vector + np.dot(U,sigmas)

  # Check for boundaries problems
  i	= 0
  flag  = 0
  for elem in parameter_names:
    value = data.mcmc_parameters[elem]['initial']
    if( (str(value[1])!=str(-1) and value[1] is not None) and vector_new[i]<value[1]):
      flag+=1 # if a boundary value is reached, increment
    elif( (str(value[2])!=str(-1) and value[1] is not None) and vector_new[i]>value[2]):
      flag+=1 # same
    i+=1

  # At this point, if a boundary condition is not fullfilled, ie, if flag is
  # different from zero, return False
  if flag!=0:
    return False
  
  # If it is not the case, proceed with normal computationThe value of
  # new_vector is then put into the 'current' point in parameter space.
  i=0
  for elem in parameter_names:
    data.mcmc_parameters[elem]['current'] = vector_new[i]
    i+=1
    
  # Propagate the information towards the cosmo arguments
  data.update_cosmo_arguments()

  # Return 
  return True

# Transfer the 'current' point in the varying parameters to the last accepted
# one.
def accept_step(data):
  for elem in data.get_mcmc_parameters(['varying']):
    data.mcmc_parameters[elem]['last_accepted'] = data.mcmc_parameters[elem]['current']
  for elem in data.get_mcmc_parameters(['derived']):
    data.mcmc_parameters[elem]['last_accepted'] = data.mcmc_parameters[elem]['current']


######################
# MCMC CHAIN
######################
def chain(_cosmo,data,command_line):

  ## Initialisation
  num_failure=1000   # Default number of accepted failure
  loglike=0
  failure=False    # Failure flag

  # Recover the covariance matrix according to the input, if the varying set of
  # parameters is non-zero
  if data.get_mcmc_parameters(['varying']) != []:
    sigma_eig,U=get_cov(data,command_line)
  # In case of a fiducial run (all parameters fixed), simply run once and print
  # out the likelihood
  else:
    print(' /|\  You are running with no varying parameters...')
    print('/_o_\ Computing model for only one point')
    failure,loglike = compute_lkl(_cosmo,data)
    io.print_vector([data.out,sys.stdout],1,loglike,data)
    return 1,loglike

  # Counter for the number of failures
  failed=0

  # If restart wanted, pick initial value for arguments
  if command_line.restart is not None:
    read_args_from_chain(data,command_line.restart)

  # Pick a position (from last accepted point if restart, from the mean value
  # else)
  for i in range(100):
    if get_new_pos(data,sigma_eig,U,failed) is True:
      break
  # Compute the starting Likelihood
  failure,loglike=compute_lkl(_cosmo,data)

  # Failure check of initialization
  while ((failure is True) and (failed<=num_failure)):
    failed +=1
    for i in range(100):
      if get_new_pos(data,sigma_eig,U,failed) is True:
	break
    failure,loglike=compute_lkl(_cosmo,data)

  if failure is True:
    print(' /|\   {0} tried {1} times to initialize with given parameters, and failed...'.format(data.cosmological_module_name,num_failure+2))
    print('/_o_\  You might want to change your starting values, or pick default ones!')
    exit()

  # If the first step was finally computed, pick it as the last accepted value
  # (accept_step), and modify accordingly the max_loglike
  accept_step(data)
  max_loglike = loglike

  acc,rej=0.0,0.0	# acceptance and rejection number count
  N=1			# number of time the system stayed in the current position

  # Print on screen the computed parameters
  io.print_parameters(sys.stdout,data)

  k = 1
  # Main loop, that goes on while the maximum number of failure is not reached,
  # and while the expected amount of steps (N) is not taken.
  while (k <= command_line.N and failed < num_failure):

    # Pick a new position ('current' flag in mcmc_parameters), and compute its
    # likelihood. If get_new_pos returns True, it means it did not encounter
    # any boundary problem. Otherwise, just increase the multiplicity of the
    # point and start the loop again
    if get_new_pos(data,sigma_eig,U,k) is True:
      failure,newloglike=compute_lkl(_cosmo,data)
    else: #reject step
      rej+=1
      N+=1
      k+=1
      continue
    
    # In case of failure in the last step, print out the faulty
    # cosmo_arguments, and start a new incrementation of the while loop. Note
    # that k was not incremented due to the continue statement: this failed
    # point will not count towards the total number of steps asked.
    if(failure==True):
      failed += 1
      print('Warning: %s failed due to choice of parameters, picking up new values' % data.cosmological_module_name)
      print(data.cosmo_arguments)
      continue

    # Harmless trick to avoid exponentiating large numbers. This decides
    # whether or not the system should move.
    if (newloglike != data.boundary_loglike):
      if (newloglike >= loglike):
	alpha = 1.
      else:
	alpha=np.exp(newloglike-loglike)
    else:
      alpha = -1

    if ((alpha == 1.) or (rd.uniform(0,1) < alpha)): # accept step

      # Print out the last accepted step (WARNING: this is NOT the one we just
      # computed ('current' flag), but really the previous one.) with its
      # proper multiplicity (number of times the system stayed there).
      io.print_vector([data.out,sys.stdout],N,loglike,data)

      # Report the 'current' point to the 'last_accepted'
      accept_step(data)
      loglike=newloglike
      if loglike > max_loglike:
	max_loglike = loglike
      acc+=1.0
      N=1 # Reset the multiplicity
      
    else: # reject step
      rej+=1.0
      N+=1 # Increase multiplicity of last accepted point
      
    # Regularly (option to set in parameter file), close and reopen the buffer
    # to force to write on file.
    if acc % data.write_step ==0:
      io.refresh_file(data)
    k += 1 # One iteration done
  # END OF WHILE LOOP

  # If at this moment, the multiplicity is higher than 1, it means the current
  # point is not yet accepted, but it also mean that we did not print out the
  # last_accepted one yet. So we do.
  if N>1:
    io.print_vector([data.out,sys.stdout],N-1,loglike,data)


  # Warn the user that the code finished because of class failures.
  if (failed == num_failure):
    sys.stdout.write('\n\n /|\   The computation failed {0} times, \n'.format(failed))
    sys.stdout.write('/_o_\  Please check the values of your parameters\n')

  # Print out some information on the finished chain
  rate=acc/(acc+rej)
  sys.stdout.write('\n#  {0} steps done, acceptance rate: {1}\n'.format(command_line.N,rate))
  
  # For a restart, and if the code did not fail too much, erase the starting
  # point to keep only the new, longer chain.
  if ((command_line.restart is not None) and (failed < num_failure)):
    os.remove(command_line.restart)
    sys.stdout.write('  deleting starting point of the chain {0}\n'.format(command_line.restart))
  return rate,max_loglike
