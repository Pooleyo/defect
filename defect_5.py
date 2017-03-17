# create realistic s_range
# calculate R_e from M value
# plot wilkens function
# plot integrand vs n and store plots in directory
# make peak intensity profile
# make comparison of peak intensity with Balogh figures
# must log everything

import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import os
from scipy import integrate

log_name = "log.defect"
program_name = "defect_3.py"



def initialisation(G, B,a_lattice,s_num,n_range,del_s_pernm,M,average_defect_density, contrast_factor):

	start_time = time.localtime()
	
	date_string = str(start_time[2]) + "/" + str(start_time[1]) + "/" + str(start_time[0])
	time_string = str(start_time[3]) + "h " + str(start_time[4]) + "m " + str(start_time[5]) + "s"


	f = open(log_name, "w")
	f.write(program_name + "\n" + date_string + "\n" + time_string)
	f.write("\n\nG = (" + str(G[0]) + ", " + str(G[1]) + ", " + str(G[2]) + ")\n"
	"B = (" + str(B[0]) + ", " + str(B[1]) + ", " + str(B[2]) + ")\n"
	"a_lattice = " + str(a_lattice) + "\n"
	"s_num = " + str(s_num) + "\n"
	"n_lower = " + str(n_range[0]) + "\n"
	"n_upper = " + str(n_range[-1]) + "\n"
	"n_num = " + str(len(n_range)) + "\n"
	"del_s_pernm = " + str(del_s_pernm) + "\n"
	"M = " + str(M) + "\n"
	"average_defect_density = " + str(average_defect_density) + "\n"
	"contrast_factor = " + str(contrast_factor) + "\n"
	)
	f.close()
	return;
	
	
def build_boundary_profile(M, average_defect_density, s_num, del_s_pernm, exponent_cutoff, samples_per_oscillation, G, B, contrast_factor, n_range_guess, a_lattice):
 
	R_e = R_e_from_M(M, average_defect_density)
	s_range = s_range_boundary_defined(G, B, contrast_factor, s_num, del_s_pernm, M, average_defect_density, n_range_guess, a_lattice)	
	
	n_range = create_n_range(exponent_cutoff, samples_per_oscillation, G, B, contrast_factor, average_defect_density, R_e, s_range)
	
	print n_range
	

	return s_range, n_range

	
	
def s_range_boundary_defined(G, B, contrast_factor, s_num, del_s_pernm, M, average_defect_density, n_range_guess, a_lattice):

	print "\nCreating intensity profile for the broadest peak...\n"

	gsqr = [0] * len(G)

	for i in range(len(G)):
	
		gsqr[i] = (G[i][0] ** 2) + (G[i][1] ** 2) + (G[i][2] ** 2)

	widest_peak_G = G[np.argmax(gsqr)]
	
	widest_peak_C = max(contrast_factor)


	calculate_intensities = True

	run(calculate_intensities, s_num, del_s_pernm, M, average_defect_density, n_range = n_range_guess, G = widest_peak_G, B = B, a_lattice = a_lattice, contrast_factor = widest_peak_C, comparison_with_balogh = False, wilkens_plot = True, cos_vs_n_plot = False, intensity_vs_n_plot = False)

	
	s_proposed, intensity_proposed = np.loadtxt("intensity_vs_S.dat", delimiter="	", skiprows=1, usecols=(0,1), unpack=True)

	
	intensity_ratio_proposed = min(intensity_proposed)/max(intensity_proposed)
	
	print "Max_intensity/min_intensity for broadest peak = " + str(intensity_ratio_proposed)
	
	
	files = ["log.defect", "intensity_vs_S.dat", "intensity_vs_S.png", "intensity_vs_n", "nsqr_by_wilkens_vs_n.png", "wilkens_plot"]
	
	new_directory_name = "I_vs_S_at_max_width"
	location = os.getcwd() + "/" + str(new_directory_name)
	subprocess.call("rm -r " + new_directory_name, shell=True)
	subprocess.call("mkdir " + str(new_directory_name), shell=True)
	for i in range(len(files)):
		subprocess.call("mv " + os.getcwd() + "/" + files[i] + " " + str(location), shell=True)
	
	
	s_range = s_range_user(s_num, del_s_pernm)
	
	return s_range;
	
	
	
def s_range_user(s_num, del_s_pernm):
	
	s_range_upper = del_s_pernm * 1.0e9
	
	s_range = np.linspace(0, s_range_upper, s_num)
	
	
	f = open(log_name, "a")
	f.write("\n\nFUNCTION CALL: s_range_method_3"
	"\nOUTPUT: "
	"\nS")
	for i in range(len(s_range)):
		f.write( "\n" + str(s_range[i]) )
	
	f.close()

	return s_range
	

def contrast_factor(G, gl_angle):

	individual_C = [0] * len(gl_angle)
	dgl = abs(gl_angle[0] - gl_angle[1])
	
	for i in range(len(gl_angle)):
	
		individual_C[i] = (np.sin(gl_angle[i]) ** 2) * (np.cos(gl_angle[i]) ** 2)
	
	average_C = np.mean(individual_C)
	print "Average contrast factor = " + str(average_C)

	return average_C



def R_e_from_M(M, average_defect_density):

	R_e = M/np.sqrt(average_defect_density)

	print "R_e = " + str(R_e)	

	f = open(log_name, "a")
	f.write("\n\nFUNCTION CALL: R_e_from_M"
	"\nOUTPUT: "
	"\nR_e"
	"\n" + str(R_e) + "")
	f.close()

	return R_e 
	
	
	
def eta(R_e, n):

	eta = (1.0/2.0)*np.exp(-0.25)*n*(1.0/R_e)

	return eta
	

def eta_2(R_e, n):

	eta = n/R_e
	
	return eta
	
	
	
def wilkens_function(R_e, n):

	wilkens_eta = eta_2(R_e, n)

	if wilkens_eta < 0.5:
	
		wilkens_function = -np.log(wilkens_eta) + 2.0 - 0.25 - np.log(2)
	
	elif wilkens_eta >= 0.5:

		wilkens_function = 512.0/(90.0 * np.pi * wilkens_eta) - ((11.0/24.0) + (np.log(2 * wilkens_eta)/4.0) ) * ( 1.0/(wilkens_eta ** 2) )
	
	#Below is the correct, full Wilkens function when eta < 1.0. Unfortunately there is a volume integral involved and it is unclear to what volume it is referring (Wilkens 1970). In the absence of a properly working estimate, the above conditions are instead used outside of their intended parameter space (for eta outside their parameter space, that is).
	"""
	else:
	
		def integrand(x): 
			return np.arcsin((4.0/3.0) * np.pi * (R_e**3))/((4.0/3.0) * np.pi * (R_e**3))
		
		vol_integration = integrate.quad(integrand, 0.0, R_e) 

		
		 #Note that this part of the function is not the function from Wilkens 1970; instead this is an evaluation made by this codes author, included in .
		wilkens_function = -np.log(wilkens_eta) + 2.0 - 0.25 - np.log(2) + 512.0/(90.0 * np.pi * wilkens_eta) - ( (np.sqrt(1 - (wilkens_eta ** 2))/np.pi) * ( (769.0/(180.0*wilkens_eta)) + (wilkens_eta * 41.0/90.0) + ((2.0/90.0) * (wilkens_eta ** 3) ) ) ) - ( (np.arcsin(wilkens_eta)/np.pi) * ( (11.0/(12.0 * (wilkens_eta ** 2))) + (7.0/2.0) + ((wilkens_eta ** 2)/3.0) ) ) + ((wilkens_eta ** 2)/6.0) + ( (2.0/np.pi) * ( 1.0 - (1.0/(4.0 * (wilkens_eta ** 2) * vol_integration[0] 
	"""
	
	return wilkens_function
	
	
	
def create_n_range(exponent_cutoff, samples_per_oscillation, G, B, contrast_factor, defect_density, R_e, s_range):

	#Since our integral	is an exponent multiplied by a cosine, to estimate convergence we must set some arbitrary cutoff, at which point we say that the oscillations around intensity=0 are negligible and no longer contribute to the intensity. We set this with the exponent_cutoff var. This means n_upper will increase with lower contrast factor, lower density, and lower G^2. Once n_upper is determined, it is necessary to determine the period of oscillation, to ensure that points are being taken at regular enough intervals so as to not cause aliasing.
	
	print "\nCreating range of n values for integration..."
	
	Q_value = np.log(exponent_cutoff)
	
	n_prelim = np.linspace(0.0, 5e-6, 1e4+1)
	
	n_prelim = np.delete(n_prelim, 0)
	
	
	nsqr_by_wilkens = [0] * len(n_prelim)
	
	for i in range(len(n_prelim)):

		nsqr_by_wilkens[i] = (n_prelim[i] ** 2) * wilkens_function(R_e, n_prelim[i])

	plotname_nsqr_by_wilkens = "nsqr_by_wilkens_vs_n.png"
	plt.plot(n_prelim, nsqr_by_wilkens)
	plt.ylim(min(nsqr_by_wilkens), max(nsqr_by_wilkens))
	plt.xlim(min(n_prelim), max(n_prelim))
	plt.savefig(plotname_nsqr_by_wilkens)
	plt.close()
	
	gsqr_all = [0] * len(G)
	
	for i in range(len(G)):
		gsqr_all[i] =  (G[i][0] ** 2) + (G[i][1] ** 2) + (G[i][2] ** 2)
		
	gsqr = min(gsqr_all)

	bsqr = (B[0] ** 2) + (B[1] ** 2) + (B[2] ** 2)
	
	if min(contrast_factor) == 0:
	
		smallest_contrast_factor = 0.0 + max(contrast_factor)/(len(contrast_factor)-1)
		
	else:
		
		smallest_contrast_factor = 0.0


	max_nsqr_by_wilkens = Q_value * 2.0 /(np.pi * smallest_contrast_factor * gsqr * bsqr * defect_density)
	
	diff_nsqr_by_wilkens = [0] * len(n_prelim)
	
	for i in range(len(n_prelim)):

		diff_nsqr_by_wilkens[i] = abs(nsqr_by_wilkens[i] - max_nsqr_by_wilkens)
			
				
	n_overshot_fraction = 1.1

	
	n_upper = n_prelim[np.argmin(diff_nsqr_by_wilkens)] * n_overshot_fraction
	
	s_per_oscillation = 1.0/n_upper
	
	s_upper = max(s_range)

	#Since we are taking the largest value of s in s_range, this var represents the maximum possible oscillations of our integral in our range of n. For values of s that are closer to 0, there will be fewer oscillations over n_range.
	max_oscillations_in_s_range = s_upper/s_per_oscillation
	
	n_num = int(samples_per_oscillation * max_oscillations_in_s_range) + 1

	n_range = np.linspace(0.0, n_upper, n_num)
	n_range = np.delete(n_range, 0)
	
	print "\nLength of n_range = " + str(len(n_range))

	return (n_range);
	
	
	
def Q(n, G, B, a_lattice, contrast_factor, average_defect_density, R_e):

	gsqr = ( (G[0] ** 2) + (G[1] ** 2) + (G[2] ** 2) ) / (a_lattice ** 2)

	bsqr = ( (B[0] ** 2) + (B[1] ** 2) + (B[2] ** 2) ) * (a_lattice ** 2) * 0.5

	Q = 0.5 * np.pi * gsqr * bsqr * contrast_factor * average_defect_density * (n ** 2) * wilkens_function(R_e, n)
	
	return Q
	

	
	
def A(n, G, B, a_lattice, contrast_factor, average_defect_density, R_e):

	A = np.exp(-Q(n, G, B, a_lattice, contrast_factor, average_defect_density, R_e))
	
	return A
	
	
	
def intensity(n_range, s_range,  G, B, a_lattice, contrast_factor, average_defect_density, R_e):

	print "\nCreating intensity profile for G = " + str(int(G[0])) + str(int(G[1])) + str(int(G[2])) + ", C = " + str(contrast_factor)

	dn = n_range[1] - n_range[0]
	
	print "Low eta = " + str(eta_2(R_e, n_range[0])) + ""
	print "High eta = " + str(eta_2(R_e, n_range[-1])) + "\n"
	

	intensity_at_n = [0] * len(s_range)
	cos = [0] * len(s_range)
	for i in range(len(s_range)):	
		cos[i] = [0] * len(n_range)
		intensity_at_n[i] = [0] * len(n_range)		

	intensity_at_s = [0] * len(s_range)			

	new_directory_name = "intensity_vs_n"
	subprocess.call("rm -r " + str(new_directory_name), shell=True)
	subprocess.call("mkdir " + str(new_directory_name), shell=True)

	for i in range(len(s_range)):
	
		print str(i+1) + " of " + str(len(s_range))
	
		for j in range(len(n_range)):
	
			cos[i][j] = np.cos(2 * np.pi * n_range[j] * s_range[i])
			
			evaluation = A(n_range[j], G, B, a_lattice, contrast_factor, average_defect_density, R_e) * cos[i][j]

			intensity_at_n[i][j] = 2 * dn * evaluation

			
		datafile_intensity_vs_n = "intensity_vs_n@S=" + str(s_range[i]) + ".dat"
		g = open(datafile_intensity_vs_n, "w")
		g.write("n\tintensity\tintegrand\tA\tcos")
		for j in range(len(n_range)):			
			g.write("\n" + str(n_range[j]) + "\t" + str(intensity_at_n[i][j]))

		g.close()		
		location = os.getcwd() + "/" + str(new_directory_name)
		subprocess.call("mv " + os.getcwd() + "/" + datafile_intensity_vs_n + " " + str(location), shell=True)
				
		intensity_at_s[i] = sum(intensity_at_n[i])
	
	
	f = open(log_name, "a")
	f.write("\n\nFUNCTION CALL: intensity"
	"\nLow eta = " + str(eta_2(R_e, n_range[0])) + ""
	"\nHigh  eta = " + str(eta_2(R_e, n_range[-1])) + ""
	"\nOUTPUT: "
	"\nS\t\tintensity_at_s\t\tdn\n")
	for i in range(len(s_range)):	
		if i == 0:
			f.write( str(s_range[i]) + "\t" + str(intensity_at_s[i]) + "\t" + str(dn) + "\n" )
			
		else:
			f.write( str(s_range[i]) + "\t" + str(intensity_at_s[i]) + "\n" )			
				
		

	f.close()		


	
	datafile_intensity_vs_s = "intensity_vs_S.dat"
	g = open(datafile_intensity_vs_s, "w")
	g.write("S\tintensity")
	for j in range(len(s_range)):			
		g.write("\n" + str(s_range[j]) + "\t" + str(intensity_at_s[j]))

	g.close()


	plotname_intensity_vs_s = "intensity_vs_S.png"
	plt.scatter(s_range, intensity_at_s)
	plt.ylim([min(intensity_at_s), max(intensity_at_s)])
	plt.savefig(plotname_intensity_vs_s)
	plt.close()


# Creates an approximation of the possible effect of the intensity profile on the Debye-Waller effect. This is done by reading in the intensity vs S data for each plot, and simulating 
def estimate_consequences_debye_waller(G, contrast_factor, plot_normI_vs_C, plot_I_vs_C, plot_I_vs_C_no_zero, plot_normI_vs_C_no_zero):

	print "\nEstimating consequences of defects on Debye-Waller Temperature measurement..."
	new_directory = os.getcwd() + "/effect_on_debye_waller/"
	subprocess.call("rm -r " + new_directory, shell=True)
	subprocess.call("mkdir " + new_directory, shell=True)
	


	label = [0] * len(contrast_factor)
	files = [0] * len(G)
	gsqr = [0] * len(G)

	for i in range(len(files)):

		files[i] = [0] * len(contrast_factor)


	for j in range(len(G)):
	
		gsqr[j] = (int(G[j][0]) ** 2) + (int(G[j][1]) ** 2) + (int(G[j][2]) ** 2)
		 
		print gsqr
		
		for i in range(len(contrast_factor)):

			files[j][i] = str(os.getcwd()) + "/" + str(int(G[j][0]))  + str(int(G[j][1])) + str(int(G[j][2])) + "/contrast_" + str(contrast_factor[i]) + "/intensity_vs_S.dat"
			label[i] = "Contrast " + str(contrast_factor[i])


	intensity = [0] * len(G)
	full_peak_intensity = [0] * len(G)
	full_peak_s = [0] * len(G)
	half_intensity = [0] * len(G)
	del_intensity = [0] * len(G)
	hwhm = [0] * len(G)
	norm_integrated_intensity = [0] * len(G)
	integrated_intensity = [0] * len(G)
	intensity_over_hwhm = [0] * len(G)
	norm_intensity_over_hwhm = [0] * len(G)
	s = [0] * len(G)

	for i in range(len(files)):
		intensity[i] = [0] * len(contrast_factor)
		full_peak_intensity[i] = [0] * len(contrast_factor)
		full_peak_s[i] = [0] * len(contrast_factor)
		half_intensity[i] = [0] * len(contrast_factor)
		del_intensity[i] = [0] * len(contrast_factor)
		hwhm[i] = [0] * len(contrast_factor)
		norm_integrated_intensity[i] = [0] * len(contrast_factor)
		integrated_intensity[i] = [0] * len(contrast_factor)
		intensity_over_hwhm[i] = [0] * len(contrast_factor)
		#norm_intensity_over_fwhm = [0] * len(contrast_factor)
		s[i] = [0] * len(contrast_factor)


	
	for j in range(len(G)):

		print G[j]

		for i in range(len(contrast_factor)):
		
	
			s[j][i], intensity[j][i] = np.loadtxt(files[j][i], skiprows=1, unpack=True, usecols=(0,1))


			intermediate = np.delete(intensity[j][i], 0)


			full_peak_intensity[j][i] = np.concatenate((intensity[j][i][::-1], intermediate), axis=0)


			intermediate_2 = np.delete(s[j][i], 0)

			full_peak_s[j][i] = np.concatenate((s[j][i][::-1], -intermediate_2), axis=0)

		
			central_s = 0.0
		
			sigma_guess = s[j][i][int(len(s[j][i])/2.0)]

		
			#popt = fit_gaussian(full_peak_s[j][i], full_peak_intensity[j][i], max(intensity[j][i]), central_s, sigma_guess)
		
			#hwhm[j][i] = np.sqrt(2 * np.log(2.0)) * popt[2]
			"""
			y = gaussian(s[j][i], popt[0], popt[1], popt[2])
			plt.plot(s[j][i], y)
			plt.scatter(s[j][i], intensity[j][i])
			intersection = gaussian(hwhm[j][i], popt[0], popt[1], popt[2])
			plt.scatter(hwhm[j][i], intersection, color="r")
			plt.show()
			plt.close()
			"""
			integrated_intensity[j][i] = sum(intensity[j][i])
			norm_integrated_intensity[j][i] = integrated_intensity[j][i]/max(integrated_intensity[j])
			#intensity_over_hwhm[j][i] = abs(integrated_intensity[j][i]/hwhm[j][i])

		
		#norm_intensity_over_hwhm[j] = intensity_over_hwhm[j]/max(intensity_over_hwhm[j])

	max_gsqr_ind = np.argmax(gsqr)
	min_gsqr_ind = np.argmin(gsqr)
	
	
	if plot_normI_vs_C == True:
	
		plt.figure(figsize=(16,9))
	
		for j in range(len(G)):

			colour = np.random.rand(1,3)
			plt.scatter(contrast_factor, norm_integrated_intensity[j], color=colour, label=G[j])
			plt.ylabel("Normalised integrated intensity (arb.)")
			plt.xlabel("Contrast Factor")
			
			print "Min(norm_intensity)/max(norm_intensity) for " + str(int(G[j][0])) +str(int(G[j][1])) + str(int(G[j][2])) + " = " + str(min(norm_integrated_intensity[j])/max(norm_integrated_intensity[j]))	


		plt.title( "Min(norm_intensity)/max(norm_intensity) for " + str(int(G[max_gsqr_ind][0])) +str(int(G[max_gsqr_ind][1])) + str(int(G[max_gsqr_ind][2])) + " = " + str(min(norm_integrated_intensity[max_gsqr_ind])/max(norm_integrated_intensity[max_gsqr_ind]))	)
		plt.legend(loc="upper right")
		plt.ylim(0.9*min(norm_integrated_intensity[max_gsqr_ind]), 1.1*max(norm_integrated_intensity[max_gsqr_ind]))
		
		plt.savefig(new_directory + "normI_vs_C.png", dpi=200)
		plt.close()
		

		
	if plot_I_vs_C == True:
	
		plt.figure(figsize=(16,9))
	
		for j in range(len(G)):

			colour = np.random.rand(1,3)
			plt.scatter(contrast_factor, integrated_intensity[j], color=colour, label=G[j])
			plt.ylabel("Integrated intensity (arb.)")
			plt.xlabel("Contrast Factor")
			
			print "Min(intensity)/max(intensity) for " + str(int(G[j][0])) +str(int(G[j][1])) + str(int(G[j][2])) + " = " + str(min(integrated_intensity[j])/max(integrated_intensity[j]))	


		plt.title( "Min(intensity)/max(intensity) for " + str(int(G[max_gsqr_ind][0])) +str(int(G[max_gsqr_ind][1])) + str(int(G[max_gsqr_ind][2])) + " = " + str(min(integrated_intensity[max_gsqr_ind])/max(integrated_intensity[max_gsqr_ind]))	)
		plt.legend(loc="upper right")
		plt.ylim(0.9*min(integrated_intensity[max_gsqr_ind]), 1.1*max(integrated_intensity[max_gsqr_ind]))
		plt.savefig(new_directory + "I_vs_C.png", dpi=80)
		plt.close()
		
		
		
		
	if plot_I_vs_C_no_zero == True:
	
		plt.figure(figsize=(16,9))
		contrast_factor_no_zero = np.delete(contrast_factor, 0)

		integrated_intensity_no_zero = [0] * len(G)
	
		for j in range(len(G)):

			colour = np.random.rand(1,3)

			if np.count_nonzero(contrast_factor) == len(contrast_factor):

				integrated_intensity_no_zero[j] = integrated_intensity[j]
				
			else:

				integrated_intensity_no_zero[j] = np.delete(integrated_intensity[j], max(integrated_intensity[j]) )



			plt.scatter(contrast_factor_no_zero, integrated_intensity_no_zero[j], color=colour, label=G[j])
			plt.ylabel("Integrated intensity (arb.)")
			plt.xlabel("Contrast Factor")
			
			print "Min(intensity)/max(intensity) for " + str(int(G[j][0])) +str(int(G[j][1])) + str(int(G[j][2])) + " = " + str(min(integrated_intensity_no_zero[j])/max(integrated_intensity_no_zero[j]))	


		plt.title( "Min(intensity)/max(intensity) for " + str(int(G[max_gsqr_ind][0])) +str(int(G[max_gsqr_ind][1])) + str(int(G[max_gsqr_ind][2])) + " = " + str(min(integrated_intensity_no_zero[max_gsqr_ind])/max(integrated_intensity_no_zero[max_gsqr_ind]))	)
		plt.legend(loc="upper right")
		plt.ylim(0.9*min(integrated_intensity_no_zero[max_gsqr_ind]), 1.1*max(integrated_intensity_no_zero[min_gsqr_ind]))
		plt.savefig(new_directory + "I_vs_C_no_zero.png", dpi=80)
		plt.close()
	


	if plot_normI_vs_C_no_zero == True:
	
		plt.figure(figsize=(16,9))
		contrast_factor_no_zero = np.delete(contrast_factor, 0)

		integrated_intensity_no_zero = [0] * len(G)
		norm_integrated_intensity_no_zero = [0] * len(G)
	
		for j in range(len(G)):

			colour = np.random.rand(1,3)

			if np.count_nonzero(contrast_factor) == len(contrast_factor):

				norm_integrated_intensity_no_zero[j] = integrated_intensity[j]
				
			else:

				norm_integrated_intensity_no_zero[j] = np.delete(integrated_intensity[j], max(integrated_intensity[j]) )

			print norm_integrated_intensity_no_zero
			
			normalise_constant = max(norm_integrated_intensity_no_zero[j])

			for i in range(len(norm_integrated_intensity_no_zero[j])):

				norm_integrated_intensity_no_zero[j][i] = norm_integrated_intensity_no_zero[j][i]/normalise_constant
				
			print norm_integrated_intensity_no_zero

			plt.scatter(contrast_factor_no_zero, norm_integrated_intensity_no_zero[j], color=colour, label=G[j])
			plt.ylabel("Normalised integrated intensity (arb.)")
			plt.xlabel("Contrast Factor")
			
			print "Min(intensity)/max(intensity) for " + str(int(G[j][0])) +str(int(G[j][1])) + str(int(G[j][2])) + " = " + str(min(norm_integrated_intensity_no_zero[j])/max(norm_integrated_intensity_no_zero[j]))	


		plt.title( "Min(intensity)/max(intensity) for " + str(int(G[max_gsqr_ind][0])) +str(int(G[max_gsqr_ind][1])) + str(int(G[max_gsqr_ind][2])) + " = " + str(min(norm_integrated_intensity_no_zero[max_gsqr_ind])/max(norm_integrated_intensity_no_zero[max_gsqr_ind])) )
		plt.legend(loc="upper right")
		plt.ylim(0.9*min(norm_integrated_intensity_no_zero[min_gsqr_ind]), 1.1*max(norm_integrated_intensity_no_zero[min_gsqr_ind]))
		plt.savefig(new_directory + "normI_vs_C_no_zero.png", dpi=80)
		plt.close()	
		

	if plot_I_vs_theta == True:

		plt.figure(figsize=(16,9))
	
		for j in range(len(G)):

			colour = np.random.rand(1,3)
			plt.scatter(contrast_factor, integrated_intensity[j], color=colour, label=G[j])
			plt.ylabel("Integrated intensity (arb.)")
			plt.xlabel("Contrast Factor")
			
			print "Min(intensity)/max(intensity) for " + str(int(G[j][0])) +str(int(G[j][1])) + str(int(G[j][2])) + " = " + str(min(integrated_intensity[j])/max(integrated_intensity[j]))	


		plt.title( "Min(intensity)/max(intensity) for " + str(int(G[max_gsqr_ind][0])) +str(int(G[max_gsqr_ind][1])) + str(int(G[max_gsqr_ind][2])) + " = " + str(min(integrated_intensity[max_gsqr_ind])/max(integrated_intensity[max_gsqr_ind]))	)
		plt.legend(loc="upper right")
		plt.ylim(0.9*min(integrated_intensity[max_gsqr_ind]), 1.1*max(integrated_intensity[max_gsqr_ind]))
		plt.savefig(new_directory + "I_vs_C.png", dpi=80)
		plt.close()	
	
	
	return
	
	
def transform_G_to_angle(G):

	# When theta = pi/2, g and l are perpindicular, and C = 0. When theta = 0, g and l are parallel and C = 0. 
	# In order to transform contrast to angle, one must first assign a direction with which the sample is facing. In this case, the sample is a single crystal containing only screw dislocations. The sample will be oriented normal to the 001 direction as is conventional. All reflections may now be assigned an angle relative to the 001 vector.
	
	reference_direction = [0.0, 0.0, 1.0] 
	reference_length = np.sqrt( (reference_direction[0] ** 2) + (reference_direction[1] ** 2) + (reference_direction[2] ** 2) )

	print "\nReference direction = (" + str(reference_direction[0]) + ", " + str(reference_direction[1]) + ", " + str(reference_direction[2]) + ")\n"
	
	dot_product = [0] * len(G)
	G_length = [0] * len(G)
	theta_range = [0] * len(G)
	
	for i in range(len(G)):

		G_length[i] = np.sqrt( (G[i][0] ** 2) + (G[i][1] ** 2) + (G[i][2] ** 2) ) 
		dot_product[i] = np.dot(G[i], reference_direction)
		theta_range[i] = np.arccos(dot_product[i]/(G_length[i] * reference_length))
		print "Angle for " + str(int(G[i][0])) + str(int(G[i][1])) + str(int(G[i][2])) + " = " + str(theta_range[i]) + "\n"

	return theta_range
	

def make_contrast_from_theta_screw(G, theta_range):

	contrast_factor = [0] * len(theta_range)
	
	for i in range(len(theta_range)):
	
		contrast_factor[i] = (np.cos(theta_range[i]) ** 2) * (np.sin(theta_range[i]) ** 2)
		print "Contrast factor for " + str(int(G[i][0])) + str(int(G[i][1])) + str(int(G[i][2])) + " = " + str(contrast_factor[i]) + "\n"


	return contrast_factor
	
	
def apply_to_debye_waller(file_intensity_profile):

	 #= np.loadtxt(file_intensity_profile)

	return
	
	
	
def make_plots(comparison_with_balogh, wilkens_plot, intensity_vs_n_plot, cos_vs_n_plot, n_range, R_e, s_num, G, s_range):


	if comparison_with_balogh == True:
	
		print "Plotting comparison of calculated intensities with Balogh..."
	
		balogh_data = ["balogh_accurate_5e12.csv", "balogh_data_density=4.csv", "balogh_data_density=3.5.csv"]	
		balogh_label = ["Balogh Accurate 5e12", "Balogh 4e12", "Balogh 3.5e12"]	
		balogh_intensity_raw = [0] * len(balogh_data)
		balogh_intensity_renorm = [0] * len(balogh_data)
		balogh_s = [0] * len(balogh_data)
		balogh_colours = ["b", "c", "k"]
		
		plotname_balogh = "compare_with_balogh.png"
		plotname_balogh_log = "log_compare_with_balogh.png"
		
		for i in range(len(balogh_data)):

			csv_data = np.genfromtxt(balogh_data[i], delimiter="    ")
	
			balogh_intensity_raw[i] = csv_data[:,1]
			balogh_intensity_renorm[i] = balogh_intensity_raw[i]/max(balogh_intensity_raw[i])
			balogh_s[i] = csv_data[:,0]
			
			plt.scatter(balogh_s[i], balogh_intensity_renorm[i], color=balogh_colours[i], label=balogh_label[i])
			plt.scatter(balogh_s[i], balogh_intensity_raw[i], color=balogh_colours[i], label=balogh_label[i])
	
			
		calculated_data = ["intensity_vs_S.dat"]
		calculated_label = ["Calculated Peak"]
		calc_colours = ["r", "y", "g"]

		for i in range(len(calculated_data)):
	
			calc_s, calc_intensity = np.loadtxt(calculated_data[i], skiprows=1, unpack=True, usecols=(0,1))

			norm_calc_intensity = calc_intensity/max(calc_intensity)

			G_length_pernm = np.sqrt( (G[0] ** 2) + (G[1] ** 2) + (G[2] ** 2) ) 		
			calc_s_nm_1 = G_length_pernm - (calc_s / 1.0e9) 
			

							
			plt.scatter(calc_s_nm_1, norm_calc_intensity, color=calc_colours[i], label=calculated_label[i])			
		
		xlim = [min(calc_s_nm_1), G_length_pernm]
		plt.legend(loc=2)
		plt.ylim(0.0, 1.1)	
		plt.xlim(xlim[0], xlim[1])
		plt.ylabel("Intensity arb.")
		plt.xlabel("S nm^-1")
		plt.suptitle("Comparison of peak intensity vs. S for calculated values and Balogh's data")
		plt.grid(True)
		plt.savefig(plotname_balogh)
		plt.yscale("log")
		plt.ylim(0.001, 1.1)
		plt.savefig(plotname_balogh_log)

		plt.close()
		

	if wilkens_plot == True:
		
		print "Plotting Wilkens function..."
		
		y_range = [0] * len(n_range)
		
		for i in range(len(n_range)):
			
			y_range[i] = wilkens_function(R_e, n_range[i])
				
		wilkens_plotname = "wilkens_plot.png"	
		plt.scatter(n_range, y_range)
		plt.xlim(n_range[0], n_range[-1])
		plt.ylim(0.0, max(y_range))
		plt.savefig(wilkens_plotname)
		plt.close()
		
		new_directory_name = "wilkens_plot"
		location = os.getcwd() + "/" + str(new_directory_name)
		subprocess.call("mkdir " + str(new_directory_name), shell=True)
		subprocess.call("mv " + os.getcwd() + "/" + wilkens_plotname + " " + str(location), shell=True)
	
		
	if intensity_vs_n_plot == True:
	
		print "Plotting intensity vs n..."
	
		for i in range(len(s_range)):
			
			desired_directory = "intensity_vs_n"
			datafile_intensity_vs_n = os.getcwd() + "/" + desired_directory + "/intensity_vs_n@S=" + str(s_range[i]) + ".dat"
			
			plotname_intensity_vs_n = "intensity_vs_n@S=" + str(s_range[i]) + ".png"
			
			n_range, intensity_at_n = np.loadtxt(datafile_intensity_vs_n, skiprows=1, usecols=(0,1), unpack=True)

	
			gnuplot_input_filename = "in.intensity_vs_n@S=" + str(s_range[i])
			
			f = open(gnuplot_input_filename, "w")
			
			f.write("set terminal png size 1600,1200 enhanced font 'Helvetica,20'"
			"\nset output '" + str(plotname_intensity_vs_n)  + "'"
			"\nplot '" + str(datafile_intensity_vs_n) + "' using 1:2"
			)
			
			f.close()
			
			subprocess.call("gnuplot '" + str(gnuplot_input_filename) + "'", shell=True)
			
			
			
			subprocess.call("mv " + os.getcwd() + "/" + plotname_intensity_vs_n + " " + os.getcwd() + "/" + str(desired_directory), shell=True)
			
			
			subprocess.call("mv " + os.getcwd() + "/" + gnuplot_input_filename + " " + os.getcwd() + "/" + str(desired_directory), shell=True)
			

			
	if cos_vs_n_plot == True:

		print "Plotting cosine function vs n..."

		
		new_directory_name = "cos_vs_n"
		location = os.getcwd() + "/" + str(new_directory_name)
		subprocess.call("rm -r " + str(new_directory_name), shell=True)
		subprocess.call("mkdir " + str(new_directory_name), shell=True)

		
		for j in range(len(s_range)):
	
			cos = [0] * len(n_range)
	
			for i in range(len(n_range)):
			
				cos[i] = np.cos(n_range[i] * s_range[j])
			
			plotname_cos = "cosine_vs_n@S=" + str(s_range[j]) + ".png"
			plt.scatter(n_range, cos)
			plt.xlim(n_range[0], n_range[-1])
			plt.savefig(plotname_cos)
			plt.close()
			
			subprocess.call("mv " + os.getcwd() + "/" + plotname_cos + " " + str(location), shell=True)
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def run(calculate_intensities, s_num, del_s_pernm, M, average_defect_density, n_range, G, B, a_lattice, contrast_factor, comparison_with_balogh, wilkens_plot, cos_vs_n_plot, intensity_vs_n_plot):

	initialisation(G, B,a_lattice,s_num,n_range,del_s_pernm,M,average_defect_density,contrast_factor)
	s_range = s_range_user(s_num, del_s_pernm)
	R_e = R_e_from_M(M, average_defect_density)
	
	if calculate_intensities == True:
		intensity(n_range, s_range,  G, B, a_lattice, contrast_factor, average_defect_density, R_e)
		
	make_plots(comparison_with_balogh, wilkens_plot, intensity_vs_n_plot, cos_vs_n_plot, n_range, R_e, s_num, G, s_range)

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Loops through ranges of input values. 

def runner():

	t_start = time.clock()

	calculate_intensities = True
	
	#contrast_factor = np.linspace(0.0, 0.25, 11)
	G = [[2.0, 0, -2.0], [1.0,1.0,1.0], [0.0, 0.0, 2.0], [2.0,0.0,0.0], [2.0, 2.0,2.0], [1.0,1.0,3.0], [1.0,3.0,3.0], [3.0, 3.0, 3.0], [2.0, 2.0,4.0], [2.0, 4.0, 4.0], [4.0, 4.0, 4.0], [3.0, 3.0, 5.0], [3.0, 5.0, 5.0], [5.0, 5.0, 5.0]]
	B = [0.5, 0, -0.5]
	a_lattice = 3.615e-10
	s_num = 11
	del_s_pernm = 0.6
	M = 13.5
	average_defect_density = 1.0e16
	
	n_range_guess = np.linspace(1e-20, 5e-7, 1e3) # Used to estimate the range of s values.
	exponent_cutoff = 1e3 # Controls convergence of integral vs n.
	samples_per_oscillation = 10.0 # Controls number of points recorded of integral vs n.
	
	
	# Determines which plots are made at the end (the "make_plots" function).
	comparison_with_balogh = False
	wilkens_plot = False
	cos_vs_n_plot = False
	intensity_vs_n_plot = False

	# These plots are made in the "estimate_consequences_debye_waller" function.
	plot_normI_vs_C = True
	plot_I_vs_C = True
	plot_I_vs_C_no_zero = True
	plot_normI_vs_C_no_zero = True
	
	theta_range = transform_G_to_angle(G)
	contrast_factor = make_contrast_from_theta_screw(G, theta_range)
	
	exit()
	
	s_range, n_range = build_boundary_profile(M, average_defect_density, s_num, del_s_pernm, exponent_cutoff, samples_per_oscillation, G, B, contrast_factor, n_range_guess, a_lattice) 
	
	
	for j in range(len(G)):
	
		G_directory = os.getcwd() + "/" + str(int(G[j][0])) + str(int(G[j][1])) + str(int(G[j][2]))
		subprocess.call("rm -r " + G_directory, shell=True)
		subprocess.call("mkdir " + G_directory, shell=True)	
	
		for i in range(len(contrast_factor)):
	
			contrast_directory = "contrast_" + str(contrast_factor[i])
			subprocess.call("rm -r " + contrast_directory, shell=True)
			subprocess.call("mkdir " + contrast_directory, shell=True)
		
		
			run(calculate_intensities, s_num, del_s_pernm, M, average_defect_density, n_range, G[j], B, a_lattice, contrast_factor[i], comparison_with_balogh, wilkens_plot, cos_vs_n_plot, intensity_vs_n_plot)
		
		
			final_directory = str(os.getcwd()) + "/" + contrast_directory
			files_to_move = ["compare_with_balogh.png", "intensity_vs_S.dat", "intensity_vs_S.png", "log.defect", "log_compare_with_balogh.png"]
			folders_to_move = ["intensity_vs_n"]
		
			for a in range(len(files_to_move)):
		
				subprocess.call("mv " + str(os.getcwd()) + "/" + files_to_move[a] + " " + final_directory, shell=True)
			
			for a in range(len(folders_to_move)):
	
				subprocess.call("mv " + str(os.getcwd()) + "/" + folders_to_move[a] + "/ " + final_directory, shell=True)
		
			
			subprocess.call("mv " + str(os.getcwd()) + "/" + contrast_directory + "/ " + G_directory, shell=True)
	
	
	estimate_consequences_debye_waller(G, contrast_factor, plot_normI_vs_C, plot_I_vs_C, plot_I_vs_C_no_zero, plot_normI_vs_C_no_zero)
	
			
	t_end = time.clock()		
			
	t_total = t_end - t_start

	f = open(log_name, "a")
	f.write("\nDone!\nThe run took " + str(t_total) + " s to complete.")
	f.close()
			
	print "\nDone!\nThe run took " + str(t_total) + " s to complete.\n"
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

runner()
