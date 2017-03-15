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
import subprocess
import os

log_name = "log.defect"
program_name = "defect_3.py"



def initialisation():

	start_time = time.localtime()
	
	date_string = str(start_time[2]) + "/" + str(start_time[1]) + "/" + str(start_time[0])
	time_string = str(start_time[3]) + "h " + str(start_time[4]) + "m " + str(start_time[5]) + "s"


	f = open(log_name, "w")
	f.write(program_name + "\n" + date_string + "\n" + time_string)
	f.close()
	return;
	
	

	
def s_range(s_num, del_s_pernm):
	
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
	
	wilkens_function = -np.log(wilkens_eta) + 2.0 - 0.25 - np.log(2)
	
	return wilkens_function
	
	
	
def Q(n, G, B, a_lattice, contrast_factor, average_defect_density, R_e):

	gsqr = ( (G[0] ** 2) + (G[1] ** 2) + (G[2] ** 2) ) / (a_lattice ** 2)

	bsqr = ( (B[0] ** 2) + (B[1] ** 2) + (B[2] ** 2) ) * (a_lattice ** 2) * 0.5

	Q = 0.5 * np.pi * gsqr * bsqr * contrast_factor * average_defect_density * (n ** 2) * wilkens_function(R_e, n)
	
	return Q
	
	
	
def A(n, G, B, a_lattice, contrast_factor, average_defect_density, R_e):

	A = np.exp(-Q(n, G, B, a_lattice, contrast_factor, average_defect_density, R_e))
	
	return A
	
	
	
def intensity(n_lower, n_upper, n_num, s_range,  G, B, a_lattice, contrast_factor, average_defect_density, R_e):

	n_range = np.linspace(n_lower, n_upper, n_num)
	dn = n_range[1] - n_range[0]
	
	print "Low eta = " + str(eta_2(R_e, n_lower)) + ""
	print "High eta = " + str(eta_2(R_e, n_upper)) + "\n"
	

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
	"\nLow eta = " + str(eta_2(R_e, n_lower)) + ""
	"\nHigh  eta = " + str(eta_2(R_e, n_upper)) + ""
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

	
	
def make_plots(comparison_with_balogh, wilkens_plot, intensity_vs_n_plot, cos_vs_n_plot, n_lower, n_upper, n_num, R_e, s_num, G):


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
			print balogh_intensity_raw
			print balogh_intensity_renorm
			
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
		
		n_range = np.linspace(n_lower, n_upper, n_num)
		y_range = [0] * len(n_range)
		
		for i in range(len(n_range)):
			
			y_range[i] = wilkens_function(R_e, n_range[i])
				
		wilkens_plotname = "wilkens_plot.png"	
		plt.scatter(n_range, y_range)
		plt.xlim(n_lower, n_upper)
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
			
			n_range, intensity_at_n = np.loadtxt(datafile_intensity_vs_n, skiprows=1, usecols=(0,1), unpack=True)
	
			plotname_intensity_vs_n = "intensity_vs_n@S=" + str(s_range[i]) + ".png"		
			plt.scatter(n_range, intensity_at_n)
			plt.xlim([0, n_upper])
			plt.ylim([min(intensity_at_n), max(intensity_at_n)])
			plt.savefig(plotname_intensity_vs_n)	
			plt.close()
			
			subprocess.call("mv " + os.getcwd() + "/" + plotname_intensity_vs_n + " " + os.getcwd() + "/" + str(desired_directory), shell=True)
			
			
	if cos_vs_n_plot == True:

		print "Plotting cosine function vs n..."

		n_range = np.linspace(n_lower, n_upper, n_num)
		
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
# Input parameters.

G = [2.0, 0, -2.0]
B = [0.5, 0, -0.5]
gl_angle = 0
a_lattice = 3.615e-10
s_num = 21
n_lower = 1.0/(1.0e25)
n_upper = 2.0e-6
n_num = 1.0e6 + 1
del_s_pernm = 0.003
M = 13.5
average_defect_density = 3.5e12

# Determines which plots are made at the end (the "make_plots" function).
comparison_with_balogh = True
wilkens_plot = False
cos_vs_n_plot = False
intensity_vs_n_plot = False

# Runs the code.
initialisation()
s_range = s_range(s_num, del_s_pernm)
R_e = R_e_from_M(M, average_defect_density)
contrast_factor = 0.25#contrast_factor(G, gl_angle)
intensity(n_lower, n_upper, n_num, s_range,  G, B, a_lattice, contrast_factor, average_defect_density, R_e)
make_plots(comparison_with_balogh, wilkens_plot, intensity_vs_n_plot, cos_vs_n_plot, n_lower, n_upper, n_num, R_e, s_num, G)
