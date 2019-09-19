import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def plot_sci(M, S_Mpu, S_Mpr, S_Msu, S_Msr, S_Mdu, S_Mdr,\
				C_Mpu, C_Mpr, C_Msu, C_Msr, C_Mdu, C_Mdr,\
				I_Mpu, I_Mpr, I_Msu, I_Msr, I_Mdu, I_Mdr, t):

	print('Equilibrium population in ' + str(M))
	print(S_Mpu[-1] + S_Mpr[-1] + S_Msu[-1] + S_Msr[-1] + S_Mdu[-1] + S_Mdr[-1] + C_Mpu[-1] + C_Mpr[-1] + C_Msu[-1] + C_Msr[-1] + C_Mdu[-1] + C_Mdr[-1] + I_Mpu[-1] + I_Mpr[-1] + I_Msu[-1] + I_Msr[-1] + I_Mdu[-1] + I_Mdr[-1])

	S = [S_Mpu, S_Mpr, S_Msu, S_Msr, S_Mdu, S_Mdr]
	C = [C_Mpu, C_Mpr, C_Msu, C_Msr, C_Mdu, C_Mdr]
	I = [I_Mpu, I_Mpr, I_Msu, I_Msr, I_Mdu, I_Mdr]
	
	if M == 'D':
		title_str = ['Deliverty ward (PHCs & Urban)',
					 'Deliverty ward (PHCs & Rural)',
					 'Deliverty ward (CHCs & Urban)',
					 'Deliverty ward (CHCs & Rural)',
					 'Deliverty ward (DHs & Urban)',
					 'Deliverty ward (DHs & Rural)']

	if M == 'B':
		 title_str = ['Neonatal Care (PHCs & Urban)',
					 'Neonatal Care (PHCs & Rural)',
					 'Neonatal Care (CHCs & Urban)',
					 'Neonatal Care (CHCs & Rural)',
					 'Neonatal Care (DHs & Urban)',
					 'Neonatal Care (DHs & Rural)']

	# ploting...
	pos = [321, 322, 323, 324, 325, 326]
	plt.figure(facecolor='w', figsize=(15, 15))
	for i in range(0,5+1):
		plt.subplot(pos[i])
		plt.plot(t, S[i], 'b', alpha=0.5, lw=2, label='Susceptible (S)')
		plt.plot(t, C[i], 'g', alpha=0.5, lw=2, label='Colonized (C)')
		plt.plot(t, I[i], 'r', alpha=0.5, lw=2, label='Infected (I)')
		plt.xlabel('Time / Days')
		plt.ylabel('Number of population S/C/I (1000s)')
		# plt.ylim(-.1,1.0)
		plt.grid(b=True, which='major', c='w', lw=2, ls='-')
		plt.title(title_str[i],fontsize=18)
		plt.legend(loc='upper right')
		
	plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)
	plt.show()


def cal_mean(alist, nparts, num_samples = 100):
	l = []
	for j in range(num_samples): #50 is for the plotting purpose
		cal = 0
		#nparts = 10
		for i in range(nparts):
			cal += alist[i][j]  
		mean = cal/nparts
		l.append(mean)
		
	return l
	
def cal_sci_means(t, nparts, S_Dpu_list, C_Dpu_list, I_Dpu_list, S_Bpu_list, C_Bpu_list, I_Bpu_list, H_Dpu_list, H_Bpu_list, T_Dpu_list, F_Dpu_list, F_Bpu_list, X_Dpu_list, X_Bpu_list, M_Dpu_list, M_Bpu_list, S_Dsu_list, C_Dsu_list, I_Dsu_list, S_Bsu_list, C_Bsu_list, I_Bsu_list, H_Dsu_list, H_Bsu_list, T_Dsu_list, F_Dsu_list, F_Bsu_list, X_Dsu_list, X_Bsu_list, M_Dsu_list, M_Bsu_list, S_Ddu_list, C_Ddu_list, I_Ddu_list, S_Bdu_list, C_Bdu_list, I_Bdu_list, H_Ddu_list, H_Bdu_list, T_Ddu_list, F_Ddu_list, F_Bdu_list, X_Ddu_list, X_Bdu_list, M_Ddu_list, M_Bdu_list,\
					S_Dpr_list, C_Dpr_list, I_Dpr_list, S_Bpr_list, C_Bpr_list, I_Bpr_list, H_Dpr_list, H_Bpr_list, T_Dpr_list, F_Dpr_list, F_Bpr_list, X_Dpr_list, X_Bpr_list, M_Dpr_list, M_Bpr_list, S_Dsr_list, C_Dsr_list, I_Dsr_list, S_Bsr_list, C_Bsr_list, I_Bsr_list, H_Dsr_list, H_Bsr_list, T_Dsr_list, F_Dsr_list, F_Bsr_list, X_Dsr_list, X_Bsr_list, M_Dsr_list, M_Bsr_list, S_Ddr_list, C_Ddr_list, I_Ddr_list, S_Bdr_list, C_Bdr_list, I_Bdr_list, H_Ddr_list, H_Bdr_list, T_Ddr_list, F_Ddr_list, F_Bdr_list, X_Ddr_list, X_Bdr_list, M_Ddr_list, M_Bdr_list):
	# urban
	S_Dpu_mean = np.array(cal_mean(S_Dpu_list, nparts)) 
	C_Dpu_mean = np.array(cal_mean(C_Dpu_list, nparts)) 
	I_Dpu_mean = np.array(cal_mean(I_Dpu_list, nparts)) 
	S_Bpu_mean = np.array(cal_mean(S_Bpu_list, nparts)) 
	C_Bpu_mean = np.array(cal_mean(C_Bpu_list, nparts)) 
	I_Bpu_mean = np.array(cal_mean(I_Bpu_list, nparts)) 
	H_Dpu_mean = np.array(cal_mean(H_Dpu_list, nparts)) 
	H_Bpu_mean = np.array(cal_mean(H_Bpu_list, nparts)) 
	T_Dpu_mean = np.array(cal_mean(T_Dpu_list, nparts)) 
	F_Dpu_mean = np.array(cal_mean(F_Dpu_list, nparts)) 
	F_Bpu_mean = np.array(cal_mean(F_Bpu_list, nparts)) 
	X_Dpu_mean = np.array(cal_mean(X_Dpu_list, nparts)) 
	X_Bpu_mean = np.array(cal_mean(X_Bpu_list, nparts)) 
	M_Dpu_mean = np.array(cal_mean(M_Dpu_list, nparts)) 
	M_Bpu_mean = np.array(cal_mean(M_Bpu_list, nparts))
	S_Dsu_mean = np.array(cal_mean(S_Dsu_list, nparts)) 
	C_Dsu_mean = np.array(cal_mean(C_Dsu_list, nparts)) 
	I_Dsu_mean = np.array(cal_mean(I_Dsu_list, nparts)) 
	S_Bsu_mean = np.array(cal_mean(S_Bsu_list, nparts)) 
	C_Bsu_mean = np.array(cal_mean(C_Bsu_list, nparts)) 
	I_Bsu_mean = np.array(cal_mean(I_Bsu_list, nparts)) 
	H_Dsu_mean = np.array(cal_mean(H_Dsu_list, nparts)) 
	H_Bsu_mean = np.array(cal_mean(H_Bsu_list, nparts)) 
	T_Dsu_mean = np.array(cal_mean(T_Dsu_list, nparts)) 
	F_Dsu_mean = np.array(cal_mean(F_Dsu_list, nparts)) 
	F_Bsu_mean = np.array(cal_mean(F_Bsu_list, nparts)) 
	X_Dsu_mean = np.array(cal_mean(X_Dsu_list, nparts)) 
	X_Bsu_mean = np.array(cal_mean(X_Bsu_list, nparts)) 
	M_Dsu_mean = np.array(cal_mean(M_Dsu_list, nparts)) 
	M_Bsu_mean = np.array(cal_mean(M_Bsu_list, nparts))
	S_Ddu_mean = np.array(cal_mean(S_Ddu_list, nparts)) 
	C_Ddu_mean = np.array(cal_mean(C_Ddu_list, nparts)) 
	I_Ddu_mean = np.array(cal_mean(I_Ddu_list, nparts)) 
	S_Bdu_mean = np.array(cal_mean(S_Bdu_list, nparts)) 
	C_Bdu_mean = np.array(cal_mean(C_Bdu_list, nparts)) 
	I_Bdu_mean = np.array(cal_mean(I_Bdu_list, nparts)) 
	H_Ddu_mean = np.array(cal_mean(H_Ddu_list, nparts)) 
	H_Bdu_mean = np.array(cal_mean(H_Bdu_list, nparts)) 
	T_Ddu_mean = np.array(cal_mean(T_Ddu_list, nparts)) 
	F_Ddu_mean = np.array(cal_mean(F_Ddu_list, nparts)) 
	F_Bdu_mean = np.array(cal_mean(F_Bdu_list, nparts)) 
	X_Ddu_mean = np.array(cal_mean(X_Ddu_list, nparts)) 
	X_Bdu_mean = np.array(cal_mean(X_Bdu_list, nparts)) 
	M_Ddu_mean = np.array(cal_mean(M_Ddu_list, nparts)) 
	M_Bdu_mean = np.array(cal_mean(M_Bdu_list, nparts)) 

	# rural
	S_Dpr_mean = np.array(cal_mean(S_Dpr_list, nparts)) 
	C_Dpr_mean = np.array(cal_mean(C_Dpr_list, nparts)) 
	I_Dpr_mean = np.array(cal_mean(I_Dpr_list, nparts)) 
	S_Bpr_mean = np.array(cal_mean(S_Bpr_list, nparts)) 
	C_Bpr_mean = np.array(cal_mean(C_Bpr_list, nparts)) 
	I_Bpr_mean = np.array(cal_mean(I_Bpr_list, nparts)) 
	H_Dpr_mean = np.array(cal_mean(H_Dpr_list, nparts)) 
	H_Bpr_mean = np.array(cal_mean(H_Bpr_list, nparts)) 
	T_Dpr_mean = np.array(cal_mean(T_Dpr_list, nparts)) 
	F_Dpr_mean = np.array(cal_mean(F_Dpr_list, nparts)) 
	F_Bpr_mean = np.array(cal_mean(F_Bpr_list, nparts)) 
	X_Dpr_mean = np.array(cal_mean(X_Dpr_list, nparts)) 
	X_Bpr_mean = np.array(cal_mean(X_Bpr_list, nparts)) 
	M_Dpr_mean = np.array(cal_mean(M_Dpr_list, nparts)) 
	M_Bpr_mean = np.array(cal_mean(M_Bpr_list, nparts))
	S_Dsr_mean = np.array(cal_mean(S_Dsr_list, nparts)) 
	C_Dsr_mean = np.array(cal_mean(C_Dsr_list, nparts)) 
	I_Dsr_mean = np.array(cal_mean(I_Dsr_list, nparts)) 
	S_Bsr_mean = np.array(cal_mean(S_Bsr_list, nparts)) 
	C_Bsr_mean = np.array(cal_mean(C_Bsr_list, nparts)) 
	I_Bsr_mean = np.array(cal_mean(I_Bsr_list, nparts)) 
	H_Dsr_mean = np.array(cal_mean(H_Dsr_list, nparts)) 
	H_Bsr_mean = np.array(cal_mean(H_Bsr_list, nparts)) 
	T_Dsr_mean = np.array(cal_mean(T_Dsr_list, nparts)) 
	F_Dsr_mean = np.array(cal_mean(F_Dsr_list, nparts)) 
	F_Bsr_mean = np.array(cal_mean(F_Bsr_list, nparts)) 
	X_Dsr_mean = np.array(cal_mean(X_Dsr_list, nparts)) 
	X_Bsr_mean = np.array(cal_mean(X_Bsr_list, nparts)) 
	M_Dsr_mean = np.array(cal_mean(M_Dsr_list, nparts)) 
	M_Bsr_mean = np.array(cal_mean(M_Bsr_list, nparts))
	S_Ddr_mean = np.array(cal_mean(S_Ddr_list, nparts)) 
	C_Ddr_mean = np.array(cal_mean(C_Ddr_list, nparts)) 
	I_Ddr_mean = np.array(cal_mean(I_Ddr_list, nparts)) 
	S_Bdr_mean = np.array(cal_mean(S_Bdr_list, nparts)) 
	C_Bdr_mean = np.array(cal_mean(C_Bdr_list, nparts)) 
	I_Bdr_mean = np.array(cal_mean(I_Bdr_list, nparts)) 
	H_Ddr_mean = np.array(cal_mean(H_Ddr_list, nparts)) 
	H_Bdr_mean = np.array(cal_mean(H_Bdr_list, nparts)) 
	T_Ddr_mean = np.array(cal_mean(T_Ddr_list, nparts)) 
	F_Ddr_mean = np.array(cal_mean(F_Ddr_list, nparts)) 
	F_Bdr_mean = np.array(cal_mean(F_Bdr_list, nparts)) 
	X_Ddr_mean = np.array(cal_mean(X_Ddr_list, nparts)) 
	X_Bdr_mean = np.array(cal_mean(X_Bdr_list, nparts)) 
	M_Ddr_mean = np.array(cal_mean(M_Ddr_list, nparts)) 
	M_Bdr_mean = np.array(cal_mean(M_Bdr_list, nparts)) 

	return S_Dpu_mean,C_Dpu_mean,I_Dpu_mean,S_Bpu_mean,C_Bpu_mean,I_Bpu_mean,H_Dpu_mean,H_Bpu_mean,T_Dpu_mean,F_Dpu_mean,F_Bpu_mean,X_Dpu_mean,X_Bpu_mean,M_Dpu_mean,M_Bpu_mean,S_Dsu_mean,C_Dsu_mean,I_Dsu_mean,S_Bsu_mean,C_Bsu_mean,I_Bsu_mean,H_Dsu_mean,H_Bsu_mean,T_Dsu_mean,F_Dsu_mean,F_Bsu_mean,X_Dsu_mean,X_Bsu_mean,M_Dsu_mean,M_Bsu_mean,S_Ddu_mean,C_Ddu_mean,I_Ddu_mean,S_Bdu_mean,C_Bdu_mean,I_Bdu_mean,H_Ddu_mean,H_Bdu_mean,T_Ddu_mean,F_Ddu_mean,F_Bdu_mean,X_Ddu_mean,X_Bdu_mean,M_Ddu_mean,M_Bdu_mean,S_Dpr_mean,C_Dpr_mean,I_Dpr_mean,S_Bpr_mean,C_Bpr_mean,I_Bpr_mean,H_Dpr_mean,H_Bpr_mean,T_Dpr_mean,F_Dpr_mean,F_Bpr_mean,X_Dpr_mean,X_Bpr_mean,M_Dpr_mean,M_Bpr_mean,S_Dsr_mean,C_Dsr_mean,I_Dsr_mean,S_Bsr_mean,C_Bsr_mean,I_Bsr_mean,H_Dsr_mean,H_Bsr_mean,T_Dsr_mean,F_Dsr_mean,F_Bsr_mean,X_Dsr_mean,X_Bsr_mean,M_Dsr_mean,M_Bsr_mean,S_Ddr_mean,C_Ddr_mean,I_Ddr_mean,S_Bdr_mean,C_Bdr_mean,I_Bdr_mean,H_Ddr_mean,H_Bdr_mean,T_Ddr_mean,F_Ddr_mean,F_Bdr_mean,X_Ddr_mean,X_Bdr_mean,M_Ddr_mean,M_Bdr_mean

