import numpy as np
import pandas as pd

def ranperm(N):
	l=np.random.uniform(0,1,N)
	return np.array([sorted(l).index(v) for v in l])
	
def randomLHS(n, k):
	P = pd.DataFrame([ranperm(N=n) for i in range(k)]).T
	P = P + pd.DataFrame([np.random.uniform(0,1,n) for i in range(k)]).T

	return(P/n)

def triangleLHS(n, a, b, c):
	# a = min, b = max, c = mode
	p_list = randomLHS(n, 1)[0]
	
	l = []
	for p in p_list:
		if p <= (c - a) / (b - a):
			tri = a + np.sqrt(p * (c - a) * (b - a))
		else:
			tri = b - np.sqrt((1 - p) * (b - c) * (b - a))
		l.append(tri)
		
	return l

def uniformLHS(n, a, b):
	p_list = randomLHS(n, 1)[0]
	
	l = []
	for p in p_list:
		uni = a + (b - a) * p
		l.append(uni)
		
	return l

def read_parms_input_lhs(file_name, nparts):

	input_df = pd.read_excel(file_name, sheet_name=0, header=0, index_col=0)
	input_df = input_df.where((pd.notnull(input_df)), None)
	l = input_df.index.tolist()

	for var in l:
		if (input_df.loc[var]['min'] and input_df.loc[var]['max'] and input_df.loc[var]['mode']) != None:
			globals()[var] = triangleLHS(n = nparts, a = input_df.loc[var]['min'], b = input_df.loc[var]['max'], c = input_df.loc[var]['mode'])
			
		# Specials
		elif input_df.loc[var]['const'] != None:
			globals()[var] = [input_df.loc[var]['const']] * nparts
			
		elif input_df.loc[var]['min'] != None and input_df.loc[var]['max'] != None and input_df.loc[var]['mode'] == None:
			globals()[var] = uniformLHS(n = nparts, a = input_df.loc[var]['min'], b = input_df.loc[var]['max'])

		# Initial populations
		elif input_df.loc[var]['initial'] != None:
			globals()[var] = input_df.loc[var]['initial']

	# Specials
	theta_n_Dig = []
	for i in range(len(theta_j_Dig)):
		theta_n_Dig.append(1 - theta_j_Dig[i])

	aS_j_Dig = []
	for i in range(len(aC_j_Dig)):
		aS_j_Dig.append(1 - aC_j_Dig[i] - aI_j_Dig[i])

	aS_n_Dig = []
	for i in range(len(aC_n_Dig)):
		aS_n_Dig.append(1 - aC_n_Dig[i] - aI_n_Dig[i])


	# Initial environment parms
	A_Dpu = A * N_Dig0
	A_Bpu = A * N_Dig0
	A_Dpr = A * N_Dig0
	A_Bpr = A * N_Dig0
	A_Dsu = A * N_Dig0
	A_Bsu = A * N_Dig0
	A_Dsr = A * N_Dig0
	A_Bsr = A * N_Dig0
	A_Ddu = A * N_Dig0
	A_Bdu = A * N_Dig0
	A_Ddr = A * N_Dig0
	A_Bdr = A * N_Dig0

	# Initial of the Populations in each compartment
	y0_u = S_Dpu0, C_Dpu0, I_Dpu0, S_Bpu0, C_Bpu0, I_Bpu0, H_Dpu0, H_Bpu0, T_Dpu0, F_Dpu0, F_Bpu0, X_Dpu0, X_Bpu0, M_Dpu0, M_Bpu0,\
		S_Dsu0, C_Dsu0, I_Dsu0, S_Bsu0, C_Bsu0, I_Bsu0, H_Dsu0, H_Bsu0, T_Dsu0, F_Dsu0, F_Bsu0, X_Dsu0, X_Bsu0, M_Dsu0, M_Bsu0,\
		S_Ddu0, C_Ddu0, I_Ddu0, S_Bdu0, C_Bdu0, I_Bdu0, H_Ddu0, H_Bdu0, T_Ddu0, F_Ddu0, F_Bdu0, X_Ddu0, X_Bdu0, M_Ddu0, M_Bdu0

	y0_r = S_Dpr0, C_Dpr0, I_Dpr0, S_Bpr0, C_Bpr0, I_Bpr0, H_Dpr0, H_Bpr0, T_Dpr0, F_Dpr0, F_Bpr0, X_Dpr0, X_Bpr0, M_Dpr0, M_Bpr0,\
		S_Dsr0, C_Dsr0, I_Dsr0, S_Bsr0, C_Bsr0, I_Bsr0, H_Dsr0, H_Bsr0, T_Dsr0, F_Dsr0, F_Bsr0, X_Dsr0, X_Bsr0, M_Dsr0, M_Bsr0,\
		S_Ddr0, C_Ddr0, I_Ddr0, S_Bdr0, C_Bdr0, I_Bdr0, H_Ddr0, H_Bdr0, T_Ddr0, F_Ddr0, F_Bdr0, X_Ddr0, X_Bdr0, M_Ddr0, M_Bdr0
	
	return beta_B,beta_D,beta_BD,b_C,b_I,phi_C,phi_I,gamma_C,gamma_I,\
			delta_D,delta_B,m,q,k_n,k_i,x_D,x_B,y,z,rho_D,rho_B,\
			alpha_D,alpha_D_F,alpha_D_X,alpha_D_M,alpha_D_H,alpha_B_F,alpha_B,alpha_B_X,alpha_B_M,alpha_B_H,\
			alpha_H_F,alpha_H_X,alpha_H_M,alpha_H_T,alpha_F_H,alpha_X_H,alpha_M_H,\
			a_B,a_D,b_B,b_D,epsilon_F,epsilon_T,epsilon_H,epsilon_M,epsilon_X,\
			sigma_H,sigma_F_D_D,sigma_X_D_D,sigma_M_D_D,sigma_T_D_D,sigma_F_B_B,sigma_X_B_B,sigma_M_B_B,\
			d_H_D,d_H_B,d_YDF_D,d_YDX_D,d_YDM_D,d_YDT_D,d_ZBF_B,d_ZBX_B,d_ZBM_B,\
			pi_Dpu,pi_Bpu,mu_Dpu,mu_Bpu,eta_Dpu,eta_Bpu,pi_Dpr,pi_Bpr,mu_Dpr,mu_Bpr,eta_Dpr,eta_Bpr,\
			pi_Dsu,pi_Bsu,mu_Dsu,mu_Bsu,eta_Dsu,eta_Bsu,pi_Dsr,pi_Bsr,mu_Dsr,mu_Bsr,eta_Dsr,eta_Bsr,\
			pi_Ddu,pi_Bdu,mu_Ddu,mu_Bdu,eta_Ddu,eta_Bdu,pi_Ddr,pi_Bdr,mu_Ddr,mu_Bdr,eta_Ddr,eta_Bdr,\
			zeta_Dsu,zeta_Dpu,zeta_Dsr,zeta_Dpr,theta_j_Dig,theta_n_Dig,aS_j_Dig,aS_n_Dig,aC_j_Dig,aC_n_Dig,aI_j_Dig,aI_n_Dig,\
			N_Dig0, S_Dpu0, C_Dpu0, I_Dpu0, S_Bpu0, C_Bpu0, I_Bpu0, H_Dpu0, H_Bpu0, T_Dpu0, F_Dpu0, F_Bpu0, X_Dpu0, X_Bpu0, M_Dpu0, M_Bpu0,\
			S_Dsu0, C_Dsu0, I_Dsu0, S_Bsu0, C_Bsu0, I_Bsu0, H_Dsu0, H_Bsu0, T_Dsu0, F_Dsu0, F_Bsu0, X_Dsu0, X_Bsu0, M_Dsu0, M_Bsu0,\
			S_Ddu0, C_Ddu0, I_Ddu0, S_Bdu0, C_Bdu0, I_Bdu0, H_Ddu0, H_Bdu0, T_Ddu0, F_Ddu0, F_Bdu0, X_Ddu0, X_Bdu0, M_Ddu0, M_Bdu0,\
			S_Dpr0, C_Dpr0, I_Dpr0, S_Bpr0, C_Bpr0, I_Bpr0, H_Dpr0, H_Bpr0, T_Dpr0, F_Dpr0, F_Bpr0, X_Dpr0, X_Bpr0, M_Dpr0, M_Bpr0,\
			S_Dsr0, C_Dsr0, I_Dsr0, S_Bsr0, C_Bsr0, I_Bsr0, H_Dsr0, H_Bsr0, T_Dsr0, F_Dsr0, F_Bsr0, X_Dsr0, X_Bsr0, M_Dsr0, M_Bsr0,\
			S_Ddr0, C_Ddr0, I_Ddr0, S_Bdr0, C_Bdr0, I_Bdr0, H_Ddr0, H_Bdr0, T_Ddr0, F_Ddr0, F_Bdr0, X_Ddr0, X_Bdr0, M_Ddr0, M_Bdr0,\
			A,A_Dpu,A_Bpu,A_Dpr,A_Bpr,A_Dsu,A_Bsu,A_Dsr,A_Bsr,A_Ddu,A_Bdu,A_Ddr,A_Bdr,\
			y0_u, y0_r
