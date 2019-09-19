import numpy as np
import pandas as pd
# from numba import jit
# @jit(nopython=True)

def find_end(alist, nparts):
	l = []
	for i in range(nparts):
		l.append(alist[i][-1])

	return l

def record_equilibrium(output_file_name,nparts,beta_B,beta_D,beta_BD,b_C,b_I,phi_C,phi_I,gamma_C,gamma_I,\
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
				A,A_Dpu,A_Bpu,A_Dpr,A_Bpr,A_Dsu,A_Bsu,A_Dsr,A_Bsr,A_Ddu,A_Bdu,A_Ddr,A_Bdr,\
				y0_u, y0_r,\
				S_Dpu_list, C_Dpu_list, I_Dpu_list, S_Bpu_list, C_Bpu_list, I_Bpu_list, H_Dpu_list, H_Bpu_list, T_Dpu_list, F_Dpu_list, F_Bpu_list, X_Dpu_list, X_Bpu_list, M_Dpu_list, M_Bpu_list, S_Dsu_list, C_Dsu_list, I_Dsu_list, S_Bsu_list, C_Bsu_list, I_Bsu_list, H_Dsu_list, H_Bsu_list, T_Dsu_list, F_Dsu_list, F_Bsu_list, X_Dsu_list, X_Bsu_list, M_Dsu_list, M_Bsu_list, S_Ddu_list, C_Ddu_list, I_Ddu_list, S_Bdu_list, C_Bdu_list, I_Bdu_list, H_Ddu_list, H_Bdu_list, T_Ddu_list, F_Ddu_list, F_Bdu_list, X_Ddu_list, X_Bdu_list, M_Ddu_list, M_Bdu_list,\
				S_Dpr_list, C_Dpr_list, I_Dpr_list, S_Bpr_list, C_Bpr_list, I_Bpr_list, H_Dpr_list, H_Bpr_list, T_Dpr_list, F_Dpr_list, F_Bpr_list, X_Dpr_list, X_Bpr_list, M_Dpr_list, M_Bpr_list, S_Dsr_list, C_Dsr_list, I_Dsr_list, S_Bsr_list, C_Bsr_list, I_Bsr_list, H_Dsr_list, H_Bsr_list, T_Dsr_list, F_Dsr_list, F_Bsr_list, X_Dsr_list, X_Bsr_list, M_Dsr_list, M_Bsr_list, S_Ddr_list, C_Ddr_list, I_Ddr_list, S_Bdr_list, C_Bdr_list, I_Bdr_list, H_Ddr_list, H_Bdr_list, T_Ddr_list, F_Ddr_list, F_Bdr_list, X_Ddr_list, X_Bdr_list, M_Ddr_list, M_Bdr_list):
	# urban
	S_Dpu_End = find_end(S_Dpu_list, nparts) 
	C_Dpu_End = find_end(C_Dpu_list, nparts) 
	I_Dpu_End = find_end(I_Dpu_list, nparts) 
	S_Bpu_End = find_end(S_Bpu_list, nparts) 
	C_Bpu_End = find_end(C_Bpu_list, nparts) 
	I_Bpu_End = find_end(I_Bpu_list, nparts) 
	H_Dpu_End = find_end(H_Dpu_list, nparts) 
	H_Bpu_End = find_end(H_Bpu_list, nparts) 
	T_Dpu_End = find_end(T_Dpu_list, nparts) 
	F_Dpu_End = find_end(F_Dpu_list, nparts) 
	F_Bpu_End = find_end(F_Bpu_list, nparts) 
	X_Dpu_End = find_end(X_Dpu_list, nparts) 
	X_Bpu_End = find_end(X_Bpu_list, nparts) 
	M_Dpu_End = find_end(M_Dpu_list, nparts) 
	M_Bpu_End = find_end(M_Bpu_list, nparts)
	S_Dsu_End = find_end(S_Dsu_list, nparts) 
	C_Dsu_End = find_end(C_Dsu_list, nparts) 
	I_Dsu_End = find_end(I_Dsu_list, nparts) 
	S_Bsu_End = find_end(S_Bsu_list, nparts) 
	C_Bsu_End = find_end(C_Bsu_list, nparts) 
	I_Bsu_End = find_end(I_Bsu_list, nparts) 
	H_Dsu_End = find_end(H_Dsu_list, nparts) 
	H_Bsu_End = find_end(H_Bsu_list, nparts) 
	T_Dsu_End = find_end(T_Dsu_list, nparts) 
	F_Dsu_End = find_end(F_Dsu_list, nparts) 
	F_Bsu_End = find_end(F_Bsu_list, nparts) 
	X_Dsu_End = find_end(X_Dsu_list, nparts) 
	X_Bsu_End = find_end(X_Bsu_list, nparts) 
	M_Dsu_End = find_end(M_Dsu_list, nparts) 
	M_Bsu_End = find_end(M_Bsu_list, nparts)
	S_Ddu_End = find_end(S_Ddu_list, nparts) 
	C_Ddu_End = find_end(C_Ddu_list, nparts) 
	I_Ddu_End = find_end(I_Ddu_list, nparts) 
	S_Bdu_End = find_end(S_Bdu_list, nparts) 
	C_Bdu_End = find_end(C_Bdu_list, nparts) 
	I_Bdu_End = find_end(I_Bdu_list, nparts) 
	H_Ddu_End = find_end(H_Ddu_list, nparts) 
	H_Bdu_End = find_end(H_Bdu_list, nparts) 
	T_Ddu_End = find_end(T_Ddu_list, nparts) 
	F_Ddu_End = find_end(F_Ddu_list, nparts) 
	F_Bdu_End = find_end(F_Bdu_list, nparts) 
	X_Ddu_End = find_end(X_Ddu_list, nparts) 
	X_Bdu_End = find_end(X_Bdu_list, nparts) 
	M_Ddu_End = find_end(M_Ddu_list, nparts) 
	M_Bdu_End = find_end(M_Bdu_list, nparts) 

	# rural
	S_Dpr_End = find_end(S_Dpr_list, nparts) 
	C_Dpr_End = find_end(C_Dpr_list, nparts) 
	I_Dpr_End = find_end(I_Dpr_list, nparts) 
	S_Bpr_End = find_end(S_Bpr_list, nparts) 
	C_Bpr_End = find_end(C_Bpr_list, nparts) 
	I_Bpr_End = find_end(I_Bpr_list, nparts) 
	H_Dpr_End = find_end(H_Dpr_list, nparts) 
	H_Bpr_End = find_end(H_Bpr_list, nparts) 
	T_Dpr_End = find_end(T_Dpr_list, nparts) 
	F_Dpr_End = find_end(F_Dpr_list, nparts) 
	F_Bpr_End = find_end(F_Bpr_list, nparts) 
	X_Dpr_End = find_end(X_Dpr_list, nparts) 
	X_Bpr_End = find_end(X_Bpr_list, nparts) 
	M_Dpr_End = find_end(M_Dpr_list, nparts) 
	M_Bpr_End = find_end(M_Bpr_list, nparts)
	S_Dsr_End = find_end(S_Dsr_list, nparts) 
	C_Dsr_End = find_end(C_Dsr_list, nparts) 
	I_Dsr_End = find_end(I_Dsr_list, nparts) 
	S_Bsr_End = find_end(S_Bsr_list, nparts) 
	C_Bsr_End = find_end(C_Bsr_list, nparts) 
	I_Bsr_End = find_end(I_Bsr_list, nparts) 
	H_Dsr_End = find_end(H_Dsr_list, nparts) 
	H_Bsr_End = find_end(H_Bsr_list, nparts) 
	T_Dsr_End = find_end(T_Dsr_list, nparts) 
	F_Dsr_End = find_end(F_Dsr_list, nparts) 
	F_Bsr_End = find_end(F_Bsr_list, nparts) 
	X_Dsr_End = find_end(X_Dsr_list, nparts) 
	X_Bsr_End = find_end(X_Bsr_list, nparts) 
	M_Dsr_End = find_end(M_Dsr_list, nparts) 
	M_Bsr_End = find_end(M_Bsr_list, nparts)
	S_Ddr_End = find_end(S_Ddr_list, nparts) 
	C_Ddr_End = find_end(C_Ddr_list, nparts) 
	I_Ddr_End = find_end(I_Ddr_list, nparts) 
	S_Bdr_End = find_end(S_Bdr_list, nparts) 
	C_Bdr_End = find_end(C_Bdr_list, nparts) 
	I_Bdr_End = find_end(I_Bdr_list, nparts) 
	H_Ddr_End = find_end(H_Ddr_list, nparts) 
	H_Bdr_End = find_end(H_Bdr_list, nparts) 
	T_Ddr_End = find_end(T_Ddr_list, nparts) 
	F_Ddr_End = find_end(F_Ddr_list, nparts) 
	F_Bdr_End = find_end(F_Bdr_list, nparts) 
	X_Ddr_End = find_end(X_Ddr_list, nparts) 
	X_Bdr_End = find_end(X_Bdr_list, nparts) 
	M_Ddr_End = find_end(M_Ddr_list, nparts) 
	M_Bdr_End = find_end(M_Bdr_list, nparts)


	df = pd.DataFrame({'a_D': a_D, 'a_B': a_B, 'b_D': b_D, 'b_B': b_B, 'beta_B': beta_B, 'beta_D': beta_D, 'beta_BD': beta_BD, 'b_C': b_C, 'b_I': b_I,\
	'phi_C': phi_C, 'phi_I': phi_I, 'gamma_C': gamma_C, 'gamma_I': gamma_I, 'delta_D': delta_D, 'delta_B': delta_B, 'm': m, 'q': q,\
	'k_n': k_n, 'k_i': k_i, 'x_D': x_D, 'x_B': x_B, 'y': y, 'z': z, 'rho_D': rho_D, 'rho_B': rho_B, 'pi_Dpu': pi_Dpu, 'pi_Dsu': pi_Dsu,\
	'pi_Ddu': pi_Ddu, 'pi_Dpr': pi_Dpr, 'pi_Dsr': pi_Dsr, 'pi_Ddr': pi_Ddr, 'pi_Bpu': pi_Bpu, 'pi_Bpr': pi_Bpr, 'pi_Bsu': pi_Bsu,\
	'pi_Bsr': pi_Bsr, 'pi_Bdu': pi_Bdu, 'pi_Bdr': pi_Bdr, 'mu_Dpu': mu_Dpu, 'mu_Dsu': mu_Dsu, 'mu_Ddu': mu_Ddu, 'mu_Dpr': mu_Dpr,\
	'mu_Dsr': mu_Dsr, 'mu_Ddr': mu_Ddr, 'mu_Bpu': mu_Bpu, 'mu_Bpr': mu_Bpr, 'mu_Bsu': mu_Bsu, 'mu_Bsr': mu_Bsr, 'mu_Bdu': mu_Bdu,\
	'mu_Bdr': mu_Bdr, 'eta_Dpu': eta_Dpu, 'eta_Dpr': eta_Dpr, 'eta_Dsu': eta_Dsu, 'eta_Dsr': eta_Dsr, 'eta_Ddu': eta_Ddu, 'eta_Ddr': eta_Ddr,\
	'eta_Bpu': eta_Bpu, 'eta_Bpr': eta_Bpr, 'eta_Bsu': eta_Bsu, 'eta_Bsr': eta_Bsr, 'eta_Bdu': eta_Bdu, 'eta_Bdr': eta_Bdr,\
	'zeta_Dpu': zeta_Dpu, 'zeta_Dsu': zeta_Dsu, 'zeta_Dpr': zeta_Dpr, 'zeta_Dsr': zeta_Dsr, 'alpha_D': alpha_D, 'alpha_D_F': alpha_D_F,\
	'alpha_D_X': alpha_D_X, 'alpha_D_M': alpha_D_M, 'alpha_D_H': alpha_D_H, 'alpha_B': alpha_B, 'alpha_B_F': alpha_B_F,\
	'alpha_B_X': alpha_B_X, 'alpha_B_M': alpha_B_M, 'alpha_B_H': alpha_B_H, 'alpha_H_F': alpha_H_F, 'alpha_H_X': alpha_H_X, 'alpha_H_M': alpha_H_M,\
	'alpha_H_T': alpha_H_T, 'alpha_F_H': alpha_F_H, 'alpha_X_H': alpha_X_H, 'alpha_M_H': alpha_M_H, 'epsilon_F': epsilon_F,\
	'epsilon_T': epsilon_T, 'epsilon_H': epsilon_H, 'epsilon_M': epsilon_M, 'epsilon_X': epsilon_X, 'sigma_H': sigma_H, 'sigma_F_D_D': sigma_F_D_D,\
	'sigma_X_D_D': sigma_X_D_D, 'sigma_M_D_D': sigma_M_D_D, 'sigma_T_D_D': sigma_T_D_D, 'sigma_F_B_B': sigma_F_B_B,\
	'sigma_X_B_B': sigma_X_B_B, 'sigma_M_B_B': sigma_M_B_B, 'd_H_D': d_H_D, 'd_H_B': d_H_B, 'd_YDF_D': d_YDF_D, 'd_YDX_D': d_YDX_D, 'd_YDM_D': d_YDM_D,\
	'd_YDT_D': d_YDT_D, 'd_ZBF_B': d_ZBF_B, 'd_ZBX_B': d_ZBX_B, 'd_ZBM_B': d_ZBM_B, 'theta_j_Dig': theta_j_Dig, 'aC_j_Dig': aC_j_Dig,\
	'aC_n_Dig': aC_n_Dig, 'aI_j_Dig': aI_j_Dig, 'aI_n_Dig': aI_n_Dig,\
	'zeta_Dsu': zeta_Dsu, 'zeta_Dpu': zeta_Dpu, 'zeta_Dsr': zeta_Dsr, 'zeta_Dpr': zeta_Dpr, 'theta_j_Dig': theta_j_Dig,'theta_n_Dig': theta_n_Dig,'aS_j_Dig': aS_j_Dig,'aS_n_Dig': aS_n_Dig,'aC_j_Dig': aC_j_Dig,'aC_n_Dig': aC_n_Dig,'aI_j_Dig': aI_j_Dig,'aI_n_Dig': aI_n_Dig,\
	'S_Dpu_End':S_Dpu_End, 'C_Dpu_End':C_Dpu_End, 'I_Dpu_End':I_Dpu_End, 'S_Bpu_End':S_Bpu_End, 'C_Bpu_End':C_Bpu_End, 'I_Bpu_End':I_Bpu_End, 'H_Dpu_End':H_Dpu_End, 'H_Bpu_End':H_Bpu_End, 'T_Dpu_End':T_Dpu_End, 'F_Dpu_End':F_Dpu_End, 'F_Bpu_End':F_Bpu_End, 'X_Dpu_End':X_Dpu_End, 'X_Bpu_End':X_Bpu_End, 'M_Dpu_End':M_Dpu_End, 'M_Bpu_End':M_Bpu_End,\
	'S_Dsu_End':S_Dsu_End, 'C_Dsu_End':C_Dsu_End, 'I_Dsu_End':I_Dsu_End, 'S_Bsu_End':S_Bsu_End, 'C_Bsu_End':C_Bsu_End, 'I_Bsu_End':I_Bsu_End, 'H_Dsu_End':H_Dsu_End, 'H_Bsu_End':H_Bsu_End, 'T_Dsu_End':T_Dsu_End, 'F_Dsu_End':F_Dsu_End, 'F_Bsu_End':F_Bsu_End, 'X_Dsu_End':X_Dsu_End, 'X_Bsu_End':X_Bsu_End, 'M_Dsu_End':M_Dsu_End, 'M_Bsu_End':M_Bsu_End,\
	'S_Ddu_End':S_Ddu_End, 'C_Ddu_End':C_Ddu_End, 'I_Ddu_End':I_Ddu_End, 'S_Bdu_End':S_Bdu_End, 'C_Bdu_End':C_Bdu_End, 'I_Bdu_End':I_Bdu_End, 'H_Ddu_End':H_Ddu_End, 'H_Bdu_End':H_Bdu_End, 'T_Ddu_End':T_Ddu_End, 'F_Ddu_End':F_Ddu_End, 'F_Bdu_End':F_Bdu_End, 'X_Ddu_End':X_Ddu_End, 'X_Bdu_End':X_Bdu_End, 'M_Ddu_End':M_Ddu_End, 'M_Bdu_End':M_Bdu_End,\
	'S_Dpr_End':S_Dpr_End, 'C_Dpr_End':C_Dpr_End, 'I_Dpr_End':I_Dpr_End, 'S_Bpr_End':S_Bpr_End, 'C_Bpr_End':C_Bpr_End, 'I_Bpr_End':I_Bpr_End, 'H_Dpr_End':H_Dpr_End, 'H_Bpr_End':H_Bpr_End, 'T_Dpr_End':T_Dpr_End, 'F_Dpr_End':F_Dpr_End, 'F_Bpr_End':F_Bpr_End, 'X_Dpr_End':X_Dpr_End, 'X_Bpr_End':X_Bpr_End, 'M_Dpr_End':M_Dpr_End, 'M_Bpr_End':M_Bpr_End,\
	'S_Dsr_End':S_Dsr_End, 'C_Dsr_End':C_Dsr_End, 'I_Dsr_End':I_Dsr_End, 'S_Bsr_End':S_Bsr_End, 'C_Bsr_End':C_Bsr_End, 'I_Bsr_End':I_Bsr_End, 'H_Dsr_End':H_Dsr_End, 'H_Bsr_End':H_Bsr_End, 'T_Dsr_End':T_Dsr_End, 'F_Dsr_End':F_Dsr_End, 'F_Bsr_End':F_Bsr_End, 'X_Dsr_End':X_Dsr_End, 'X_Bsr_End':X_Bsr_End, 'M_Dsr_End':M_Dsr_End, 'M_Bsr_End':M_Bsr_End,\
	'S_Ddr_End':S_Ddr_End, 'C_Ddr_End':C_Ddr_End, 'I_Ddr_End':I_Ddr_End, 'S_Bdr_End':S_Bdr_End, 'C_Bdr_End':C_Bdr_End, 'I_Bdr_End':I_Bdr_End, 'H_Ddr_End':H_Ddr_End, 'H_Bdr_End':H_Bdr_End, 'T_Ddr_End':T_Ddr_End, 'F_Ddr_End':F_Ddr_End, 'F_Bdr_End':F_Bdr_End, 'X_Ddr_End':X_Ddr_End, 'X_Bdr_End':X_Bdr_End, 'M_Ddr_End':M_Ddr_End, 'M_Bdr_End':M_Bdr_End})

	# Export to spreadsheet
	df.to_excel(output_file_name, sheet_name='sheet1', index=False)

