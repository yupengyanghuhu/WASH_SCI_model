import numpy as np
import pandas as pd
from scipy.integrate import odeint

# the SUM in the Lambda
def betaPoisson(x, a, b):
	#x = {d_YDLprime * Lprime_Dig, d_ZBLprime * Lprime_Dig}
	return 1-(1+x/b)**(-a*b)
	# return x


def sum_FOI_D_Lprime_combination(sigma_F_D_D, sigma_X_D_D, sigma_M_D_D, sigma_T_D_D, \
								d_YDF_D, d_YDX_D, d_YDM_D, d_YDT_D,\
								F_Dig, X_Dig, M_Dig, T_Dig, a, b):

	# lambda_F_D_YD * F_Dig = betaPoisson(d_YDF_D * F_Dig)

	return sigma_F_D_D * betaPoisson(d_YDF_D * F_Dig, a, b) +\
			sigma_X_D_D * betaPoisson(d_YDX_D * X_Dig, a, b) +\
			sigma_M_D_D * betaPoisson(d_YDM_D * M_Dig, a, b) +\
			sigma_T_D_D * betaPoisson(d_YDT_D * T_Dig, a, b)


def sum_FOI_B_Lprime_combination(sigma_F_B_B, sigma_X_B_B, sigma_M_B_B,\
								d_ZBF_B, d_ZBX_B, d_ZBM_B,\
								F_Big, X_Big, M_Big, a, b):

	# lambda_F_B_ZB * F_Big = betaPoisson(d_ZBF_B * F_Big)

	return sigma_F_B_B * betaPoisson(d_ZBF_B * F_Big, a, b) +\
			sigma_X_B_B * betaPoisson(d_ZBX_B * X_Big, a, b) +\
			sigma_M_B_B * betaPoisson(d_ZBM_B * M_Big, a, b)


# The SCI ode model differential equations.
def deriv(Y, t, beta_B,beta_D,beta_BD,b_C,b_I,phi_C,phi_I,gamma_C,gamma_I,\
				delta_D,delta_B,m,q,k_n,k_i,x_D,x_B,y,z,rho_D,rho_B,\
				alpha_D,alpha_D_F,alpha_D_X,alpha_D_M,alpha_D_H,alpha_B_F,alpha_B,alpha_B_X,alpha_B_M,alpha_B_H,\
				alpha_H_F,alpha_H_X,alpha_H_M,alpha_H_T,alpha_F_H,alpha_X_H,alpha_M_H,\
				a_B,a_D,b_B,b_D,epsilon_F,epsilon_T,epsilon_H,epsilon_M,epsilon_X,\
				sigma_H,sigma_F_D_D,sigma_X_D_D,sigma_M_D_D,sigma_T_D_D,sigma_F_B_B,sigma_X_B_B,sigma_M_B_B,\
				d_H_D,d_H_B,d_YDF_D,d_YDX_D,d_YDM_D,d_YDT_D,d_ZBF_B,d_ZBX_B,d_ZBM_B,\
				A,A_Dpg,A_Bpg,pi_Dpg,pi_Bpg,mu_Dpg,mu_Bpg,eta_Dpg,eta_Bpg,\
				A_Dsg,A_Bsg,pi_Dsg,pi_Bsg,mu_Dsg,mu_Bsg,eta_Dsg,eta_Bsg,\
				A_Ddg,A_Bdg,pi_Ddg,pi_Bdg,mu_Ddg,mu_Bdg,eta_Ddg,eta_Bdg,\
				zeta_Dsg,zeta_Dpg,theta_j_Dig,theta_n_Dig,aS_j_Dig,aS_n_Dig,aC_j_Dig,aC_n_Dig,aI_j_Dig,aI_n_Dig):
	
	# Initial S/C/I
	S_Dpg, C_Dpg, I_Dpg, S_Bpg, C_Bpg, I_Bpg, H_Dpg, H_Bpg, T_Dpg, F_Dpg, F_Bpg, X_Dpg, X_Bpg, M_Dpg, M_Bpg,\
	S_Dsg, C_Dsg, I_Dsg, S_Bsg, C_Bsg, I_Bsg, H_Dsg, H_Bsg, T_Dsg, F_Dsg, F_Bsg, X_Dsg, X_Bsg, M_Dsg, M_Bsg,\
	S_Ddg, C_Ddg, I_Ddg, S_Bdg, C_Bdg, I_Bdg, H_Ddg, H_Bdg, T_Ddg, F_Ddg, F_Bdg, X_Ddg, X_Bdg, M_Ddg, M_Bdg = Y
	
	# psi_Dig
	psi_Ddg = (pi_Ddg + mu_Ddg) * (S_Ddg + C_Ddg) + (mu_Ddg + delta_D) * I_Ddg - zeta_Dsg * (S_Dsg + C_Dsg + I_Dsg) - zeta_Dpg * (S_Dpg + C_Dpg + I_Dpg)
	psi_Dsg = (pi_Dsg + mu_Dsg + zeta_Dsg) * (S_Dsg + C_Dsg) + (mu_Dsg + delta_D + zeta_Dsg) * I_Dsg
	psi_Dpg = (pi_Dpg + mu_Dpg + zeta_Dpg) * (S_Dpg + C_Dpg) + (mu_Dpg + delta_D + zeta_Dpg) * I_Dpg

	# Lambda_Dig
	Lambda_Dpg = beta_D * (C_Dpg + I_Dpg) + beta_BD * (C_Bpg + I_Bpg) + \
					sigma_H * betaPoisson(d_H_D * H_Dpg, a_D, b_D) +\
					sum_FOI_D_Lprime_combination(sigma_F_D_D, sigma_X_D_D, sigma_M_D_D, sigma_T_D_D, \
									d_YDF_D, d_YDX_D, d_YDM_D, d_YDT_D,\
									F_Dpg, X_Dpg, M_Dpg, T_Dpg, a_D, b_D)
	
	Lambda_Dsg = beta_D * (C_Dsg + I_Dsg) + beta_BD * (C_Bsg + I_Bsg) + \
					sigma_H * betaPoisson(d_H_D * H_Dsg, a_D, b_D) +\
					sum_FOI_D_Lprime_combination(sigma_F_D_D, sigma_X_D_D, sigma_M_D_D, sigma_T_D_D, \
									d_YDF_D, d_YDX_D, d_YDM_D, d_YDT_D,\
									F_Dsg, X_Dsg, M_Dsg, T_Dsg, a_D, b_D)
	
	Lambda_Ddg = beta_D * (C_Ddg + I_Ddg) + beta_BD * (C_Bdg + I_Bdg) + \
					sigma_H * betaPoisson(d_H_D * H_Ddg, a_D, b_D) +\
					sum_FOI_D_Lprime_combination(sigma_F_D_D, sigma_X_D_D, sigma_M_D_D, sigma_T_D_D, \
									d_YDF_D, d_YDX_D, d_YDM_D, d_YDT_D,\
									F_Ddg, X_Ddg, M_Ddg, T_Ddg, a_D, b_D)
	
	# Lambda_Big
	Lambda_Bpg = beta_B * (C_Bpg + I_Bpg) + beta_BD * (C_Dpg + I_Dpg) + \
					sigma_H * betaPoisson(d_H_B * H_Bpg, a_B, b_B) +\
					sum_FOI_B_Lprime_combination(sigma_F_B_B, sigma_X_B_B, sigma_M_B_B,\
									d_ZBF_B, d_ZBX_B, d_ZBM_B,\
									F_Bpg, X_Bpg, M_Bpg, a_B, b_B)

	Lambda_Bsg = beta_B * (C_Bsg + I_Bsg) + beta_BD * (C_Dsg + I_Dsg) + \
					sigma_H * betaPoisson(d_H_B * H_Bsg, a_B, b_B) +\
					sum_FOI_B_Lprime_combination(sigma_F_B_B, sigma_X_B_B, sigma_M_B_B,\
									d_ZBF_B, d_ZBX_B, d_ZBM_B,\
									F_Bsg, X_Bsg, M_Bsg, a_B, b_B)
	
	Lambda_Bdg = beta_B * (C_Bdg + I_Bdg) + beta_BD * (C_Ddg + I_Ddg) + \
					sigma_H * betaPoisson(d_H_B * H_Bdg, a_B, b_B) +\
					sum_FOI_B_Lprime_combination(sigma_F_B_B, sigma_X_B_B, sigma_M_B_B,\
									d_ZBF_B, d_ZBX_B, d_ZBM_B,\
									F_Bdg, X_Bdg, M_Bdg, a_B, b_B)
	
	
	############# DELIVERY WARD ############
	# For S_Dpg:
	dS_Dpgdt = (gamma_C + b_C * phi_C) * C_Dpg + (gamma_I + b_I * phi_I) * I_Dpg + psi_Dpg * (theta_j_Dig * aS_j_Dig + theta_n_Dig * aS_n_Dig)\
	- (Lambda_Dpg * ((1 - y) + y * rho_D) + pi_Dpg + mu_Dpg + zeta_Dpg) * S_Dpg
	
	# For C_Dpg:
	dC_Dpgdt = (1 - x_D) * (Lambda_Dpg * ((1 - y) + y * rho_D)) * S_Dpg + psi_Dpg * (theta_j_Dig * aC_j_Dig + theta_n_Dig * aC_n_Dig)\
	- (x_D * (Lambda_Dpg * ((1 - y) + y * rho_D)) + gamma_C + b_C * phi_C + eta_Dpg + pi_Dpg + mu_Dpg + zeta_Dpg) * C_Dpg

	# For I_Dpg:
	dI_Dpgdt = x_D * (Lambda_Dpg * ((1 - y) + y * rho_D)) * (S_Dpg + C_Dpg) + eta_Dpg * C_Dpg + psi_Dpg * (theta_j_Dig * aI_j_Dig + theta_n_Dig * aI_n_Dig)\
	- (gamma_I + b_I * phi_I + mu_Dpg + delta_D + zeta_Dpg) * I_Dpg
	
	# For S_Dsg:
	dS_Dsgdt = (gamma_C + b_C * phi_C) * C_Dsg + (gamma_I + b_I * phi_I) * I_Dsg + psi_Dsg * (theta_j_Dig * aS_j_Dig + theta_n_Dig * aS_n_Dig)\
	- (Lambda_Dsg * ((1 - y) + y * rho_D) + pi_Dsg + mu_Dsg + zeta_Dsg) * S_Dsg
	
	# For C_Dsg:
	dC_Dsgdt = (1 - x_D) * (Lambda_Dsg * ((1 - y) + y * rho_D)) * S_Dsg + psi_Dsg * (theta_j_Dig * aC_j_Dig + theta_n_Dig * aC_n_Dig)\
	- (x_D * (Lambda_Dsg * ((1 - y) + y * rho_D)) + gamma_C + b_C * phi_C + eta_Dsg + pi_Dsg + mu_Dsg + zeta_Dsg) * C_Dsg

	# For I_Dsg:
	dI_Dsgdt = x_D * (Lambda_Dsg * ((1 - y) + y * rho_D)) * (S_Dsg + C_Dsg) + eta_Dsg * C_Dsg + psi_Dsg * (theta_j_Dig * aI_j_Dig + theta_n_Dig * aI_n_Dig)\
	- (gamma_I + b_I * phi_I + mu_Dsg + delta_D + zeta_Dsg) * I_Dsg
	
	# For S_Ddg:
	dS_Ddgdt = (gamma_C + b_C * phi_C) * C_Ddg + (gamma_I + b_I * phi_I) * I_Ddg + psi_Ddg * (theta_j_Dig * aS_j_Dig + theta_n_Dig * aS_n_Dig) + zeta_Dsg * S_Dsg + zeta_Dpg * S_Dpg\
	- (Lambda_Ddg * ((1 - y) + y * rho_D) + pi_Ddg + mu_Ddg) * S_Ddg
	
	# For C_Ddg:
	dC_Ddgdt = (1 - x_D) * (Lambda_Ddg * ((1 - y) + y * rho_D)) * S_Ddg + psi_Ddg * (theta_j_Dig * aC_j_Dig + theta_n_Dig * aC_n_Dig) + zeta_Dsg * C_Dsg + zeta_Dpg * C_Dpg\
	- (x_D * (Lambda_Ddg * ((1 - y) + y * rho_D)) + gamma_C + b_C * phi_C + eta_Ddg + pi_Ddg + mu_Ddg) * C_Ddg

	# For I_Ddg:
	dI_Ddgdt = x_D * (Lambda_Ddg * ((1 - y) + y * rho_D)) * (S_Ddg + C_Ddg) + eta_Ddg * C_Ddg + psi_Ddg * (theta_j_Dig * aI_j_Dig + theta_n_Dig * aI_n_Dig) + zeta_Dsg * I_Dsg + zeta_Dpg * I_Dpg\
	- (gamma_I + b_I * phi_I + mu_Ddg + delta_D) * I_Ddg
	
	
	############# NEONATAL WARD ############
	# For S_Bpg:
	dS_Bpgdt = (gamma_C + b_C * phi_C) * C_Bpg + (gamma_I + b_I * phi_I) * I_Bpg\
	+ pi_Dpg * m * q * ((1 - k_n) * S_Dpg + (1 - k_i) * (C_Dpg + I_Dpg))\
	- (Lambda_Bpg * ((1 - z) + z * rho_B) + pi_Bpg + mu_Bpg) * S_Bpg

	# For C_Bpg:
	dC_Bpgdt = (1 - x_B) * (Lambda_Bpg * ((1 - z) + z * rho_B)) * S_Bpg\
	+ pi_Dpg * m * q * (k_n * S_Dpg + k_i * (C_Dpg + I_Dpg)) * (1 - x_B)\
	- (x_B * (Lambda_Bpg * ((1 - z) + z * rho_B)) + gamma_C + b_C * phi_C + eta_Bpg + pi_Bpg + mu_Bpg) * C_Bpg

	# For I_Bpg:
	dI_Bpgdt = x_B * (Lambda_Bpg * ((1 - z) + z * rho_B)) * (S_Dpg + C_Dpg) + eta_Bpg * C_Bpg\
	+ pi_Dpg * m * (k_n * S_Dpg + k_i * (C_Dpg + I_Dpg)) * x_B\
	- (gamma_I + b_I * phi_I + mu_Bpg + delta_B) * I_Bpg
	
	# For S_Bsg:
	dS_Bsgdt = (gamma_C + b_C * phi_C) * C_Bsg + (gamma_I + b_I * phi_I) * I_Bsg\
	+ pi_Dsg * m * q * ((1 - k_n) * S_Dsg + (1 - k_i) * (C_Dsg + I_Dsg))\
	- (Lambda_Bsg * ((1 - z) + z * rho_B) + pi_Bsg + mu_Bsg) * S_Bsg

	# For C_Bsg:
	dC_Bsgdt = (1 - x_B) * (Lambda_Bsg * ((1 - z) + z * rho_B)) * S_Bsg\
	+ pi_Dsg * m * q * (k_n * S_Dsg + k_i * (C_Dsg + I_Dsg)) * (1 - x_B)\
	- (x_B * (Lambda_Bsg * ((1 - z) + z * rho_B)) + gamma_C + b_C * phi_C + eta_Bsg + pi_Bsg + mu_Bsg) * C_Bsg

	# For I_Bsg:
	dI_Bsgdt = x_B * (Lambda_Bsg * ((1 - z) + z * rho_B)) * (S_Dsg + C_Dsg) + eta_Bsg * C_Bsg\
	+ pi_Dsg * m * (k_n * S_Dsg + k_i * (C_Dsg + I_Dsg)) * x_B\
	- (gamma_I + b_I * phi_I + mu_Bsg + delta_B) * I_Bsg
	
	# For S_Bdg:
	dS_Bdgdt = (gamma_C + b_C * phi_C) * C_Bdg + (gamma_I + b_I * phi_I) * I_Bdg\
	+ pi_Ddg * m * q * ((1 - k_n) * S_Ddg + (1 - k_i) * (C_Ddg + I_Ddg))\
	- (Lambda_Bdg * ((1 - z) + z * rho_B) + pi_Bdg + mu_Bdg) * S_Bdg

	# For C_Bdg:
	dC_Bdgdt = (1 - x_B) * (Lambda_Bdg * ((1 - z) + z * rho_B)) * S_Bdg\
	+ pi_Ddg * m * q * (k_n * S_Ddg + k_i * (C_Ddg + I_Ddg)) * (1 - x_B)\
	- (x_B * (Lambda_Bdg * ((1 - z) + z * rho_B)) + gamma_C + b_C * phi_C + eta_Bdg + pi_Bdg + mu_Bdg) * C_Bdg

	# For I_Bdg:
	dI_Bdgdt = x_B * (Lambda_Bdg * ((1 - z) + z * rho_B)) * (S_Ddg + C_Ddg) + eta_Bdg * C_Bdg\
	+ pi_Ddg * m * (k_n * S_Ddg + k_i * (C_Ddg + I_Ddg)) * x_B\
	- (gamma_I + b_I * phi_I + mu_Bdg + delta_B) * I_Bdg


	############ Environmental compartments {F,X,M,H(nsep),T} ############
	
	# i = p
	## F_Dpg
	dF_Dpgdt = alpha_D_F * (C_Dpg + I_Dpg)/A_Dpg + alpha_H_F * H_Dpg/A_Dpg\
				- epsilon_F * F_Dpg

	## X_Dpg
	dX_Dpgdt = alpha_D_X * (C_Dpg + I_Dpg)/A_Dpg + alpha_H_X * H_Dpg/A_Dpg\
				- epsilon_X * X_Dpg

	## M_Dpg
	dM_Dpgdt = alpha_D_M * (C_Dpg + I_Dpg)/A_Dpg + alpha_H_M * H_Dpg/A_Dpg\
				- epsilon_M * M_Dpg

	## F_Bpg
	dF_Bpgdt = alpha_B_F * (C_Bpg + I_Bpg)/A_Bpg + alpha_H_F * H_Bpg/A_Bpg\
				- epsilon_F * F_Bpg

	## X_Bpg
	dX_Bpgdt = alpha_B_X * (C_Bpg + I_Bpg)/A_Bpg + alpha_H_X * H_Bpg/A_Bpg\
				- epsilon_X * X_Bpg

	## M_Bpg
	dM_Bpgdt = alpha_B_M * (C_Bpg + I_Bpg)/A_Bpg + alpha_H_M * H_Bpg/A_Bpg\
				- epsilon_M * M_Bpg

	## H_Dpg
	dH_Dpgdt = alpha_D_H * (C_Dpg + I_Dpg)/A_Dpg\
				+ (alpha_F_H * F_Dpg + alpha_M_H * M_Dpg + alpha_X_H * X_Dpg)/A_Dpg\
				- epsilon_H * H_Dpg

	## H_Bpg
	dH_Bpgdt = alpha_B_H * (C_Bpg + I_Bpg)/A_Dpg\
				+ (alpha_F_H * F_Bpg + alpha_M_H * M_Bpg + alpha_X_H * X_Bpg)/A_Dpg\
				- epsilon_H * H_Bpg

	## T_Dpg
	dT_Dpgdt = alpha_D * (C_Dpg + I_Dpg)/A_Dpg - epsilon_T * T_Dpg
	
	
	# i = s
	## F_Dsg
	dF_Dsgdt = alpha_D_F * (C_Dsg + I_Dsg)/A_Dsg + alpha_H_F * H_Dsg/A_Dsg\
				- epsilon_F * F_Dsg

	## X_Dsg
	dX_Dsgdt = alpha_D_X * (C_Dsg + I_Dsg)/A_Dsg + alpha_H_X * H_Dsg/A_Dsg\
				- epsilon_X * X_Dsg

	## M_Dsg
	dM_Dsgdt = alpha_D_M * (C_Dsg + I_Dsg)/A_Dsg + alpha_H_M * H_Dsg/A_Dsg\
				- epsilon_M * M_Dsg

	## F_Bsg
	dF_Bsgdt = alpha_B_F * (C_Bsg + I_Bsg)/A_Bsg + alpha_H_F * H_Bsg/A_Bsg\
				- epsilon_F * F_Bsg

	## X_Bsg
	dX_Bsgdt = alpha_B_X * (C_Bsg + I_Bsg)/A_Bsg + alpha_H_X * H_Bsg/A_Bsg\
				- epsilon_X * X_Bsg

	## M_Bsg
	dM_Bsgdt = alpha_B_M * (C_Bsg + I_Bsg)/A_Bsg + alpha_H_M * H_Bsg/A_Bsg\
				- epsilon_M * M_Bsg

	## H_Dsg
	dH_Dsgdt = alpha_D_H * (C_Dsg + I_Dsg)/A_Dsg\
				+ (alpha_F_H * F_Dsg + alpha_M_H * M_Dsg + alpha_X_H * X_Dsg)/A_Dsg\
				- epsilon_H * H_Dsg

	## H_Bsg
	dH_Bsgdt = alpha_B_H * (C_Bsg + I_Bsg)/A_Dsg\
				+ (alpha_F_H * F_Bsg + alpha_M_H * M_Bsg + alpha_X_H * X_Bsg)/A_Dsg\
				- epsilon_H * H_Bsg

	## T_Dsg
	dT_Dsgdt = alpha_D * (C_Dsg + I_Dsg)/A_Dsg - epsilon_T * T_Dsg

	
	# i = d
	## F_Ddg
	dF_Ddgdt = alpha_D_F * (C_Ddg + I_Ddg)/A_Ddg + alpha_H_F * H_Ddg/A_Ddg\
				- epsilon_F * F_Ddg

	## X_Ddg
	dX_Ddgdt = alpha_D_X * (C_Ddg + I_Ddg)/A_Ddg + alpha_H_X * H_Ddg/A_Ddg\
				- epsilon_X * X_Ddg

	## M_Ddg
	dM_Ddgdt = alpha_D_M * (C_Ddg + I_Ddg)/A_Ddg + alpha_H_M * H_Ddg/A_Ddg\
				- epsilon_M * M_Ddg

	## F_Bdg
	dF_Bdgdt = alpha_B_F * (C_Bdg + I_Bdg)/A_Bdg + alpha_H_F * H_Bdg/A_Bdg\
				- epsilon_F * F_Bdg

	## X_Bdg
	dX_Bdgdt = alpha_B_X * (C_Bdg + I_Bdg)/A_Bdg + alpha_H_X * H_Bdg/A_Bdg\
				- epsilon_X * X_Bdg

	## M_Bdg
	dM_Bdgdt = alpha_B_M * (C_Bdg + I_Bdg)/A_Bdg + alpha_H_M * H_Bdg/A_Bdg\
				- epsilon_M * M_Bdg

	## H_Ddg
	dH_Ddgdt = alpha_D_H * (C_Ddg + I_Ddg)/A_Ddg\
				+ (alpha_F_H * F_Ddg + alpha_M_H * M_Ddg + alpha_X_H * X_Ddg)/A_Ddg\
				- epsilon_H * H_Ddg

	## H_Bdg
	dH_Bdgdt = alpha_B_H * (C_Bdg + I_Bdg)/A_Ddg\
				+ (alpha_F_H * F_Bdg + alpha_M_H * M_Bdg + alpha_X_H * X_Bdg)/A_Ddg\
				- epsilon_H * H_Bdg

	## T_Ddg
	dT_Ddgdt = alpha_D * (C_Ddg + I_Ddg)/A_Ddg - epsilon_T * T_Ddg
	
	
	return dS_Dpgdt, dC_Dpgdt, dI_Dpgdt, dS_Bpgdt, dC_Bpgdt, dI_Bpgdt, dH_Dpgdt, dH_Bpgdt, dT_Dpgdt, dF_Dpgdt, dF_Bpgdt, dX_Dpgdt, dX_Bpgdt, dM_Dpgdt, dM_Bpgdt,\
			dS_Dsgdt, dC_Dsgdt, dI_Dsgdt, dS_Bsgdt, dC_Bsgdt, dI_Bsgdt, dH_Dsgdt, dH_Bsgdt, dT_Dsgdt, dF_Dsgdt, dF_Bsgdt, dX_Dsgdt, dX_Bsgdt, dM_Dsgdt, dM_Bsgdt,\
			dS_Ddgdt, dC_Ddgdt, dI_Ddgdt, dS_Bdgdt, dC_Bdgdt, dI_Bdgdt, dH_Ddgdt, dH_Bdgdt, dT_Ddgdt, dF_Ddgdt, dF_Bdgdt, dX_Ddgdt, dX_Bdgdt, dM_Ddgdt, dM_Bdgdt


def simulation_sci(nparts, start, end, num_samples, beta_B,beta_D,beta_BD,b_C,b_I,phi_C,phi_I,gamma_C,gamma_I,\
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
				y0_u, y0_r):

	S_Dpu_list = []
	C_Dpu_list = []
	I_Dpu_list = []
	S_Bpu_list = []
	C_Bpu_list = []
	I_Bpu_list = []
	H_Dpu_list = []
	H_Bpu_list = []
	T_Dpu_list = []
	F_Dpu_list = []
	F_Bpu_list = []
	X_Dpu_list = []
	X_Bpu_list = []
	M_Dpu_list = []
	M_Bpu_list = []
	S_Dsu_list = []
	C_Dsu_list = []
	I_Dsu_list = []
	S_Bsu_list = []
	C_Bsu_list = []
	I_Bsu_list = []
	H_Dsu_list = []
	H_Bsu_list = []
	T_Dsu_list = []
	F_Dsu_list = []
	F_Bsu_list = []
	X_Dsu_list = []
	X_Bsu_list = []
	M_Dsu_list = []
	M_Bsu_list = []
	S_Ddu_list = []
	C_Ddu_list = []
	I_Ddu_list = []
	S_Bdu_list = []
	C_Bdu_list = []
	I_Bdu_list = []
	H_Ddu_list = []
	H_Bdu_list = []
	T_Ddu_list = []
	F_Ddu_list = []
	F_Bdu_list = []
	X_Ddu_list = []
	X_Bdu_list = []
	M_Ddu_list = []
	M_Bdu_list = []
	S_Dpr_list = []
	C_Dpr_list = []
	I_Dpr_list = []
	S_Bpr_list = []
	C_Bpr_list = []
	I_Bpr_list = []
	H_Dpr_list = []
	H_Bpr_list = []
	T_Dpr_list = []
	F_Dpr_list = []
	F_Bpr_list = []
	X_Dpr_list = []
	X_Bpr_list = []
	M_Dpr_list = []
	M_Bpr_list = []
	S_Dsr_list = []
	C_Dsr_list = []
	I_Dsr_list = []
	S_Bsr_list = []
	C_Bsr_list = []
	I_Bsr_list = []
	H_Dsr_list = []
	H_Bsr_list = []
	T_Dsr_list = []
	F_Dsr_list = []
	F_Bsr_list = []
	X_Dsr_list = []
	X_Bsr_list = []
	M_Dsr_list = []
	M_Bsr_list = []
	S_Ddr_list = []
	C_Ddr_list = []
	I_Ddr_list = []
	S_Bdr_list = []
	C_Bdr_list = []
	I_Bdr_list = []
	H_Ddr_list = []
	H_Bdr_list = []
	T_Ddr_list = []
	F_Ddr_list = []
	F_Bdr_list = []
	X_Ddr_list = []
	X_Bdr_list = []
	M_Ddr_list = []
	M_Bdr_list = []

	t = np.linspace(start, end, num_samples) #(start, end, number of samples to generate -- turning points)

	for i in range(nparts):
		if i == int(round(nparts * 0.2)) or i == int(round(nparts * 0.4)) or i == int(round(nparts * 0.6)) or i == int(round(nparts * 0.8)):
			print('Time Process: ' + str(i) + '/' + str(nparts))
		if i == nparts - 1:
			print("COMPLETE!")
		
		# Parameters assignment
		beta_B_,beta_D_,beta_BD_,b_C_,b_I_,phi_C_,phi_I_,gamma_C_,gamma_I_,\
		delta_D_,delta_B_,m_,q_,k_n_,k_i_,x_D_,x_B_,y_,z_,rho_D_,rho_B_,\
		alpha_D_,alpha_D_F_,alpha_D_X_,alpha_D_M_,alpha_D_H_,alpha_B_F_,alpha_B_,alpha_B_X_,alpha_B_M_,alpha_B_H_,\
		alpha_H_F_,alpha_H_X_,alpha_H_M_,alpha_H_T_,alpha_F_H_,alpha_X_H_,alpha_M_H_,\
		a_B_,a_D_,b_B_,b_D_,epsilon_F_,epsilon_T_,epsilon_H_,epsilon_M_,epsilon_X_,\
		sigma_H_,sigma_F_D_D_,sigma_X_D_D_,sigma_M_D_D_,sigma_T_D_D_,sigma_F_B_B_,sigma_X_B_B_,sigma_M_B_B_,\
		d_H_D_,d_H_B_,d_YDF_D_,d_YDX_D_,d_YDM_D_,d_YDT_D_,d_ZBF_B_,d_ZBX_B_,d_ZBM_B_,\
		A_,A_Dpu_,A_Bpu_,pi_Dpu_,pi_Bpu_,mu_Dpu_,mu_Bpu_,eta_Dpu_,eta_Bpu_,\
		A_Dsu_,A_Bsu_,pi_Dsu_,pi_Bsu_,mu_Dsu_,mu_Bsu_,eta_Dsu_,eta_Bsu_,\
		A_Ddu_,A_Bdu_,pi_Ddu_,pi_Bdu_,mu_Ddu_,mu_Bdu_,eta_Ddu_,eta_Bdu_,\
		zeta_Dsu_,zeta_Dpu_,theta_j_Dig_,theta_n_Dig_,aS_j_Dig_,aS_n_Dig_,aC_j_Dig_,aC_n_Dig_,aI_j_Dig_,aI_n_Dig_\
		= beta_B[i],beta_D[i],beta_BD[i],b_C[i],b_I[i],phi_C[i],phi_I[i],gamma_C[i],gamma_I[i],\
		delta_D[i],delta_B[i],m[i],q[i],k_n[i],k_i[i],x_D[i],x_B[i],y[i],z[i],rho_D[i],rho_B[i],\
		alpha_D[i],alpha_D_F[i],alpha_D_X[i],alpha_D_M[i],alpha_D_H[i],alpha_B_F[i],alpha_B[i],alpha_B_X[i],alpha_B_M[i],alpha_B_H[i],\
		alpha_H_F[i],alpha_H_X[i],alpha_H_M[i],alpha_H_T[i],alpha_F_H[i],alpha_X_H[i],alpha_M_H[i],\
		a_B[i],a_D[i],b_B[i],b_D[i],epsilon_F[i],epsilon_T[i],epsilon_H[i],epsilon_M[i],epsilon_X[i],\
		sigma_H[i],sigma_F_D_D[i],sigma_X_D_D[i],sigma_M_D_D[i],sigma_T_D_D[i],sigma_F_B_B[i],sigma_X_B_B[i],sigma_M_B_B[i],\
		d_H_D[i],d_H_B[i],d_YDF_D[i],d_YDX_D[i],d_YDM_D[i],d_YDT_D[i],d_ZBF_B[i],d_ZBX_B[i],d_ZBM_B[i],\
		A,A_Dpu,A_Bpu,pi_Dpu[i],pi_Bpu[i],mu_Dpu[i],mu_Bpu[i],eta_Dpu[i],eta_Bpu[i],\
		A_Dsu,A_Bsu,pi_Dsu[i],pi_Bsu[i],mu_Dsu[i],mu_Bsu[i],eta_Dsu[i],eta_Bsu[i],\
		A_Ddu,A_Bdu,pi_Ddu[i],pi_Bdu[i],mu_Ddu[i],mu_Bdu[i],eta_Ddu[i],eta_Bdu[i],\
		zeta_Dsu[i],zeta_Dpu[i],theta_j_Dig[i],theta_n_Dig[i],aS_j_Dig[i],aS_n_Dig[i],aC_j_Dig[i],aC_n_Dig[i],aI_j_Dig[i],aI_n_Dig[i]

		A_Dpr_,A_Bpr_,pi_Dpr_,pi_Bpr_,mu_Dpr_,mu_Bpr_,eta_Dpr_,eta_Bpr_,\
		A_Dsr_,A_Bsr_,pi_Dsr_,pi_Bsr_,mu_Dsr_,mu_Bsr_,eta_Dsr_,eta_Bsr_,\
		A_Ddr_,A_Bdr_,pi_Ddr_,pi_Bdr_,mu_Ddr_,mu_Bdr_,eta_Ddr_,eta_Bdr_,\
		zeta_Dsr_,zeta_Dpr_\
		= A_Dpr,A_Bpr,pi_Dpr[i],pi_Bpr[i],mu_Dpr[i],mu_Bpr[i],eta_Dpr[i],eta_Bpr[i],\
		A_Dsr,A_Bsr,pi_Dsr[i],pi_Bsr[i],mu_Dsr[i],mu_Bsr[i],eta_Dsr[i],eta_Bsr[i],\
		A_Ddr,A_Bdr,pi_Ddr[i],pi_Bdr[i],mu_Ddr[i],mu_Bdr[i],eta_Ddr[i],eta_Bdr[i],\
		zeta_Dsr[i],zeta_Dpr[i]

		
		# Call deriv for integration -- URBAN
		ret_u = odeint(deriv, y0_u, t, args=(beta_B_,beta_D_,beta_BD_,b_C_,b_I_,phi_C_,phi_I_,gamma_C_,gamma_I_,\
						delta_D_,delta_B_,m_,q_,k_n_,k_i_,x_D_,x_B_,y_,z_,rho_D_,rho_B_,\
						alpha_D_,alpha_D_F_,alpha_D_X_,alpha_D_M_,alpha_D_H_,alpha_B_F_,alpha_B_,alpha_B_X_,alpha_B_M_,alpha_B_H_,\
						alpha_H_F_,alpha_H_X_,alpha_H_M_,alpha_H_T_,alpha_F_H_,alpha_X_H_,alpha_M_H_,\
						a_B_,a_D_,b_B_,b_D_,epsilon_F_,epsilon_T_,epsilon_H_,epsilon_M_,epsilon_X_,\
						sigma_H_,sigma_F_D_D_,sigma_X_D_D_,sigma_M_D_D_,sigma_T_D_D_,sigma_F_B_B_,sigma_X_B_B_,sigma_M_B_B_,\
						d_H_D_,d_H_B_,d_YDF_D_,d_YDX_D_,d_YDM_D_,d_YDT_D_,d_ZBF_B_,d_ZBX_B_,d_ZBM_B_,\
						A_,A_Dpu_,A_Bpu_,pi_Dpu_,pi_Bpu_,mu_Dpu_,mu_Bpu_,eta_Dpu_,eta_Bpu_,\
						A_Dsu_,A_Bsu_,pi_Dsu_,pi_Bsu_,mu_Dsu_,mu_Bsu_,eta_Dsu_,eta_Bsu_,\
						A_Ddu_,A_Bdu_,pi_Ddu_,pi_Bdu_,mu_Ddu_,mu_Bdu_,eta_Ddu_,eta_Bdu_,\
						zeta_Dsu_,zeta_Dpu_,theta_j_Dig_,theta_n_Dig_,aS_j_Dig_,aS_n_Dig_,aC_j_Dig_,aC_n_Dig_,aI_j_Dig_,aI_n_Dig_))

		S_Dpu, C_Dpu, I_Dpu, S_Bpu, C_Bpu, I_Bpu, H_Dpu, H_Bpu, T_Dpu, F_Dpu, F_Bpu, X_Dpu, X_Bpu, M_Dpu, M_Bpu,\
		S_Dsu, C_Dsu, I_Dsu, S_Bsu, C_Bsu, I_Bsu, H_Dsu, H_Bsu, T_Dsu, F_Dsu, F_Bsu, X_Dsu, X_Bsu, M_Dsu, M_Bsu,\
		S_Ddu, C_Ddu, I_Ddu, S_Bdu, C_Bdu, I_Bdu, H_Ddu, H_Bdu, T_Ddu, F_Ddu, F_Bdu, X_Ddu, X_Bdu, M_Ddu, M_Bdu = ret_u.T

		S_Dpu_list.append(S_Dpu) 
		C_Dpu_list.append(C_Dpu) 
		I_Dpu_list.append(I_Dpu) 
		S_Bpu_list.append(S_Bpu) 
		C_Bpu_list.append(C_Bpu) 
		I_Bpu_list.append(I_Bpu) 
		H_Dpu_list.append(H_Dpu) 
		H_Bpu_list.append(H_Bpu) 
		T_Dpu_list.append(T_Dpu) 
		F_Dpu_list.append(F_Dpu) 
		F_Bpu_list.append(F_Bpu) 
		X_Dpu_list.append(X_Dpu) 
		X_Bpu_list.append(X_Bpu) 
		M_Dpu_list.append(M_Dpu) 
		M_Bpu_list.append(M_Bpu)
		S_Dsu_list.append(S_Dsu) 
		C_Dsu_list.append(C_Dsu) 
		I_Dsu_list.append(I_Dsu) 
		S_Bsu_list.append(S_Bsu) 
		C_Bsu_list.append(C_Bsu) 
		I_Bsu_list.append(I_Bsu) 
		H_Dsu_list.append(H_Dsu) 
		H_Bsu_list.append(H_Bsu) 
		T_Dsu_list.append(T_Dsu) 
		F_Dsu_list.append(F_Dsu) 
		F_Bsu_list.append(F_Bsu) 
		X_Dsu_list.append(X_Dsu) 
		X_Bsu_list.append(X_Bsu) 
		M_Dsu_list.append(M_Dsu) 
		M_Bsu_list.append(M_Bsu)
		S_Ddu_list.append(S_Ddu) 
		C_Ddu_list.append(C_Ddu) 
		I_Ddu_list.append(I_Ddu) 
		S_Bdu_list.append(S_Bdu) 
		C_Bdu_list.append(C_Bdu) 
		I_Bdu_list.append(I_Bdu) 
		H_Ddu_list.append(H_Ddu) 
		H_Bdu_list.append(H_Bdu) 
		T_Ddu_list.append(T_Ddu) 
		F_Ddu_list.append(F_Ddu) 
		F_Bdu_list.append(F_Bdu) 
		X_Ddu_list.append(X_Ddu) 
		X_Bdu_list.append(X_Bdu) 
		M_Ddu_list.append(M_Ddu) 
		M_Bdu_list.append(M_Bdu)

		
		# Call deriv for integration -- RURAL
		ret_r = odeint(deriv, y0_r, t, args=(beta_B_,beta_D_,beta_BD_,b_C_,b_I_,phi_C_,phi_I_,gamma_C_,gamma_I_,\
						delta_D_,delta_B_,m_,q_,k_n_,k_i_,x_D_,x_B_,y_,z_,rho_D_,rho_B_,\
						alpha_D_,alpha_D_F_,alpha_D_X_,alpha_D_M_,alpha_D_H_,alpha_B_F_,alpha_B_,alpha_B_X_,alpha_B_M_,alpha_B_H_,\
						alpha_H_F_,alpha_H_X_,alpha_H_M_,alpha_H_T_,alpha_F_H_,alpha_X_H_,alpha_M_H_,\
						a_B_,a_D_,b_B_,b_D_,epsilon_F_,epsilon_T_,epsilon_H_,epsilon_M_,epsilon_X_,\
						sigma_H_,sigma_F_D_D_,sigma_X_D_D_,sigma_M_D_D_,sigma_T_D_D_,sigma_F_B_B_,sigma_X_B_B_,sigma_M_B_B_,\
						d_H_D_,d_H_B_,d_YDF_D_,d_YDX_D_,d_YDM_D_,d_YDT_D_,d_ZBF_B_,d_ZBX_B_,d_ZBM_B_,\
						A_,A_Dpr_,A_Bpr_,pi_Dpr_,pi_Bpr_,mu_Dpr_,mu_Bpr_,eta_Dpr_,eta_Bpr_,\
						A_Dsr_,A_Bsr_,pi_Dsr_,pi_Bsr_,mu_Dsr_,mu_Bsr_,eta_Dsr_,eta_Bsr_,\
						A_Ddr_,A_Bdr_,pi_Ddr_,pi_Bdr_,mu_Ddr_,mu_Bdr_,eta_Ddr_,eta_Bdr_,\
						zeta_Dsr_,zeta_Dpr_,theta_j_Dig_,theta_n_Dig_,aS_j_Dig_,aS_n_Dig_,aC_j_Dig_,aC_n_Dig_,aI_j_Dig_,aI_n_Dig_))

		S_Dpr, C_Dpr, I_Dpr, S_Bpr, C_Bpr, I_Bpr, H_Dpr, H_Bpr, T_Dpr, F_Dpr, F_Bpr, X_Dpr, X_Bpr, M_Dpr, M_Bpr,\
		S_Dsr, C_Dsr, I_Dsr, S_Bsr, C_Bsr, I_Bsr, H_Dsr, H_Bsr, T_Dsr, F_Dsr, F_Bsr, X_Dsr, X_Bsr, M_Dsr, M_Bsr,\
		S_Ddr, C_Ddr, I_Ddr, S_Bdr, C_Bdr, I_Bdr, H_Ddr, H_Bdr, T_Ddr, F_Ddr, F_Bdr, X_Ddr, X_Bdr, M_Ddr, M_Bdr = ret_r.T

		S_Dpr_list.append(S_Dpr) 
		C_Dpr_list.append(C_Dpr) 
		I_Dpr_list.append(I_Dpr) 
		S_Bpr_list.append(S_Bpr) 
		C_Bpr_list.append(C_Bpr) 
		I_Bpr_list.append(I_Bpr) 
		H_Dpr_list.append(H_Dpr) 
		H_Bpr_list.append(H_Bpr) 
		T_Dpr_list.append(T_Dpr) 
		F_Dpr_list.append(F_Dpr) 
		F_Bpr_list.append(F_Bpr) 
		X_Dpr_list.append(X_Dpr) 
		X_Bpr_list.append(X_Bpr) 
		M_Dpr_list.append(M_Dpr) 
		M_Bpr_list.append(M_Bpr)
		S_Dsr_list.append(S_Dsr) 
		C_Dsr_list.append(C_Dsr) 
		I_Dsr_list.append(I_Dsr) 
		S_Bsr_list.append(S_Bsr) 
		C_Bsr_list.append(C_Bsr) 
		I_Bsr_list.append(I_Bsr) 
		H_Dsr_list.append(H_Dsr) 
		H_Bsr_list.append(H_Bsr) 
		T_Dsr_list.append(T_Dsr) 
		F_Dsr_list.append(F_Dsr) 
		F_Bsr_list.append(F_Bsr) 
		X_Dsr_list.append(X_Dsr) 
		X_Bsr_list.append(X_Bsr) 
		M_Dsr_list.append(M_Dsr) 
		M_Bsr_list.append(M_Bsr)
		S_Ddr_list.append(S_Ddr) 
		C_Ddr_list.append(C_Ddr) 
		I_Ddr_list.append(I_Ddr) 
		S_Bdr_list.append(S_Bdr) 
		C_Bdr_list.append(C_Bdr) 
		I_Bdr_list.append(I_Bdr) 
		H_Ddr_list.append(H_Ddr) 
		H_Bdr_list.append(H_Bdr) 
		T_Ddr_list.append(T_Ddr) 
		F_Ddr_list.append(F_Ddr) 
		F_Bdr_list.append(F_Bdr) 
		X_Ddr_list.append(X_Ddr) 
		X_Bdr_list.append(X_Bdr) 
		M_Ddr_list.append(M_Ddr) 
		M_Bdr_list.append(M_Bdr)

	return t, S_Dpu_list, C_Dpu_list, I_Dpu_list, S_Bpu_list, C_Bpu_list, I_Bpu_list, H_Dpu_list, H_Bpu_list, T_Dpu_list, F_Dpu_list, F_Bpu_list, X_Dpu_list, X_Bpu_list, M_Dpu_list, M_Bpu_list, S_Dsu_list, C_Dsu_list, I_Dsu_list, S_Bsu_list, C_Bsu_list, I_Bsu_list, H_Dsu_list, H_Bsu_list, T_Dsu_list, F_Dsu_list, F_Bsu_list, X_Dsu_list, X_Bsu_list, M_Dsu_list, M_Bsu_list, S_Ddu_list, C_Ddu_list, I_Ddu_list, S_Bdu_list, C_Bdu_list, I_Bdu_list, H_Ddu_list, H_Bdu_list, T_Ddu_list, F_Ddu_list, F_Bdu_list, X_Ddu_list, X_Bdu_list, M_Ddu_list, M_Bdu_list,\
			S_Dpr_list, C_Dpr_list, I_Dpr_list, S_Bpr_list, C_Bpr_list, I_Bpr_list, H_Dpr_list, H_Bpr_list, T_Dpr_list, F_Dpr_list, F_Bpr_list, X_Dpr_list, X_Bpr_list, M_Dpr_list, M_Bpr_list, S_Dsr_list, C_Dsr_list, I_Dsr_list, S_Bsr_list, C_Bsr_list, I_Bsr_list, H_Dsr_list, H_Bsr_list, T_Dsr_list, F_Dsr_list, F_Bsr_list, X_Dsr_list, X_Bsr_list, M_Dsr_list, M_Bsr_list, S_Ddr_list, C_Ddr_list, I_Ddr_list, S_Bdr_list, C_Bdr_list, I_Bdr_list, H_Ddr_list, H_Bdr_list, T_Ddr_list, F_Ddr_list, F_Bdr_list, X_Ddr_list, X_Bdr_list, M_Ddr_list, M_Bdr_list 
