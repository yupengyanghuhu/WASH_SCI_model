from read_parms_input_lhs import *
from simulation_odemodel_sci import *
from plot_sci import *
from record_equilibrium import *

# Partitions for LHS, can be changed
nparts = 1000
# A grid of time points (in days)
num_samples = 100
# Range (in days)
start = 0
end = 1000

# Read Parameters
beta_B,beta_D,beta_BD,b_C,b_I,phi_C,phi_I,gamma_C,gamma_I,\
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
y0_u, y0_r = read_parms_input_lhs('parms_input_final.xlsx', nparts) 


# Simulation
t, S_Dpu_list, C_Dpu_list, I_Dpu_list, S_Bpu_list, C_Bpu_list, I_Bpu_list, H_Dpu_list, H_Bpu_list, T_Dpu_list, F_Dpu_list, F_Bpu_list, X_Dpu_list, X_Bpu_list, M_Dpu_list, M_Bpu_list, S_Dsu_list, C_Dsu_list, I_Dsu_list, S_Bsu_list, C_Bsu_list, I_Bsu_list, H_Dsu_list, H_Bsu_list, T_Dsu_list, F_Dsu_list, F_Bsu_list, X_Dsu_list, X_Bsu_list, M_Dsu_list, M_Bsu_list, S_Ddu_list, C_Ddu_list, I_Ddu_list, S_Bdu_list, C_Bdu_list, I_Bdu_list, H_Ddu_list, H_Bdu_list, T_Ddu_list, F_Ddu_list, F_Bdu_list, X_Ddu_list, X_Bdu_list, M_Ddu_list, M_Bdu_list,\
S_Dpr_list, C_Dpr_list, I_Dpr_list, S_Bpr_list, C_Bpr_list, I_Bpr_list, H_Dpr_list, H_Bpr_list, T_Dpr_list, F_Dpr_list, F_Bpr_list, X_Dpr_list, X_Bpr_list, M_Dpr_list, M_Bpr_list, S_Dsr_list, C_Dsr_list, I_Dsr_list, S_Bsr_list, C_Bsr_list, I_Bsr_list, H_Dsr_list, H_Bsr_list, T_Dsr_list, F_Dsr_list, F_Bsr_list, X_Dsr_list, X_Bsr_list, M_Dsr_list, M_Bsr_list, S_Ddr_list, C_Ddr_list, I_Ddr_list, S_Bdr_list, C_Bdr_list, I_Bdr_list, H_Ddr_list, H_Bdr_list, T_Ddr_list, F_Ddr_list, F_Bdr_list, X_Ddr_list, X_Bdr_list, M_Ddr_list, M_Bdr_list =\
simulation_sci(nparts, start, end, num_samples, beta_B,beta_D,beta_BD,b_C,b_I,phi_C,phi_I,gamma_C,gamma_I,\
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
				y0_u, y0_r)


# Calculate & Plot the mean values of compartments:
S_Dpu_mean, C_Dpu_mean, I_Dpu_mean, S_Bpu_mean, C_Bpu_mean, I_Bpu_mean, H_Dpu_mean, H_Bpu_mean, T_Dpu_mean, F_Dpu_mean, F_Bpu_mean, X_Dpu_mean, X_Bpu_mean, M_Dpu_mean, M_Bpu_mean, S_Dsu_mean, C_Dsu_mean, I_Dsu_mean, S_Bsu_mean, C_Bsu_mean, I_Bsu_mean, H_Dsu_mean, H_Bsu_mean, T_Dsu_mean, F_Dsu_mean, F_Bsu_mean, X_Dsu_mean, X_Bsu_mean, M_Dsu_mean, M_Bsu_mean, S_Ddu_mean, C_Ddu_mean, I_Ddu_mean, S_Bdu_mean, C_Bdu_mean, I_Bdu_mean, H_Ddu_mean, H_Bdu_mean, T_Ddu_mean, F_Ddu_mean, F_Bdu_mean, X_Ddu_mean, X_Bdu_mean, M_Ddu_mean, M_Bdu_mean,\
S_Dpr_mean, C_Dpr_mean, I_Dpr_mean, S_Bpr_mean, C_Bpr_mean, I_Bpr_mean, H_Dpr_mean, H_Bpr_mean, T_Dpr_mean, F_Dpr_mean, F_Bpr_mean, X_Dpr_mean, X_Bpr_mean, M_Dpr_mean, M_Bpr_mean, S_Dsr_mean, C_Dsr_mean, I_Dsr_mean, S_Bsr_mean, C_Bsr_mean, I_Bsr_mean, H_Dsr_mean, H_Bsr_mean, T_Dsr_mean, F_Dsr_mean, F_Bsr_mean, X_Dsr_mean, X_Bsr_mean, M_Dsr_mean, M_Bsr_mean, S_Ddr_mean, C_Ddr_mean, I_Ddr_mean, S_Bdr_mean, C_Bdr_mean, I_Bdr_mean, H_Ddr_mean, H_Bdr_mean, T_Ddr_mean, F_Ddr_mean, F_Bdr_mean, X_Ddr_mean, X_Bdr_mean, M_Ddr_mean, M_Bdr_mean = \
cal_sci_means(t, nparts, S_Dpu_list, C_Dpu_list, I_Dpu_list, S_Bpu_list, C_Bpu_list, I_Bpu_list, H_Dpu_list, H_Bpu_list, T_Dpu_list, F_Dpu_list, F_Bpu_list, X_Dpu_list, X_Bpu_list, M_Dpu_list, M_Bpu_list, S_Dsu_list, C_Dsu_list, I_Dsu_list, S_Bsu_list, C_Bsu_list, I_Bsu_list, H_Dsu_list, H_Bsu_list, T_Dsu_list, F_Dsu_list, F_Bsu_list, X_Dsu_list, X_Bsu_list, M_Dsu_list, M_Bsu_list, S_Ddu_list, C_Ddu_list, I_Ddu_list, S_Bdu_list, C_Bdu_list, I_Bdu_list, H_Ddu_list, H_Bdu_list, T_Ddu_list, F_Ddu_list, F_Bdu_list, X_Ddu_list, X_Bdu_list, M_Ddu_list, M_Bdu_list,\
				S_Dpr_list, C_Dpr_list, I_Dpr_list, S_Bpr_list, C_Bpr_list, I_Bpr_list, H_Dpr_list, H_Bpr_list, T_Dpr_list, F_Dpr_list, F_Bpr_list, X_Dpr_list, X_Bpr_list, M_Dpr_list, M_Bpr_list, S_Dsr_list, C_Dsr_list, I_Dsr_list, S_Bsr_list, C_Bsr_list, I_Bsr_list, H_Dsr_list, H_Bsr_list, T_Dsr_list, F_Dsr_list, F_Bsr_list, X_Dsr_list, X_Bsr_list, M_Dsr_list, M_Bsr_list, S_Ddr_list, C_Ddr_list, I_Ddr_list, S_Bdr_list, C_Bdr_list, I_Bdr_list, H_Ddr_list, H_Bdr_list, T_Ddr_list, F_Ddr_list, F_Bdr_list, X_Ddr_list, X_Bdr_list, M_Ddr_list, M_Bdr_list)

## Delivery ward
plot_sci('D', S_Dpu_mean, S_Dpr_mean, S_Dsu_mean, S_Dsr_mean, S_Ddu_mean, S_Ddr_mean,\
				C_Dpu_mean, C_Dpr_mean, C_Dsu_mean, C_Dsr_mean, C_Ddu_mean, C_Ddr_mean,\
				I_Dpu_mean, I_Dpr_mean, I_Dsu_mean, I_Dsr_mean, I_Ddu_mean, I_Ddr_mean, t)

# Neonatal care
# plot_sci('B', S_Bpu_mean, S_Bpr_mean, S_Bsu_mean, S_Bsr_mean, S_Bdu_mean, S_Bdr_mean,\
# 				 C_Bpu_mean, C_Bpr_mean, C_Bsu_mean, C_Bsr_mean, C_Bdu_mean, C_Bdr_mean,\
# 				 I_Bpu_mean, I_Bpr_mean, I_Bsu_mean, I_Bsr_mean, I_Bdu_mean, I_Bdr_mean, t)


# find the final Equilibrium values & parms of compartments:
record_equilibrium('parms_output.xlsx', nparts, beta_B,beta_D,beta_BD,b_C,b_I,phi_C,phi_I,gamma_C,gamma_I,\
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
				S_Dpr_list, C_Dpr_list, I_Dpr_list, S_Bpr_list, C_Bpr_list, I_Bpr_list, H_Dpr_list, H_Bpr_list, T_Dpr_list, F_Dpr_list, F_Bpr_list, X_Dpr_list, X_Bpr_list, M_Dpr_list, M_Bpr_list, S_Dsr_list, C_Dsr_list, I_Dsr_list, S_Bsr_list, C_Bsr_list, I_Bsr_list, H_Dsr_list, H_Bsr_list, T_Dsr_list, F_Dsr_list, F_Bsr_list, X_Dsr_list, X_Bsr_list, M_Dsr_list, M_Bsr_list, S_Ddr_list, C_Ddr_list, I_Ddr_list, S_Bdr_list, C_Bdr_list, I_Bdr_list, H_Ddr_list, H_Bdr_list, T_Ddr_list, F_Ddr_list, F_Bdr_list, X_Ddr_list, X_Bdr_list, M_Ddr_list, M_Bdr_list)

