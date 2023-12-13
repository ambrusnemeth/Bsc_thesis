import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

def rho_APC(t):
    if t < 720:
        return 1
    elif t >= 1440:
        return 0
    else:
        return 0.5

def Wnt(t):
    if t<720:
        return 0
    else:
        return 2

# kezdeti értékek beállítása
y0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 30.0, 25.0, 9.0, 30.0, 27.0,
      0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 30.0, 25.0, 9.0, 30.0, 27.0,
      25, 10, 10, 10, 10, 10, 10, 35,
      25, 10, 10, 10, 10, 10, 10, 35]
#y0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 30.0, 25.0, 9.0, 30.0, 27.0, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 0.51, 30.0, 25.0, 9.0, 30.0, 27.0]
#y0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 30.0, 25.0, 9.0, 30.0, 27.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 30.0, 25.0, 9.0, 30.0, 27.0]
# N = 0.5, F = 0.5, I_1 = 0.5, H_1 = 0.5, P = 0.5, D = 0.5, H2 = 0.5, G = 30.0, C = 25.0, B = 9.0, I_2 = 30.0, A = 27.0

t1 = np.linspace(0, 2160, num=2000)
#t1 = np.linspace(0, 4320, num=2000)



# paraméterek megadása

Dbar = 0                # the mean Delta level expressed by neighbouring cells
m_1 = 3                 # Hill coefficient for the Hill function for Delta
kappa_1 = 0.6           # dissociation constant for the Hill function for Delta
mu_N = 0.017            # decay rate of Notch

mu_F = 0.00385          # decay rate of NICD
alpha_frag = 0.8        # the rate of Notch fragmentation
alpha_1 = 6.8           # the rate at which it binds with β-catenin

mu_I1 = 0.04            # decay rate of I_1

mu_H1 = 0.065           # decay rate of Hes1
sigma_K = 1             # inhibition constant in the Dsh function
# Psi_W=sigma_K/(sigma_k+W)
theta_2 = 0.75          # relative contribution of I1 to Hes1 transcription upregulation
m_2 = 3                 # Hill coefficient for the Hill function for I_2
kappa_2 = 0.9           # dissociation constant for the Hill function for I_2
theta_7 = 0.25             # relative contribution of B to Hes1 transcription upregulation
m_7 = 3                 # Hill coefficient for the Hill function for I_2
kappa_7 = 10            # dissociation constant for the Hill function for B
sigma_2 = 3.5           # inhibition constant for the Hes1 autorepression
n_2 = 3                 # exponent for the Hes1 autorepression

mu_P = 0.035            # decay rate of Ngn3
sigma_3 = 1.21          # inhibition constant for Ngn3 inhibition by Hes1
n_3 = 3                 # exponent for Ngn3 inhibition by Hes1

mu_D = 0.049            # decay rate of the Delta ligand
theta_4 = 0.4           # maximal rate of Delta synthesis
m_4 = 3                 # Hill coefficient for the Hill function for Delta synthesis
kappa_4 = 14.7          # dissociation constant for the Hill function for Delta synthesis

mu_H2 = 0.0311          # decay rate of Hath1
xi_5 = 0.9              # maximal rate Hes1 influence upon the production of Hath1
sigma_5 = 1.7           # inhibition constant for the Hes1 influence upon the production of Hath1
n_5 = 3                 # exponent for the Hes1 influence upon the production of Hath1

mu_G = 9.36e-4          # decay rate of GSK3β
W = 0                   # Wnt stimulus
alpha_2 = 0.0174        # basal rate of GSK3β production
alpha_5 = 0.1044        # rate of formation of the destruction complex
# Psi_WA = 1.4*A**2/(1+(1+W)**4)

mu_C = 1.061            # decay rate of the destruction complex
alpha_3 = 1.465e-4      # the rate at which the destruction complex binds with β-catenin

mu_B = 0.00636          # decay rate of β-catenin
alpha_4 = 0.472         # production rate of β-catenin

mu_I2 = 4.204           # decay rate of I_2

mu_A = 6.23e-4          # decay rate of Axin
theta_6 = 1.64e-4       # maximal rate of Axin synthesis
m_6 = 3                 # Hill coefficient for Axin synthesis
kappa_6 = 0.026         # dissociation constant for axin synthesis

# és ami kimaradt
xi_2 = 0.5              # maximal rate of Hes1 transcription
theta_1 = 0.06          # maximum rate of Notch synthesis
xi_3 = 0.9              # maximal rate of Ngn3 inhibition by Hes1

K_d = 0.4989        # PER:CLOCK dissociation constant, 1 < K_d < 2
K_m = 7.88          # Michaelis constant for the degradation of nuclear PER [c]
K_A = 20            # Dissociation constant of the BMAL1:Ebox complex
A_T = 33.25         # Concentration of total CLOCK:BMAL1 (bound and unbound) in the nucleus [c], "AT=30 is not unreasonable" (cikk: 20)
alpha = 103.7575/60    # ! Rate of Per mRNA transcription [c/t], "α<45" (cikk: 50)
beta = 0.92249/60      # ! Rate of Per phosphorylation in the cytoplasm [1/t]
delta = 0.56249/60     # ! Rate of Per degradation [1/t]
mu = 8.6561/60         # ! Rate of Per translation [1/t]
gamma = 0.98999/60     # ! Rate of Per transport from the cytoplasm to the nucleus [1/t]
delta_N = 6.0074/60    # ! Michaelis-Menten limiting rate of degradation [c/t] Szerintem kisebb, mint 1...

params = [Dbar, m_1, kappa_1, mu_N, mu_F, alpha_frag, alpha_1, mu_I1, mu_H1, sigma_K, theta_2, m_2, kappa_2, theta_7,
          m_7, kappa_7, sigma_2, n_2, mu_P, sigma_3, n_3, mu_D, theta_4, m_4, kappa_4, mu_H2, xi_5, sigma_5, n_5, mu_G,
          W, alpha_2, alpha_5, mu_C, alpha_3, mu_B, alpha_4, mu_I2, mu_A, theta_6, m_6, kappa_6, xi_2, theta_1, xi_3,
          alpha, beta, delta, mu, A_T, gamma, K_d, K_m, delta_N]

# ode-solver
def sim(variables, t, params):
    N = variables[0]            # Membrane-bound Notch receptor
    F = variables[1]            # Notch Intracellular Domain (NICD)
    I_1 = variables[2]          # Intermediate 1 (NICD/β-catenin)
    H_1 = variables[3]          # Hes1
    P = variables[4]            # Ngn3
    D = variables[5]            # Delta ligand
    H_2 = variables[6]          # Hath1
    G = variables[7]            # GSK3β
    C = variables[8]            # Destruction complex
    B = variables[9]            # Active β-catenin
    I_2 = variables[10]         # Intermediate 2 (GSK3β/β-catenin)
    A = variables[11]           # Axin

    N2 = variables[12]  # Membrane-bound Notch receptor
    F2 = variables[13]  # Notch Intracellular Domain (NICD)
    I_12 = variables[14]  # Intermediate 1 (NICD/β-catenin)
    H_12 = variables[15]  # Hes1
    P2 = variables[16]  # Ngn3
    D2 = variables[17]  # Delta ligand
    H_22 = variables[18]  # Hath1
    G2 = variables[19]  # GSK3β
    C2 = variables[20]  # Destruction complex
    B2 = variables[21]  # Active β-catenin
    I_22 = variables[22]  # Intermediate 2 (GSK3β/β-catenin)
    A2 = variables[23]  # Axin

    M1 = variables[24]  # PER mRNA
    P_11 = variables[25]
    P_21 = variables[26]
    P_31 = variables[27]
    P_41 = variables[28]
    P_51 = variables[29]
    P_61 = variables[30]
    P_N1 = variables[31]  # Nuclear PER protein

    M2 = variables[32]  # PER mRNA
    P_12 = variables[33]
    P_22 = variables[34]
    P_32 = variables[35]
    P_42 = variables[36]
    P_52 = variables[37]
    P_62 = variables[38]
    P_N2 = variables[39]  # Nuclear PER protein

    Dbar = params[0]
    m_1 = params[1]
    kappa_1 = params[2]
    mu_N = params[3]
    mu_F = params[4]
    alpha_frag = params[5]
    alpha_1 = params[6]
    mu_I1 = params[7]
    mu_H1 = params[8]
    sigma_K = params[9]
    theta_2 = params[10]
    m_2 = params[11]
    kappa_2 = params[12]
    theta_7 = params[13]
    m_7 = params[14]
    kappa_7 = params[15]
    sigma_2 = params[16]
    n_2 = params[17]
    mu_P = params[18]
    sigma_3 = params[19]
    n_3 = params[20]
    mu_D = params[21]
    theta_4 = params[22]
    m_4 = params[23]
    kappa_4 = params[24]
    mu_H2 = params[25]
    xi_5 = params[26]
    sigma_5 = params[27]
    n_5 = params[28]
    mu_G = params[29]
    W = params[30]
    alpha_2 = params[31]
    alpha_5 = params[32]
    mu_C = params[33]
    alpha_3 = params[34]
    mu_B = params[35]
    alpha_4 = params[36]
    mu_I2 = params[37]
    mu_A = params[38]
    theta_6 = params[39]
    m_6 = params[40]
    kappa_6 = params[41]
    xi_2 = params[42]
    theta_1= params[43]
    xi_3 = params[44]

    alpha = params[45]
    beta = params[46]
    delta = params[47]
    mu = params[48]
    A_T = params[49]
    gamma = params[50]
    K_d = params[51]
    K_m = params[52]
    delta_N = params[53]

    dNdt = -mu_N*N + theta_1*D2**m_1/(kappa_1**m_1 + D2**m_1)
    dFdt = -mu_F*F + alpha_frag*mu_N*N - alpha_1*B*F
    dI_1dt = -mu_I1*I_1 + alpha_1*B*F
    dH_1dt = -mu_H1*H_1 + sigma_K/(sigma_K + W)*(theta_2*I_1**m_2/(kappa_2**m_2 + I_1**m_2) + theta_7*B**m_7/(kappa_7**m_7 + B**m_7))*xi_2*sigma_2**n_2/(sigma_2**n_2 + H_1**n_2)*0.9*0.5*(A_T - P_N1 - K_d + ((A_T - P_N1 - K_d)**2 + 4*K_d*A_T)**0.5)
    dPdt = -mu_P*P + xi_3*sigma_3**n_3/(sigma_3**n_3 + H_1**n_3)
    dDdt = -mu_D*D + theta_4*P**m_4/(kappa_4**m_4 + P**m_4)
    dH_2dt = -mu_H2*H_2 + xi_5*sigma_5**n_5/(sigma_5**n_5 + H_1**n_5)
    dGdt = -mu_G*G + (1 + W)*alpha_2 + mu_C*C - alpha_5*1.4*(A**2)/(1 + (1 + W)**4)*G # alpha_5*1.4*(A**2)/(1 + (1 + W)**4)*G
    dCdt = -mu_C*C + alpha_5*1.4*(A**2)/(1 + (1 + W)**4)*G + mu_I2*I_2 - alpha_3*B*C
    dBdt = -mu_B*B + (1 + W)*alpha_4 - alpha_1*B*F - alpha_3*B*C
    dI_2dt = -mu_I2*I_2 + alpha_3*B*C
    dAdt = -mu_A*A + theta_6*B**m_6/(kappa_6**m_6 + B**m_6)     # 11

    dN2dt = -mu_N*N2 + theta_1*D**m_1/(kappa_1**m_1 + D**m_1)
    dF2dt = -mu_F*F2 + alpha_frag*mu_N*N2 - alpha_1*B2*F2
    dI_12dt = -mu_I1*I_12 + alpha_1*B2*F2
    dH_12dt = -mu_H1*H_12 + sigma_K/(sigma_K + W)*(theta_2*I_12**m_2/(kappa_2**m_2 + I_12**m_2) + theta_7*B2**m_7/(kappa_7**m_7 + B2**m_7))*xi_2*sigma_2**n_2/(sigma_2**n_2 + H_12**n_2)*0.9*0.5*(A_T - P_N2 - K_d + ((A_T - P_N2 - K_d)**2 + 4*K_d*A_T)**0.5)
    dP2dt = -mu_P*P2 + xi_3*sigma_3**n_3/(sigma_3**n_3 + H_12**n_3)
    dD2dt = -mu_D*D2 + theta_4*P2**m_4/(kappa_4**m_4 + P2**m_4)
    dH_22dt = -mu_H2*H_22 + xi_5*sigma_5**n_5/(sigma_5**n_5 + H_12**n_5)
    dG2dt = -mu_G*G2 + (1 + W)*alpha_2 + mu_C*C2 - alpha_5*1.4*(A2**2)/(1 + (1 + W)**4)*G2
    dC2dt = -mu_C*C2 + alpha_5*1.4*(A2**2)/(1 + (1 + W)**4)*G2 + mu_I2*I_22 - alpha_3*B2*C2
    dB2dt = -mu_B*B2 + (1 + W)*alpha_4 - alpha_1*B2*F2 - alpha_3*B2*C2
    dI_22dt = -mu_I2*I_22 + alpha_3*B2*C2
    dA2dt = -mu_A*A2 + theta_6*B2**m_6/(kappa_6**m_6 + B2**m_6)

    dM1dt = alpha*0.5*(A_T - P_N1 - K_d + ((A_T - P_N1 - K_d)**2 + 4*K_d*A_T)**0.5)/(A_T) - delta*M1
    # a segédábrán P volt, de az szerintem a P_N-t akarta jelenteni...
    dP_11dt = mu * M1 - (beta + delta) * P_11
    dP_21dt = beta * P_11 - (beta + delta) * P_21
    dP_31dt = beta * P_21 - (beta + delta) * P_31
    dP_41dt = beta * P_31 - (beta + delta) * P_41
    dP_51dt = beta * P_41 - (beta + delta) * P_51
    dP_61dt = beta * P_51 - (gamma + delta) * P_61
    dP_N1dt = gamma*P_61 - delta_N * P_N1/(K_m + P_N1)

    dM2dt = alpha*0.5*(A_T - P_N2 - K_d + ((A_T - P_N2 - K_d)**2 + 4*K_d*A_T)**0.5)/(A_T) - delta * M2
    # a segédábrán P volt, de az szerintem a P_N-t akarta jelenteni...
    dP_12dt = mu * M2 - (beta + delta) * P_12
    dP_22dt = beta * P_12 - (beta + delta) * P_22
    dP_32dt = beta * P_22 - (beta + delta) * P_32
    dP_42dt = beta * P_32 - (beta + delta) * P_42
    dP_52dt = beta * P_42 - (beta + delta) * P_52
    dP_62dt = beta * P_52 - (gamma + delta) * P_62
    dP_N2dt = gamma * P_62 - delta_N * P_N2/(K_m + P_N2)

    return ([dNdt, dFdt, dI_1dt, dH_1dt, dPdt, dDdt, dH_2dt, dGdt, dCdt, dBdt, dI_2dt, dAdt,
             dN2dt, dF2dt, dI_12dt, dH_12dt, dP2dt, dD2dt, dH_22dt, dG2dt, dC2dt, dB2dt, dI_22dt, dA2dt,
             dM1dt, dP_11dt, dP_21dt, dP_31dt, dP_41dt, dP_51dt, dP_61dt, dP_N1dt,
             dM2dt, dP_12dt, dP_22dt, dP_32dt, dP_42dt, dP_52dt, dP_62dt, dP_N2dt])

#y = odeint(sim,y0,t,args=(params,))

tt1 = t1/60

fig = plt.figure()

fig.text(0.04, 0.5, 'Concentration (nM)', va='center', rotation='vertical')
fig.text(0.5, 0.04, 'Time (h)', ha='center')

y = odeint(sim,y0,t1,args=(params,))
fig.add_subplot(4,3,1)
plt.title("Healthy, two coupled cells")
plt.ylabel("β-cat")
plt.plot(tt1, y[:,9], 'k')
plt.plot(tt1, y[:,21], 'k:')
plt.ylim([0, 200])

fig.add_subplot(4,3,4)
plt.ylabel("Hes1")
plt.plot(tt1, y[:,3], 'k')
plt.plot(tt1, y[:,15], 'k:')
plt.ylim([0, 5])

params[30] = 1
y = odeint(sim,y0,t1,args=(params,))
fig.add_subplot(4,3,7)
plt.ylabel("β-cat")
plt.plot(tt1, y[:,9], 'r')
plt.plot(tt1, y[:,21], 'r:')
plt.ylim([0, 200])

fig.add_subplot(4,3,10)
plt.ylabel("Hes1")
plt.plot(tt1, y[:,3], 'r')
plt.plot(tt1, y[:,15], 'r:')
plt.ylim([0, 5])

def sim(variables, t, params):
    N = variables[0]  # Membrane-bound Notch receptor
    F = variables[1]  # Notch Intracellular Domain (NICD)
    I_1 = variables[2]  # Intermediate 1 (NICD/β-catenin)
    H_1 = variables[3]  # Hes1
    P = variables[4]  # Ngn3
    D = variables[5]  # Delta ligand
    H_2 = variables[6]  # Hath1
    G = variables[7]  # GSK3β
    C = variables[8]  # Destruction complex
    B = variables[9]  # Active β-catenin
    I_2 = variables[10]  # Intermediate 2 (GSK3β/β-catenin)
    A = variables[11]  # Axin

    N2 = variables[12]  # Membrane-bound Notch receptor
    F2 = variables[13]  # Notch Intracellular Domain (NICD)
    I_12 = variables[14]  # Intermediate 1 (NICD/β-catenin)
    H_12 = variables[15]  # Hes1
    P2 = variables[16]  # Ngn3
    D2 = variables[17]  # Delta ligand
    H_22 = variables[18]  # Hath1
    G2 = variables[19]  # GSK3β
    C2 = variables[20]  # Destruction complex
    B2 = variables[21]  # Active β-catenin
    I_22 = variables[22]  # Intermediate 2 (GSK3β/β-catenin)
    A2 = variables[23]  # Axin

    M1 = variables[24]  # PER mRNA
    P_11 = variables[25]
    P_21 = variables[26]
    P_31 = variables[27]
    P_41 = variables[28]
    P_51 = variables[29]
    P_61 = variables[30]
    P_N1 = variables[31]  # Nuclear PER protein

    M2 = variables[32]  # PER mRNA
    P_12 = variables[33]
    P_22 = variables[34]
    P_32 = variables[35]
    P_42 = variables[36]
    P_52 = variables[37]
    P_62 = variables[38]
    P_N2 = variables[39]  # Nuclear PER protein

    Dbar = params[0]
    m_1 = params[1]
    kappa_1 = params[2]
    mu_N = params[3]
    mu_F = params[4]
    alpha_frag = params[5]
    alpha_1 = params[6]
    mu_I1 = params[7]
    mu_H1 = params[8]
    sigma_K = params[9]
    theta_2 = params[10]
    m_2 = params[11]
    kappa_2 = params[12]
    theta_7 = params[13]
    m_7 = params[14]
    kappa_7 = params[15]
    sigma_2 = params[16]
    n_2 = params[17]
    mu_P = params[18]
    sigma_3 = params[19]
    n_3 = params[20]
    mu_D = params[21]
    theta_4 = params[22]
    m_4 = params[23]
    kappa_4 = params[24]
    mu_H2 = params[25]
    xi_5 = params[26]
    sigma_5 = params[27]
    n_5 = params[28]
    mu_G = params[29]
    W = params[30]
    alpha_2 = params[31]
    alpha_5 = params[32]
    mu_C = params[33]
    alpha_3 = params[34]
    mu_B = params[35]
    alpha_4 = params[36]
    mu_I2 = params[37]
    mu_A = params[38]
    theta_6 = params[39]
    m_6 = params[40]
    kappa_6 = params[41]
    xi_2 = params[42]
    theta_1 = params[43]
    xi_3 = params[44]

    alpha = params[45]
    beta = params[46]
    delta = params[47]
    mu = params[48]
    A_T = params[49]
    gamma = params[50]
    K_d = params[51]
    K_m = params[52]
    delta_N = params[53]

    dNdt = -mu_N*N + theta_1*D2**m_1/(kappa_1**m_1 + D2**m_1)
    dFdt = -mu_F*F + alpha_frag*mu_N*N - alpha_1*B*F
    dI_1dt = -mu_I1*I_1 + alpha_1*B*F
    dH_1dt = -mu_H1*H_1 + sigma_K/(sigma_K + W)*(theta_2*I_1**m_2/(kappa_2**m_2 + I_1**m_2) + theta_7*B**m_7/(kappa_7**m_7 + B**m_7))*xi_2*sigma_2**n_2/(sigma_2**n_2 + H_1**n_2)*0.9*0.5*(A_T - P_N1 - K_d + ((A_T - P_N1 - K_d)**2 + 4*K_d*A_T)**0.5)
    dPdt = -mu_P*P + xi_3*sigma_3**n_3/(sigma_3**n_3 + H_1**n_3)
    dDdt = -mu_D*D + theta_4*P**m_4/(kappa_4**m_4 + P**m_4)
    dH_2dt = -mu_H2*H_2 + xi_5*sigma_5**n_5/(sigma_5**n_5 + H_1**n_5)
    dGdt = -mu_G*G + (1 + W)*alpha_2 + mu_C*C - rho_APC(t)*alpha_5*1.4*(A**2)/(1 + (1 + W)**4)*G
    dCdt = -mu_C*C + rho_APC(t)*alpha_5*1.4*(A**2)/(1 + (1 + W)**4)*G + mu_I2*I_2 - alpha_3*B*C
    dBdt = -mu_B*B + (1 + W)*alpha_4 - alpha_1*B*F - alpha_3*B*C
    dI_2dt = -mu_I2*I_2 + alpha_3*B*C
    dAdt = -mu_A*A + theta_6*B**m_6/(kappa_6**m_6 + B**m_6)  # 11

    dN2dt = -mu_N*N2 + theta_1*D**m_1/(kappa_1**m_1 + D**m_1)
    dF2dt = -mu_F*F2 + alpha_frag*mu_N*N2 - alpha_1*B2*F2
    dI_12dt = -mu_I1*I_12 + alpha_1*B2*F2
    dH_12dt = -mu_H1*H_12 + sigma_K/(sigma_K + W)*(theta_2*I_12**m_2/(kappa_2**m_2 + I_12**m_2) + theta_7*B2**m_7/(kappa_7**m_7 + B2**m_7))*xi_2*sigma_2**n_2/(sigma_2**n_2 + H_12**n_2)*0.9*0.5*(A_T - P_N2 - K_d + ((A_T - P_N2 - K_d)**2 + 4*K_d*A_T)**0.5)
    dP2dt = -mu_P*P2 + xi_3*sigma_3**n_3/(sigma_3**n_3 + H_12**n_3)
    dD2dt = -mu_D*D2 + theta_4*P2**m_4/(kappa_4**m_4 + P2**m_4)
    dH_22dt = -mu_H2*H_22 + xi_5*sigma_5**n_5/(sigma_5**n_5 + H_12**n_5)
    dG2dt = -mu_G*G2 + (1 + W)*alpha_2 + mu_C*C2 - alpha_5*1.4*(A2**2)/(1 + (1 + W)**4)*G2
    dC2dt = -mu_C*C2 + alpha_5*1.4*(A2**2)/(1 + (1 + W)**4)*G2 + mu_I2*I_22 - alpha_3*B2*C2
    dB2dt = -mu_B*B2 + (1 + W)*alpha_4 - alpha_1*B2*F2 - alpha_3*B2*C2
    dI_22dt = -mu_I2*I_22 + alpha_3*B2*C2
    dA2dt = -mu_A*A2 + theta_6*B2**m_6/(kappa_6**m_6 + B2**m_6)

    dM1dt = alpha * 0.5 * (A_T - P_N1 - K_d + ((A_T - P_N1 - K_d) ** 2 + 4 * K_d * A_T) ** 0.5) / (A_T) - delta * M1
    # a segédábrán P volt, de az szerintem a P_N-t akarta jelenteni...
    dP_11dt = mu * M1 - (beta + delta) * P_11
    dP_21dt = beta * P_11 - (beta + delta) * P_21
    dP_31dt = beta * P_21 - (beta + delta) * P_31
    dP_41dt = beta * P_31 - (beta + delta) * P_41
    dP_51dt = beta * P_41 - (beta + delta) * P_51
    dP_61dt = beta * P_51 - (gamma + delta) * P_61
    dP_N1dt = gamma * P_61 - delta_N * P_N1 / (K_m + P_N1)

    dM2dt = alpha * 0.5 * (A_T - P_N2 - K_d + ((A_T - P_N2 - K_d) ** 2 + 4 * K_d * A_T) ** 0.5) / (A_T) - delta * M2
    # a segédábrán P volt, de az szerintem a P_N-t akarta jelenteni...
    dP_12dt = mu * M2 - (beta + delta) * P_12
    dP_22dt = beta * P_12 - (beta + delta) * P_22
    dP_32dt = beta * P_22 - (beta + delta) * P_32
    dP_42dt = beta * P_32 - (beta + delta) * P_42
    dP_52dt = beta * P_42 - (beta + delta) * P_52
    dP_62dt = beta * P_52 - (gamma + delta) * P_62
    dP_N2dt = gamma * P_62 - delta_N * P_N2 / (K_m + P_N2)

    return ([dNdt, dFdt, dI_1dt, dH_1dt, dPdt, dDdt, dH_2dt, dGdt, dCdt, dBdt, dI_2dt, dAdt,
             dN2dt, dF2dt, dI_12dt, dH_12dt, dP2dt, dD2dt, dH_22dt, dG2dt, dC2dt, dB2dt, dI_22dt, dA2dt,
             dM1dt, dP_11dt, dP_21dt, dP_31dt, dP_41dt, dP_51dt, dP_61dt, dP_N1dt,
             dM2dt, dP_12dt, dP_22dt, dP_32dt, dP_42dt, dP_52dt, dP_62dt, dP_N2dt])

'''
f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12) = plt.subplots(12, sharex=True, sharey=False)

line1, = ax1.plot(t , y[:,0], color="b",label="N")
line2, = ax2.plot(t , y[:,1], color="r",label="Na")
line3, = ax3.plot(t , y[:,2], color="g",label="Na_n")
line4, = ax4.plot(t , y[:,3], color="c",label="M_F")
line5, = ax5.plot(t , y[:,4], color="m",label="F")
line6, = ax6.plot(t , y[:,5], color="y",label="M")
line7, = ax7.plot(t , y[:,6], color="k",label="CP")
line8, = ax8.plot(t , y[:,7], color="y",label="CP2")
line9, = ax9.plot(t , y[:,8], color="k",label="TF")
line10, = ax10.plot(t , y[:,9], color="b",label="A")
line11, = ax11.plot(t , y[:,10], color="r",label="M_Ax")
line12, = ax12.plot(t , y[:,11], color="g",label="B")

ax1.set_ylabel('Number')
ax1.set_xlabel('Time')

ax2.legend(handles=[line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12])

plt.show()
'''

params[30] = 0
fig.add_subplot(4,3,2)
plt.title("APC mutant")
y = odeint(sim,y0,t1,args=(params,))
plt.plot(tt1, y[:,9], 'k:')
plt.plot(tt1, y[:,21], 'k')
plt.ylim([0, 200])

fig.add_subplot(4,3,5)
plt.plot(tt1, y[:,3], 'k:')
plt.plot(tt1, y[:,15], 'k')
plt.ylim([0, 5])

params[30] = 1
y = odeint(sim,y0,t1,args=(params,))
fig.add_subplot(4,3,8)
plt.plot(tt1, y[:,9], 'r:')
plt.plot(tt1, y[:,21], 'r')
plt.ylim([0, 200])

fig.add_subplot(4,3,11)
plt.plot(tt1, y[:,3], 'r:')
plt.plot(tt1, y[:,15], 'r')
plt.ylim([0, 5])

def sim2(variables, t, params):
    N = variables[0]  # Membrane-bound Notch receptor
    F = variables[1]  # Notch Intracellular Domain (NICD)
    I_1 = variables[2]  # Intermediate 1 (NICD/β-catenin)
    H_1 = variables[3]  # Hes1
    P = variables[4]  # Ngn3
    D = variables[5]  # Delta ligand
    H_2 = variables[6]  # Hath1
    G = variables[7]  # GSK3β
    C = variables[8]  # Destruction complex
    B = variables[9]  # Active β-catenin
    I_2 = variables[10]  # Intermediate 2 (GSK3β/β-catenin)
    A = variables[11]  # Axin

    N2 = variables[12]  # Membrane-bound Notch receptor
    F2 = variables[13]  # Notch Intracellular Domain (NICD)
    I_12 = variables[14]  # Intermediate 1 (NICD/β-catenin)
    H_12 = variables[15]  # Hes1
    P2 = variables[16]  # Ngn3
    D2 = variables[17]  # Delta ligand
    H_22 = variables[18]  # Hath1
    G2 = variables[19]  # GSK3β
    C2 = variables[20]  # Destruction complex
    B2 = variables[21]  # Active β-catenin
    I_22 = variables[22]  # Intermediate 2 (GSK3β/β-catenin)
    A2 = variables[23]  # Axin

    M1 = variables[24]  # PER mRNA
    P_11 = variables[25]
    P_21 = variables[26]
    P_31 = variables[27]
    P_41 = variables[28]
    P_51 = variables[29]
    P_61 = variables[30]
    P_N1 = variables[31]  # Nuclear PER protein

    M2 = variables[32]  # PER mRNA
    P_12 = variables[33]
    P_22 = variables[34]
    P_32 = variables[35]
    P_42 = variables[36]
    P_52 = variables[37]
    P_62 = variables[38]
    P_N2 = variables[39]  # Nuclear PER protein

    Dbar = params[0]
    m_1 = params[1]
    kappa_1 = params[2]
    mu_N = params[3]
    mu_F = params[4]
    alpha_frag = params[5]
    alpha_1 = params[6]
    mu_I1 = params[7]
    mu_H1 = params[8]
    sigma_K = params[9]
    theta_2 = params[10]
    m_2 = params[11]
    kappa_2 = params[12]
    theta_7 = params[13]
    m_7 = params[14]
    kappa_7 = params[15]
    sigma_2 = params[16]
    n_2 = params[17]
    mu_P = params[18]
    sigma_3 = params[19]
    n_3 = params[20]
    mu_D = params[21]
    theta_4 = params[22]
    m_4 = params[23]
    kappa_4 = params[24]
    mu_H2 = params[25]
    xi_5 = params[26]
    sigma_5 = params[27]
    n_5 = params[28]
    mu_G = params[29]
    W = params[30]
    alpha_2 = params[31]
    alpha_5 = params[32]
    mu_C = params[33]
    alpha_3 = params[34]
    mu_B = params[35]
    alpha_4 = params[36]
    mu_I2 = params[37]
    mu_A = params[38]
    theta_6 = params[39]
    m_6 = params[40]
    kappa_6 = params[41]
    xi_2 = params[42]
    theta_1 = params[43]
    xi_3 = params[44]

    alpha = params[45]
    beta = params[46]
    delta = params[47]
    mu = params[48]
    A_T = params[49]
    gamma = params[50]
    K_d = params[51]
    K_m = params[52]
    delta_N = params[53]

    dNdt = -mu_N * N + theta_1 * D2 ** m_1 / (kappa_1 ** m_1 + D2 ** m_1)
    dFdt = -mu_F * F + alpha_frag * mu_N * N - alpha_1 * B * F
    dI_1dt = -mu_I1 * I_1 + alpha_1 * B * F
    dH_1dt = -mu_H1 * H_1 + sigma_K / (sigma_K + Wnt(t)) * (
                theta_2 * I_1 ** m_2 / (kappa_2 ** m_2 + I_1 ** m_2) + theta_7 * B ** m_7 / (
                    kappa_7 ** m_7 + B ** m_7)) * xi_2 * sigma_2 ** n_2 / (sigma_2 ** n_2 + H_1 ** n_2) * 0.9 * 0.5 * (
                         A_T - P_N1 - K_d + ((A_T - P_N1 - K_d) ** 2 + 4 * K_d * A_T) ** 0.5)
    dPdt = -mu_P * P + xi_3 * sigma_3 ** n_3 / (sigma_3 ** n_3 + H_1 ** n_3)
    dDdt = -mu_D * D + theta_4 * P ** m_4 / (kappa_4 ** m_4 + P ** m_4)
    dH_2dt = -mu_H2 * H_2 + xi_5 * sigma_5 ** n_5 / (sigma_5 ** n_5 + H_1 ** n_5)
    dGdt = -mu_G * G + (1 + Wnt(t)) * alpha_2 + mu_C * C - alpha_5 * 1.4 * A ** 2 / (1 + (1 + Wnt(t)) ** 4) * G
    dCdt = -mu_C * C + alpha_5 * 1.4 * A ** 2 / (1 + (1 + Wnt(t)) ** 4) * G + mu_I2 * I_2 - alpha_3 * B * C
    dBdt = -mu_B * B + (1 + Wnt(t)) * alpha_4 - alpha_1 * B * F - alpha_3 * B * C
    dI_2dt = -mu_I2 * I_2 + alpha_3 * B * C
    dAdt = -mu_A * A + theta_6 * B ** m_6 / (kappa_6 ** m_6 + B ** m_6)  # 11

    dN2dt = -mu_N * N2 + theta_1 * D ** m_1 / (kappa_1 ** m_1 + D ** m_1)
    dF2dt = -mu_F * F2 + alpha_frag * mu_N * N2 - alpha_1 * B2 * F2
    dI_12dt = -mu_I1 * I_12 + alpha_1 * B2 * F2
    dH_12dt = -mu_H1 * H_12 + sigma_K / (sigma_K + W) * (
                theta_2 * I_12 ** m_2 / (kappa_2 ** m_2 + I_12 ** m_2) + theta_7 * B2 ** m_7 / (
                    kappa_7 ** m_7 + B2 ** m_7)) * xi_2 * sigma_2 ** n_2 / (
                          sigma_2 ** n_2 + H_12 ** n_2) * 0.9 * 0.5 * (
                          A_T - P_N2 - K_d + ((A_T - P_N2 - K_d) ** 2 + 4 * K_d * A_T) ** 0.5)
    dP2dt = -mu_P * P2 + xi_3 * sigma_3 ** n_3 / (sigma_3 ** n_3 + H_12 ** n_3)
    dD2dt = -mu_D * D2 + theta_4 * P2 ** m_4 / (kappa_4 ** m_4 + P2 ** m_4)
    dH_22dt = -mu_H2 * H_22 + xi_5 * sigma_5 ** n_5 / (sigma_5 ** n_5 + H_12 ** n_5)
    dG2dt = -mu_G * G2 + (1 + W) * alpha_2 + mu_C * C2 - alpha_5 * 1.4 * (A2 ** 2) / (1 + (1 + W) ** 4) * G2
    dC2dt = -mu_C * C2 + alpha_5 * 1.4 * (A2 ** 2) / (1 + (1 + W) ** 4) * G2 + mu_I2 * I_22 - alpha_3 * B2 * C2
    dB2dt = -mu_B * B2 + (1 + W) * alpha_4 - alpha_1 * B2 * F2 - alpha_3 * B2 * C2
    dI_22dt = -mu_I2 * I_22 + alpha_3 * B2 * C2
    dA2dt = -mu_A * A2 + theta_6 * B2 ** m_6 / (kappa_6 ** m_6 + B2 ** m_6)

    dM1dt = alpha * 0.5 * (A_T - P_N1 - K_d + ((A_T - P_N1 - K_d) ** 2 + 4 * K_d * A_T) ** 0.5) / (A_T) - delta * M1
    # a segédábrán P volt, de az szerintem a P_N-t akarta jelenteni...
    dP_11dt = mu * M1 - (beta + delta) * P_11
    dP_21dt = beta * P_11 - (beta + delta) * P_21
    dP_31dt = beta * P_21 - (beta + delta) * P_31
    dP_41dt = beta * P_31 - (beta + delta) * P_41
    dP_51dt = beta * P_41 - (beta + delta) * P_51
    dP_61dt = beta * P_51 - (gamma + delta) * P_61
    dP_N1dt = gamma * P_61 - delta_N * P_N1 / (K_m + P_N1)

    dM2dt = alpha * 0.5 * (A_T - P_N2 - K_d + ((A_T - P_N2 - K_d) ** 2 + 4 * K_d * A_T) ** 0.5) / (A_T) - delta * M2
    # a segédábrán P volt, de az szerintem a P_N-t akarta jelenteni...
    dP_12dt = mu * M2 - (beta + delta) * P_12
    dP_22dt = beta * P_12 - (beta + delta) * P_22
    dP_32dt = beta * P_22 - (beta + delta) * P_32
    dP_42dt = beta * P_32 - (beta + delta) * P_42
    dP_52dt = beta * P_42 - (beta + delta) * P_52
    dP_62dt = beta * P_52 - (gamma + delta) * P_62
    dP_N2dt = gamma * P_62 - delta_N * P_N2 / (K_m + P_N2)

    return ([dNdt, dFdt, dI_1dt, dH_1dt, dPdt, dDdt, dH_2dt, dGdt, dCdt, dBdt, dI_2dt, dAdt,
             dN2dt, dF2dt, dI_12dt, dH_12dt, dP2dt, dD2dt, dH_22dt, dG2dt, dC2dt, dB2dt, dI_22dt, dA2dt,
             dM1dt, dP_11dt, dP_21dt, dP_31dt, dP_41dt, dP_51dt, dP_61dt, dP_N1dt,
             dM2dt, dP_12dt, dP_22dt, dP_32dt, dP_42dt, dP_52dt, dP_62dt, dP_N2dt])

params[30] = 0
fig.add_subplot(4,3,3)
plt.title("Wnt mutant")
y = odeint(sim2,y0,t1,args=(params,))
plt.plot(tt1, y[:,9], 'k:')
plt.plot(tt1, y[:,21], 'k')
plt.ylim([0, 200])

fig.add_subplot(4,3,6)
plt.plot(tt1, y[:,3], 'k:')
plt.plot(tt1, y[:,15], 'k')
plt.ylim([0, 5])

def Wnt(t):
    if t<720:
        return 1
    else:
        return 2

params[30] = 1
fig.add_subplot(4,3,9)
y = odeint(sim2,y0,t1,args=(params,))
plt.plot(tt1, y[:,9], 'r:')
plt.plot(tt1, y[:,21], 'r')
plt.ylim([0, 200])

fig.add_subplot(4,3,12)
plt.plot(tt1, y[:,3], 'r:')
plt.plot(tt1, y[:,15], 'r')
plt.ylim([0, 5])

plt.show()