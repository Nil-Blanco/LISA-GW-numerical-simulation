import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. CONFIGURACIÓN Y CONSTANTES
# ==========================================
c = 2.998e8
L_km = 2.5e6
L_sec = (L_km * 1000) / c
H0 = 2.26e-18
f_gw = 1e-1
omega = 2 * np.pi * f_gw
Z_A_Gpc = 0.1
Z_A_sec = (Z_A_Gpc * 1000 * 3.086e22) / c
epsilon = 1.0

def get_coefficients(alpha, H_val, w, L_s, Z_s):
    cos_a = np.cos(alpha)
    T_s = Z_s 
    A = (w * H_val * (L_s**2) / 2) * (cos_a - 2) * cos_a
    B = L_s * w * (1 - H_val * Z_s - cos_a)
    C = (w / 2) * (-2 * H_val * T_s * Z_s + H_val * (Z_s**2) + 2 * T_s - 2 * Z_s)
    return A, B, C

def calcular_residuo_integral(alpha, H_val):
    if np.sin(alpha) == 0:
        return 0.0
    A, B, C = get_coefficients(alpha, H_val, omega, L_sec, Z_A_sec)
    x = np.linspace(-1, 0, 1000)
    Theta = A * x**2 + B * x + C
    R_x = Z_A_sec + x * L_sec * np.cos(alpha)
    integrand = np.cos(Theta) / R_x
    integral_val = np.trapezoid(integrand, x)
    prefactor = - (L_sec * epsilon / 2) * (np.sin(alpha)**2)
    return np.abs(prefactor * integral_val)

# ==========================================
# 2. GENERACIÓN ALEATORIA
# ==========================================
u1 = np.array([1, 0, 0])
deg60 = np.radians(60)
u2 = np.array([np.cos(deg60), np.sin(deg60), 0])

N_sources = 2000 
print(f"Generando {N_sources} fuentes...")

phi_vals = np.random.uniform(0, 2*np.pi, N_sources)
costheta_vals = np.random.uniform(-1, 1, N_sources)
theta_vals = np.arccos(costheta_vals)

alpha1_vals = np.arccos(np.dot(np.array([np.sin(theta_vals)*np.cos(phi_vals), 
                                         np.sin(theta_vals)*np.sin(phi_vals), 
                                         np.cos(theta_vals)]).T, u1).clip(-1, 1))
alpha2_vals = np.arccos(np.dot(np.array([np.sin(theta_vals)*np.cos(phi_vals), 
                                         np.sin(theta_vals)*np.sin(phi_vals), 
                                         np.cos(theta_vals)]).T, u2).clip(-1, 1))

# ==========================================
# 3. BARRIDO FINO DE H Y CÁLCULO DE PROPORCIÓN
# ==========================================
# Evaluamos 30 puntos entre H=0 y H=2.0*H0
h_factors = np.linspace(0.7, 1.3, 333)
proporciones = []

# DEFINIMOS LA ZONA OSCURA (El cuadrado inferior izquierdo)
# Ajusta este valor si consideras que la zona es más grande o más pequeña
umbral = 0.55e-16

print("Calculando la concentración en la zona de baja señal para distintos valores de H...")

for factor in h_factors:
    H_val = factor * H0
    
    # Contadores para este H
    en_zona = 0
    
    for i in range(N_sources):
        tau1 = calcular_residuo_integral(alpha1_vals[i], H_val)
        tau2 = calcular_residuo_integral(alpha2_vals[i], H_val)
        
        # Condición para estar en la zona baja
        if tau1 > umbral and tau2 > umbral:
            en_zona += 1
            
    # Calculamos el porcentaje
    proporcion = (en_zona / N_sources) * 100
    proporciones.append(proporcion)

print("Cálculos finalizados. Graficando...")

# ==========================================
# 4. GRAFICAR RESULTADOS
# ==========================================
plt.figure(figsize=(10, 6))
plt.plot(h_factors, proporciones, marker='o', linestyle='-', color='purple', linewidth=2)

plt.axvline(x=1.0, color='gray', linestyle='--', alpha=0.7, label=r'Valor Real ($H_0$)')

plt.title(rf'Concentración de eventos en zona de baja señal vs $H$ ($Z_A = {Z_A_Gpc}$ Gpc)', fontsize=14)
plt.xlabel(r'Factor multiplicativo de la Constante de Hubble ($H / H_0$)', fontsize=12)
plt.ylabel(r'% de fuentes en la zona ($\tau_1, \tau_2 < 0.4 \times 10^{-16}$)', fontsize=12)

plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()