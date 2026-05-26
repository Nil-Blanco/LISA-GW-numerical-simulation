import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. CONFIGURACIÓN Y CONSTANTES
# ==========================================
c = 2.998e8             # m/s
L_km = 2.5e6            # Longitud brazo LISA en km
L_sec = (L_km * 1000) / c

H0 = 2.26e-18           # s^-1 (Valor real de Hubble)
H_zero = 0.0            # Caso de control (H=0)

f_gw = 1e-1             # Frecuencia 0.1 Hz
omega = 2 * np.pi * f_gw

# Distancia a la fuente
Z_A_Gpc = 0.1           # Distancia en Gpc
Z_A_sec = (Z_A_Gpc * 1000 * 3.086e22) / c # Distancia en segundos-luz

epsilon = 1.0           # Amplitud relativa

# ==========================================
# 2. FUNCIONES FÍSICAS (Ecuaciones 10, 11 y 12 del TFG)
# ==========================================

def get_coefficients(alpha, H_val, w, L_s, Z_s):
    """
    Calcula los coeficientes A, B y C exactos según el TFG.
    """
    cos_a = np.cos(alpha)
    
    # Asumimos T_A approx Z_A (tiempo de viaje dominado por la distancia)
    # como se indica en la sección IV.C del TFG.
    T_s = Z_s 
    
    # Ec. 11: Coeficiente A
    A = (w * H_val * (L_s**2) / 2) * (cos_a - 2) * cos_a
    
    # Ec. 11: Coeficiente B
    # Simplificando la expresión con la aproximación T_A ~ Z_A
    B = L_s * w * (1 - H_val * Z_s - cos_a)
    
    # Ec. 12: Coeficiente C (VALOR REAL)
    # Fórmula: C = (w/2) * (-2*H*T*Z + H*Z^2 + 2*T - 2*Z)
    # Si T_s = Z_s, los términos 2*T - 2*Z se cancelan.
    # Queda el término cuadrático proporcional a H.
    C = (w / 2) * (-2 * H_val * T_s * Z_s + H_val * (Z_s**2) + 2 * T_s - 2 * Z_s)
    
    return A, B, C

def calcular_residuo_integral(alpha, H_val):
    """
    Calcula la integral del residuo temporal usando la regla del trapecio.
    """
    # Evitar división por cero o geometría nula
    if np.sin(alpha) == 0:
        return 0.0

    # Obtenemos coeficientes A, B, C con el valor de H correspondiente
    A, B, C = get_coefficients(alpha, H_val, omega, L_sec, Z_A_sec)
    
    # Malla de integración (x va de -1 a 0)
    # Usamos 1000 pasos para buena precisión (trapecios)
    x = np.linspace(-1, 0, 3000)
    
    # Fase Theta(x) = Ax^2 + Bx + C
    Theta = A * x**2 + B * x + C
    
    # Distancia R(x) aprox Z_A + x*L*cos(alpha)
    R_x = Z_A_sec + x * L_sec * np.cos(alpha)
    
    # Integrando: cos(Theta) / R(x)
    integrand = np.cos(Theta) / R_x
    
    # Resolver integral con trapecios
    integral_val = np.trapezoid(integrand, x)
    
    # Prefactor geométrico: - (L * epsilon / 2) * sin^2(alpha)
    prefactor = - (L_sec * epsilon / 2) * (np.sin(alpha)**2)
    
    # Devolvemos el valor absoluto (magnitud)
    return np.abs(prefactor * integral_val)

# ==========================================
# 3. BARRIDO SISTEMÁTICO (GRID DE LA ESFERA)
# ==========================================
# Geometría de los brazos (60 grados)
u1 = np.array([1, 0, 0])
deg60 = np.radians(60)
u2 = np.array([np.cos(deg60), np.sin(deg60), 0])

# ==========================================
# 3. GENERACIÓN ALEATORIA (MONTECARLO)
# ==========================================
# Vectores de los brazos (60 grados)
u1 = np.array([1, 0, 0])
deg60 = np.radians(60)
u2 = np.array([np.cos(deg60), np.sin(deg60), 0])

# Definir número total de fuentes aleatorias
N_sources = 100  # Ponemos 10000 para probar rápido, luego puedes subirlo a 45000

print(f"Generando {N_sources} fuentes de forma aleatoria...")

# Para distribuir aleatoriamente de forma uniforme sobre una ESFERA:
phi_vals = np.random.uniform(0, 2*np.pi, N_sources)
costheta_vals = np.random.uniform(-1, 1, N_sources)
theta_vals = np.arccos(costheta_vals)

# Listas para almacenar los resultados
res1_H0, res2_H0 = [], []      
res1_Hzero, res2_Hzero = [], [] 
alpha1_deg, alpha2_deg = [], [] 

# UN SOLO BUCLE para recorrer las fuentes aleatorias
for i in range(N_sources):
    th = theta_vals[i]
    ph = phi_vals[i]
    
    # Vector dirección de la fuente
    nx = np.sin(th) * np.cos(ph)
    ny = np.sin(th) * np.sin(ph)
    nz = np.cos(th)
    n_vec = np.array([nx, ny, nz])
    
    # Calcular ángulos de incidencia (.clip para evitar errores numéricos)
    cos_a1 = np.dot(n_vec, u1).clip(-1, 1)
    cos_a2 = np.dot(n_vec, u2).clip(-1, 1)
    
    alpha1 = np.arccos(cos_a1)
    alpha2 = np.arccos(cos_a2)
    
    # Guardamos los ángulos
    alpha1_deg.append(np.degrees(alpha1))
    alpha2_deg.append(np.degrees(alpha2))
    
    # --- CÁLCULO 1: Con Constante Cosmológica (H = H0) ---
    res1_H0.append(calcular_residuo_integral(alpha1, H0))
    res2_H0.append(calcular_residuo_integral(alpha2, H0))
    
    # --- CÁLCULO 2: Sin Constante Cosmológica (H = 0) ---
    res1_Hzero.append(calcular_residuo_integral(alpha1, H_zero))
    res2_Hzero.append(calcular_residuo_integral(alpha2, H_zero))

print("Cálculos terminados. Generando gráficos...")
# ==========================================
# 4. GRAFICAR RESULTADOS (Comparativa 2x2)
# ==========================================
# Hacemos la figura un poco más "cuadrada" de base (12x12)
fig, axs = plt.subplots(2, 2, figsize=(12, 12)) 

# Calcular la señal combinada (suma de magnitudes)
suma_H0 = np.array(res1_H0) + np.array(res2_H0)
suma_Hzero = np.array(res1_Hzero) + np.array(res2_Hzero)

# --- FILA 1: H = H0 ---
sc1 = axs[0, 0].scatter(alpha1_deg, alpha2_deg, c=suma_H0, cmap='viridis', s=5)
axs[0, 0].set_title(r'Geometría ($H=H_0$): Intensidad Combinada ($\tau_1 + \tau_2$)')
axs[0, 0].set_ylabel(r'$\alpha_2$ [deg]')
axs[0, 0].grid(True, alpha=0.3)
axs[0, 0].set_aspect('equal', adjustable='box') # <--- FUERZA ESCALA 1:1
plt.colorbar(sc1, ax=axs[0, 0], label=r'$\tau_{GW} (Total)$', shrink=0.8)

axs[0, 1].scatter(res1_H0, res2_H0, c='darkblue', s=2, alpha=0.5)
axs[0, 1].set_title(r'Consistencia de Señal ($H=H_0$)')
axs[0, 1].set_ylabel(r'$\tau_2$')
axs[0, 1].grid(True, alpha=0.3)
axs[0, 1].set_aspect('equal', adjustable='box') # <--- FUERZA ESCALA 1:1

# --- FILA 2: H = 0 ---
sc2 = axs[1, 0].scatter(alpha1_deg, alpha2_deg, c=suma_Hzero, cmap='inferno', s=5)
axs[1, 0].set_title(r'Geometría CONTROL ($H=0$): Intensidad Combinada ($\tau_1 + \tau_2$)')
axs[1, 0].set_xlabel(r'$\alpha_1$ [deg]')
axs[1, 0].set_ylabel(r'$\alpha_2$ [deg]')
axs[1, 0].grid(True, alpha=0.3)
axs[1, 0].set_aspect('equal', adjustable='box') # <--- FUERZA ESCALA 1:1
plt.colorbar(sc2, ax=axs[1, 0], label=r'$\tau_{GW} (Total)$', shrink=0.8)

axs[1, 1].scatter(res1_Hzero, res2_Hzero, c='darkred', s=2, alpha=0.5)
axs[1, 1].set_title(r'Consistencia de Señal CONTROL ($H=0$)')
axs[1, 1].set_xlabel(r'$\tau_1$')
axs[1, 1].set_ylabel(r'$\tau_2$')
axs[1, 1].grid(True, alpha=0.3)
axs[1, 1].set_aspect('equal', adjustable='box') # <--- FUERZA ESCALA 1:1

plt.tight_layout()
plt.show()