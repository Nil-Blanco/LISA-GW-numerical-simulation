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
# 3. GENERACIÓN ALEATORIA (MONTECARLO) Y BARRIDO DE H
# ==========================================
# Vectores de los brazos (60 grados)
u1 = np.array([1, 0, 0])
deg60 = np.radians(60)
u2 = np.array([np.cos(deg60), np.sin(deg60), 0])

N_sources = 100 # Ajustado para que 6 graficos no tarden una eternidad

print(f"Generando {N_sources} fuentes de forma aleatoria...")

phi_vals = np.random.uniform(0, 2*np.pi, N_sources)
costheta_vals = np.random.uniform(-1, 1, N_sources)
theta_vals = np.arccos(costheta_vals)

alpha1_deg, alpha2_deg = [], [] 

# Definir los factores multiplicativos para H0
h_factors = [0.0, 0.75, 0.9, 1.0, 1.1, 1.25]

# Usaremos un diccionario para guardar los resultados de cada caso
# Clave: factor de H (ej: 0.5). Valor: [lista_tau1, lista_tau2, lista_suma]
resultados = {factor: [[], [], []] for factor in h_factors}

# UN SOLO BUCLE para recorrer las fuentes aleatorias
for i in range(N_sources):
    th = theta_vals[i]
    ph = phi_vals[i]
    
    nx = np.sin(th) * np.cos(ph)
    ny = np.sin(th) * np.sin(ph)
    nz = np.cos(th)
    n_vec = np.array([nx, ny, nz])
    
    cos_a1 = np.dot(n_vec, u1).clip(-1, 1)
    cos_a2 = np.dot(n_vec, u2).clip(-1, 1)
    
    alpha1 = np.arccos(cos_a1)
    alpha2 = np.arccos(cos_a2)
    
    alpha1_deg.append(np.degrees(alpha1))
    alpha2_deg.append(np.degrees(alpha2))
    
    # Bucle interno para calcular los diferentes casos de H
    for factor in h_factors:
        H_val = factor * H0
        tau1 = calcular_residuo_integral(alpha1, H_val)
        tau2 = calcular_residuo_integral(alpha2, H_val)
        
        resultados[factor][0].append(tau1)
        resultados[factor][1].append(tau2)
        resultados[factor][2].append(tau1 + tau2) # Suma total

print("Cálculos terminados. Generando gráficos...")

# ==========================================
# 4. GRAFICAR RESULTADOS (Evolución de la Consistencia)
# ==========================================
# Creamos una figura con 2 filas y 3 columnas para los 6 casos de consistencia
fig, axs = plt.subplots(2, 3, figsize=(15, 10))

# Aplanamos el array de ejes para iterar fácilmente
axs = axs.flatten()

# Colores diferentes para cada caso para que quede bonito
colores = ['darkred', 'purple', 'indigo', 'darkblue', 'teal', 'darkgreen']

for idx, factor in enumerate(h_factors):
    tau1_data = resultados[factor][0]
    tau2_data = resultados[factor][1]
    
    # Seleccionar el subplot correspondiente
    ax = axs[idx]
    
    ax.scatter(tau1_data, tau2_data, c=colores[idx], s=1, alpha=0.5)
    
    # Título indicando el factor
    if factor == 0.0:
        ax.set_title(r'Control ($H = 0$)')
    else:
        ax.set_title(rf'$H = {factor} H_0$')
        
    ax.set_xlabel(r'$\tau_1$')
    ax.set_ylabel(r'$\tau_2$')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    # Ajustar límites ligeramente para que todas se vean a la misma escala si es posible
    # Calculamos el maximo global aproximado para fijar los ejes (opcional)
    # ax.set_xlim(0, max_global)
    # ax.set_ylim(0, max_global)

plt.tight_layout()
plt.suptitle(rf'Evolución de la Consistencia de Señal ($Z_A = {Z_A_Gpc}$ Gpc)', fontsize=16, y=1.02)
plt.show()

# ==========================================
# 5. GRAFICAR GEOMETRÍAS (Opcional, guardado aparte)
# ==========================================
# Si también quieres ver cómo cambia la "mancha", puedes hacer otra figura solo para H=0 vs H=1.5H0
# (Te lo dejo omitido para no saturar, pero los datos están en `resultados[factor][2]`)