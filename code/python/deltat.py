import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm

# --- 1. CARGAR DATOS ---
nombre_archivo = 'lisa_evolucion_T.dat' 

if not os.path.exists(nombre_archivo):
    print(f"ERROR: No encuentro '{nombre_archivo}'. Ejecuta el Fortran primero.")
    exit()

data = np.loadtxt(nombre_archivo)
alpha = data[:, 0]
# Ahora taus contiene 101 columnas (desde T=Z+0s hasta T=Z+1000s cada 10s)
taus = data[:, 1:] 
num_casos = taus.shape[1]

# --- 2. PREPARAR DATOS (Normalización a 10^-16) ---
escala_display = 1e16  
taus_norm = taus * escala_display

# Seleccionamos casos para la tabla (mostramos 5 hitos: 0, 250, 500, 750, 1000s)
indices_tabla = [0, 25, 50, 75, 100] 
table_data = []
for i in indices_tabla:
    tau_actual = np.abs(taus_norm[:, i])
    idx_max = np.argmax(tau_actual)
    dt = i * 10
    table_data.append([f"{dt}s", f"{np.max(tau_actual):.2f}", f"{alpha[idx_max]:.1f}°"])

# --- 3. CONFIGURACIÓN DE LA GRÁFICA ---
fig, ax = plt.subplots(figsize=(12, 7))

# Mapa de colores para las 101 líneas (degradado suave)
colores = cm.viridis(np.linspace(0, 0.8, num_casos))

for i in range(num_casos):
    delta_t = i * 10
    # Solo ponemos etiqueta en la leyenda a los hitos principales para no saturar
    label_t = f'$\Delta T = {delta_t}$s' if i in indices_tabla else None
    
    ax.plot(alpha, np.abs(taus_norm[:, i]), color=colores[i], linewidth=0.8, 
            alpha=0.4, label=label_t)

# --- 4. LÍNEA DE LA MEDIA (La "Mitjana" pedida por el tutor) ---
# Calculamos el promedio de todos los instantes temporales
tau_media = np.mean(np.abs(taus_norm), axis=1)
ax.plot(alpha, tau_media, color='black', linewidth=3, linestyle='-', 
        label='Average (Mean Response)')

# --- 5. TABLA (Arriba Izquierda) ---
col_labels = ['Time', r'$|\tau|$ ($10^{-16}$)', r'$\alpha$ (deg)']
the_table = ax.table(cellText=table_data, colLabels=col_labels,
                     loc='upper left', bbox=[0.05, 0.68, 0.40, 0.25])
the_table.auto_set_font_size(False)
the_table.set_fontsize(9)

# --- 6. LEYENDA (Arriba Derecha) ---
# Mostramos solo los casos con etiqueta (los de indices_tabla + la media)
ax.legend(fontsize=9, loc='upper right', framealpha=0.9, shadow=True, title="Representative Times")

# --- 7. TÍTULO Y EJES ---
ax.set_title(r'LISA: $|\tau(\alpha)|$ Temporal Evolution (101 Cases, $\Delta T = 10$s) Z=0.05GPC', fontsize=14)
ax.set_xlabel(r'Angle $\alpha$ (deg)', fontsize=12)
ax.set_ylabel(r'$|\tau|$ Amplitude ($10^{-16}$)', fontsize=12)

# Formato de ejes profesional
ax.ticklabel_format(style='plain', axis='y') 
ax.grid(True, linestyle=':', alpha=0.5)
ax.set_xlim(0, 90)

plt.tight_layout()
plt.savefig('lisa_evolucion_completa.png', dpi=300)
print(f"Gráfica generada con {num_casos} casos y línea de promedio.")
plt.show()