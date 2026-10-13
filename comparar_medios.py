"""
Compara cómo se enfría un eje de acero de 50 mm de diámetro
templado desde 850 °C en distintos medios.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RAIZ = Path(__file__).resolve().parent
sys.path.insert(0, str(RAIZ / "src"))
from temple import MEDIOS, simular, t85

DIAMETRO = 0.050
R = DIAMETRO / 2
colores = ["#1f5fa8", "#d35400", "#27ae60", "#7f8c8d"]

fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
print(f"Eje de {DIAMETRO * 1000:.0f} mm, austenizado a 850 °C\n")
print(f"{'Medio':<16}{'h [W/m²K]':>10}{'t8/5 sup [s]':>14}{'t8/5 centro [s]':>17}{'ΔT máx [°C]':>13}")

resultados = {}
for (medio, h), col in zip(MEDIOS.items(), colores):
    t_fin = 300 if h > 1000 else 3000
    t, r, T = simular(R, h, t_final=t_fin, n=41)
    resultados[medio] = (t, r, T)
    sup, cen = T[:, -1], T[:, 0]
    dT = np.max(cen - sup)
    print(f"{medio:<16}{h:>10.0f}{t85(t, sup):>14.1f}{t85(t, cen):>17.1f}{dT:>13.0f}")

    axs[0].semilogx(t[1:], cen[1:], color=col, lw=2, label=f"{medio} (centro)")
    axs[0].semilogx(t[1:], sup[1:], color=col, lw=1.2, ls="--")

axs[0].axhspan(500, 800, color="0.85", zorder=0)
axs[0].text(0.15, 650, "rango t8/5", fontsize=8, color="0.35")
axs[0].set_xlim(0.1, 3000)
axs[0].set_xlabel("Tiempo [s]")
axs[0].set_ylabel("Temperatura [°C]")
axs[0].set_title("Curvas de enfriamiento (— centro, -- superficie)")
axs[0].legend(fontsize=8)
axs[0].grid(True, which="both", alpha=0.3)

# perfiles radiales en el temple en agua
t, r, T = resultados["agua (agitada)"]
for ts in (1, 5, 15, 30, 60, 120):
    k = np.argmin(np.abs(t - ts))
    axs[1].plot(r * 1000, T[k], lw=2, label=f"t = {ts} s")
axs[1].set_xlabel("Radio [mm]")
axs[1].set_ylabel("Temperatura [°C]")
axs[1].set_title("Perfil radial durante el temple en agua")
axs[1].legend(fontsize=8)
axs[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(RAIZ / "figuras" / "curvas_enfriamiento.png", dpi=150)
print("\nFigura guardada en figuras/curvas_enfriamiento.png")
