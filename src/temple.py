"""
Enfriamiento de un cilindro largo de acero durante el temple.

Ecuación de calor en coordenadas cilíndricas (simetría radial):

    ∂T/∂t = α [ ∂²T/∂r² + (1/r) ∂T/∂r ]

con convección en la superficie:  -k ∂T/∂r = h (T - T∞)

Se resuelve con diferencias finitas explícitas. En el centro (r = 0)
se usa el límite de L'Hôpital: el laplaciano vale 2·∂²T/∂r².
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import j0, j1


@dataclass
class Acero:
    k: float = 30.0       # W/m·K  (promedio entre 200 y 850 °C)
    rho: float = 7850.0   # kg/m³
    cp: float = 600.0     # J/kg·K

    @property
    def alpha(self):
        return self.k / (self.rho * self.cp)


# Coeficientes de convección efectivos (valores medios de referencia, W/m²·K)
MEDIOS = {
    "agua (agitada)": 5000.0,
    "aceite": 1200.0,
    "aire forzado": 150.0,
    "aire calmo": 25.0,
}


def simular(radio, h, T0=850.0, Tinf=25.0, t_final=600.0, n=41,
            material=Acero(), guardar_cada=None):
    """Devuelve (t, r, T) con T[k, i] = temperatura en el instante t[k] y radio r[i]."""
    a, k = material.alpha, material.k
    r = np.linspace(0, radio, n)
    dr = r[1] - r[0]
    bi_celda = h * dr / k
    dt = 0.9 * dr**2 / (4 * a * (1 + bi_celda))   # criterio de estabilidad
    pasos = int(np.ceil(t_final / dt))
    dt = t_final / pasos
    guardar_cada = guardar_cada or max(1, pasos // 2000)

    T = np.full(n, T0, dtype=float)
    hist_t, hist_T = [0.0], [T.copy()]
    c = a * dt / dr**2

    for p in range(1, pasos + 1):
        Tn = T.copy()
        # centro
        Tn[0] = T[0] + 4 * c * (T[1] - T[0])
        # nodos interiores
        i = np.arange(1, n - 1)
        Tn[i] = T[i] + c * ((T[i + 1] - 2 * T[i] + T[i - 1])
                            + (T[i + 1] - T[i - 1]) * dr / (2 * r[i]))
        # superficie con nodo fantasma
        fantasma = T[-2] - 2 * dr * h / k * (T[-1] - Tinf)
        Tn[-1] = T[-1] + c * ((fantasma - 2 * T[-1] + T[-2])
                              + (fantasma - T[-2]) * dr / (2 * r[-1]))
        T = Tn
        if p % guardar_cada == 0 or p == pasos:
            hist_t.append(p * dt)
            hist_T.append(T.copy())

    return np.array(hist_t), r, np.array(hist_T)


def solucion_analitica(r, t, radio, h, T0, Tinf, material=Acero(), terminos=30):
    """Serie exacta para un cilindro infinito con convección (Incropera, cap. 5)."""
    a, k = material.alpha, material.k
    Bi = h * radio / k
    f = lambda z: z * j1(z) - Bi * j0(z)
    raices, z = [], 1e-6
    while len(raices) < terminos:      # una raíz por intervalo entre ceros de J0
        z2 = z + 0.1
        if f(z) * f(z2) < 0:
            raices.append(brentq(f, z, z2))
        z = z2
    Fo = a * t / radio**2
    theta = np.zeros_like(np.asarray(r, dtype=float))
    for zn in raices:
        Cn = 2 / zn * j1(zn) / (j0(zn)**2 + j1(zn)**2)
        theta = theta + Cn * np.exp(-zn**2 * Fo) * j0(zn * np.asarray(r) / radio)
    return Tinf + (T0 - Tinf) * theta


def tiempo_cruce(t, curva, temperatura):
    """Instante en que una curva de enfriamiento baja de cierta temperatura."""
    idx = np.argmax(curva <= temperatura)
    if curva[idx] > temperatura:
        return np.nan          # nunca llega a esa temperatura
    if idx == 0:
        return t[0]
    t1, t2 = t[idx - 1], t[idx]
    T1, T2 = curva[idx - 1], curva[idx]
    return t1 + (temperatura - T1) * (t2 - t1) / (T2 - T1)


def t85(t, curva):
    """Tiempo de enfriamiento entre 800 y 500 °C (t8/5), indicador clásico
    de la microestructura que se forma al enfriar un acero."""
    return tiempo_cruce(t, curva, 500.0) - tiempo_cruce(t, curva, 800.0)
