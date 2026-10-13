import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from temple import Acero, simular, solucion_analitica, t85, tiempo_cruce


@pytest.mark.parametrize("h", [25.0, 1200.0, 5000.0])
def test_coincide_con_solucion_analitica(h):
    R = 0.025
    t, r, T = simular(R, h, T0=850, Tinf=25, t_final=120, n=61)
    exacta_centro = solucion_analitica(0.0, t[-1], R, h, 850, 25)
    exacta_sup = solucion_analitica(R, t[-1], R, h, 850, 25)
    assert T[-1, 0] == pytest.approx(exacta_centro, abs=2.0)
    assert T[-1, -1] == pytest.approx(exacta_sup, abs=2.0)


def test_limite_capacitancia_concentrada():
    # Con Bi << 0.1 el cilindro se enfría casi uniforme: T = T∞ + ΔT·exp(-t/τ)
    R, h, m = 0.01, 10.0, Acero()
    t, r, T = simular(R, h, t_final=2000, n=21)
    tau = m.rho * m.cp * R / (2 * h)
    esperado = 25 + 825 * np.exp(-t[-1] / tau)
    assert T[-1].mean() == pytest.approx(esperado, rel=0.01)
    assert T[-1].max() - T[-1].min() < 1.0


def test_superficie_enfria_antes_que_el_centro():
    t, r, T = simular(0.03, 5000, t_final=60)
    assert np.all(T[1:, -1] <= T[1:, 0] + 1e-9)


def test_tiempo_cruce_interpola():
    t = np.array([0.0, 10.0, 20.0])
    c = np.array([900.0, 700.0, 500.0])
    assert tiempo_cruce(t, c, 800.0) == pytest.approx(5.0)
    assert t85(t, c) == pytest.approx(15.0)
