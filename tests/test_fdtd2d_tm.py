# HUELLA-EJECUTOR-FDTD2D-TM-20261021
"""Ejecutor directo de `fdtd2d_tm.py` — fotonica-monotile.

`fdtd2d_tm.py` es el solver FDTD Yee con PML de campo partido que produce la
cifra C1 (piso de ruido del instrumento, medio homogeneo eps_r = 1) y el
contrato `eps` que `inyector_diel.py` consume. Hasta hoy ningun test lo abria:
el censo `censo_numpy_dependencia_20261020.py` lo listaba `SIN_EJECUTOR`
(junto con `inyector_diel.py`), y `malla_hat.py` era el unico portador
ejecutado. Un import que nadie ejecuta no puede fallar.

Lo que se mide aqui, medido antes de escribirse (tres sondas read-only):

  * contrato: `correr(L_a, modo, eps=None) -> dict` con 31 claves, `nx =
    round(L_a/DX)`, `ny = 1` en C1a y `ny = nx` en C1b, `nt = round(T_TOTAL/DT)`,
    29 bins espectrales en [NU_MIN, NU_MAX];
  * fallos: la sonda corta levanta `SystemExit` con `.code` **str** (no hay RC
    numerico); la matriz `eps` de forma ajena y la que toca el borde de la PML
    levantan `AssertionError`;
  * imports: stdlib + numpy. Cero dependencia cross-repo.

El oraculo barato es la corrida C1a de L = 24a: 384x1 celdas, 2304 pasos,
0.03..0.06 s. La corrida C1b de la misma L (384x384, 5.14 s) se hace UNA vez,
porque es la unica forma de medir la rama 2D (`ny = nx`, reparto cilindrico).

Desviaciones declaradas (no ocultas):

  1. `main()` NO se corre: encadena cuatro corridas (C1b en L = 36a son 576^2,
     ~11 s) y escribe dos ficheros en `crudos/`. `main()` solo orquesta;
     el contrato fisico vive en `correr()`, y es `correr()` lo que se asierta.
  2. `correr()` NO falla cerrado ante un modo desconocido: `correr(24, "ZZZ")`
     devuelve `modo="ZZZ"` con `ny = nx` y `reparto_teorico = sqrt(1/2)` en
     3.93 s, es decir un C1b silencioso. Es un fail-open del objeto (misma
     familia que `sim_acustica2d.correr`), medido y declarado aqui, y
     deliberadamente NO asertado: asertarlo seria bendecir el agujero.
  3. `TOL_CIFRA = 1e-3` (relativa) no es una promesa de bit-exactitud: dos
     corridas en esta maquina son identicas, pero la leccion de `8117005` es
     que el bit depende de la maquina y el numero no. La tolerancia mide
     deriva entre corridas; el umbral que manda es el que el propio objeto
     declara (UMBRAL, FACTOR_Q4).

Falsador declarado (lo que hace rojo a este test):

  1. numpy ausente -> `import fdtd2d_tm` levanta `ModuleNotFoundError`.
  2. la dispersion de C1a deja de caber en `FACTOR_Q4 * yee_max` -> cae el
     test de Q4 (mutante `S` o malla alterada).
  3. el reparto deja de cumplir su teorico -> caen Q3 (C1a y C1b).
  4. la guarda de sonda desaparece -> el `SystemExit` no llega y cae su test.
  5. la validacion de `eps` desaparece -> caen los dos tests de `eps`.

Mutantes medidos en copia (`/tmp/fot-fals`, restaurado despues): 5 de 5 muertos.
  M1 `S 0.5->0.45`            -> cae la regresion de cifras (1 failed)
  M2 `reparto_teorico` sin la rama C1a -> caen Q3-C1a y la regresion (2 failed)
  M3 guarda de sonda `if False:` -> cae el `SystemExit` esperado
  M4 `assert eps.shape` anulado  -> cae el rechazo de forma ajena
  M5 `assert` del borde anulado  -> cae el rechazo de eps en la PML
M1 solo lo muerde la regresion de cifras: Q4 compara la fase medida contra el
compañero analitico de Yee, y los dos usan la misma `S`, asi que un cambio de
Courant los mueve juntos. El umbral auto-consistente no sustituye a la cifra.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

RAIZ = Path(__file__).resolve().parents[1]
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

import fdtd2d_tm as solver  # noqa: E402

L_NULO = 24.0
TOL_CIFRA = 1e-3

CLAVES = frozenset(
    {
        "modo", "L_a", "nx", "ny", "nt", "segundos_corrida", "nu",
        "A_d_1", "A_d_2", "lnA_d_1", "lnA_d_2", "riz_d_1", "riz_d_2",
        "reparto_media", "reparto_desv_max", "reparto_teorico", "reparto_riz",
        "reparto_min", "reparto_max", "fase_error", "fase_error_max", "rama_m",
        "cola_sonda", "cola_dominio", "cola_dominio_celda", "cola_dominio_paso",
        "cola_dominio_t", "cola_final_celda", "cola_final", "residuo_media_Ez",
        "pico_sonda",
    }
)

CIFRAS_C1A = {
    "reparto_media": 0.9999071789731048,
    "reparto_desv_max": 0.005474581720536387,
    "reparto_riz": 0.005394204639334759,
    "fase_error_max": 0.0014045811674006092,
}


@pytest.fixture(scope="module")
def c1a() -> dict:
    return solver.correr(L_NULO, "C1a")


@pytest.fixture(scope="module")
def c1b() -> dict:
    return solver.correr(L_NULO, "C1b")


def test_el_solver_corre_y_devuelve_su_contrato(c1a):
    assert set(c1a) == set(CLAVES), set(c1a) ^ set(CLAVES)
    assert c1a["modo"] == "C1a"
    assert c1a["L_a"] == L_NULO
    n = len(c1a["nu"])
    assert n == len(c1a["A_d_1"]) == len(c1a["A_d_2"]) == len(c1a["fase_error"])
    assert n == len(c1a["lnA_d_1"]) == len(c1a["lnA_d_2"])


def test_la_malla_de_C1a_es_cuasi_1d_y_la_de_C1b_es_2d(c1a, c1b):
    nx = int(round(L_NULO / solver.DX))
    assert c1a["nx"] == nx == c1b["nx"]
    assert c1a["ny"] == 1
    assert c1b["ny"] == c1b["nx"] == nx
    assert c1a["nt"] == int(round(solver.T_TOTAL / solver.DT))


def test_la_ventana_espectral_cumple_sus_bordes_declarados(c1a):
    nu = np.asarray(c1a["nu"], dtype=np.float64)
    assert len(nu) == 29
    assert nu.min() >= solver.NU_MIN and nu.max() <= solver.NU_MAX
    paso = float(np.diff(nu).max())
    esperado = 1.0 / ((c1a["nt"] + 1) * solver.DT)  # rfftfreq usa Nt+1 muestras, no T_total
    assert math.isclose(paso, esperado, rel_tol=1e-9), (paso, esperado)


def test_la_dispersion_de_C1a_cabe_en_el_factor_Q4_del_error_de_Yee(c1a):
    """Q4 del instrumento: fase medida <= FACTOR_Q4 x error analitico de Yee."""
    nus = np.linspace(solver.NU_MIN, solver.NU_MAX, 512)
    yee = float(np.abs(solver.yee_error_fase(nus)).max())
    tope = solver.FACTOR_Q4 * yee
    assert yee > 0.0
    assert 0.0 < c1a["fase_error_max"] <= tope, (c1a["fase_error_max"], tope)


def test_el_reparto_plano_de_C1a_cumple_su_teorico(c1a):
    """Q3 en cuasi-1D: el reparto entre sondas es 1.0 (onda plana)."""
    assert c1a["reparto_teorico"] == 1.0
    assert c1a["reparto_desv_max"] <= solver.UMBRAL, c1a["reparto_desv_max"]
    assert c1a["reparto_riz"] <= solver.UMBRAL, c1a["reparto_riz"]


def test_el_reparto_cilindrico_de_C1b_cumple_su_teorico(c1b):
    """Q3 en 2D: el reparto decae como sqrt(d1/d2) = 0.7071."""
    teorico = math.sqrt(solver.D_SONDAS[0] / solver.D_SONDAS[1])
    assert math.isclose(c1b["reparto_teorico"], teorico, rel_tol=0.0, abs_tol=1e-15)
    assert c1b["reparto_desv_max"] <= solver.UMBRAL, c1b["reparto_desv_max"]
    assert c1b["reparto_riz"] <= solver.UMBRAL, c1b["reparto_riz"]


def test_la_regresion_C1a_reproduce_las_cifras_medidas(c1a):
    """Deriva entre corridas, no bit-exactitud (leccion de `8117005`)."""
    for clave, medido in CIFRAS_C1A.items():
        assert math.isclose(c1a[clave], medido, rel_tol=TOL_CIFRA), (clave, c1a[clave])


def test_el_pulso_no_deja_una_cola_por_encima_del_umbral(c1a):
    assert c1a["cola_sonda"] <= solver.UMBRAL, c1a["cola_sonda"]
    assert c1a["cola_dominio"] <= solver.UMBRAL, c1a["cola_dominio"]
    assert c1a["pico_sonda"] > 0.0


def test_una_malla_mas_corta_que_la_sonda_falla_cerrado():
    for L in (16.0, 19.0):
        with pytest.raises(SystemExit) as exc:
            solver.correr(L, "C1a")
        assert isinstance(exc.value.code, str), exc.value.code
        assert "demasiado corta para la sonda" in exc.value.code


def test_una_eps_de_forma_ajena_a_la_malla_no_se_acepta():
    nx = int(round(L_NULO / solver.DX))
    with pytest.raises(AssertionError, match="malla activa"):
        solver.correr(L_NULO, "C1a", eps=np.zeros((8, 8), dtype=np.float64))
    assert nx > 0


def test_una_eps_que_toca_la_pml_se_rechaza():
    """La conductividad nominal de la PML no admite dielectrico en el borde."""
    nx = int(round(L_NULO / solver.DX))
    eps = np.ones((nx, 1), dtype=np.float64)
    eps[0, 0] = 2.0
    with pytest.raises(AssertionError, match="PML debe conservar"):
        solver.correr(L_NULO, "C1a", eps=eps)


def test_una_eps_homogenea_reproduce_el_camino_certificado(c1a):
    """El hook `eps` con eps_r = 1 no puede mover una sola cifra de C1."""
    nx = int(round(L_NULO / solver.DX))
    con_hook = solver.correr(L_NULO, "C1a", eps=np.ones((nx, 1), dtype=np.float64))
    distintas = [k for k in c1a if k != "segundos_corrida" and con_hook[k] != c1a[k]]
    assert distintas == [], distintas


def test_correr_no_escribe_en_disco():
    """Solo `main()` escribe crudos; la funcion que se audita es pura."""
    antes = sorted(p.name for p in solver.CRUDOS.iterdir())
    solver.correr(L_NULO, "C1a")
    despues = sorted(p.name for p in solver.CRUDOS.iterdir())
    assert antes == despues, set(despues) ^ set(antes)
