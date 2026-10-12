#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fotonica-monotile · acto 1 · instrumento FDTD 2D TM + piso de ruido C1.

UNIVERSO DECLARADO (se imprime antes de toda cifra):
  unidades normalizadas a = 1, c = 1, eps0 = mu0 = 1; frecuencia reducida
  nu = omega a / (2 pi c); malla Yee 2D en modo TM (Ez, Hx, Hy) con
  dx = dy = a/16; Courant S = c dt/dx = 0.5 (limite 2D 1/sqrt(2) = 0.7071);
  esponja de 20 celdas, perfil borde^3, alfa = 0.30 por paso (sin PML: es una
  esponja declarada y su reflexion ES lo que C1 mide); fuente = pulso gaussiano
  sigma = 0.5 a/c, t0 = 3 sigma, corriente BLANDA en una sola columna i = 2a;
  sondas en columnas de posicion ABSOLUTA d = 8a y d = 16a sobre la fila central
  (al crecer L el dominio solo se alarga hacia la derecha); ventana de medida
  nu en [0.15, 0.55] (el borde bajo lo fija que la sonda cercana quede a
  >= 1 lambda); T_total = 72 a/c, Nt = 2304 pasos, dnu = 1/T_total;
  umbral de planitud declarado ANTES de medir = 0.05 en |delta ln A|;
  doble precision; motor numpy; unidad de conteo = celda.

QUE MIDE C1 (dos sub-corridas, ambas a L = 24a y L = 36a)
  Dominio HOMOGENEO eps_r = 1 con la misma malla, paso, esponja, fuente y
  sondas. No hay objeto. Entrega:
    A(L, nu)  = |FFT(sonda)| / |S(omega)|   con |S| = exp(-sigma^2 omega^2 / 2)
    riz(L)    = max_nu |ln A - media_nu ln A|          -> criterio <= 0.05
    delta     = max_nu |ln A(36a) - ln A(24a)|         -> criterio <= 0.05
    reparto   = A(16a)/A(8a)  (ley 2D: sqrt(8/16) = 0.7071 · en 1D: 1.0000)
    cola      = max|sonda|(t > 0.75 T_total) / pico     -> piso de reflexion
  C1a = cuasi-1D (ny=1, esponja SOLO en x, onda plana): aqui vive la
  dispersion numerica, comparable EXACTA con la formula de Yee, y el piso de
  reflexion mas afilado. C1b = 2D real (ny=nx, esponja en los cuatro bordes,
  fuente de columna finita): aqui viven el reparto cilindrico y la dependencia
  espuria con L en la geometria que el objeto va a usar.
  Y el acompaneante ANALITICO de la malla de Yee
    sin(omega dt/2) = S sin(k dx/2)   ->  v_p = omega/k   -> error de fase.
  El error de fase MEDIDO sale de la diferencia de fase entre las dos sondas
  (k = fase / delta_d) con la convencion de rama declarada: m se elige una sola
  vez en el primer punto de la ventana.

REGLA DURA
  Un decaimiento solo cuenta POR ENCIMA de este piso. Si C1 no es plano se
  para y se repara el instrumento ANTES de medir el objeto (pre-registro §4.1).

Uso:
  python3 -B fdtd2d_tm.py
  python3 -B fdtd2d_tm.py --json
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

# ------------------------------------------------------------------ declarado
C = 1.0
A = 1.0
DX = A / 16.0
DY = DX
S = 0.5                      # Courant = c dt / dx
DT = S * DX / C
ESPONJA = 20
ALFA = 0.30
ORDEN = 3
SIGMA = 0.5
T0 = 3.0 * SIGMA
I_FUENTE = int(round(2.0 * A / DX))
D_SONDAS = (8.0 * A, 16.0 * A)
L_LISTA = (24.0 * A, 36.0 * A)
MODOS = ("C1a", "C1b")
NU_MIN = 0.15
NU_MAX = 0.55
T_TOTAL = 72.0
UMBRAL = 0.05
FACTOR_Q4 = 3.0
FRACCION_COLA = 0.75

AQUI = Path(__file__).resolve().parent
CRUDOS = AQUI / "crudos"

UNIVERSO = (
    "unidades normalizadas a=1, c=1, eps0=mu0=1; nu = omega a/(2 pi c); "
    f"malla Yee dx=dy=a/16 ({DX:g}); Courant S={S} (limite 2D "
    f"{1/math.sqrt(2):.4f}); esponja {ESPONJA} celdas, perfil borde^{ORDEN}, "
    f"alfa={ALFA}; fuente = pulso gaussiano sigma={SIGMA} a/c, t0={T0} a/c, "
    f"corriente blanda en la columna i={I_FUENTE}; sondas en posicion ABSOLUTA "
    f"d={D_SONDAS[0]:g}a y d={D_SONDAS[1]:g}a; L en {L_LISTA[0]:g}a y "
    f"{L_LISTA[1]:g}a; ventana nu en [{NU_MIN}, {NU_MAX}] (el borde bajo lo fija "
    f"que la sonda cercana quede a >= 1 lambda: 1/{NU_MIN} = {1.0/NU_MIN:.2f} a); "
    f"T_total={T_TOTAL} a/c, Nt={int(round(T_TOTAL/DT))}, dnu={1.0/T_TOTAL:.4f}; "
    f"umbral de planitud = {UMBRAL} en |delta ln A|; doble precision; "
    "homogeneo eps_r=1 (NO hay objeto en C1); motor numpy; unidad de conteo = celda; "
    "dos sub-corridas: C1a = cuasi-1D (ny=1, esponja solo en x, onda plana -> aqui vive "
    "la dispersion numerica, comparable EXACTA con la formula de Yee) y C1b = 2D real "
    "(ny=nx, esponja en los cuatro bordes, fuente de columna finita -> aqui viven el "
    "reparto cilindrico y la dependencia espuria con L)"
)


def t_utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# ------------------------------------------------------------------ instrumento
def mascara_esponja(nx: int, ny: int, ejes=(0, 1)) -> np.ndarray:
    m = np.ones((nx, ny), dtype=np.float64)
    for eje, n in ((0, nx), (1, ny)):
        if eje not in ejes:
            continue
        idx = np.arange(n)
        borde = np.maximum(
            np.clip(ESPONJA - idx, 0, ESPONJA),
            np.clip(idx - (n - 1 - ESPONJA), 0, ESPONJA),
        ) / float(ESPONJA)
        s = 1.0 - ALFA * borde ** ORDEN
        m *= s.reshape((-1, 1) if eje == 0 else (1, -1))
    return m


def fuente(t: float) -> float:
    return math.exp(-((t - T0) ** 2) / (2.0 * SIGMA ** 2))


def espectro_fuente(nu: np.ndarray) -> np.ndarray:
    om = 2.0 * math.pi * nu
    return np.exp(-(SIGMA ** 2) * om ** 2 / 2.0)


def yee_error_fase(nu: np.ndarray) -> np.ndarray:
    """(1/dt)^2 sin^2(omega dt/2) = (1/dx)^2 [sin^2(kx dx/2)+sin^2(ky dy/2)]
    sobre el eje: sin(k dx/2) = sin(omega dt/2)/S ; v_p = omega/k."""
    om = 2.0 * math.pi * nu
    arg = np.sin(om * DT / 2.0) / S
    arg = np.clip(arg, -1.0, 1.0)
    k = 2.0 * np.arcsin(arg) / DX
    vp = np.divide(om, k, out=np.ones_like(om), where=(k > 0))
    return vp / C - 1.0


def correr(L_a: float, modo: str) -> dict:
    nx = int(round(L_a / DX))
    ny = 1 if modo == "C1a" else nx
    ejes = (0,) if modo == "C1a" else (0, 1)
    eps = np.ones((nx, ny), dtype=np.float64)      # C1: sin objeto
    inv_eps = 1.0 / eps
    m = mascara_esponja(nx, ny, ejes)
    nt = int(round(T_TOTAL / DT))
    jc = ny // 2
    idx = [I_FUENTE + int(round(d / DX)) for d in D_SONDAS]
    if idx[-1] >= nx - ESPONJA:
        raise SystemExit(
            f"malla L={L_a:g}a demasiado corta para la sonda d={D_SONDAS[-1]:g}a"
        )
    ez = np.zeros((nx, ny))
    hx = np.zeros((nx, ny))
    hy = np.zeros((nx, ny))
    s = [np.zeros(nt + 1), np.zeros(nt + 1)]
    cola_dominio = 0.0
    n_cola = int(round(FRACCION_COLA * nt))
    kh = DT / DY
    t0 = time.time()
    for n in range(nt + 1):
        t = n * DT
        hx[:, :-1] -= kh * (ez[:, 1:] - ez[:, :-1])
        hy[:-1, :] += kh * (ez[1:, :] - ez[:-1, :])
        ez[1:, 1:] += DT * inv_eps[1:, 1:] * (
            (hy[1:, 1:] - hy[:-1, 1:]) / DX - (hx[1:, 1:] - hx[1:, :-1]) / DY
        )
        ez[I_FUENTE, :] += fuente(t)
        ez *= m
        hx *= m
        hy *= m
        s[0][n] = ez[idx[0], jc]
        s[1][n] = ez[idx[1], jc]
        if n >= n_cola:
            v = float(np.abs(ez).max())
            if v > cola_dominio:
                cola_dominio = v
    transcurrido = time.time() - t0

    frec = np.fft.rfftfreq(nt + 1, DT)                    # unidades c/a = nu
    fft = [np.fft.rfft(x) * DT for x in s]
    sel = (frec >= NU_MIN) & (frec <= NU_MAX)
    nu = frec[sel]
    sf = espectro_fuente(nu)
    A = [np.abs(fft[i][sel]) / sf for i in (0, 1)]
    lnA = [np.log(np.maximum(A[i], 1e-300)) for i in (0, 1)]

    # --- fase medida entre las dos sondas (convencion de rama declarada) ---
    dd = D_SONDAS[1] - D_SONDAS[0]
    om = 2.0 * math.pi * nu
    ang = np.unwrap(np.angle(fft[1][sel] * np.conj(fft[0][sel])))
    mm = float(np.round((om[0] / C * dd - ang[0]) / (2.0 * math.pi)))
    k = (ang + 2.0 * math.pi * mm) / dd
    fase_err = np.divide(
        om, k * C, out=np.ones_like(om), where=(np.abs(k) > 1e-12)
    ) - 1.0

    pico2 = float(np.abs(s[1]).max())
    cola_sonda = float(np.abs(s[1][n_cola:]).max()) / pico2 if pico2 > 0 else float("nan")
    cola_dom = cola_dominio / pico2 if pico2 > 0 else float("nan")
    reparto_teorico = 1.0 if modo == "C1a" else math.sqrt(D_SONDAS[0] / D_SONDAS[1])

    return {
        "modo": modo,
        "L_a": L_a,
        "nx": nx,
        "ny": ny,
        "nt": nt,
        "segundos_corrida": round(transcurrido, 2),
        "nu": nu.tolist(),
        "A_d_1": A[0].tolist(),
        "A_d_2": A[1].tolist(),
        "lnA_d_1": lnA[0].tolist(),
        "lnA_d_2": lnA[1].tolist(),
        "riz_d_1": float(np.abs(lnA[0] - lnA[0].mean()).max()),
        "riz_d_2": float(np.abs(lnA[1] - lnA[1].mean()).max()),
        "reparto_media": float(np.mean(A[1] / A[0])),
        "reparto_desv_max": float(np.abs(A[1] / A[0] - reparto_teorico).max()),
        "reparto_teorico": reparto_teorico,
        "fase_error": fase_err.tolist(),
        "fase_error_max": float(np.abs(fase_err).max()),
        "rama_m": mm,
        "cola_sonda": cola_sonda,
        "cola_dominio": cola_dom,
        "pico_sonda": pico2,
    }


def sha16(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def main() -> int:
    nus = np.linspace(NU_MIN, NU_MAX, 512)
    yee_max = float(np.abs(yee_error_fase(nus)).max())
    salida = {
        "instrumento": "fotonica-monotile/fdtd2d_tm.py",
        "acto": "C1 · piso de ruido del instrumento (homogeneo eps_r=1)",
        "UNIVERSO DECLARADO": UNIVERSO,
        "t_utc_inicio": t_utc(),
        "umbral_planitud": UMBRAL,
        "ventana_nu": [NU_MIN, NU_MAX],
        "sondas_a": [d / A for d in D_SONDAS],
        "modos": list(MODOS),
        "L_a": [l / A for l in L_LISTA],
        "yee_error_fase_max_ventana": yee_max,
        "factor_q4": FACTOR_Q4,
        "corridas": [],
    }
    for modo in MODOS:
        for L in L_LISTA:
            salida["corridas"].append(correr(L, modo))

    por = {(c["modo"], c["L_a"]): c for c in salida["corridas"]}
    deltas = {}
    for modo in MODOS:
        a = np.array(por[(modo, L_LISTA[0])]["lnA_d_2"])
        b = np.array(por[(modo, L_LISTA[1])]["lnA_d_2"])
        deltas[modo] = float(np.abs(b - a).max())
    salida["delta_lnA_36_vs_24"] = deltas
    riz_max = max(max(c["riz_d_1"], c["riz_d_2"]) for c in salida["corridas"])
    salida["riz_max"] = riz_max
    salida["Q1_riz_ok"] = bool(riz_max <= UMBRAL)
    salida["Q2_delta_ok"] = bool(max(deltas.values()) <= UMBRAL)
    desv = {m: max(por[(m, L)]["reparto_desv_max"] for L in L_LISTA) for m in MODOS}
    salida["reparto_desv_max_por_modo"] = desv
    salida["Q3_reparto_ok"] = bool(max(desv.values()) <= UMBRAL)
    fase_max = max(por[("C1a", L)]["fase_error_max"] for L in L_LISTA)
    salida["fase_error_max_C1a"] = fase_max
    salida["Q4_fase_ok"] = bool(fase_max <= FACTOR_Q4 * yee_max)
    salida["t_utc_fin"] = t_utc()
    salida["sha16_script"] = sha16(Path(__file__).resolve())

    etiqueta = salida["t_utc_inicio"].replace("-", "").replace(":", "")
    CRUDOS.mkdir(parents=True, exist_ok=True)
    pj = CRUDOS / f"crudo_C1_{etiqueta}.json"
    pj.write_text(json.dumps(salida, indent=1, ensure_ascii=False), encoding="utf-8")
    salida["ruta_crudo_json"] = str(pj)
    salida["sha16_crudo_json"] = sha16(pj)

    Ln = []
    Ln.append("UNIVERSO DECLARADO: " + salida["UNIVERSO DECLARADO"])
    Ln.append(f"instrumento: {salida['instrumento']} · sha16 {salida['sha16_script']}")
    Ln.append(f"t_utc: {salida['t_utc_inicio']} -> {salida['t_utc_fin']}")
    Ln.append(f"umbral de planitud declarado: {UMBRAL} · factor Q4: {FACTOR_Q4}")
    Ln.append("")
    Ln.append("---- C1 · homogeneo eps_r = 1 (NO hay objeto) ----")
    for c in salida["corridas"]:
        Ln.append(
            f"{c['modo']} · L = {c['L_a']:g}a · malla {c['nx']}x{c['ny']} · "
            f"{c['nt']} pasos · {c['segundos_corrida']} s"
        )
        Ln.append(
            f"   riz(d={D_SONDAS[0]/A:g}a) = {c['riz_d_1']:.6f}   "
            f"riz(d={D_SONDAS[1]/A:g}a) = {c['riz_d_2']:.6f}"
        )
        Ln.append(
            f"   reparto A({D_SONDAS[1]/A:g}a)/A({D_SONDAS[0]/A:g}a) = "
            f"{c['reparto_media']:.6f}  (teorico {c['reparto_teorico']:.6f}, "
            f"desv max {c['reparto_desv_max']:.6f})"
        )
        Ln.append(
            f"   fase: error max = {c['fase_error_max']:.6e}  ·  rama m = {c['rama_m']:.0f}"
        )
        Ln.append(
            f"   cola de sonda (t > {FRACCION_COLA:g} T) = {c['cola_sonda']:.3e}  ·  "
            f"cola de dominio = {c['cola_dominio']:.3e}"
        )
    Ln.append("")
    for modo in MODOS:
        Ln.append(
            f"delta({modo}) = max_nu |ln A(36a) - ln A(24a)| = {deltas[modo]:.6f}"
        )
    Ln.append(
        f"acompanante analitico Yee: error de fase max en la ventana = {yee_max:.6e}"
    )
    Ln.append("")
    Ln.append("---- predicciones selladas en el pre-registro (§6.b) ----")
    Ln.append(f"Q1 riz <= {UMBRAL}: {'VERDE' if salida['Q1_riz_ok'] else 'ROJO'}")
    Ln.append(f"Q2 delta <= {UMBRAL}: {'VERDE' if salida['Q2_delta_ok'] else 'ROJO'}")
    Ln.append(
        f"Q3 reparto (2D 0.7071 · 1D 1.0000) +- {UMBRAL}: "
        f"{'VERDE' if salida['Q3_reparto_ok'] else 'ROJO'}"
    )
    Ln.append(
        f"Q4 fase medida {fase_max:.3e} <= {FACTOR_Q4:g} x Yee {yee_max:.3e}"
        f" = {FACTOR_Q4 * yee_max:.3e}: {'VERDE' if salida['Q4_fase_ok'] else 'ROJO'}"
    )
    Ln.append("")
    Ln.append("---- ventana: nu · A_C1a(8a) · A_C1a(16a) · A_C1b(8a) · A_C1b(16a) ----")
    ref = por[("C1a", L_LISTA[0])]
    for i in range(len(ref["nu"])):
        Ln.append(
            f"  {ref['nu'][i]:.4f}  {ref['A_d_1'][i]:.6e}  {ref['A_d_2'][i]:.6e}  "
            f"{por[('C1b', L_LISTA[0])]['A_d_1'][i]:.6e}  "
            f"{por[('C1b', L_LISTA[0])]['A_d_2'][i]:.6e}"
        )
    texto = "\n".join(Ln) + "\n"
    pt = CRUDOS / f"crudo_C1_{etiqueta}.txt"
    pt.write_text(texto, encoding="utf-8")
    print(texto, end="")
    if "--json" in sys.argv:
        print(json.dumps(salida, indent=1, ensure_ascii=False))
    print(f"crudo json: {pj} · sha16 {salida['sha16_crudo_json']}")
    print(f"crudo txt : {pt} · sha16 {sha16(pt)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
