#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""inyector_diel.py -- ACTO 3 del pre-registro de fotónica (monotilo).

# HUELLA_UNICA: inyector_diel_acto3_20261012

Este módulo NO mide transmitancia por sí solo: construye y VERIFICA el campo
dieléctrico (`eps`) del acto 3 según lo sellado en
`PREREGISTRO_INSTRUMENTO_FDTD_20261012.md` (§9 · §9.8) y, sólo con la bandera
`--correr`, lanza las corridas FDTD por el hook `eps` del solver.

ACTO 3.0 = construir + verificar (por defecto). No lanza ninguna corrida FDTD
           y no escribe en disco: re-mide y coteja cifras YA selladas.
ACTO 3.1 = hook `eps` en `fdtd2d_tm.correr` (acto propio, no vive aquí: §9.7).
ACTO 3.2 = `--correr`: U1 + 5 semillas y el criterio sellado de §9.4.

UNIVERSO DECLARADO (se imprime antes de cualquier cifra):
  geometría heredada de C1 (dx = a/16, PML 20 celdas, fuente i = 32, sondas
  8a/16a, T_total = 72 a/c, Nt = 2304, ventana 0.15..0.55 c/a = 29 bins)
  objeto = parche de nivel 2 rasterizado a 24a con 1 celda de máscara =
  1 celda de Yee; lado 175.803156968241 u; celda 0.457820721 u; 3603 teselas
  eps_mat = 4.0 · eps_0 = 1.0 · L = 24 a (espesor) · 5 semillas 101..105
  criterio: |T_U1 - Tbar_U2| > 3 sigma_U2 + riz(T)_C1, con riz(T)_C1 = 0.05
  (§9.5 lectura (a), la mayor) en al menos 3 bins espectrales consecutivos
  unidad de conteo: celda (1 celda = a/16 = DX)

Convención de ejes (DECLARADA, no elegida después de ver la cifra): la máscara
de `malla_hat.rasterizar` tiene ejes (i_y, j_x) por `np.meshgrid(xs, ys)`; el
solver usa `eps(i_x, j_y)` ⇒ se TRANSPONE (`union.T`) para que el eje x de la
máscara sea el eje x del dominio (dirección de propagación). Sin la
transposición la muestra queda reflejada respecto de y = x: otro objeto.

Lo que este módulo NO hace (§9.7): no construye el nivel 3, ni U3 (hexagonal
periódica), ni C2 (Penrose); no toca el repo local de sfqed-paper; no instala ni
deposita tarballs; no toca el registro de Zenodo. No introduce ningún parámetro
nuevo: todo lo que usa está sellado en §9/§9.8 y se cita con su llave.
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
import sys
from pathlib import Path

import numpy as np

import fdtd2d_tm
import malla_hat

AQUI = Path(__file__).resolve().parent
CRUDOS = AQUI / "crudos"
CRUDO_REF = CRUDOS / "crudo_C1_20261012T062648Z.json"
REF_CORRIDA_DECLARADA = 2

# ---------------- parámetros sellados (§9.3 · §9.6) ----------------
NIVEL = 2
N = 384
L_A = 24.0
MODO = "C1b"
EPS_MAT = 4.0
EPS_0 = 1.0
BARRIDO = 24
SEMILLAS = (101, 102, 103, 104, 105)

# ---------------- cifras selladas que este módulo RE-MIDE (§9.8) ----------------
MALLA_SHA16_CERTIFICADO = "914a8666f1d13ddf"
FDTD_SHA16_DECLARADO = "ef29596d00b873b5"
FDTD_BLOB_DECLARADO = "0d7820a8ea566e403e9528f30945a4a04ff9cbc4"
LADO_DECLARADO = 175.803156968241
N_TESELAS_DECLARADO = 3603
TOL_LLENADO_DECLARADA = 8.077e-03
FILL_SELLADA = 0.403832758842046
COBERTURA_U1_DECLARADA = 0.403679741753472
SOLAPE_U1_DECLARADO = 2.0345052083e-05
ARISTA_MIN_DECLARADA_U = 0.5
ARISTA_MIN_DECLARADA_CELDAS = 1.0921
ARISTA_MEDIA_DECLARADA_CELDAS = 1.5451
ARISTA_MAX_DECLARADA_CELDAS = 2.1843
ANCHO_MIN_DECLARADO_U = 2.598076211
ANCHO_MIN_DECLARADO_CELDAS = 5.6749
ANCHO_MEDIA_DECLARADA_CELDAS = 5.9680
ANCHO_MAX_DECLARADA_CELDAS = 6.5528
COBERTURA_U2_DECLARADA = (
    0.382473416,
    0.385084364,
    0.383477105,
    0.384636773,
    0.382995605,
)
CELDAS_DIF_DECLARADAS = (69807, 69700, 69985, 70190, 70026)
D_ELEGIDO_DECLARADO = 3.045985824
DESCARTES_DECLARADOS = 225
SOLAPE_U2_DECLARADO = 0.016994900

# ---------------- criterio sellado (§9.4 · §9.5) ----------------
K_SIGMA = 3.0
RIZ_T_C1 = 0.05
BINS_CONSECUTIVOS = 3


# ------------------------------------------------------------------------------
# utilidades
# ------------------------------------------------------------------------------

def sha16_de(datos: bytes) -> str:
    return hashlib.sha256(datos).hexdigest()[:16]


def sha16_de_archivo(p: Path) -> str:
    return sha16_de(p.read_bytes())


def t_utc() -> str:
    import time

    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def head_git() -> str:
    try:
        h = (AQUI / ".git" / "HEAD").read_text(encoding="utf-8").strip()
        if h.startswith("ref: "):
            ref = h[5:].strip()
            p = AQUI / ".git" / ref
            if p.exists():
                return p.read_text(encoding="utf-8").strip()[:7]
            pk = AQUI / ".git" / "packed-refs"
            if pk.exists():
                for ln in pk.read_text(encoding="utf-8").splitlines():
                    if ln.endswith(ref):
                        return ln.split()[0][:7]
            return "ref-sin-archivo:" + ref
        return h[:7]
    except OSError:
        return "desconocido"


def linea_universo() -> str:
    dx = float(fdtd2d_tm.DX)
    pml = int(fdtd2d_tm.PML_CELDAS)
    return (
        "UNIVERSO DECLARADO: geometría de C1 (dx = a/16 = "
        f"{dx:g}, PML {pml} celdas por borde, fuente i = {int(fdtd2d_tm.I_FUENTE)}, "
        f"sondas {fdtd2d_tm.D_SONDAS[0] / fdtd2d_tm.A:g}a/{fdtd2d_tm.D_SONDAS[1] / fdtd2d_tm.A:g}a, "
        f"T_total = {fdtd2d_tm.T_TOTAL:g} a/c, Nt = {int(round(fdtd2d_tm.T_TOTAL / fdtd2d_tm.DT))}, "
        f"ventana {fdtd2d_tm.NU_MIN:g}..{fdtd2d_tm.NU_MAX:g} c/a) · objeto = parche de nivel "
        f"{NIVEL} rasterizado a {L_A:g}a con 1 celda de máscara = 1 celda de Yee · "
        f"N = {N} celdas · modo {MODO} · eps_mat {EPS_MAT:g} · eps_0 {EPS_0:g} · "
        f"L = {L_A:g} a (espesor) · semillas {SEMILLAS} · "
        f"criterio 3 sigma + riz(T)_C1 = {RIZ_T_C1:g} en >= {BINS_CONSECUTIVOS} bins "
        "consecutivos · unidad de conteo: celda"
    )


class Fallos:
    """Acumula los cotejos de la corrida. Fail-closed: si hay ROJO, RC != 0."""

    def __init__(self) -> None:
        self.lista: list[str] = []

    def chk(self, nombre: str, ok: bool, detalle: str = "") -> bool:
        if not ok:
            self.lista.append(nombre)
        etiqueta = "OK" if ok else "ROJO"
        cola = f" · {detalle}" if detalle else ""
        print(f"   [{etiqueta}] {nombre}{cola}")
        return ok

    def cerca(self, nombre: str, medido: float, declarado: float, rtol: float) -> bool:
        ok = math.isclose(medido, declarado, rel_tol=rtol, abs_tol=0.0)
        return self.chk(
            nombre,
            ok,
            f"medido {medido:.12g} · declarado {declarado:.12g} · rtol {rtol:g}",
        )


def rachas(mask: list[bool]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    ini = None
    for i, v in enumerate(mask):
        if v and ini is None:
            ini = i
        if not v and ini is not None:
            out.append((ini, i - 1))
            ini = None
    if ini is not None:
        out.append((ini, len(mask) - 1))
    return out


def aristas_y_anchos(poligonos: list) -> tuple[list[float], list[float]]:
    """Unidad de salida: unidad u (la misma del polígono)."""
    aristas: list[float] = []
    anchos: list[float] = []
    for poly in poligonos:
        m = len(poly)
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
        anchos.append(max(max(xs) - min(xs), max(ys) - min(ys)))
        for k in range(m):
            xa, ya = poly[k]
            xb, yb = poly[(k + 1) % m]
            aristas.append(math.hypot(xb - xa, yb - ya))
    return aristas, anchos


# ------------------------------------------------------------------------------
# ACTO 3.0 · construir y verificar la máscara y el campo eps
# ------------------------------------------------------------------------------

def construir(verbose: bool = True) -> dict:
    f = Fallos()
    dx = float(fdtd2d_tm.DX)
    pml = int(fdtd2d_tm.PML_CELDAS)
    nx = int(round(L_A / dx))
    ny = nx if MODO == "C1b" else 1

    print("---- ACTO 3.0 · construir y verificar (cero corridas FDTD) ----")
    print(
        f"T_UTC {t_utc()} · git HEAD {head_git()} · malla_hat sha16 "
        f"{sha16_de_archivo(AQUI / 'malla_hat.py')} · fdtd2d_tm sha16 "
        f"{sha16_de_archivo(AQUI / 'fdtd2d_tm.py')}"
    )
    print(linea_universo())
    print("")

    # --- procedencia y geometría del dominio (hereda de C1, sin grados nuevos)
    print("-- procedencia y dominio (unidad: celda; longitud: u) --")
    f.chk(
        "malla_hat.py es el motor de geometría certificado",
        sha16_de_archivo(AQUI / "malla_hat.py") == MALLA_SHA16_CERTIFICADO,
        f"sha16 {sha16_de_archivo(AQUI / 'malla_hat.py')} · certificado {MALLA_SHA16_CERTIFICADO}",
    )
    mismo_fdtd = sha16_de_archivo(AQUI / "fdtd2d_tm.py") == FDTD_SHA16_DECLARADO
    if not mismo_fdtd:
        print(
            "   [nota] fdtd2d_tm.py NO es el blob certificado "
            f"({FDTD_SHA16_DECLARADO} / {FDTD_BLOB_DECLARADO}): esperado si el ACTO 3.1 "
            "ya añadió el hook eps. La referencia sellada se CITA, no se re-mide."
        )
    f.chk("N = nx = ny (1 celda de máscara = 1 celda de Yee)", nx == N and ny == N,
          f"nx {nx} · ny {ny} · N {N}")
    f.chk("PML por borde = 20 celdas", pml == 20, f"pml {pml}")

    # --- U1 (la geometría sellada del acto 2)
    print("")
    print("-- U1: parche de nivel 2, ventana llevada a 24a (unidad: celda / u) --")
    cif1, teselas1, polis1, ven1 = malla_hat.u1(NIVEL, N)
    x0, x1, y0, y1 = ven1
    lado = x1 - x0
    celda_u = lado / N
    f.cerca("lado de la ventana = 175.803156968241 u", lado, LADO_DECLARADO, 1e-9)
    f.cerca("celda = lado / N (a/16 = 0.0625 en unidades de malla)",
            celda_u, LADO_DECLARADO / N, 1e-12)
    f.chk("n_teselas de U1 = 3603", int(cif1["n_teselas"]) == N_TESELAS_DECLARADO,
          f"n_teselas {int(cif1['n_teselas'])} · declarado {N_TESELAS_DECLARADO}")
    f.cerca("cobertura de U1 = 0.403679741753472", float(cif1["cobertura"]),
            COBERTURA_U1_DECLARADA, 1e-9)
    f.cerca("solape de U1 (2.03e-05)", float(cif1["solape"]),
            SOLAPE_U1_DECLARADO, 1e-2)
    f.cerca("suma_areas de U1 (fracción de llenado sellada 0.403833)",
            float(cif1["suma_areas"]), FILL_SELLADA, 1e-9)

    union1, conteo1 = malla_hat.rasterizar(polis1, *ven1, N)
    f.chk("máscara de U1 = máscara medida por medir() (sha16)",
          sha16_de(union1.tobytes()) == cif1["sha16_mascara"],
          f"sha16 {sha16_de(union1.tobytes())}")
    print(
        f"   celdas cubiertas U1 = {int(union1.sum())} de {N * N}"
        f" · celdas con solape = {int((conteo1 >= 2).sum())}"
        f" · sha16_mascara {cif1['sha16_mascara']}"
    )

    # --- resolución realizada (declaración de §9.8.3, re-medida aquí)
    print("")
    print("-- resolución realizada (re-mide §9.8.3; divergencia se imprime, "
          "no aborta: no toca el criterio) --")
    aristas_u, anchos_u = aristas_y_anchos(polis1)
    a_min = min(aristas_u)
    a_med = sum(aristas_u) / len(aristas_u)
    a_max = max(aristas_u)
    f.cerca("arista min = 0.500000000 u", a_min, ARISTA_MIN_DECLARADA_U, 1e-9)
    print(
        f"   aristas (u): min {a_min:.9f} · media {a_med:.9f} · max {a_max:.9f}"
        f" · n = {len(aristas_u)}"
    )
    print(
        f"   aristas (celdas): min {a_min / celda_u:.4f} (declarado "
        f"{ARISTA_MIN_DECLARADA_CELDAS:g}) · media {a_med / celda_u:.4f} (declarado "
        f"{ARISTA_MEDIA_DECLARADA_CELDAS:g}) · max {a_max / celda_u:.4f} (declarado "
        f"{ARISTA_MAX_DECLARADA_CELDAS:g})"
    )
    an_min = min(anchos_u)
    an_med = sum(anchos_u) / len(anchos_u)
    an_max = max(anchos_u)
    f.cerca("ancho min = 2.598076211 u", an_min, ANCHO_MIN_DECLARADO_U, 1e-9)
    print(
        f"   anchos (celdas): min {an_min / celda_u:.4f} (declarado "
        f"{ANCHO_MIN_DECLARADO_CELDAS:g}) · media {an_med / celda_u:.4f} (declarado "
        f"{ANCHO_MEDIA_DECLARADA_CELDAS:g}) · max {an_max / celda_u:.4f} (declarado "
        f"{ANCHO_MAX_DECLARADA_CELDAS:g})"
    )
    print(
        f"   nota de atribución (§9.8.3): un POSITIVO es seguro; un NEGATIVO no es "
        f"atribuible del todo a la geometría (arista min {a_min / celda_u:.4f} celdas)"
    )

    # --- U2: el ensemble del nulo (5 semillas, §9.2)
    print("")
    print("-- U2: ensemble de 5 semillas 101..105 (§9.2) (unidad: celda) --")
    semillas: list[dict] = []
    for k, s in enumerate(SEMILLAS):
        cif2, polis2 = malla_hat.u2(polis1, teselas1, ven1, N,
                                    float(cif1["cobertura"]), s, BARRIDO)
        union2, conteo2 = malla_hat.rasterizar(polis2, *ven1, N)
        dif = int((union1 != union2).sum())
        inter = int((union1 & union2).sum())
        union_c = int((union1 | union2).sum())
        jac = inter / union_c if union_c else 0.0
        fill = float(cif2["fraccion_de_llenado"])
        d_fill = abs(fill - float(cif1["suma_areas"]))
        print(f"   semilla {s}")
        # u2 NO devuelve la cifra "n_teselas": su cif sale de medir() y el conjunto
        # de teselas se mide con len(polis2). Se cuenta el objeto, no se lee una llave.
        n_teselas2 = len(polis2)
        f.chk(f"   n_teselas = 3603 (semilla {s}, contado en len(polis2))",
              n_teselas2 == N_TESELAS_DECLARADO,
              f"n_teselas {n_teselas2}")
        f.chk(f"   congruencia: n_vertices_igual (semilla {s})",
              bool(cif2["n_vertices_igual"]), f"{cif2['n_vertices_igual']}")
        f.chk(f"   |Delta f| <= 8.077e-03 (semilla {s})",
              d_fill <= TOL_LLENADO_DECLARADA,
              f"Delta f {d_fill:.3e} · tope {TOL_LLENADO_DECLARADA:.3e}")
        print(
            f"   [nota] Delta f es IDENTIDAD por construcción (§9.8.2: teselas "
            f"congruentes, rotación = isometría) => NO es evidencia; margen medido "
            f"{TOL_LLENADO_DECLARADA / d_fill if d_fill else float('inf'):.3e}x"
        )
        f.chk(f"   la semilla mueve la máscara (celdas distintas > 0) (semilla {s})",
              dif > 0, f"celdas distintas {dif} · declarado {CELDAS_DIF_DECLARADAS[k]}")
        semillas.append({
            "semilla": int(s),
            "n_teselas": n_teselas2,
            "sha16_mascara": sha16_de(union2.tobytes()),
            "cobertura": float(cif2["cobertura"]),
            "solape": float(cif2["solape"]),
            "fraccion_de_llenado": fill,
            "delta_f": d_fill,
            "delta_cobertura_vs_U1": float(cif2["delta_cobertura_vs_U1"]),
            "celdas_distintas_vs_U1": dif,
            "jaccard_vs_U1": jac,
            "d_elegido": float(cif2["d_elegido"]),
            "descartes_de_reticulo": int(cif2["descartes_de_reticulo"]),
            "areas_max_dif_vs_U1": float(cif2["areas_max_dif_vs_U1"]),
            "n_puntos_reticulo": int(cif2["n_puntos_reticulo"]),
            "union": union2,
            "cif": cif2,
        })

    print("")
    print("-- ensemble: cobertura de la unión y reparto espacial (unidad: celda) --")
    f.chk("las 5 máscaras de U2 son distintas entre sí (sha16)",
          len({x["sha16_mascara"] for x in semillas}) == len(SEMILLAS),
          f"distintas {len({x['sha16_mascara'] for x in semillas})} de {len(SEMILLAS)}")
    cob2 = [x["cobertura"] for x in semillas]
    cobs = " · ".join(f"{c:.9f}" for c in cob2)
    print(f"   cobertura de la unión U2 (declaradas {COBERTURA_U2_DECLARADA}): {cobs}")
    print(f"   cobertura de la unión U1 = {float(cif1['cobertura']):.9f}")
    d = [abs(float(cif1["cobertura"]) - c) for c in cob2]
    print(
        "   delta cobertura U1 - U2 (unidad: fracción): "
        + " · ".join(f"{x:.4e}" for x in d)
        + f"  => relativo {max(d) / float(cif1['cobertura']) * 100:.2f} %"
    )
    print(
        "   [nota §9.8.4] en la rejilla U1 y U2 NO cargan la misma cantidad de "
        f"dieléctrico (solape U1 {float(cif1['solape']):.3e} vs solape U2 por semilla "
        + " · ".join(f"{x['semilla']}:{x['solape']:.6f}" for x in semillas)
        + "): toda lectura positiva se reporta con las DOS lecturas nombradas"
    )
    print(
        f"   [nota de procedencia] el valor sellado {SOLAPE_U2_DECLARADO:g} es de la "
        "corrida de calibración de malla_hat (semilla 20261012, otro ensemble): NO es "
        "la expectativa de estas 5 semillas y NO gatea ningún cotejo (se imprime para "
        "que la divergencia quede a la vista, no para taparla)"
    )
    print(
        "   celdas distintas vs U1 (declaradas "
        f"{CELDAS_DIF_DECLARADAS}): "
        + " · ".join(f"{x['celdas_distintas_vs_U1']}" for x in semillas)
        + f" de {N * N} · Jaccard "
        + " · ".join(f"{x['jaccard_vs_U1']:.6f}" for x in semillas)
    )
    print(
        "   d_elegido (u): " + " · ".join(f"{x['d_elegido']:.9f}" for x in semillas)
        + f" · declarado {D_ELEGIDO_DECLARADO:g} · descartes de retículo "
        + " · ".join(f"{x['descartes_de_reticulo']}" for x in semillas)
        + f" · declarado {DESCARTES_DECLARADOS}"
    )

    # --- campo eps (ACTO 3.0: se construye y se verifica, no se corre)
    print("")
    print("-- campo eps (convención de ejes DECLARADA: eps(i_x, j_y) = union.T) --")
    mascara = np.ascontiguousarray(union1.T)
    f.chk("máscara transpuesta: forma (nx, ny)", mascara.shape == (nx, ny),
          f"forma {mascara.shape} · (nx, ny) = {(nx, ny)}")
    eps = EPS_0 + (EPS_MAT - EPS_0) * mascara.astype(np.float64)
    interior = np.zeros_like(eps, dtype=bool)
    interior[pml:nx - pml, pml:ny - pml] = True
    eps[:pml, :] = EPS_0
    eps[-pml:, :] = EPS_0
    eps[:, :pml] = EPS_0
    eps[:, -pml:] = EPS_0
    f.chk("borde de PML con eps_0 = 1.0 (absorbente certificado)",
          bool(np.all(eps[~interior] == EPS_0)),
          f"celdas de borde {int((~interior).sum())}")
    f.chk("interior sólo toma eps_0 o eps_mat",
          bool(np.all((eps[interior] == EPS_0) | (eps[interior] == EPS_MAT))),
          f"valores distintos {sorted(set(np.unique(eps).tolist()))}")
    frac_int = float(mascara[interior].mean())
    print(
        f"   fracción de máscara en el interior = {frac_int:.9f}"
        f" · cobertura de la unión (rejilla completa) = {float(cif1['cobertura']):.9f}"
        f" · celdas en eps_mat = {int((eps == EPS_MAT).sum())}"
    )
    print(
        f"   sha16 máscara (i_y, j_x) {sha16_de(union1.tobytes())} · "
        f"sha16 máscara transpuesta {sha16_de(mascara.tobytes())} · "
        f"sha16 eps (bytes float64) {sha16_de(np.ascontiguousarray(eps).tobytes())}"
    )

    # --- criterio sellado, recordado tal cual (sin tocarlo)
    print("")
    print("-- criterio sellado (§9.4 · §9.5), recordado sin tocar --")
    print(f"   |T_U1(nu) - Tbar_U2(nu)| > {K_SIGMA:g} sigma_U2(nu) + riz(T)_C1, "
          f"riz(T)_C1 = {RIZ_T_C1:g}, en >= {BINS_CONSECUTIVOS} bins consecutivos")
    print(
        f"   lectura descartada de riz(T)_C1 (§9.5 (b), registrada): "
        f"reparto_riz C1b = 0.00552322568652841 => {RIZ_T_C1 / 0.00552322568652841:.4f}x "
        f"más fácil que la adoptada ({RIZ_T_C1:g})"
    )
    print(
        "   referencia sellada (se CITA, no se re-mide): "
        f"{CRUDO_REF.name} · corridas[{REF_CORRIDA_DECLARADA}] · A_d_2 (sonda 16a)"
    )

    resumen = {
        "acto": "3.0 (construir y verificar; cero corridas FDTD)",
        "t_utc": t_utc(),
        "git_head": head_git(),
        "universo": linea_universo(),
        "malla_hat_sha16": sha16_de_archivo(AQUI / "malla_hat.py"),
        "fdtd2d_tm_sha16": sha16_de_archivo(AQUI / "fdtd2d_tm.py"),
        "nivel": NIVEL,
        "N": N,
        "nx": nx,
        "ny": ny,
        "pml_celdas": pml,
        "lado_u": lado,
        "celda_u": celda_u,
        "modo": MODO,
        "eps_mat": EPS_MAT,
        "eps_0": EPS_0,
        "L_a": L_A,
        "barrido": BARRIDO,
        "semillas": list(SEMILLAS),
        "criterio": {
            "k_sigma": K_SIGMA,
            "riz_T_C1": RIZ_T_C1,
            "riz_T_C1_descartada": 0.00552322568652841,
            "bins_consecutivos": BINS_CONSECUTIVOS,
            "crudo_referencia": CRUDO_REF.name,
            "corrida_referencia": REF_CORRIDA_DECLARADA,
        },
        "u1": {
            "n_teselas": int(cif1["n_teselas"]),
            "cobertura": float(cif1["cobertura"]),
            "solape": float(cif1["solape"]),
            "suma_areas": float(cif1["suma_areas"]),
            "sha16_mascara": sha16_de(union1.tobytes()),
            "celdas_cubiertas": int(union1.sum()),
            "arista_min_u": a_min,
            "arista_media_u": a_med,
            "arista_max_u": a_max,
            "arista_min_celdas": a_min / celda_u,
            "arista_media_celdas": a_med / celda_u,
            "arista_max_celdas": a_max / celda_u,
            "ancho_min_u": an_min,
            "ancho_media_u": an_med,
            "ancho_max_u": an_max,
        },
        "u2": [
            {k: v for k, v in x.items() if k not in ("union", "cif")} for x in semillas
        ],
        "eps": {
            "sha16_mascara_transpuesta": sha16_de(mascara.tobytes()),
            "sha16_eps_bytes": sha16_de(np.ascontiguousarray(eps).tobytes()),
            "celdas_en_eps_mat": int((eps == EPS_MAT).sum()),
            "fraccion_mascara_interior": frac_int,
            "borde_pml_con_eps_0": True,
        },
        "cotejos_rojos": list(f.lista),
        "rc": 2 if f.lista else 0,
    }
    print("")
    rojos = ", ".join(f.lista) if f.lista else "ninguno"
    print(f"-- ACTO 3.0: cotejos ROJOS: {rojos} --")
    return resumen


# ------------------------------------------------------------------------------
# ACTO 3.2 · el hook eps y las corridas (sólo con --correr)
# ------------------------------------------------------------------------------

MENSAJE_SIN_HOOK = (
    "ERROR: fdtd2d_tm.correr no acepta un campo eps externo (ACTO 3.1 pendiente).\n"
    "       El solver certificado construye eps = np.ones((nx, ny)) dentro de\n"
    "       correr() (§9.6): sin el hook, cualquier 'inyección' sería un parche.\n"
    "       Orden para el ACTO 3.1 (acto propio, §9.7): pre-vuelo A2.1 con\n"
    f"       BLOB_PRE = {FDTD_BLOB_DECLARADO} · LINEAS_PRE = 472 · huella = la línea\n"
    "       nueva y distintiva; después, re-corrida de C1 que reproduzca Q1..Q4\n"
    "       VERDE y las amplitudes de gen 6.\n"
    "       Este módulo NO parchea el instrumento certificado."
)


def sonda_hook():
    """Devuelve el nombre del parámetro del hook. Fail-closed si no existe."""
    firma = inspect.signature(fdtd2d_tm.correr)
    if "eps" not in firma.parameters:
        raise SystemExit(MENSAJE_SIN_HOOK)
    p = firma.parameters["eps"]
    if p.default is inspect.Parameter.empty:
        raise SystemExit(
            MENSAJE_SIN_HOOK + "       (el parámetro 'eps' existe pero NO tiene default:\n"
            "       rompería la ruta certificada.)"
        )
    if p.kind not in (
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        inspect.Parameter.KEYWORD_ONLY,
    ):
        raise SystemExit(MENSAJE_SIN_HOOK + "       (el parámetro 'eps' no es nombrable.)")
    return "eps"


def campo_eps(union: np.ndarray) -> np.ndarray:
    nx = int(round(L_A / float(fdtd2d_tm.DX)))
    ny = nx if MODO == "C1b" else 1
    if union.shape != (nx, ny):
        raise SystemExit(
            f"ERROR: la máscara tiene forma {union.shape} y el dominio pide {(nx, ny)}: "
            "1 celda de máscara = 1 celda de Yee exige igualdad exacta de forma."
        )
    eps = EPS_0 + (EPS_MAT - EPS_0) * np.ascontiguousarray(union.T).astype(np.float64)
    pml = int(fdtd2d_tm.PML_CELDAS)
    eps[:pml, :] = EPS_0
    eps[-pml:, :] = EPS_0
    eps[:, :pml] = EPS_0
    eps[:, -pml:] = EPS_0
    return eps


def referencia_sellada() -> tuple[np.ndarray, np.ndarray, dict]:
    if not CRUDO_REF.exists():
        raise SystemExit(f"ERROR: falta la referencia sellada {CRUDO_REF} (§9.4).")
    raw = json.loads(CRUDO_REF.read_text(encoding="utf-8"))
    corridas = raw["corridas"]
    if not (0 <= REF_CORRIDA_DECLARADA < len(corridas)):
        raise SystemExit(
            f"ERROR: la referencia declara corridas[{REF_CORRIDA_DECLARADA}] y el crudo "
            f"tiene {len(corridas)} corridas."
        )
    c = corridas[REF_CORRIDA_DECLARADA]
    if c.get("modo") != MODO or abs(float(c.get("L_a", -1)) - L_A) > 1e-12:
        raise SystemExit(
            f"ERROR: corridas[{REF_CORRIDA_DECLARADA}] es modo {c.get('modo')} "
            f"L = {c.get('L_a')} y §9.4 fija {MODO} con L = {L_A:g}a: el índice "
            "declarado NO es la referencia sellada."
        )
    if "A_d_2" not in c or "nu" not in c:
        raise SystemExit("ERROR: la referencia sellada no trae A_d_2 / nu.")
    return (
        np.asarray(c["nu"], dtype=np.float64),
        np.asarray(c["A_d_2"], dtype=np.float64),
        c,
    )


def correr_acto3(escribir_crudo: bool = True) -> int:
    nombre_param = sonda_hook()
    print("---- ACTO 3.2 · corridas FDTD por el hook eps ----")
    print(f"T_UTC {t_utc()} · git HEAD {head_git()}")
    print(linea_universo())
    print("")
    nu_ref, a_ref, c_ref = referencia_sellada()
    print(
        f"referencia: {CRUDO_REF.name} · corridas[{REF_CORRIDA_DECLARADA}] · "
        f"modo {c_ref['modo']} · L_a {c_ref['L_a']:g} · "
        f"riz_d_2 {c_ref.get('riz_d_2')} · {len(nu_ref)} bins"
    )

    res_u1 = construir(verbose=False)
    if res_u1["rc"] != 0:
        print("ACTO 3.0 con ROJOS: no se corre nada (fail-closed).")
        return 2

    cif1, teselas1, polis1, ven1 = malla_hat.u1(NIVEL, N)
    union1, _ = malla_hat.rasterizar(polis1, *ven1, N)
    eps_u1 = campo_eps(union1)
    print("")
    print("-- U1 --")
    r_u1 = fdtd2d_tm.correr(**{"L_a": L_A, "modo": MODO, nombre_param: eps_u1})
    nu = np.asarray(r_u1["nu"], dtype=np.float64)
    if nu.shape != nu_ref.shape or not np.allclose(nu, nu_ref):
        raise SystemExit(
            "ERROR: la ventana de bins de la corrida no coincide con la referencia "
            "sellada: el cociente T no está definido en la misma rejilla."
        )
    t_u1 = np.asarray(r_u1["A_d_2"], dtype=np.float64) / a_ref
    print(f"   A_d_2 (16a) · T = A/U_ref · riz_d_2 {r_u1['riz_d_2']:.6f}")

    t_sem: list[np.ndarray] = []
    for s in SEMILLAS:
        cif2, polis2 = malla_hat.u2(polis1, teselas1, ven1, N,
                                    float(cif1["cobertura"]), s, BARRIDO)
        union2, _ = malla_hat.rasterizar(polis2, *ven1, N)
        eps_s = campo_eps(union2)
        print(f"-- semilla {s} --")
        r = fdtd2d_tm.correr(**{"L_a": L_A, "modo": MODO, nombre_param: eps_s})
        if not np.allclose(np.asarray(r["nu"], dtype=np.float64), nu_ref):
            raise SystemExit(f"ERROR: ventana de bins distinta en la semilla {s}.")
        t_sem.append(np.asarray(r["A_d_2"], dtype=np.float64) / a_ref)

    mat = np.vstack(t_sem)
    t_bar = mat.mean(axis=0)
    sigma = mat.std(axis=0, ddof=1)
    dif = np.abs(t_u1 - t_bar)
    banda = K_SIGMA * sigma + RIZ_T_C1
    cumple = (dif > banda).tolist()
    rr = rachas(cumple)
    larga = max((b - a + 1) for a, b in rr) if rr else 0
    positivo = larga >= BINS_CONSECUTIVOS

    print("")
    print("---- tabla cruda (unidad: fracción; nu en c/a) ----")
    cab = ("     nu     T_U1" + "".join(f"    T_s{s}" for s in SEMILLAS)
           + "    Tbar_U2     sigma   |dif|     banda  cumple")
    print(cab)
    lineas = [cab]
    for i in range(len(nu)):
        ln = (
            f"  {nu[i]:.4f}  {t_u1[i]:.6f}"
            + "".join(f"  {mat[k][i]:.6f}" for k in range(len(SEMILLAS)))
            + f"  {t_bar[i]:.6f}  {sigma[i]:.6f}  {dif[i]:.6f}  {banda[i]:.6f}  "
            + ("SI" if cumple[i] else "no")
        )
        print(ln)
        lineas.append(ln)

    print("")
    print(f"rachas que cumplen: {rr if rr else 'ninguna'} · racha más larga {larga} bin(s)")
    if rr:
        for a, b in rr:
            print(f"   racha de {b - a + 1} bin(s): nu {nu[a]:.4f} .. {nu[b]:.4f} c/a")
    if positivo:
        veredicto = (
            f"VEREDICTO POSITIVO: la desigualdad sellada se cumple en {larga} bins "
            f"consecutivos (>= {BINS_CONSECUTIVOS}). La correlación espacial del "
            "monotilo se distingue en transmisión a este contraste."
        )
    else:
        veredicto = (
            "VEREDICTO NEGATIVO / NULO SOSTENIDO: en toda la ventana vale "
            "|T_U1 - Tbar_U2| <= 3 sigma_U2 + riz(T)_C1 => la hipótesis (que la "
            "correlación espacial del monotilo se distinga en transmisión a este "
            "contraste) queda formalmente refutada en este régimen de contraste. "
            "Se publica como resultado, no como fracaso (§9.4)."
        )
    print(veredicto)

    salida = {
        "acto": "3.2",
        "t_utc": t_utc(),
        "git_head": head_git(),
        "universo": linea_universo(),
        "malla_hat_sha16": sha16_de_archivo(AQUI / "malla_hat.py"),
        "fdtd2d_tm_sha16": sha16_de_archivo(AQUI / "fdtd2d_tm.py"),
        "referencia": {
            "crudo": CRUDO_REF.name,
            "corrida": REF_CORRIDA_DECLARADA,
            "modo": c_ref["modo"],
            "L_a": c_ref["L_a"],
        },
        "nu": nu.tolist(),
        "A_d_2_u1": np.asarray(r_u1["A_d_2"], dtype=np.float64).tolist(),
        "A_d_2_semillas": [x.tolist() for x in
                           [np.asarray(t, dtype=np.float64) for t in mat]],
        "T_u1": t_u1.tolist(),
        "T_semillas": [x.tolist() for x in mat],
        "T_bar_u2": t_bar.tolist(),
        "sigma_u2": sigma.tolist(),
        "dif": dif.tolist(),
        "banda": banda.tolist(),
        "cumple": cumple,
        "rachas": rr,
        "racha_mas_larga": int(larga),
        "criterio": {
            "k_sigma": K_SIGMA,
            "riz_T_C1": RIZ_T_C1,
            "bins_consecutivos": BINS_CONSECUTIVOS,
        },
        "veredicto": "POSITIVO" if positivo else "NEGATIVO / NULO SOSTENIDO",
    }
    salida["rc"] = 0

    if escribir_crudo:
        etiqueta = t_utc().replace("-", "").replace(":", "")
        pj = CRUDOS / f"crudo_acto3_{etiqueta}.json"
        textos = salida.copy()
        pj.write_text(json.dumps(textos, indent=1, ensure_ascii=False), encoding="utf-8")
        salida["sha16_crudo_json"] = sha16_de_archivo(pj)
        pt = CRUDOS / f"crudo_acto3_{etiqueta}.txt"
        cuerpo = "\n".join(
            [linea_universo(), "", lineas[0], *lineas[1:], "", veredicto, ""]
        )
        pt.write_text(cuerpo, encoding="utf-8")
        print(f"crudo json: {pj} · sha16 {salida['sha16_crudo_json']}")
        print(f"crudo txt : {pt} · sha16 {sha16_de_archivo(pt)}")

    print("")
    print(f"-- ACTO 3.2: veredicto {salida['veredicto']} · rc {salida['rc']} --")
    return salida["rc"]


# ------------------------------------------------------------------------------
# main
# ------------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="ACTO 3 (dieléctrico en la malla certificada). Por defecto: "
        "construye y verifica; no corre FDTD ni escribe en disco."
    )
    ap.add_argument("--correr", action="store_true",
                    help="ACTO 3.2: lanza U1 + las 5 semillas (requiere el hook eps)")
    ap.add_argument("--sin-crudo", action="store_true",
                    help="con --correr: no escribe el crudo en crudos/")
    ap.add_argument("--json", action="store_true", help="imprime el resumen JSON")
    args = ap.parse_args(argv)

    if args.correr:
        rc = correr_acto3(escribir_crudo=not args.sin_crudo)
    else:
        res = construir(verbose=True)
        if args.json:
            print(json.dumps(res, indent=1, ensure_ascii=False))
        rc = res["rc"]
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
