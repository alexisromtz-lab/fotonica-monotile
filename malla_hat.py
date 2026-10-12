# -*- coding: utf-8 -*-
"""malla_hat.py — generador de la red dielectrica (malla) del monotilo Hat: U1 y U2.

LINEA: fotonica-monotile (acto 2, entrada). No depende del paper: no lo importa ni lo toca.

UNIVERSO DECLARADO (se imprime en cada corrida; ninguna cifra sale sin el):
  OBJETO     teselas del monotilo Hat (polikite de 13 vertices, area 8*sqrt(3) unidades^2) y
             su jerarquia de sustitucion (metatiles H/T/P/F, reglas de inflacion de hatviz).
  UNIDAD     tesela (hoja Hat) · vertice (vertice de contorno) · celda (celda de la rejilla
             raster de la ventana). Toda cifra dice su unidad.
  MOTOR      python3 + math (doble precision) para la geometria afin; numpy para rasterizar
             por rejilla regular; random.Random(semilla) para el nulo U2. Cero azar sin semilla.
  VENTANA    UTC de la corrida (hora del reloj, impresa en el crudo).
  FUENTE     reglas, contornos iniciales y matrices de los cuatro metatiles iniciales COPIADOS
             de la salida que la fuente misma imprime al ejecutarse bajo node v24.15.0
             (volcado verbatim), no de transcripcion a mano:
               hat.js       sha16 577deec32bf67d16 · 688 lineas · 15021 bytes
               geometry.js  sha16 a6d4ab9b8baf5513 ·  96 lineas ·  2122 bytes
  ALCANCE    U1 = jerarquia determinista: mismas reglas de sustitucion, jerarquia intacta.
             U2 = nulo decisivo: mismas teselas (mismas formas, mismo numero) y misma fraccion
             de llenado de la ventana, con ORIENTACIONES barajadas y POSICIONES barajadas.
             Declarado: en U2 las teselas pueden solaparse; la estructura es la UNION.
  CERTIFICADO por MATRIZ (aritmetica entera, exacta): h_{k+1} = transpuesta(M) * h_k,
             h_0 = (4,1,2,2); autovalor de Perron esperado phi^4 = 6.854101966249685.
             Los contornos de las cajas NO son certificado: son esquematicos (cortan las
             teselas). La evidencia medida de eso va impresa, no afirmada.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import sys
import time

import numpy as np

# --------------------------------------------------------------------------------------
# 0. Constantes del universo
# --------------------------------------------------------------------------------------

R3 = 1.7320508075688772
HR3 = 0.8660254037844386
PHI4_ESPERADO = 6.854101966249685          # phi^4, factor de inflacion por nivel (area)
IDENT = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0)

FUENTE_HAT = "hat.js"
FUENTE_GEO = "geometry.js"
SHA16_HAT = "577deec32bf67d16"
SHA16_GEO = "a6d4ab9b8baf5513"
VENTANA_FUENTE = "2026-10-12T07:19Z (volcado verbatim de la fuente bajo node)"

UNIVERSO = (
    "teselas del monotilo Hat (13 vertices, area 8*sqrt(3)) + jerarquia de sustitucion "
    "H/T/P/F de hatviz | UNIDAD: tesela / vertice / celda | MOTOR: python3 + math + numpy, "
    "random.Random(semilla) | FUENTE: " + FUENTE_HAT + " sha16 " + SHA16_HAT + " y "
    + FUENTE_GEO + " sha16 " + SHA16_GEO
)

M_MADRE = ((3, 1, 3, 3), (1, 0, 0, 0), (2, 0, 1, 2), (2, 0, 1, 3))  # columnas H,T,P,F
H0_TESELAS = (4, 1, 2, 2)                                           # h_0


def t_utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha16_de(datos: bytes) -> str:
    return hashlib.sha256(datos).hexdigest()[:16]


# --------------------------------------------------------------------------------------
# 1. Capa afin (verbatim de geometry.js, sha16 a6d4ab9b8baf5513)
# --------------------------------------------------------------------------------------

def pt(x, y):
    """geometry.js:5  function pt( x, y ) { return { x : x, y : y }; }"""
    return (float(x), float(y))


def hex_pt(x, y):
    """geometry.js  hexPt(x,y) = pt(x + 0.5*y, hr3*y)"""
    return pt(x + 0.5 * y, HR3 * y)


def inv(T):
    """geometry.js  inv(T): inversa afin (det = T0*T4 - T1*T3)."""
    det = T[0] * T[4] - T[1] * T[3]
    return (T[4] / det, -T[1] / det, (T[1] * T[5] - T[4] * T[2]) / det,
            -T[3] / det, T[0] / det, (T[3] * T[2] - T[0] * T[5]) / det)


def mul(A, B):
    """geometry.js  mul(A,B): composicion afin."""
    return (A[0] * B[0] + A[1] * B[3],
            A[0] * B[1] + A[1] * B[4],
            A[0] * B[2] + A[1] * B[5] + A[2],
            A[3] * B[0] + A[4] * B[3],
            A[3] * B[1] + A[4] * B[4],
            A[3] * B[2] + A[4] * B[5] + A[5])


def padd(P, Q):
    return (P[0] + Q[0], P[1] + Q[1])


def psub(P, Q):
    return (P[0] - Q[0], P[1] - Q[1])


def trot(ang):
    return (math.cos(ang), -math.sin(ang), 0.0, math.sin(ang), math.cos(ang), 0.0)


def ttrans(tx, ty):
    return (1.0, 0.0, float(tx), 0.0, 1.0, float(ty))


def rot_about(p, ang):
    return mul(ttrans(p[0], p[1]), mul(trot(ang), ttrans(-p[0], -p[1])))


def trans_pt(M, P):
    return (M[0] * P[0] + M[1] * P[1] + M[2], M[3] * P[0] + M[4] * P[1] + M[5])


def match_seg(p, q):
    """geometry.js  matchSeg: lleva el segmento (0,1)->(0,0)? forma verbatim de la fuente."""
    return (q[0] - p[0], p[1] - q[1], p[0], q[1] - p[1], q[0] - p[0], p[1])


def match_two(p1, q1, p2, q2):
    """geometry.js  matchTwo(p1,q1,p2,q2) = mul(matchSeg(p2,q2), inv(matchSeg(p1,q1)))."""
    return mul(match_seg(p2, q2), inv(match_seg(p1, q1)))


def intersect(p1, q1, p2, q2):
    """geometry.js  intersect: cruce de dos RECTAS infinitas (p1,q1) y (p2,q2)."""
    d = (q2[1] - p2[1]) * (q1[0] - p1[0]) - (q2[0] - p2[0]) * (q1[1] - p1[1])
    uA = ((q2[0] - p2[0]) * (p1[1] - p2[1]) - (q2[1] - p2[1]) * (p1[0] - p2[0])) / d
    return (p1[0] + uA * (q1[0] - p1[0]), p1[1] + uA * (q1[1] - p1[1]))


def area(poligono):
    """Area por la formula del cordon (valor absoluto). Unidad: unidades^2."""
    s = 0.0
    n = len(poligono)
    for i in range(n):
        x1, y1 = poligono[i]
        x2, y2 = poligono[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2.0


def punto_en_poligono(p, poligono):
    """Pertenencia escalar por paridad (mismo criterio que el raster)."""
    x, y = p
    dentro = False
    n = len(poligono)
    for i in range(n):
        xa, ya = poligono[i]
        xb, yb = poligono[(i + 1) % n]
        if (ya > y) != (yb > y):
            xint = (xb - xa) * (y - ya) / (yb - ya) + xa
            if x < xint:
                dentro = not dentro
    return dentro


# --------------------------------------------------------------------------------------
# 2. Tesela y metatiles (verbatim de hat.js, sha16 577deec32bf67d16)
# --------------------------------------------------------------------------------------

# hat_outline (13 vertices) impreso por la fuente bajo node. Unidad: unidades.
HAT_OUTLINE = (
    (0.0, 0.0),
    (-1.5, -HR3),
    (-1.0, -2.0 * HR3),
    (1.0, -2.0 * HR3),
    (1.5, -HR3),
    (3.0, -2.0 * HR3),
    (4.5, -HR3),
    (4.0, 0.0),
    (3.0, 0.0),
    (3.0, 2.0 * HR3),
    (1.5, 3.0 * HR3),
    (1.0, 2.0 * HR3),
    (0.0, 2.0 * HR3),
)

# Contornos de las cuatro cajas iniciales (locales a cada *_init en la fuente; aqui, datos).
H_OUTLINE = ((0.0, 0.0), (4.0, 0.0), (4.5, HR3), (2.5, 5.0 * HR3), (1.5, 5.0 * HR3), (-0.5, HR3))
T_OUTLINE = ((0.0, 0.0), (3.0, 0.0), (1.5, 3.0 * HR3))
P_OUTLINE = ((0.0, 0.0), (4.0, 0.0), (3.0, 2.0 * HR3), (-1.0, 2.0 * HR3))
F_OUTLINE = ((0.0, 0.0), (3.0, 0.0), (3.5, HR3), (3.0, 2.0 * HR3), (-1.0, 2.0 * HR3))

# Matrices de los hijos de cada metatile inicial: impresas por la fuente bajo node.
# Cada 6-tupla es [T0..T5]; el comentario dice la hoja (H/H1/T/P/F).
H_INIT = (
    (-0.24999999999999997, 0.4330127018922194, 1.0, -0.4330127018922194,
     -0.24999999999999997, 1.7320508075688774),                     # hoja H
    (-0.25, 0.4330127018922193, 4.0, -0.4330127018922193,
     -0.25, 1.7320508075688772),                                    # hoja H
    (-0.25000000000000006, -0.43301270189221935, 2.5, 0.43301270189221935,
     -0.25000000000000006, 2.598076211353315),                      # hoja H
    (-0.25, 0.4330127018922193, 2.5, 0.4330127018922193,
     0.25, 0.8660254037844386),                                     # hoja H1
)
H_INIT_ETIQUETAS = ("H", "H", "H", "H1")
T_INIT = ((0.5, 0.0, 0.5, 0.0, 0.5, HR3),)
T_INIT_ETIQUETAS = ("T",)
P_INIT = ((0.5, 0.0, 1.5, 0.0, 0.5, HR3),
          (0.25, 0.4330127018922193, 0.0, -0.4330127018922193, 0.25, 1.7320508075688772))
P_INIT_ETIQUETAS = ("P", "P")
F_INIT = P_INIT
F_INIT_ETIQUETAS = ("F", "F")


class Hoja:
    """HatTile: una tesela Hat. Unidad: tesela."""
    __slots__ = ("etiqueta", "contorno")

    def __init__(self, etiqueta):
        self.etiqueta = etiqueta
        self.contorno = HAT_OUTLINE


class Meta:
    """MetaTile: contorno (caja) + hijos (cada uno: transformada + geometria)."""
    __slots__ = ("tipo", "contorno", "hijos")

    def __init__(self, tipo, contorno):
        self.tipo = tipo
        self.contorno = [pt(*p) for p in contorno]
        self.hijos = []          # lista de (T, geometria)

    def add_child(self, T, geom):
        self.hijos.append((tuple(T), geom))

    def eval_child(self, n, i):
        """hat.js  evalChild(n,i) = transPt(children[n].T, children[n].geom.shape[i])"""
        T, geom = self.hijos[n]
        return trans_pt(T, geom.contorno[i])

    def recentre(self):
        """hat.js  recentre(): promedio de vertices; traslada contorno e hijos (solo traslada)."""
        cx = sum(p[0] for p in self.contorno) / len(self.contorno)
        cy = sum(p[1] for p in self.contorno) / len(self.contorno)
        M = ttrans(-cx, -cy)
        self.contorno = [trans_pt(M, p) for p in self.contorno]
        self.hijos = [(mul(M, T), geom) for (T, geom) in self.hijos]


def h_init():
    return _init("H", H_OUTLINE, H_INIT, H_INIT_ETIQUETAS)


def t_init():
    return _init("T", T_OUTLINE, T_INIT, T_INIT_ETIQUETAS)


def p_init():
    return _init("P", P_OUTLINE, P_INIT, P_INIT_ETIQUETAS)


def f_init():
    return _init("F", F_OUTLINE, F_INIT, F_INIT_ETIQUETAS)


def _init(tipo, contorno, matrices, etiquetas):
    m = Meta(tipo, contorno)
    for T, et in zip(matrices, etiquetas):
        m.add_child(T, Hoja(et))
    return m


# --------------------------------------------------------------------------------------
# 3. Reglas de sustitucion (verbatim de constructPatch en hat.js)
# --------------------------------------------------------------------------------------

REGLAS = (
    ("H",),
    (0, 0, "P", 2),
    (1, 0, "H", 2),
    (2, 0, "P", 2),
    (3, 0, "H", 2),
    (4, 4, "P", 2),
    (0, 4, "F", 3),
    (2, 4, "F", 3),
    (4, 1, 3, 2, "F", 0),
    (8, 3, "H", 0),
    (9, 2, "P", 0),
    (10, 2, "H", 0),
    (11, 4, "P", 2),
    (12, 0, "H", 2),
    (13, 0, "F", 3),
    (14, 2, "F", 1),
    (15, 3, "H", 4),
    (8, 2, "F", 1),
    (17, 3, "H", 0),
    (18, 2, "P", 0),
    (19, 2, "H", 2),
    (20, 4, "F", 3),
    (20, 0, "P", 2),
    (22, 0, "H", 2),
    (23, 4, "F", 3),
    (23, 0, "F", 3),
    (16, 0, "P", 2),
    (9, 4, 0, 2, "T", 2),
    (4, 0, "F", 3),
)

def etiqueta_regla(r):
    """Etiqueta de la regla: ('H',) -> r[0]; (i,v,'X',n) -> r[2]; (i,v,j,w,'X',n) -> r[4]."""
    return r[0] if len(r) == 1 else r[-2]


TIPOS_PATCH = "".join(etiqueta_regla(r) for r in REGLAS)


def construct_patch(H, T, P, F):
    """hat.js  constructPatch(H,T,P,F): las 29 reglas, en orden. Devuelve el parche."""
    shapes = {"H": H, "T": T, "P": P, "F": F}
    ret = Meta("parche", ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)))
    for r in REGLAS:
        if len(r) == 1:
            ret.add_child(IDENT, shapes[r[0]])
        elif len(r) == 4:
            i, v, lab, ni = r
            Tc, ch = ret.hijos[i]
            poly = ch.contorno
            P_ = trans_pt(Tc, poly[(v + 1) % len(poly)])
            Q_ = trans_pt(Tc, poly[v])
            npoly = shapes[lab].contorno
            TT = match_two(npoly[ni], npoly[(ni + 1) % len(npoly)], P_, Q_)
            ret.add_child(TT, shapes[lab])
        else:
            i0, v0, i2, v2, lab, ni = r
            TcP, chP = ret.hijos[i0]
            TcQ, chQ = ret.hijos[i2]
            P_ = trans_pt(TcQ, chQ.contorno[v2])
            Q_ = trans_pt(TcP, chP.contorno[v0])
            npoly = shapes[lab].contorno
            TT = match_two(npoly[ni], npoly[(ni + 1) % len(npoly)], P_, Q_)
            ret.add_child(TT, shapes[lab])
    return ret


HIJOS_NUEVO_H = (0, 9, 16, 27, 26, 6, 1, 8, 10, 15)
HIJOS_NUEVO_P = (7, 2, 3, 4, 28)
HIJOS_NUEVO_F = (21, 20, 22, 23, 24, 25)
HIJOS_NUEVO_T = (11,)


def construct_metatiles(patch):
    """hat.js  constructMetatiles(patch): verbatim (indices, angulos, hijos, recentre)."""
    bps1 = patch.eval_child(8, 2)
    bps2 = patch.eval_child(21, 2)
    rbps = trans_pt(rot_about(bps1, -2.0 * math.pi / 3.0), bps2)
    p72 = patch.eval_child(7, 2)
    p252 = patch.eval_child(25, 2)
    llc = intersect(bps1, rbps, patch.eval_child(6, 2), p72)

    w = psub(patch.eval_child(6, 2), llc)
    oH = [llc, bps1]
    w = trans_pt(trot(-math.pi / 3.0), w)
    oH.append(padd(oH[1], w))
    oH.append(patch.eval_child(14, 2))
    w = trans_pt(trot(-math.pi / 3.0), w)
    oH.append(psub(oH[3], w))
    oH.append(patch.eval_child(6, 2))
    new_H = Meta("H", oH)
    for i in HIJOS_NUEVO_H:
        new_H.add_child(patch.hijos[i][0], patch.hijos[i][1])

    oP = [p72, padd(p72, psub(bps1, llc)), bps1, llc]
    new_P = Meta("P", oP)
    for i in HIJOS_NUEVO_P:
        new_P.add_child(patch.hijos[i][0], patch.hijos[i][1])

    oF = [bps2, patch.eval_child(24, 2), patch.eval_child(25, 0), p252,
          padd(p252, psub(llc, bps1))]
    new_F = Meta("F", oF)
    for i in HIJOS_NUEVO_F:
        new_F.add_child(patch.hijos[i][0], patch.hijos[i][1])

    AAA = oH[2]
    BBB = padd(oH[1], psub(oH[4], oH[5]))
    CCC = trans_pt(rot_about(BBB, -math.pi / 3.0), AAA)
    new_T = Meta("T", [BBB, CCC, AAA])
    for i in HIJOS_NUEVO_T:
        new_T.add_child(patch.hijos[i][0], patch.hijos[i][1])

    for m in (new_H, new_P, new_F, new_T):
        m.recentre()
    return [new_H, new_T, new_P, new_F]


# --------------------------------------------------------------------------------------
# 4. Niveles, hojas y certificado por matriz
# --------------------------------------------------------------------------------------

def nivel_k(k):
    """Aplica k veces el bucle de la app: parche sobre las cuatro metatiles, luego sustitucion."""
    H, T, P, F = h_init(), t_init(), p_init(), f_init()
    for _ in range(k):
        H, T, P, F = construct_metatiles(construct_patch(H, T, P, F))
    return H, T, P, F


def hojas(meta, T=IDENT, salida=None):
    """Hojas (teselas) con su transformada de mundo. Unidad: tesela."""
    if salida is None:
        salida = []
    for Tc, geom in meta.hijos:
        Tw = mul(T, Tc)
        if isinstance(geom, Hoja):
            salida.append((Tw, geom.etiqueta))
        else:
            hojas(geom, Tw, salida)
    return salida


def n_hojas(meta):
    return len(hojas(meta))


def contorno_tesela(T):
    return [trans_pt(T, p) for p in HAT_OUTLINE]


def tipo_de(geom):
    """Especie del hijo: Meta -> su .tipo; Hoja -> 'H' si la etiqueta empieza por H, si no su etiqueta."""
    if isinstance(geom, Hoja):
        return "H" if geom.etiqueta.startswith("H") else geom.etiqueta
    return geom.tipo


def cuenta_tipos(meta):
    """Cuenta hijos por especie (unidad: metatile en nivel>=1; tesela en nivel 0)."""
    c = {"H": 0, "T": 0, "P": 0, "F": 0}
    for _, geom in meta.hijos:
        c[tipo_de(geom)] += 1
    return (c["H"], c["T"], c["P"], c["F"])


def cuenta_hojas_por_especie(meta):
    """Cuenta hojas por especie (unidad: tesela)."""
    c = {"H": 0, "T": 0, "P": 0, "F": 0}
    for _, etiqueta in hojas(meta):
        c["H" if etiqueta.startswith("H") else etiqueta] += 1
    return (c["H"], c["T"], c["P"], c["F"])


def matvec(M, v):
    return tuple(sum(M[i][j] * v[j] for j in range(4)) for i in range(4))


def nivel_matriz(k):
    """h_k por aritmetica entera. Convencion declarada: M[i][j] = numero de hijos de especie j
    dentro del metatil de especie i (filas y columnas 0..3 = H,T,P,F), luego
    h_{k+1}[s] = sum_j M[s][j] * h_k[j], es decir h_{k+1} = M * h_k (SIN transponer).
    La traspuesta daria (21, 4, 16, 22) en el nivel 1; el arbol mide (25, 4, 14, 16)."""
    h = H0_TESELAS
    for _ in range(k):
        h = matvec(M_MADRE, h)
    return h


def perron(M):
    """Autovalor de mayor modulo (doble precision, numpy)."""
    return max(abs(complejo) for complejo in np.linalg.eigvals(np.array(M, dtype=float)))


def certificado_matriz(nivel_max, verbose=True):
    """Mide, nivel por nivel: hijos por tipo, hojas por metatile y total del parche."""
    filas = []
    for k in range(nivel_max + 1):
        H, T, P, F = nivel_k(k)
        c = [cuenta_tipos(m) for m in (H, T, P, F)]
        esp = [cuenta_hojas_por_especie(m) for m in (H, T, P, F)]
        hojas_met = tuple(n_hojas(m) for m in (H, T, P, F))
        h_esperado = nivel_matriz(k)
        filas.append({
            "nivel": k,
            "hijos_por_tipo": c,
            "hojas_por_especie": esp,
            "hojas_por_metatile": hojas_met,
            "hojas_por_metatile_esperado": h_esperado,
            "hojas_coincide": tuple(hojas_met) == tuple(h_esperado),
            "tipos_H": "".join(tipo_de(g) for _, g in H.hijos),
            "tipos_T": "".join(tipo_de(g) for _, g in T.hijos),
            "tipos_P": "".join(tipo_de(g) for _, g in P.hijos),
            "tipos_F": "".join(tipo_de(g) for _, g in F.hijos),
        })
        if k == 0:
            filas[-1]["tipos_parche"] = TIPOS_PATCH
        else:
            parche = construct_patch(*nivel_k(k - 1))
            filas[-1]["tipos_parche"] = "".join(tipo_de(g) for _, g in parche.hijos)
            filas[-1]["hojas_parche"] = n_hojas(parche)
    return filas


# --------------------------------------------------------------------------------------
# 5. Rasterizado (rejilla regular en la ventana)
# --------------------------------------------------------------------------------------

def rasterizar(poligonos, x0, x1, y0, y1, n):
    """Union de poligonos en una rejilla n x n. Devuelve (union, conteo).
    Unidad: celda. union = celda cubierta al menos una vez; conteo = cuantas veces."""
    xs = x0 + (np.arange(n) + 0.5) * (x1 - x0) / n
    ys = y0 + (np.arange(n) + 0.5) * (y1 - y0) / n
    X, Y = np.meshgrid(xs, ys)
    union = np.zeros((n, n), dtype=bool)
    conteo = np.zeros((n, n), dtype=np.int16)
    for poly in poligonos:
        px = [p[0] for p in poly]
        py = [p[1] for p in poly]
        j0 = max(0, int((min(px) - x0) / (x1 - x0) * n) - 1)
        j1 = min(n, int((max(px) - x0) / (x1 - x0) * n) + 2)
        i0 = max(0, int((min(py) - y0) / (y1 - y0) * n) - 1)
        i1 = min(n, int((max(py) - y0) / (y1 - y0) * n) + 2)
        if j0 >= j1 or i0 >= i1:
            continue
        XX = X[i0:i1, j0:j1]
        YY = Y[i0:i1, j0:j1]
        dentro = np.zeros(XX.shape, dtype=bool)
        m = len(px)
        for k in range(m):
            xa, ya = px[k], py[k]
            xb, yb = px[(k + 1) % m], py[(k + 1) % m]
            if ya == yb:
                continue
            cond = (ya > YY) != (yb > YY)
            xint = (xb - xa) * (YY - ya) / (yb - ya) + xa
            dentro ^= cond & (XX < xint)
        union[i0:i1, j0:j1] |= dentro
        conteo[i0:i1, j0:j1] += dentro.astype(np.int16)
    return union, conteo


def bbox(puntos):
    xs = [p[0] for p in puntos]
    ys = [p[1] for p in puntos]
    return min(xs), max(xs), min(ys), max(ys)


def ventana_para(puntos, lado=None):
    """Ventana CUADRADA centrada en el bbox de las teselas. Unidad: unidades."""
    x0, x1, y0, y1 = bbox(puntos)
    ancho = max(x1 - x0, y1 - y0)
    lado = ancho if lado is None else float(lado)
    cx = 0.5 * (x0 + x1)
    cy = 0.5 * (y0 + y1)
    return cx - 0.5 * lado, cx + 0.5 * lado, cy - 0.5 * lado, cy + 0.5 * lado


def medir(poligonos, ven, n):
    """Devuelve cifras de la ventana. Unidad de las cifras: celda / fraccion."""
    x0, x1, y0, y1 = ven
    union, conteo = rasterizar(poligonos, x0, x1, y0, y1, n)
    total = float(n * n)
    return {
        "union": union,
        "conteo": conteo,
        "cobertura": float(union.sum()) / total,
        "solape": float((conteo >= 2).sum()) / total,
        "sha16_mascara": sha16_de(np.ascontiguousarray(union).tobytes()),
        "suma_areas": sum(area(p) for p in poligonos) / ((x1 - x0) * (y1 - y0)),
    }


# --------------------------------------------------------------------------------------
# 6. U1 y U2
# --------------------------------------------------------------------------------------

def u1(nivel, n, lado=None):
    """U1: jerarquia determinista. Devuelve diccionario con cifras y mascara."""
    H, T, P, F = nivel_k(nivel)
    parche = construct_patch(H, T, P, F)
    teselas = hojas(parche)
    polis = [contorno_tesela(T_) for T_, _ in teselas]
    ven = ventana_para([p for poly in polis for p in poly], lado=lado)
    cif = medir(polis, ven, n)
    ar = [area(p) for p in polis]
    det = [T_[0] * T_[4] - T_[1] * T_[3] for T_, _ in teselas]
    cif.update({
        "nivel": nivel,
        "n_teselas": len(teselas),
        "ventana": ven,
        "lado": ven[1] - ven[0],
        "celda": (ven[1] - ven[0]) / n,
        "orientaciones": [math.atan2(T_[3], T_[0]) for T_, _ in teselas],
        "centroides": [centroide(contorno_tesela(T_)) for T_, _ in teselas],
        "area_tesela_min": min(ar),
        "area_tesela_max": max(ar),
        "area_tesela_media": sum(ar) / len(ar),
        "det_min": min(det),
        "det_max": max(det),
        "n_vertices_por_tesela": len(polis[0]),
        "tipos_parche": "".join(tipo_de(g) for _, g in parche.hijos),
    })
    del cif["union"], cif["conteo"]
    cif["_mascara"] = None
    return cif, teselas, polis, ven


def centroide(poly):
    x = sum(p[0] for p in poly) / len(poly)
    y = sum(p[1] for p in poly) / len(poly)
    return (x, y)


def angulo_borde0(poly):
    """Orientacion MEDIDA de una tesela: angulo de su primer borde (radianes).
    Mide el poligono, no la transformada que lo coloco: la cifra no hereda la convencion."""
    return math.atan2(poly[1][1] - poly[0][1], poly[1][0] - poly[0][0])


def acoplamiento_orientacion(centros, angs, nbin=10, bloque=512):
    """C(r) = media de cos(2*dtheta) sobre los pares con distancia en el bin r. Unidad: adimensional.
    NULO DECLARADO: si las orientaciones viven en multiplos de 30 grados y el emparejamiento es
    uniforme, la media vale 1/12 = 0.083333 (sum cos(60*k), k=0..11, = 0). Un valor distinto del
    nulo mide acoplamiento orientacion-distancia; un valor igual al nulo, desacoplamiento."""
    X = np.array([c[0] for c in centros], dtype=float)
    Y = np.array([c[1] for c in centros], dtype=float)
    A = np.array(angs, dtype=float)
    N = len(X)
    rmax = 0.0
    for i0 in range(0, N, bloque):
        i1 = min(N, i0 + bloque)
        dx = X[i0:i1, None] - X[None, :]
        dy = Y[i0:i1, None] - Y[None, :]
        rmax = max(rmax, float(np.sqrt(dx * dx + dy * dy).max()))
    bordes = np.linspace(0.0, rmax, nbin + 1)
    suma = np.zeros(nbin)
    cuen = np.zeros(nbin, dtype=np.int64)
    for i0 in range(0, N, bloque):
        i1 = min(N, i0 + bloque)
        dx = X[i0:i1, None] - X[None, :]
        dy = Y[i0:i1, None] - Y[None, :]
        r = np.sqrt(dx * dx + dy * dy)
        w = np.cos(2.0 * (A[i0:i1, None] - A[None, :]))
        mimo = np.arange(N)[None, :] > np.arange(i0, i1)[:, None]
        r = r[mimo]
        w = w[mimo]
        idx = np.clip(np.digitize(r, bordes) - 1, 0, nbin - 1)
        np.add.at(suma, idx, w)
        np.add.at(cuen, idx, 1)
    filas = []
    for b in range(nbin):
        m = float(suma[b] / cuen[b]) if cuen[b] else float("nan")
        filas.append((float(bordes[b]), float(bordes[b + 1]), int(cuen[b]), m))
    return filas, rmax


def u2(polis1, teselas, ven, n, cobertura_objetivo, semilla, barrido, jitter_rel=0.35):
    """U2 (nulo decisivo): MISMAS teselas. El multiconjunto de poligonos es identico al de U1 (cada
    tesela se lleva a su centroide y se gira a orientacion 0: no cambia forma ni area); orientaciones
    barajadas reasignadas a posiciones barajadas sobre un reticulo triangular con jitter, y espaciado
    ajustado por BARRIDO medido para igualar la fraccion de llenado de U1. Declarado: solape permitido
    (lo que se mide es la union)."""
    rng = random.Random(semilla)
    N = len(polis1)
    orient = [math.atan2(T_[3], T_[0]) for T_, _ in teselas]
    barajadas = list(orient)
    rng.shuffle(barajadas)
    canon = []
    for p, th in zip(polis1, orient):
        cx, cy = centroide(p)
        ca, sa = math.cos(-th), math.sin(-th)
        q = []
        for (x, y) in p:
            bx, by = x - cx, y - cy
            q.append((ca * bx - sa * by, sa * bx + ca * by))
        canon.append(q)
    x0, x1, y0, y1 = ven
    ancho, alto = x1 - x0, y1 - y0
    a0 = math.sqrt(ancho * alto / N)
    curva = []
    mejor = None
    for paso in range(barrido):
        d = a0 * (0.40 + 0.08 * paso)
        jit = jitter_rel * d
        crudos = []
        ny = int(alto / (d * HR3)) + 4
        nx = int(ancho / d) + 4
        for fila in range(ny):
            for col in range(nx):
                x = x0 - d + col * d + (0.5 * d if fila % 2 else 0.0)
                y = y0 - d + fila * d * HR3
                if x0 <= x <= x1 and y0 <= y <= y1:
                    crudos.append((x, y))
        rng.shuffle(crudos)
        puntos = crudos[:N]
        if len(puntos) < N:
            curva.append((d, None, len(crudos)))
            continue
        polis = []
        for i, (x, y) in enumerate(puntos):
            dxp = rng.uniform(-jit, jit)
            dyp = rng.uniform(-jit, jit)
            ca, sa = math.cos(barajadas[i]), math.sin(barajadas[i])
            q = [(x + dxp + ca * bx - sa * by, y + dyp + sa * bx + ca * by)
                 for (bx, by) in canon[i]]
            polis.append(q)
        cif = medir(polis, ven, n)
        costo = abs(len(crudos) - N)
        curva.append((d, cif["cobertura"], len(crudos)))
        if mejor is None or costo < mejor[0]:
            mejor = (costo, d, cif, polis, len(crudos))
    costo, d, cif, polis, n_reticulo = mejor
    ar1 = sorted(area(p) for p in polis1)
    ar2 = sorted(area(p) for p in polis)
    cif.update({
        "n_puntos_reticulo": n_reticulo,
        "descartes_de_reticulo": n_reticulo - N,
        "d_elegido": d,
        "jitter": jitter_rel * d,
        "semilla": semilla,
        "fraccion_de_llenado": cif["suma_areas"],
        "delta_cobertura_vs_U1": abs(cif["cobertura"] - cobertura_objetivo),
        "curva_barrido": curva,
        "areas_max_dif_vs_U1": max(abs(a - b) for a, b in zip(ar1, ar2)),
        "n_vertices_igual": len(polis[0]) == len(polis1[0]),
        "solape_permitido": True,
    })
    del cif["union"], cif["conteo"]
    return cif, polis


# --------------------------------------------------------------------------------------
# 7. Evidencia de que los contornos son esquematicos
# --------------------------------------------------------------------------------------

def evidencia_esquema(nivel):
    """Mide: (a) suma de areas de hojas contra area del contorno en el nivel 0;
    (b) cuantas vertices de las hojas de H_k caen FUERA del contorno de H_k."""
    filas0 = []
    for m in (h_init(), t_init(), p_init(), f_init()):
        s = sum(area(contorno_tesela(T_)) for T_, _ in hojas(m))
        filas0.append({
            "tipo": m.tipo,
            "contorno_area": area(m.contorno),
            "suma_hojas": s,
            "suma_sobre_contorno": s / area(m.contorno),
            "n_hojas": n_hojas(m),
        })
    H, T, P, F = nivel_k(nivel)
    dentro = fuera = 0
    for T_, _ in hojas(H):
        for p in contorno_tesela(T_):
            if punto_en_poligono(p, H.contorno):
                dentro += 1
            else:
                fuera += 1
    return filas0, {"nivel": nivel, "dentro": dentro, "fuera": fuera,
                    "total": dentro + fuera}


# --------------------------------------------------------------------------------------
# 8. main
# --------------------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description="Red dielectrica del monotilo Hat: U1 y U2.")
    ap.add_argument("--nivel", type=int, default=2, help="nivel de sustitucion de U1")
    ap.add_argument("--n", type=int, default=384, help="celdas por lado de la ventana")
    ap.add_argument("--semilla", type=int, default=20261012, help="semilla de U2")
    ap.add_argument("--lado", type=float, default=None, help="lado de la ventana (unidades)")
    ap.add_argument("--tol", type=float, default=0.02, help="tolerancia de fraccion de llenado")
    ap.add_argument("--barrido", type=int, default=24, help="pasos del barrido de densidad de U2")
    ap.add_argument("--json", action="store_true", help="imprime el crudo en JSON")
    ap.add_argument("--salida", type=str, default=None,
                    help="carpeta donde guardar las mascaras .npy (por omision: NO se escribe nada)")
    args = ap.parse_args(argv)

    t0 = t_utc()
    cmd = "python3 " + " ".join(sys.argv[1:])
    filas = certificado_matriz(args.nivel)
    ev0, evk = evidencia_esquema(args.nivel)
    cif1, teselas, polis1, ven = u1(args.nivel, args.n, lado=args.lado)
    cif2, polis2 = u2(polis1, teselas, ven, args.n, cif1["cobertura"], args.semilla, args.barrido)
    aco1, rmax1 = acoplamiento_orientacion([centroide(p) for p in polis1],
                                           [angulo_borde0(p) for p in polis1])
    aco2, rmax2 = acoplamiento_orientacion([centroide(p) for p in polis2],
                                           [angulo_borde0(p) for p in polis2])
    lam = perron(M_MADRE)
    frac_delta = abs(cif2["suma_areas"] - cif1["suma_areas"])
    veredicto_frac = ("FRACCION_LLENADO_OK" if frac_delta <= args.tol * cif1["suma_areas"]
                      else "FRACCION_LLENADO_DIFIERE")

    guardado = {}
    if args.salida:
        os.makedirs(args.salida, exist_ok=True)
        for nom, polis in (("u1", polis1), ("u2", polis2)):
            u, _ = rasterizar(polis, *ven, args.n)
            ruta = os.path.join(args.salida, "malla_%s_n%d_k%d.npy" % (nom, args.n, args.nivel))
            np.save(ruta, u)
            guardado[nom] = ruta

    crudo = {
        "comando": cmd,
        "ventana_utc": t0,
        "universo_declarado": UNIVERSO,
        "fuente": {"hat_js": FUENTE_HAT, "sha16_hat_js": SHA16_HAT,
                   "geometry_js": FUENTE_GEO, "sha16_geometry_js": SHA16_GEO,
                   "volcado_verbatim_utc": VENTANA_FUENTE},
        "unidades": {"conteo": ["tesela", "vertice", "celda"],
                     "superficie": "unidades^2", "fraccion": "adimensional"},
        "matriz": {"M": M_MADRE, "h0": H0_TESELAS,
                   "perron_medido": lam, "perron_esperado_phi4": PHI4_ESPERADO,
                   "delta_perron": abs(lam - PHI4_ESPERADO),
                   "filas": filas},
        "contornos_esquematicos": {"nivel0": ev0, "nivel_k": evk},
        "u1": cif1,
        "u2": cif2,
        "acoplamiento_orientacion": {"C_r_U1": aco1, "C_r_U2": aco2,
                                     "nulo_uniforme_1_12": 1.0 / 12.0,
                                     "rmax_U1": rmax1, "rmax_U2": rmax2},
        "fraccion_de_llenado": {
            "definicion": "suma de areas de las teselas / area de la ventana (cuenta solape)",
            "U1": cif1["suma_areas"], "U2": cif2["suma_areas"],
            "delta_abs": frac_delta, "tol_abs": args.tol * cif1["suma_areas"],
            "veredicto": veredicto_frac,
            "cobertura_union_U1": cif1["cobertura"], "cobertura_union_U2": cif2["cobertura"],
            "solape_U1": cif1["solape"], "solape_U2": cif2["solape"],
            "nota": ("la cobertura de la union NO es la magnitud calibrada: barajar crea solape "
                     "donde U1 no lo tiene, y la union baja por construccion. Se reporta, no se "
                     "usa como criterio.")},
        "guardado": guardado,
    }

    if args.json:
        print(json.dumps(crudo, ensure_ascii=False, default=str, indent=1))
        return 0

    print("COMANDO: " + cmd)
    print("VENTANA UTC: " + t0)
    print("UNIVERSO DECLARADO: " + UNIVERSO)
    print("FUENTE: %s sha16 %s | %s sha16 %s | volcado verbatim %s"
          % (FUENTE_HAT, SHA16_HAT, FUENTE_GEO, SHA16_GEO, VENTANA_FUENTE))
    print()
    print("--- certificado por matriz (aritmetica entera, exacta; unidad: metatile / tesela) ---")
    print("M (filas H,T,P,F; columnas H,T,P,F) = %s" % (M_MADRE,))
    print("h_0 (teselas por metatile inicial) = %s" % (H0_TESELAS,))
    for f in filas:
        print("nivel %d: hijos_por_tipo=%s hojas_por_especie=%s (unidad: metatile) hojas_metatile=%s"
              " esperado=%s coincide=%s"
              % (f["nivel"], f["hijos_por_tipo"], f["hojas_por_especie"],
                 f["hojas_por_metatile"], f["hojas_por_metatile_esperado"], f["hojas_coincide"]))
        print("         tipos H/T/P/F = %s / %s / %s / %s"
              % (f["tipos_H"], f["tipos_T"], f["tipos_P"], f["tipos_F"]))
    print("Perron(M) medido = %.12f | phi^4 esperado = %.12f | delta = %.3e"
          % (lam, PHI4_ESPERADO, abs(lam - PHI4_ESPERADO)))
    print("tipos del parche (29 hijos) = %s" % TIPOS_PATCH)
    print()
    print("--- contornos esquematicos (evidencia medida; NO son certificado) ---")
    for f in ev0:
        print("nivel 0 %s: contorno=%.6f suma_hojas=%.6f suma/contorno=%.6f (unidad: unidades^2; %d hojas)"
              % (f["tipo"], f["contorno_area"], f["suma_hojas"], f["suma_sobre_contorno"], f["n_hojas"]))
    print("nivel %d: vertices de las hojas de H_%d fuera de su contorno = %d de %d"
          % (evk["nivel"], evk["nivel"], evk["fuera"], evk["total"]))
    print()
    print("--- U1 (jerarquia determinista) ---")
    print("nivel=%d teselas=%d (unidad: tesela) ventana=%.4f..%.4f x %.4f..%.4f lado=%.4f "
          "celdas=%dx%d celda=%.6f"
          % (cif1["nivel"], cif1["n_teselas"], ven[0], ven[1], ven[2], ven[3], cif1["lado"],
             args.n, args.n, cif1["celda"]))
    print("suma_areas/ventana=%.6f cobertura=%.6f solape=%.6f sha16_mascara=%s"
          % (cif1["suma_areas"], cif1["cobertura"], cif1["solape"], cif1["sha16_mascara"]))
    print("area por tesela (unidad: unidades^2) min=%.6f media=%.6f max=%.6f | det(T) min=%.6f max=%.6f"
          % (cif1["area_tesela_min"], cif1["area_tesela_media"], cif1["area_tesela_max"],
             cif1["det_min"], cif1["det_max"]))
    print("tipos del parche = %s" % cif1["tipos_parche"])
    print()
    print("--- U2 (nulo decisivo: mismas teselas, orientaciones y posiciones barajadas) ---")
    print("N=%d puntos de reticulo=%d descartes=%d semilla=%d jitter=%.6f d=%.6f (unidad: unidades)"
          % (len(polis2), cif2["n_puntos_reticulo"], cif2["descartes_de_reticulo"],
             cif2["semilla"], cif2["jitter"], cif2["d_elegido"]))
    print("barrido d -> cobertura: " + " ".join(
        "%.3f:%s" % (d, ("%.4f" % c) if c is not None else "n/a") for d, c, _ in cif2["curva_barrido"]))
    print("cobertura de la union U1=%.6f U2=%.6f delta=%.6f (adimensional; NO es la magnitud calibrada)"
          % (cif1["cobertura"], cif2["cobertura"], cif2["delta_cobertura_vs_U1"]))
    print("FRACCION DE LLENADO (suma de areas / area de la ventana, cuenta solape) U1=%.6f U2=%.6f"
          " delta=%.3e tol_abs=%.3e VEREDICTO=%s"
          % (cif1["suma_areas"], cif2["suma_areas"], frac_delta, args.tol * cif1["suma_areas"],
             veredicto_frac))
    print("suma_areas/ventana=%.6f solape_U2=%.6f (declarado: solape permitido) sha16_mascara=%s"
          % (cif2["suma_areas"], cif2["solape"], cif2["sha16_mascara"]))
    print("areas: max|area_U2-area_U1| sobre el multiconjunto ordenado = %.3e | n_vertices igual = %s"
          % (cif2["areas_max_dif_vs_U1"], cif2["n_vertices_igual"]))
    print()
    print("--- acoplamiento orientacion-distancia (C(r)=media de cos(2*dtheta); nulo uniforme = 1/12 = 0.083333) ---")
    for nom, aco in (("U1", aco1), ("U2", aco2)):
        print("%s: " % nom + " ".join("%.2f-%.2f:n=%d,C=%.4f" % (a, b, k, m) for a, b, k, m in aco))
    if guardado:
        print("mascaras guardadas: %s" % guardado)
    print()
    print("CIFRAS: unidad declarada en cada linea; comando y hora UTC en la cabeza del crudo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
