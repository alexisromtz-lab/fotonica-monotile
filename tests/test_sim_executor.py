# HUELLA-EJECUTOR-SIM-FOTONICA-20261020
"""Ejecutor del instrumento — fotonica-monotile.

`malla_hat.py` importa numpy para rasterizar la red dielectrica del monotilo y
produce las cifras de fraccion de llenado que el README publica, pero hasta hoy
ningun test ni paso del aparato lo abria: el censo
`censo_numpy_dependencia_20261020.py` lo clasifico `SIN_EJECUTOR`
(3 portadores, 0 ejecutados). Un import que nadie ejecuta no puede fallar, y
por tanto no era una dependencia: era una suposicion.

Este test corre el instrumento por su propia puerta (`main --json`, en proceso)
y **asierte sobre su salida**: la invariante de area de la sustitucion U1->U2 y
el veredicto de fraccion de llenado con el umbral que el propio instrumento
declara.

Falsador declarado (lo que hace rojo a este test):

  1. numpy ausente -> `import malla_hat` levanta `ModuleNotFoundError` en la
     coleccion: el test no puede ni empezar.
  2. la sustitucion U1->U2 deja de conservar el area -> el test de la
     invariante cae.
  3. el veredicto se mueve a `FRACCION_LLENADO_DIFIERE` -> cae el test del
     veredicto, que compara contra el umbral que el instrumento imprime.
  4. numpy 1.x activo -> cae el test del piso: el manifiesto declara 2.0.

La corrida usa `--nivel 1 --n 96 --barrido 4` (sub-malla de control), no la
ventana sellada: mide la conducta del instrumento, no vuelve a producir la
cifra publicada, y por eso no escribe `crudos/`.
"""

from __future__ import annotations

import contextlib
import io
import json
import re
import sys
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

import malla_hat as instrumento  # noqa: E402

ARGV = ["--nivel", "1", "--n", "96", "--barrido", "4", "--json"]


def _crudo():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = instrumento.main(ARGV)
    return rc, json.loads(buf.getvalue())


def test_el_instrumento_corre_y_emite_su_crudo():
    rc, crudo = _crudo()
    assert rc == 0, rc
    for clave in ("u1", "u2", "fraccion_de_llenado", "acoplamiento_orientacion"):
        assert clave in crudo, clave
    assert crudo["u1"]["n_teselas"] > 0
    assert re.fullmatch(r"[0-9a-f]{16}", crudo["u1"]["sha16_mascara"])


def test_la_sustitucion_conserva_el_area():
    """U2 sustituye teselas de U1 una a una: la suma de areas no se mueve."""
    _, crudo = _crudo()
    a1 = crudo["u1"]["suma_areas"]
    a2 = crudo["u2"]["suma_areas"]
    assert a1 > 0.0
    assert abs(a2 - a1) <= 1e-9 * a1, (a1, a2)


def test_el_veredicto_de_llenado_cumple_su_propio_umbral():
    _, crudo = _crudo()
    frac = crudo["fraccion_de_llenado"]
    assert frac["veredicto"] == "FRACCION_LLENADO_OK", frac["veredicto"]
    assert frac["delta_abs"] <= frac["tol_abs"], (frac["delta_abs"], frac["tol_abs"])


def test_la_cobertura_de_la_union_no_excede_la_suma_de_areas():
    """La suma cuenta solape: la union nunca puede ser mayor."""
    _, crudo = _crudo()
    for clave in ("u1", "u2"):
        bloque = crudo[clave]
        assert bloque["cobertura"] <= bloque["suma_areas"] + 1e-12, clave


def test_numpy_cargado_cumple_el_piso_declarado():
    version = instrumento.np.__version__
    assert re.fullmatch(r"\d+\.\d+\.\d+.*", version), version
    assert int(version.split(".")[0]) >= 2, version
