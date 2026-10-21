# HUELLA-EJECUTOR-INYECTOR-DIEL-20261021
"""Ejecutor directo de `inyector_diel.py` — fotonica-monotile (ACTO 3.0).

`inyector_diel.py` construye el campo `eps` del dielectrico sobre la malla
certificada y verifica con cotejos fail-closed que la inyeccion no rompe la
PML ni la geometria sellada. Hasta hoy ningun test lo abria: el censo
`censo_numpy_dependencia_20261020.py` lo listaba `SIN_EJECUTOR`. Su unica
dependencia de numpy era una suposicion.

Contrato medido antes de escribirse (tres sondas read-only):

  * `construir(verbose=False) -> dict` con 25 claves y `rc = 2 if f.lista else
    0`; `main([...])` devuelve ese rc (fail-closed real, un entero, no un
    `SystemExit` decorativo). Cuesta 15.1 s: se corre UNA vez por modulo.
  * los fallos de `campo_eps`, `sonda_hook` y `referencia_sellada` son
    `SystemExit` con `.code` **str**.
  * imports: stdlib + numpy + `malla_hat` + `fdtd2d_tm` (los dos del repo).
    Cero dependencia cross-repo.

Desviaciones declaradas (no ocultas):

  1. `correr_acto3()` NO se corre: encadena U1 + 5 semillas FDTD en 384^2
     (~46 s contando `construir`) y escribe crudos. El ACTO 3.0 es el que se
     audita aqui; el 3.2 ya tiene su propio instrumento.
  2. la clave `eps["borde_pml_con_eps_0"]` del resumen es el literal `True`,
     no la lectura de su `chk`: NO se asierta esa clave (seria asertar una
     constante). El borde se mide directamente con `campo_eps`.
  3. `FDTD_SHA16_DECLARADO` ya no casa con `fdtd2d_tm.py` (el ACTO 3.1 anadio
     el hook `eps`): el instrumento imprime una [nota] y NO se pone rojo. Se
     asierta la procedencia real (hash del fichero), no la nota.

Falsador declarado (lo que hace rojo a este test):

  1. numpy ausente -> `import inyector_diel` levanta `ModuleNotFoundError`.
  2. `malla_hat.py` mutado -> su `chk` de certificado se pone ROJO -> `rc != 0`
     y cae el test del rc.
  3. `campo_eps` deja de transponer o de limpiar el borde -> caen los dos tests
     directos de `eps`.
  4. la guarda de forma de `campo_eps` desaparece -> cae su `SystemExit`.
  5. `main()` deja de propagar el rc -> cae el test del monkeypatch.

Mutantes medidos en copia (`/tmp/fot-fals`, restaurado despues): 5 de 5 muertos.
  M6 `campo_eps` sin `.T`      -> cae la transposicion
  M7 `pml = 0` en `campo_eps`  -> caen el borde y la transposicion
  M8 `main: return 0`          -> cae la propagacion del rc
  M9 `MODO = "C1a"`            -> el instrumento firma rc 2 y cae el test del rc
  M10 `malla_hat.py` mutado    -> el chk del certificado se pone ROJO (rc 2)
M7 SOBREVIVIO a la primera version del test del borde: con `pml = 0` el slice
`eps[-0:, :]` barre el array entero y un borde vaciado seguia cumpliendo
`eps[borde] == EPS_0`. Por eso el test ahora cuenta las celdas en `eps_mat`:
limpiar el borde no puede borrar el objeto.
"""

from __future__ import annotations

import contextlib
import io
import inspect
import math
import sys
from pathlib import Path

import numpy as np
import pytest

RAIZ = Path(__file__).resolve().parents[1]
if str(RAIZ) not in sys.path:
    sys.path.insert(0, str(RAIZ))

import fdtd2d_tm  # noqa: E402
import inyector_diel as acto  # noqa: E402

NX = acto.N


@pytest.fixture(scope="module")
def crudo() -> dict:
    with contextlib.redirect_stdout(io.StringIO()):
        return acto.construir(verbose=False)


def _union_asimetrica() -> np.ndarray:
    union = np.zeros((NX, NX), dtype=bool)
    union[100:110, 200:205] = True
    return union


def test_el_acto_corre_y_firma_rc_cero(crudo):
    assert crudo["rc"] == 0
    assert crudo["cotejos_rojos"] == [], crudo["cotejos_rojos"]
    assert crudo["acto"].startswith("3.0")


def test_el_resumen_declara_la_geometria_sellada(crudo):
    assert crudo["nivel"] == 2
    assert crudo["N"] == crudo["nx"] == crudo["ny"] == NX == 384
    assert crudo["pml_celdas"] == 20
    assert crudo["modo"] == "C1b"
    assert crudo["eps_mat"] == 4.0 and crudo["eps_0"] == 1.0
    assert crudo["L_a"] == 24.0 and crudo["barrido"] == 24
    assert crudo["semillas"] == [101, 102, 103, 104, 105]


def test_la_procedencia_por_hash_es_la_del_repositorio(crudo):
    """`malla_hat.py` es el motor de geometria certificado; `fdtd2d_tm.py` se
    cita por el hash de su fichero (su certificado quedo atras al anadir el
    hook `eps`, y eso NO gatea: se mide el fichero, no la nota)."""
    assert crudo["malla_hat_sha16"] == acto.MALLA_SHA16_CERTIFICADO
    assert crudo["malla_hat_sha16"] == acto.sha16_de_archivo(RAIZ / "malla_hat.py")
    assert crudo["fdtd2d_tm_sha16"] == acto.sha16_de_archivo(RAIZ / "fdtd2d_tm.py")


def test_la_inyeccion_conserva_la_fraccion_de_llenado(crudo):
    sellada = crudo["u1"]["suma_areas"]
    assert sellada > 0.0
    for semilla in crudo["u2"]:
        assert semilla["delta_f"] <= acto.TOL_LLENADO_DECLARADA, semilla["delta_f"]
        assert math.isclose(
            semilla["fraccion_de_llenado"], sellada, abs_tol=acto.TOL_LLENADO_DECLARADA
        )


def test_las_cinco_semillas_mueven_la_mascara_y_son_distintas(crudo):
    hashes = {s["sha16_mascara"] for s in crudo["u2"]}
    assert len(hashes) == len(acto.SEMILLAS)
    assert crudo["u1"]["sha16_mascara"] not in hashes
    for semilla in crudo["u2"]:
        assert semilla["celdas_distintas_vs_U1"] > 0, semilla["semilla"]


def test_el_campo_eps_limpia_el_borde_de_la_pml():
    """El borde absorbente no admite dielectrico, ni con la union en la esquina.

    El conteo de celdas en `eps_mat` impide que el test lo firme un `eps`
    vaciado entero: limpiar el borde no puede borrar el objeto."""
    pml = int(fdtd2d_tm.PML_CELDAS)
    borde = np.ones((NX, NX), dtype=bool)
    borde[pml:-pml, pml:-pml] = False
    assert int(borde.sum()) > 0

    union = _union_asimetrica()
    eps = acto.campo_eps(union)
    assert bool(np.all(eps[borde] == acto.EPS_0)), int((eps[borde] != acto.EPS_0).sum())
    assert int((eps == acto.EPS_MAT).sum()) == int(union.sum())

    eps_esquina = acto.campo_eps(np.ones((NX, NX), dtype=bool))
    assert bool(np.all(eps_esquina[borde] == acto.EPS_0))
    assert int((eps_esquina == acto.EPS_MAT).sum()) == (NX - 2 * pml) ** 2


def test_el_campo_eps_transpone_la_mascara():
    """La convencion es eps(i_x, j_y) = union.T, y con una union asimetrica la
    transposicion es observable: sin el `.T` el test cae."""
    union = _union_asimetrica()
    eps = acto.campo_eps(union)
    pml = int(fdtd2d_tm.PML_CELDAS)
    dentro = (slice(pml, NX - pml), slice(pml, NX - pml))
    esperado = acto.EPS_0 + (acto.EPS_MAT - acto.EPS_0) * union.T
    assert bool(np.array_equal(eps[dentro], esperado[dentro]))
    assert not bool(np.array_equal(eps[dentro], esperado.T[dentro]))
    assert set(np.unique(eps).tolist()) == {acto.EPS_0, acto.EPS_MAT}


def test_el_campo_eps_rechaza_una_mascara_de_forma_ajena():
    with pytest.raises(SystemExit) as exc:
        acto.campo_eps(np.zeros((8, 8), dtype=bool))
    assert isinstance(exc.value.code, str), exc.value.code
    assert "1 celda de máscara = 1 celda de Yee" in exc.value.code


def test_el_hook_eps_del_solver_esta_nombrable():
    nombre = acto.sonda_hook()
    assert nombre == "eps"
    firma = inspect.signature(fdtd2d_tm.correr)
    param = firma.parameters[nombre]
    assert param.default is None
    assert param.kind in (
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        inspect.Parameter.KEYWORD_ONLY,
    )


def test_la_referencia_sellada_casa_con_el_modo_y_la_L_declarados():
    nu, a_dos, corrida = acto.referencia_sellada()
    assert acto.CRUDO_REF.exists()
    assert corrida["modo"] == acto.MODO
    assert math.isclose(float(corrida["L_a"]), acto.L_A, rel_tol=0.0, abs_tol=1e-12)
    assert len(nu) == len(a_dos) == 29
    assert np.all(nu >= fdtd2d_tm.NU_MIN) and np.all(nu <= fdtd2d_tm.NU_MAX)


def test_el_criterio_de_rachas_mide_la_racha_mas_larga():
    assert acto.rachas([True, True, False, True, True, True]) == [(0, 1), (3, 5)]
    assert acto.rachas([]) == []
    assert acto.rachas([False, False]) == []
    assert acto.rachas([True]) == [(0, 0)]
    rr = acto.rachas([False, True, True, True, False])
    assert max(b - a + 1 for a, b in rr) == acto.BINS_CONSECUTIVOS


def test_main_propaga_el_rc_de_construir(monkeypatch):
    """Fail-closed: si los cotejos se ponen rojos, `main` NO firma cero."""
    for rc in (0, 2, 7):
        monkeypatch.setattr(acto, "construir", lambda verbose=True, rc=rc: {"rc": rc})
        with contextlib.redirect_stdout(io.StringIO()):
            assert acto.main(["--json"]) == rc
