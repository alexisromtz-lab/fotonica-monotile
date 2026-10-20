# HUELLA-FOTONICA-TEST-CIFRAS-20261020
"""Cada cifra que este README publica tiene portador en un crudo sellado.

Este test hace exigible esa promesa en cada push. No vuelve a correr el FDTD
(72 a/c, ~40 s de pared): mide procedencia e identidad, que es exactamente lo
que un crudo sellado permite medir sin repetir la fisica.

Los falsadores, uno por cada forma de mentir de una cifra:

  1. una cifra de seis o mas decimales en el README sin hoja numerica que la
     respalde dentro de su propio redondeo ni formula analitica declarada;
  2. un crudo editado a mano: `sigma_u2`, `banda`, `dif` o `cumple` dejan de
     cumplir las identidades que el propio crudo declara;
  3. un contraste que cambia de veredicto: `28 de 29` deja de ser el numero
     derivado, o la racha `1` deja de estar por debajo de los `3` exigidos;
  4. un instrumento (`fdtd2d_tm.py`, `malla_hat.py`) tocado sin volver a sellar
     su `sha16` en el crudo que lo cita;
  5. un rojo borrado del escaparate (`ROJO x4`, commit `92834c4`).

No hay red y no hay tolerancias elegidas a mano: cada cifra se compara con su
propio redondeo (media unidad del ultimo digito publicado).
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import statistics
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
CRUDOS = RAIZ / "crudos"
README = RAIZ / "README.md"
CITATION = RAIZ / "CITATION.cff"
ZENODO = RAIZ / ".zenodo.json"

NOMBRE_C1 = "crudo_C1_20261012T062648Z.json"
NOMBRE_ACTO3 = "crudo_acto3_20261012T205230Z.json"

DOI_VERSION = "10.5281/zenodo.23329943"
DOI_CONCEPTO = "10.5281/zenodo.23329942"

# Una cifra "publicada" es la que se escribe con seis decimales o mas; las de
# menos digitos son redondeos de prosa y reclaman portador por otra via.
PATRON_CIFRA = re.compile(r"-?\d+\.\d{6,}(?:[eE][-+]?\d+)?")
PATRON_PROSA = re.compile(r"-?\d+\.\d{2,}(?:[eE][-+]?\d+)?")
PATRON_VERSION_CFF = re.compile(r'^version:\s*"?([0-9][^"\n]*?)"?\s*$', re.MULTILINE)
PATRON_VERSION_JSON = re.compile(r'"version"\s*:\s*"([^"]+)"')
PATRON_TITULO_CFF = re.compile(r'^title:\s*"(.+)"\s*$', re.MULTILINE)

REL = 1e-12


def _ruta(nombre: str) -> Path:
    ruta = CRUDOS / nombre
    assert ruta.is_file(), f"crudo ausente: {ruta.relative_to(RAIZ)}"
    return ruta


def _crudo(nombre: str) -> dict:
    ruta = _ruta(nombre)
    try:
        dato = json.loads(ruta.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:  # fail-closed: un crudo ilegible es rojo
        raise AssertionError(f"crudo ilegible: {ruta.name}: {exc}") from exc
    assert isinstance(dato, dict) and dato, f"crudo vacio: {ruta.name}"
    return dato


def _acto3() -> dict:
    return _crudo(NOMBRE_ACTO3)


def _c1() -> dict:
    return _crudo(NOMBRE_C1)


def _readme() -> str:
    return README.read_text(encoding="utf-8")


def _exige(texto: str, literal: str, donde: str) -> None:
    assert literal in texto, f"{donde} dejo de publicar `{literal}`"


def _quantum(literal: str) -> float:
    """Media unidad del ultimo digito publicado: la tolerancia de su redondeo."""
    m = re.fullmatch(r"(\d+)\.(\d+)(?:[eE]([-+]?\d+))?", literal.lstrip("-"))
    assert m, f"literal no reconocido: {literal}"
    return 0.5 * 10.0 ** (int(m.group(3) or 0) - len(m.group(2)))


def _hojas_numericas(objeto: object):
    if isinstance(objeto, dict):
        for valor in objeto.values():
            yield from _hojas_numericas(valor)
    elif isinstance(objeto, list):
        for valor in objeto:
            yield from _hojas_numericas(valor)
    elif isinstance(objeto, (int, float)) and not isinstance(objeto, bool):
        if math.isfinite(objeto):
            yield objeto


def _corrida(modo: str, largo: float) -> dict:
    """La corrida de C1 con ese modo y ese largo de dominio (L en unidades de a)."""
    c1 = _c1()
    for corrida in c1["corridas"]:
        if corrida["modo"] == modo and abs(corrida["L_a"] - largo) < 1e-12:
            return corrida
    raise AssertionError(f"el crudo de C1 no contiene la corrida {modo} con L = {largo} a")


def _cumplen_banda(acto: dict) -> int:
    """Bins que caen dentro de la banda, derivados del criterio (no de `cumple`)."""
    return sum(1 for dif, banda in zip(acto["dif"], acto["banda"]) if abs(dif) <= banda)


def _racha_mas_larga(falla: list[bool]) -> int:
    larga = actual = 0
    for incumple in falla:
        actual = actual + 1 if incumple else 0
        larga = max(larga, actual)
    return larga


def test_el_readme_publica_cifras():
    """Un README vaciado no puede pasar este test por vacuidad."""
    cifras = set(PATRON_PROSA.findall(_readme()))
    assert len(cifras) >= 20, f"el README perdio sus cifras: solo quedan {len(cifras)}"


def test_toda_cifra_de_seis_decimales_tiene_portador():
    """Ninguna cifra precisa del README vive sin hoja, texto de crudo o formula."""
    objetos = [_c1(), _acto3()]
    hojas = [valor for objeto in objetos for valor in _hojas_numericas(objeto)]
    textos = [_ruta(n).read_text(encoding="utf-8") for n in (NOMBRE_C1, NOMBRE_ACTO3)]
    assert hojas, "los crudos no tienen hojas numericas: sin portadores no hay auditoria"

    acto = _acto3()
    sigma0 = acto["sigma_u2"][0]
    riz_c1b = _c1()["riz_reparto_por_modo"]["C1b"]
    analiticos = {
        1.0 / math.sqrt(2.0): "1/sqrt(2)",
        1.0e0: "1/1",
        acto["criterio"]["k_sigma"] * sigma0 + riz_c1b: "k_sigma*sigma_u2[0]+riz_C1b",
    }

    huerfanas = []
    for literal in sorted(set(PATRON_CIFRA.findall(_readme()))):
        publicado = float(literal)
        tolerancia = _quantum(literal)
        if any(abs(valor - publicado) <= tolerancia for valor in hojas):
            continue
        if any(literal in texto for texto in textos):
            continue
        if any(abs(valor - publicado) <= tolerancia for valor in analiticos):
            continue
        huerfanas.append(literal)

    assert not huerfanas, f"cifras publicadas sin crudo que las respalde: {huerfanas}"


def test_los_crudos_sellan_el_instrumento_que_citan():
    """El crudo declara el `sha16` de su instrumento: un byte de deriva y es rojo."""
    acto = _acto3()
    for nombre, campo in (("fdtd2d_tm.py", "fdtd2d_tm_sha16"), ("malla_hat.py", "malla_hat_sha16")):
        ruta = RAIZ / nombre
        assert ruta.is_file(), f"instrumento ausente: {nombre}"
        medido = hashlib.sha256(ruta.read_bytes()).hexdigest()[:16]
        assert acto[campo] == medido, (
            f"{nombre} derivo respecto al crudo ({acto[campo]} -> {medido}): "
            "hay que re-sellar el crudo, no el test"
        )


def test_la_ventana_tiene_29_bins_y_es_la_declarada():
    acto, c1 = _acto3(), _c1()
    nu = acto["nu"]
    assert len(nu) == 29, f"la ventana tiene {len(nu)} bins, el acto midio 29"
    assert all(b > a for a, b in zip(nu, nu[1:])), "la ventana dejo de ser creciente"
    paso = nu[1] - nu[0]
    assert max(abs((b - a) - paso) for a, b in zip(nu, nu[1:])) < 1e-12, "el paso dejo de ser uniforme"
    bajo, alto = c1["ventana_nu"]
    assert bajo <= nu[0] and nu[-1] <= alto, f"los bins se salen de la ventana declarada ({bajo}, {alto})"
    _exige(_readme(), f"{len(nu)} bins", "el README")


def test_dif_es_la_resta_de_las_dos_transmisiones():
    acto = _acto3()
    peor = max(abs(dif) - abs(u1 - bar) for dif, u1, bar in zip(acto["dif"], acto["T_u1"], acto["T_bar_u2"]))
    assert abs(peor) < REL, f"`dif` no es |T_u1 - T_bar_u2|: residuo {peor}"


def test_sigma_es_la_desviacion_de_las_cinco_semillas():
    acto = _acto3()
    semillas = acto["T_semillas"]
    assert len(semillas) == 5, f"el acto declaro {len(semillas)} semillas"
    assert all(len(s) == len(acto["nu"]) for s in semillas), "una semilla no cubre los 29 bins"
    for i, sigma in enumerate(acto["sigma_u2"]):
        # Desviacion muestral (n - 1): con cinco semillas, la poblacional daria 0.237656 en el bin 0.
        medido = statistics.stdev([s[i] for s in semillas])
        assert abs(sigma - medido) <= REL * max(1.0, abs(sigma)), f"bin {i}: sigma_u2 {sigma} != {medido}"


def test_la_banda_es_k_sigma_mas_el_riz_del_crudo_de_c1():
    acto = _acto3()
    criterio = acto["criterio"]
    for i, banda in enumerate(acto["banda"]):
        medido = criterio["k_sigma"] * acto["sigma_u2"][i] + criterio["riz_T_C1"]
        assert abs(banda - medido) < 1e-12, f"bin {i}: banda {banda} != 3 sigma + riz {medido}"
    _exige(_readme(), f"k_sigma = {criterio['k_sigma']}", "el README")
    _exige(_readme(), f"riz_T_C1 = {criterio['riz_T_C1']}", "el README")


def test_cumple_es_la_negacion_del_criterio_y_no_su_afirmacion():
    """El campo `cumple` marca el bin que FALLA (README, lectura acotada)."""
    acto = _acto3()
    derivado = [abs(dif) > banda for dif, banda in zip(acto["dif"], acto["banda"])]
    assert acto["cumple"] == derivado, "`cumple` dejo de ser la negacion de |dif| <= banda"
    assert sum(acto["cumple"]) == 1, f"`cumple` marca {sum(acto['cumple'])} bins, el acto declaro 1"
    cumplen = _cumplen_banda(acto)
    assert cumplen == len(acto["nu"]) - sum(acto["cumple"]), "los conteos de dentro y fuera no casan"
    _exige(_readme(), f"{cumplen} de {len(acto['nu'])}", "el README")


def test_la_racha_mas_larga_no_alcanza_los_consecutivos_exigidos():
    acto = _acto3()
    exigidos = acto["criterio"]["bins_consecutivos"]
    falla = [abs(dif) > banda for dif, banda in zip(acto["dif"], acto["banda"])]
    assert acto["racha_mas_larga"] == _racha_mas_larga(falla), "`racha_mas_larga` no deriva de la banda"
    assert acto["racha_mas_larga"] < exigidos, (
        f"la racha ({acto['racha_mas_larga']}) alcanzo los {exigidos} exigidos: "
        "el veredicto NEGATIVO ya no se sostiene"
    )
    for inicio, fin in acto["rachas"]:
        assert fin - inicio + 1 <= acto["racha_mas_larga"] + 1, f"tramo {inicio}..{fin} fuera de la racha"
    texto = _readme()
    _exige(texto, f"racha_mas_larga={acto['racha_mas_larga']}", "el README")
    _exige(texto, f"exigidos={exigidos}", "el README")
    _exige(texto, acto["veredicto"], "el README")


def test_el_bin_que_falla_es_el_borde_inferior_y_sus_cifras():
    acto = _acto3()
    falla = [i for i, (dif, banda) in enumerate(zip(acto["dif"], acto["banda"])) if abs(dif) > banda]
    assert falla == [0], f"el unico bin fuera de banda dejo de ser el borde inferior: {falla}"
    texto = _readme()
    _exige(texto, repr(acto["nu"][0]), "el README")
    _exige(texto, f"{abs(acto['dif'][0]):.6f}", "el README")
    _exige(texto, f"{acto['banda'][0]:.6f}", "el README")
    assert abs(acto["dif"][0]) > acto["banda"][0], "el bin senalado no incumple la banda"


def test_q1_a_q4_en_verde_lo_dicen_los_numeros_no_la_etiqueta():
    """Los cuatro booleanos del crudo y las cifras que los sostienen."""
    c1 = _c1()
    for campo in ("Q1_riz_ok", "Q2_delta_ok", "Q3_reparto_ok", "Q4_fase_ok"):
        assert c1[campo] is True, f"el crudo de C1 declaro {campo} fuera de verde"
    umbral = c1["umbral_planitud"]
    assert max(c1["riz_reparto_por_modo"].values()) <= umbral, "Q1: el riz del reparto se paso del umbral"
    assert max(c1["delta_lnA_36_vs_24"].values()) <= umbral, "Q2: el borde lejano cambio la medida"
    for corrida in c1["corridas"]:
        assert corrida["reparto_desv_max"] <= umbral, f"Q3: desvio {corrida['reparto_desv_max']}"
        assert corrida["fase_error_max"] <= c1["factor_q4"] * c1["yee_error_fase_max_ventana"], (
            f"Q4: fase {corrida['fase_error_max']} contra {c1['factor_q4']} x Yee"
        )


def test_los_resumenes_del_crudo_casan_con_sus_corridas():
    """`riz_reparto_por_modo` y `delta_lnA_36_vs_24` son el peor caso de sus corridas."""
    c1 = _c1()
    for modo in ("C1a", "C1b"):
        riz = [c["reparto_riz"] for c in c1["corridas"] if c["modo"] == modo]
        assert riz, f"el crudo no contiene corridas {modo}"
        assert abs(c1["riz_reparto_por_modo"][modo] - max(riz)) < 1e-15, (
            f"{modo}: el peor riz declarado no es el de sus corridas"
        )
        corto = [c for c in c1["corridas"] if c["modo"] == modo and abs(c["L_a"] - 24.0) < 1e-12][0]
        largo = [c for c in c1["corridas"] if c["modo"] == modo and abs(c["L_a"] - 36.0) < 1e-12][0]
        medido = max(abs(b - a) for a, b in zip(largo["lnA_d_2"], corto["lnA_d_2"]))
        assert abs(c1["delta_lnA_36_vs_24"][modo] - medido) < 1e-15, (
            f"{modo}: delta_lnA_36_vs_24 no es el peor |dlnA| entre L = 24 a y L = 36 a"
        )


def test_las_cuatro_razones_de_reparto_estan_publicadas_y_tienen_portador():
    """Cada `reparto_media` de las cuatro corridas se publica con su redondeo exacto."""
    texto = _readme()
    for corrida in _c1()["corridas"]:
        _exige(texto, f"{corrida['reparto_media']:.6f}", "el README")
        assert abs(corrida["reparto_media"] - corrida["reparto_teorico"]) <= 0.05, (
            f"reparto {corrida['modo']} L = {corrida['L_a']}: fuera de 0.05 del teorico"
        )


def test_las_dos_cifras_del_piso_de_ruido_son_el_redondeo_de_su_crudo():
    c1 = _c1()
    fase = c1["fase_error_max_C1a"]
    assert fase == _corrida("C1a", 24.0)["fase_error_max"], (
        "`fase_error_max_C1a` dejo de ser el de la corrida C1a con L = 24 a"
    )
    yee = c1["yee_error_fase_max_ventana"]
    texto = _readme()
    _exige(texto, f"{fase:.3e}", "el README")
    _exige(texto, f"{c1['factor_q4'] * yee:.3e}", "el README")
    _exige(texto, f"{_acto3()['sigma_u2'][0]:.6f}", "el README")


def test_el_readme_ancla_la_identidad_del_deposito():
    texto = _readme()
    _exige(texto, DOI_VERSION, "el README")
    _exige(texto, DOI_CONCEPTO, "el README")

    cff = CITATION.read_text(encoding="utf-8")
    zenodo = json.loads(ZENODO.read_text(encoding="utf-8"))
    titulo_cff = PATRON_TITULO_CFF.search(cff)
    assert titulo_cff, "CITATION.cff perdio su titulo"
    titulo = titulo_cff.group(1)
    assert zenodo.get("title") == titulo, "el titulo de .zenodo.json no es el de CITATION.cff"
    assert titulo in texto, "el README dejo de citar el titulo verbatim de CITATION.cff"

    version_cff = PATRON_VERSION_CFF.search(cff)
    version_zenodo = PATRON_VERSION_JSON.search(ZENODO.read_text(encoding="utf-8"))
    assert version_cff and version_zenodo, "falta la version en el empaque"
    assert version_cff.group(1) == version_zenodo.group(1), (
        f"versiones desalineadas: CITATION {version_cff.group(1)} vs .zenodo.json {version_zenodo.group(1)}"
    )


def test_el_rojo_se_conserva():
    """Quitar lo rojo del escaparate es el unico cambio que este repo no admite."""
    texto = _readme()
    for marca in ("ROJO x4", "92834c4", f"{_acto3()['nu'][0]:.4f}"):
        _exige(texto, marca, "el README")
    primera = CRUDOS / "crudo_C1_20261012T055738Z.txt"
    assert primera.is_file(), "se perdio el crudo de la primera corrida (la roja)"
    rojo = primera.read_text(encoding="utf-8")
    assert rojo.count("ROJO") >= 4, "el crudo rojo perdio sus cuatro criterios en ROJO"
    _exige(rojo, f"{_c1()['yee_error_fase_max_ventana']:.3e}", "el crudo rojo")
