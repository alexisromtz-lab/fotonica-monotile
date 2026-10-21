"""Testigo de `.github/claims.json`: el contrato tiene que medir, no narrar.

HUELLA = "HUELLA-FOTONICA-CLAIMS-WITNESS-20261020"

Cada claim declara cinco cosas y este testigo las mide todas:

  · `dato`     — el portador: tiene que existir como archivo y parsear como JSON;
  · `deriva`   — la expresion que produce la cifra: tiene que leerse del portador
                 (`dato` en la raiz) y evaluar al numero que la prosa publica;
  · `afirma`   — la cifra publicada, con su `polaridad` y su `tolerancia`;
  · `total`    — (opcional) la expresion del denominador, cuando la prosa publica `N de M`;
  · `prosa`    — el documento donde `afirma` aparece literalmente.

La derivacion NO se evalua con `eval`: se interpreta un arbol `ast` con una lista
blanca de nodos (atributos sobre `dato`, `sum`, `len`, `count`, aritmetica basica).
Cualquier otra cosa es `EXPRESION_INVALIDA`, no un silencio.

El testigo sabe fallar: `test_el_testigo_de_claims_sabe_fallar` inyecta cada
mentira sobre una copia en memoria y exige su causa exacta.
"""

from __future__ import annotations

import ast
import copy
import json
import re
from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parents[1]
CLAIMS = RAIZ / ".github" / "claims.json"

CAMPOS = ("id", "afirma", "dato", "deriva", "polaridad", "prosa")
POLARIDADES = ("exacto", "complemento", "cerca_de")
NUMERO = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
FUNCIONES = {
    "sum": sum,
    "len": len,
    "count": lambda secuencia, valor: sum(1 for x in secuencia if x == valor),
}
EPS = 1e-9


class DerivaInvalida(Exception):
    """La derivacion declarada no se puede evaluar contra el portador."""


class DerivaNoLeeElDato(DerivaInvalida):
    """La derivacion no toca el portador: mide una constante disfrazada."""


def cargar_claims(ruta: Path = CLAIMS) -> dict:
    """Lee el contrato; falla cerrado si no es JSON."""
    texto = Path(ruta).read_text(encoding="utf-8")
    try:
        objeto = json.loads(texto)
    except json.JSONDecodeError as exc:
        raise ValueError(f"CLAIMS_INVALIDO: {Path(ruta).name} no parsea ({exc.lineno}:{exc.colno})") from exc
    if not isinstance(objeto, dict):
        raise ValueError(f"CLAIMS_INVALIDO: la raiz es {type(objeto).__name__}, no un objeto")
    return objeto


def _leer_campo(objeto, clave: str):
    if isinstance(objeto, dict):
        if clave not in objeto:
            raise DerivaInvalida(f"'{clave}' no existe en el portador")
        return objeto[clave]
    raise DerivaInvalida(f"'{clave}' se pide sobre {type(objeto).__name__}")


def _evalua(nodo, dato):
    if isinstance(nodo, ast.Constant):
        return nodo.value
    if isinstance(nodo, ast.Name):
        if nodo.id == "dato":
            return dato
        raise DerivaInvalida(f"nombre no permitido '{nodo.id}'")
    if isinstance(nodo, ast.Attribute):
        return _leer_campo(_evalua(nodo.value, dato), nodo.attr)
    if isinstance(nodo, ast.Subscript):
        corte = nodo.slice
        if type(corte).__name__ == "Index":
            corte = corte.value
        indice = _evalua(corte, dato)
        base = _evalua(nodo.value, dato)
        try:
            return base[indice]
        except (IndexError, KeyError, TypeError) as exc:
            raise DerivaInvalida(f"indice {indice!r} fuera del portador") from exc
    if isinstance(nodo, ast.Call):
        if not isinstance(nodo.func, ast.Name) or nodo.func.id not in FUNCIONES:
            nombre = getattr(nodo.func, "id", type(nodo.func).__name__)
            raise DerivaInvalida(f"funcion no permitida '{nombre}'")
        if nodo.keywords:
            raise DerivaInvalida("argumentos con nombre no permitidos")
        return FUNCIONES[nodo.func.id](*[_evalua(arg, dato) for arg in nodo.args])
    if isinstance(nodo, ast.BinOp) and isinstance(nodo.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
        izquierda, derecha = _evalua(nodo.left, dato), _evalua(nodo.right, dato)
        if isinstance(nodo.op, ast.Add):
            return izquierda + derecha
        if isinstance(nodo.op, ast.Sub):
            return izquierda - derecha
        if isinstance(nodo.op, ast.Mult):
            return izquierda * derecha
        try:
            return izquierda / derecha
        except ZeroDivisionError as exc:
            raise DerivaInvalida("division por cero") from exc
    if isinstance(nodo, ast.UnaryOp) and isinstance(nodo.op, (ast.USub, ast.UAdd)):
        valor = _evalua(nodo.operand, dato)
        return -valor if isinstance(nodo.op, ast.USub) else +valor
    raise DerivaInvalida(f"nodo no permitido '{type(nodo).__name__}'")


def evaluar(expresion: str, dato) -> float:
    """Evalua la derivacion declarada; exige que lea `dato` y devuelva numero."""
    try:
        arbol = ast.parse(expresion, mode="eval").body
    except SyntaxError as exc:
        raise DerivaInvalida(f"no parsea: {exc.msg}") from exc
    if not any(isinstance(nodo, ast.Name) and nodo.id == "dato" for nodo in ast.walk(arbol)):
        raise DerivaNoLeeElDato(f"'{expresion}' no lee el portador")
    valor = _evalua(arbol, dato)
    if isinstance(valor, bool) or not isinstance(valor, (int, float)):
        raise DerivaInvalida(f"devuelve {type(valor).__name__}, no un numero")
    return valor


def cifras(afirma: str) -> list[float]:
    return [float(hallado) for hallado in NUMERO.findall(afirma)]


def sostiene(polaridad: str, derivado: float, afirmado: float, total, tolerancia) -> bool:
    if polaridad == "exacto":
        return abs(derivado - afirmado) <= EPS
    if polaridad == "complemento":
        return total is not None and abs((total - derivado) - afirmado) <= EPS
    if polaridad == "cerca_de":
        return tolerancia is not None and abs(derivado - afirmado) <= tolerancia
    return False


def inconsistencias(objeto: dict, raiz: Path = RAIZ) -> list[str]:
    """Devuelve la lista de mentiras del contrato; vacia si el contrato mide."""
    fallos: list[str] = []
    if not isinstance(objeto, dict):
        return [f"CLAIMS_INVALIDO: la raiz es {type(objeto).__name__}, no un objeto"]
    if objeto.get("clase") != "A":
        fallos.append(f"CLAIMS_NO_MECANIZABLE: clase={objeto.get('clase')!r} no declara claims medibles")
    claims = objeto.get("claims")
    if not isinstance(claims, list) or not claims:
        return fallos + ["CLAIMS_VACIO: 'claims' ausente o vacio"]
    for claim in claims:
        fallos.extend(_inconsistencias_de_claim(claim, raiz))
    return fallos


def _inconsistencias_de_claim(claim, raiz: Path) -> list[str]:
    if not isinstance(claim, dict):
        return [f"CLAIM_INCOMPLETO: {type(claim).__name__} no es un objeto"]
    faltan = [campo for campo in CAMPOS if campo not in claim]
    if faltan:
        return [f"CLAIM_INCOMPLETO: '{claim.get('id', '?')}' sin {', '.join(faltan)}"]

    cid, rel = claim["id"], claim["dato"]
    fallos: list[str] = []

    expresion = claim["deriva"]
    if expresion.startswith("expr:"):
        partes = expresion.split(":", 2)
        if len(partes) != 3:
            return [f"EXPRESION_INVALIDA: {cid} '{expresion}' no es expr:<ruta>:<expresion>"]
        _, ruta_declarada, expresion = partes
        if ruta_declarada.lstrip("./") != rel.lstrip("./"):
            fallos.append(f"DERIVA_NO_APUNTA_AL_DATO: {cid} deriva lee '{ruta_declarada}' y dato declara '{rel}'")

    portador = (raiz / rel).resolve()
    if not portador.is_file():
        return fallos + [f"DATO_SIN_PORTADOR: {cid} '{rel}' no existe"]
    try:
        dato = json.loads(portador.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return fallos + [f"DATO_ILEGIBLE: {cid} '{rel}' no parsea ({exc.lineno}:{exc.colno})"]

    total = None
    try:
        derivado = evaluar(expresion, dato)
        if "total" in claim:
            total = evaluar(claim["total"], dato)
    except DerivaNoLeeElDato as exc:
        return fallos + [f"DERIVA_NO_LEE_EL_DATO: {cid} {exc}"]
    except DerivaInvalida as exc:
        return fallos + [f"EXPRESION_INVALIDA: {cid} '{expresion}': {exc}"]

    numeros = cifras(claim["afirma"])
    if not numeros:
        return fallos + [f"AFIRMA_SIN_CIFRA: {cid} '{claim['afirma']}' no publica numero"]
    afirmado = numeros[0]
    if len(numeros) > 1:
        if total is None:
            fallos.append(f"TOTAL_SIN_DERIVA: {cid} afirma un total ({numeros[1]:g}) y no declara 'total'")
        elif abs(total - numeros[1]) > EPS:
            fallos.append(f"TOTAL_DESALINEADO: {cid} deriva {total:g} y afirma {numeros[1]:g}")

    polaridad = claim["polaridad"]
    if polaridad not in POLARIDADES:
        fallos.append(f"POLARIDAD_DESCONOCIDA: {cid} '{polaridad}'")
    elif not sostiene(polaridad, derivado, afirmado, total, claim.get("tolerancia")):
        fallos.append(
            f"DERIVACION_CIEGA: {cid} deriva {derivado!r} no sostiene afirma '{claim['afirma']}' ({polaridad})"
        )

    prosa = (raiz / claim["prosa"]).resolve()
    if not prosa.is_file():
        fallos.append(f"PROSA_AUSENTE: {cid} '{claim['prosa']}' no existe")
    elif claim["afirma"] not in prosa.read_text(encoding="utf-8"):
        fallos.append(f"README_DESALINEADO: {cid} '{claim['afirma']}' no aparece en {claim['prosa']}")
    return fallos


def test_el_claims_declara_solo_lo_que_el_parser_mide():
    objeto = cargar_claims()
    fallos = inconsistencias(objeto)
    assert not fallos, "el claims.json declara lo que su aparato no mide:\n  " + "\n  ".join(fallos)
    print(f"CLAIMS_SANO: {len(objeto['claims'])} claims medidos contra su portador")


def test_el_testigo_de_claims_sabe_fallar(tmp_path):
    """El falsador del testigo: cada mentira tiene que producir su causa exacta."""
    base = cargar_claims()
    casos: list[tuple[str, dict, str]] = []

    roto = copy.deepcopy(base)
    roto["claims"][0]["dato"] = "crudos/no_existe_20261020.json"
    casos.append(("portador ausente", roto, "DATO_SIN_PORTADOR"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["deriva"] = "sum(dato.campo_que_no_existe)"
    casos.append(("atributo ausente", roto, "EXPRESION_INVALIDA"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["deriva"] = "sum([1, 2, 3])"
    casos.append(("deriva sin portador", roto, "DERIVA_NO_LEE_EL_DATO"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["deriva"] = f"{base['claims'][0]['deriva']} + open('x').read()"
    casos.append(("deriva con funcion prohibida", roto, "EXPRESION_INVALIDA"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["afirma"] = "999 de 999"
    casos.append(("cifra movida", roto, "DERIVACION_CIEGA"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["afirma"] = "sin numeros"
    casos.append(("afirma sin cifra", roto, "AFIRMA_SIN_CIFRA"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["polaridad"] = "mas_o_menos"
    casos.append(("polaridad inventada", roto, "POLARIDAD_DESCONOCIDA"))

    roto = copy.deepcopy(base)
    del roto["claims"][0]["deriva"]
    casos.append(("campo ausente", roto, "CLAIM_INCOMPLETO"))

    roto = copy.deepcopy(base)
    roto["clase"] = "B"
    casos.append(("clase no mecanizable", roto, "CLAIMS_NO_MECANIZABLE"))

    roto = copy.deepcopy(base)
    roto["claims"] = []
    casos.append(("claims vacio", roto, "CLAIMS_VACIO"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["prosa"] = "PROSA_QUE_NO_EXISTE.md"
    casos.append(("prosa ausente", roto, "PROSA_AUSENTE"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["deriva"] = f"expr:otra/ruta.json:{roto['claims'][0]['deriva']}"
    casos.append(("expr con ruta ajena", roto, "DERIVA_NO_APUNTA_AL_DATO"))

    roto = copy.deepcopy(base)
    roto["claims"][0]["total"] = "len(dato.campo_que_no_existe)"
    casos.append(("total que no evalua", roto, "EXPRESION_INVALIDA"))

    ilegible = tmp_path / "crudo_ilegible.json"
    ilegible.write_text("{ esto no es json", encoding="utf-8")
    roto = copy.deepcopy(base)
    roto["claims"][0]["dato"] = str(ilegible)
    casos.append(("portador ilegible", roto, "DATO_ILEGIBLE"))

    for nombre, objeto, causa in casos:
        fallos = inconsistencias(objeto)
        assert any(fallo.startswith(causa) for fallo in fallos), f"{nombre}: se esperaba {causa}, hubo {fallos}"

    malformado = tmp_path / "claims_malformado.json"
    malformado.write_text('{"huella": "x", "claims": [', encoding="utf-8")
    with pytest.raises(ValueError) as error:
        cargar_claims(malformado)
    assert str(error.value).startswith("CLAIMS_INVALIDO"), str(error.value)


def test_el_testigo_acepta_la_deriva_prefijada():
    """`expr:<dato>:<expresion>` es la forma explicita: se acepta sin reescribir el contrato."""
    base = cargar_claims()
    for claim in base["claims"]:
        prefijado = copy.deepcopy(base)
        prefijado["claims"] = [copy.deepcopy(claim)]
        prefijado["claims"][0]["deriva"] = f"expr:{claim['dato']}:{claim['deriva']}"
        assert not inconsistencias(prefijado), inconsistencias(prefijado)
