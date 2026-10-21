<!-- HUELLA-README-FOTONICA-MONOTILE-20261013 -->

# fotonica-monotile — Instrumento FDTD 2D TM y acto de contraste del monotilo Hat

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23329943-blue)](https://doi.org/10.5281/zenodo.23329943)
[![Licencia: MIT](https://img.shields.io/badge/licencia-MIT-green)](LICENSE)
[![Motor](https://img.shields.io/badge/motor-python%203.14.7%20%C2%B7%20numpy%202.4.6-blue)](PREREGISTRO_INSTRUMENTO_FDTD_20261012.md)
[![Criterios Q1..Q4](https://img.shields.io/badge/Q1--Q4-VERDE-brightgreen)](crudos/crudo_C1_20261012T062648Z.json)
[![Acto 3.2](https://img.shields.io/badge/acto%203.2-NEGATIVO%20%2F%20NULO%20SOSTENIDO-lightgrey)](crudos/crudo_acto3_20261012T205230Z.json)

**Gancho.** Instrumento FDTD 2D TM (esquema Yee) y su acto de contraste del monotilo Hat: la pieza publica su piso de ruido como cota (Q1..Q4 VERDE) y su contraste como **resultado negativo** — `NEGATIVO / NULO SOSTENIDO` — sostenido por `28 de 29` bins dentro de la banda `3 sigma + riz` y una racha de fallo de `1` contra `3` exigidos. Este repositorio es la **pieza fuente** del enlace que consume `acustica-monotile`.

*2D TM Yee-FDTD instrument and its Hat-monotile contrast act: the piece publishes its noise floor as a bound (Q1..Q4 GREEN) and its contrast as a **negative result** — `NEGATIVE / SUSTAINED NULL` — held by `28 of 29` bins inside the `3 sigma + riz` band and a failure run of `1` against `3` required. This repository is the **source piece** of the link consumed by `acustica-monotile`.*

**Nota sobre las insignias.** Cada insignia apunta a su artefacto. La de motor cita la firma **declarada** del entorno (`PREREGISTRO_INSTRUMENTO_FDTD_20261012.md`, §"Motor": `numpy 2.4.6` · `python 3.14.7` · unidad de conteo = celda), no una promesa de compatibilidad. **No se acuña insignia de orden de convergencia en `p`**: este repositorio no declara `p = 2` ni una escalera de convergencia medida (búsqueda de `segundo orden`, `orden 2`, `O(h^2)` en el pre-registro: **0** líneas). La insignia `Q1..Q4-VERDE` la respalda el crudo `crudos/crudo_C1_20261012T062648Z.json` (campos `Q1_riz_ok`, `Q2_delta_ok`, `Q3_reparto_ok`, `Q4_fase_ok` = `True`); la del acto 3.2, `crudos/crudo_acto3_20261012T205230Z.json` (campo `veredicto`).


Instrumento de simulación por diferencias finitas en el dominio del tiempo (FDTD, esquema Yee
2D polarización TM) y su **acto de contraste** sobre la geometría del monotilo Hat. Este
repositorio es la **pieza fuente** del enlace que consume `acustica-monotile`.

Todo lo que sigue está leído de artefactos sellados de este repositorio: el pre-registro del
instrumento y los crudos de `crudos/`. **Ninguna cifra de este README se deriva de nuevo**: se
cita con su artefacto, su unidad y su hora.

---

## Inicio rápido (15 s)

Solo biblioteca estándar: se **lee** el crudo sellado, no se corre FDTD (una corrida completa tarda minutos y vive en §6).

```sh
git clone https://github.com/alexisromtz-lab/fotonica-monotile.git
cd fotonica-monotile
python3 -B - <<'PY'
import json
d = json.load(open('crudos/crudo_acto3_20261012T205230Z.json'))
viola = [abs(x) > b for x, b in zip(d['dif'], d['banda'])]
print('UT_HORA=%s  rc=%s' % (d['t_utc'], d['rc']))
print('bins_que_cumplen=%d de %d  racha_mas_larga=%d  exigidos=%d  veredicto=%s'
      % (len(viola) - sum(viola), len(viola), d['racha_mas_larga'],
         d['criterio']['bins_consecutivos'], d['veredicto']))
PY
```

Salida (verbatim, medida el `2026-10-13`):

```text
UT_HORA=2026-10-12T20:52:30Z  rc=0
bins_que_cumplen=28 de 29  racha_mas_larga=1  exigidos=3  veredicto=NEGATIVO / NULO SOSTENIDO
```

Notas de instrumento declaradas: los scripts de esta pieza viven en la **raíz** del repositorio (`malla_hat.py`, `fdtd2d_tm.py`, `inyector_diel.py`); el directorio `sim/` está **vacío** y **`sim/run_c1.py` no existe** aquí (**0** hallazgos en el árbol). Los comandos de corrida completa, con sus banderas, están en §6.

## Resultados medidos

UNIVERSO DECLARADO: los dos crudos sellados de `crudos/` (`crudo_C1_20261012T062648Z.json` · `crudo_acto3_20261012T205230Z.json`), leídos con `json` de la biblioteca estándar; unidad de conteo = celda; hora de cada cifra = la del crudo que la contiene. Ventana espectral `nu ∈ [0.15, 0.55]` c/a, 29 bins.

| Criterio | Medido | Umbral | Clase |
|---|---|---|---|
| Q1 · riz del reparto, peor caso por modo (`riz_reparto_por_modo`) | C1a `0.005394` · C1b `0.005523` | `0.05` | VERDE |
| Q2 · `delta_lnA_36_vs_24` (el borde lejano no cambia la medida) | C1a `0.004548` · C1b `0.006001` | `0.05` | VERDE |
| Q3 · `A(16a)/A(8a)` (`reparto_media`) · C1a, `L = 24 a` | `0.999907` (teórico `1.000000`) | en `C1a` el teórico es `1.000` | VERDE |
| Q3 · `A(16a)/A(8a)` · C1a, `L = 36 a` | `1.000055` | teórico `1.000` | VERDE |
| Q3 · `A(16a)/A(8a)` · C1b, `L = 24 a` | `0.707449` (teórico `0.707107`) | `0.7071 ± 0.05` | VERDE |
| Q3 · `A(16a)/A(8a)` · C1b, `L = 36 a` | `0.707223` | `0.7071 ± 0.05` | VERDE |
| Q4 · `fase_error_max_C1a` | `1.405e-03` | `3 × yee_error_fase_max_ventana` = `4.389e-03` (rama `m = 1`) | VERDE |
| Acto 3.2 · bins que cumplen la banda `3 sigma + riz` | `28 de 29` (falla el bin `nu = 0.1527114967462039`) | `29` | 1 bin fuera |
| Acto 3.2 · `racha_mas_larga` | `1` | `>= 3` (`bins_consecutivos`) | NULO |
| Acto 3.2 · `veredicto` (campo del crudo) | `NEGATIVO / NULO SOSTENIDO` | — | `rc 0` |

Lectura acotada a lo medido: los cuatro criterios de entrada (`Q1_riz_ok`, `Q2_delta_ok`, `Q3_reparto_ok`, `Q4_fase_ok`) están en `True` en el crudo de C1 y el acto 3.2 **no** se declara verde: `28 de 29` bins caen dentro de la banda declarada antes de medir (`k_sigma = 3.0`, `riz_T_C1 = 0.05`) y la racha más larga de incumplimiento es `1`, contra `3` consecutivos exigidos. El bin que falla es `nu = 0.1527114967462039` (`|dif| = 1.584702` contra banda `0.847124`). El veredicto del crudo es `NEGATIVO / NULO SOSTENIDO` y se publica como resultado.

Dos advertencias de lectura, declaradas:

1. **El campo llamado `cumple` guarda la NEGACIÓN del criterio.** Medido: `len(campo) = 29`, `sum(campo) = 1`, índices `True` = `[0]`; derivado `abs(dif) <= banda` da `28` verdaderos con índice `False` = `[0]`; y `campo == NOT derivado` → `True`. Es decir: el array marca el bin que **falla**, y es lo que consumen `rachas`/`racha_mas_larga`. Leído literalmente, §4.2 queda en contradicción con §4.1 (`28 de 29` cumplen); leído como negación, ambas coinciden. **No se edita §4: se declara aquí.** El veredicto no cambia (`1 < 3`).
2. **Identidad de árbol del acto 3.2:** el crudo registra `git_head = 7ffbcfa` y el `HEAD` actual de este repositorio es otro; ya está declarado en §4.4 y no se re-deriva aquí.

<!-- HUELLA-A21-STANDARD-README-FOTONICA-20261013 -->

## Epistemology & Rigor (The A2.1 Standard)

> This repository operates under strict pre-registration. No parameters are tuned to fit the data post-execution. Analytical predictions (such as the 3σ bounding for stochastic noise) are sealed before running. Red/Null results are preserved and declared, guaranteeing immunity to confirmation bias and making the code entirely verifiable by both human auditors and AI agents.

Los tres axiomas del estándar, cada uno con su anclaje medido en esta casa:

1. **Predicción antes del motor.** El pre-registro (`PREREGISTRO_INSTRUMENTO_FDTD_20261012.md`, blob `23496c0592fa32c9e2ff7a9e4d5f59a6bf88a22b`) se sella a `T_UTC 2026-10-12T05:49Z` (l.4) con los cuatro criterios y sus umbrales (l.193–198, `cero figuras vistas`); la primera corrida es posterior (`055738Z`): el umbral es anterior a la salida, no se elige mirándola.
2. **El rojo se conserva.** La primera corrida salió **ROJO x4** y vive en el repositorio con su commit: `crudos/crudo_C1_20261012T055738Z.txt:33–36`, commit `92834c4`, inscrita en el libro del pre-registro (l.260, `t > 0.75 T`). Lo que la volteó fue **el instrumento** —absorbedor tipo esponja y fila `j=0` inerte, declaradas en `fdtd2d_tm.py:58`—, nunca el umbral.
3. **Anclaje sin institución.** Sello Zenodo `10.5281/zenodo.23329943` (concepto `…23329942`) y sello del instrumento `14bc745` (`v0.1.0`): identidad y re-ejecutabilidad verificables por cualquiera —persona o agente— sin institución en medio. El anclaje no certifica corrección: certifica identidad y re-ejecutabilidad.


## Cómo citar (BibTeX)

```bibtex
@software{roman_martinez_fotonica_monotile_2026,
  author   = {Roman Martinez, Jose Alexis},
  title    = {fotonica-monotile: instrumento FDTD 2D TM y acto de contraste del monotilo Hat (Paso 2)},
  year     = {2026},
  version  = {0.1.0},
  doi      = {10.5281/zenodo.23329943},
  url      = {https://github.com/alexisromtz-lab/fotonica-monotile},
  license  = {MIT},
  keywords = {monotile Hat, aperiodicidad, FDTD, Yee 2D TM, metrologia A2.1, contraste dielectrico}
}
```

El título del `@software` es el de `CITATION.cff` (verbatim). **No se declara ORCID: no existe ninguno en este ecosistema.** DOI de concepto (todas las versiones): `10.5281/zenodo.23329942`.

---
## 1. El instrumento y su pre-registro sellado

`PREREGISTRO_INSTRUMENTO_FDTD_20261012.md` — sellado en `T_UTC 2026-10-12T05:49Z`. Declara el
instrumento **antes** de correrlo:

| Pieza | Valor declarado |
|---|---|
| Esquema | Yee 2D, polarización TM (`Ez`, `Hx`, `Hy`) |
| Paso de malla | `dx = dy = a/16` |
| Número de Courant | `S = 0.5` |
| Frontera absorbente | esponja de 20 celdas por borde, perfil `borde^3`, `alfa = 0.30` |
| Fuente | gaussiana, `sigma = 0.5 a/c`, retardo `t0 = 3 sigma` |
| Sondas | absolutas en `d = 8a` y `d = 16a` |
| Ventana de frecuencias | `nu ∈ [0.15, 0.55]` c/a (§2.b del pre-registro) |
| Sub-corridas | **C1a** cuasi-1D (`ny = 1`, esponja solo en x) y **C1b** 2D real (`ny = nx`, esponja en los cuatro bordes) |

Criterios de la pieza, declarados con su umbral (l.195–198 del pre-registro):

- **Q1** — `riz <= 0.05` en C1a **y** en C1b
- **Q2** — `delta = max_nu |ln A(36a) − ln A(24a)| <= 0.05` en C1a **y** en C1b
- **Q3** — `A(16a)/A(8a) = 0.7071 ± 0.05` en C1b; en C1a el valor teórico es `1.000`
- **Q4** — `max |error de fase medido| <= 3 · max |error de Yee|` en la ventana, medido en C1a

## 2. Veredicto de la pieza (acto 1, gen 6, medido)

Del bloque `§2.c.1` del pre-registro — **medido**, no predicho:

| Criterio | Veredicto | Valores impresos |
|---|---|---|
| Q1 | **VERDE** | `[C1a 0.005394 · C1b 0.005523]` |
| Q2 | **VERDE** | `[0.004548 · 0.006001]` |
| Q3 | **VERDE** | `[C1a/24a 0.999907 · C1a/36a 1.000055 · C1b/24a 0.707449 · C1b/36a 0.707223]` |
| Q4 | **VERDE** | `fase medida 1.405e-03 <= 3 x Yee 1.463133e-03 = 4.389e-03` (rama `m = 1`) |

Con Q1..Q4 en verde, la condición de entrada al acto 2 del pre-registro queda **cumplida**
(habilitado, no corrido).

## 3. Acto 3.2 — contraste del monotilo Hat sobre la pieza C1

Crudo sellado: `crudos/crudo_acto3_20261012T205230Z.json` (y su gemelo `.txt`), tomado a
`t_utc = 2026-10-12T20:52:30Z`, `acto = 3.2`, `rc = 0`.

**Universo declarado** (tal como el crudo lo imprime): geometría de C1 — `dx = a/16 = 0.0625`,
PML de 20 celdas por borde, fuente en `i = 32`, sondas en `8a` y `16a`, `T_total = 72 a/c`,
`Nt = 2304`, ventana `0.15..0.55` c/a. Objeto: parche de nivel 2 rasterizado a `24a`, con
**1 celda de máscara = 1 celda de Yee**; `N = 384` celdas; modo **C1b**; contraste `eps_mat 4`
sobre `eps_0 1`; espesor `L = 24 a`; semillas `(101, 102, 103, 104, 105)`. **Unidad de conteo:
celda.** Referencia del contraste: `crudo_C1_20261012T062648Z.json`, corrida 2, modo C1b,
`L_a = 24.0`.

**Criterio** (declarado en el propio crudo): `k_sigma = 3.0`, `riz_T_C1 = 0.05`,
`bins_consecutivos = 3`. Por bin: `|T_U1 − T_bar_u2| <= 3 · sigma_u2 + riz(T)_C1`.

**Veredicto sellado: `NEGATIVO / NULO SOSTENIDO`**, con `rc = 0`. La hipótesis —que la
correlación espacial del monotilo se distinga en transmisión a este contraste— **queda
refutada en este régimen de contraste**. Se publica como **resultado, no como fracaso**.

## 4. Declarado, no adjudicado

Estas cuatro observaciones se declaran tal como las mide el artefacto. Ninguna se repara aquí y
ninguna cambia el veredicto sellado:

1. **28 de 29 bins cumplen la banda; 1 no.** El bin que no cumple es el **borde inferior de la
   ventana**, `nu = 0.1527114967462039`:
   `T_U1 = 1.993564` · `T_bar_u2 = 0.408862` · `sigma_u2 = 0.265708` ·
   `|dif| = 1.584702` contra `banda = 0.847124`.
2. **El crudo ya marca exactamente ese bin.** Su arreglo `cumple` vale
   `[True, False x 28]` y su único tramo marcado es `rachas = [[0, 0]]`, con
   `racha_mas_larga = 1`. Como el criterio exige **3 bins consecutivos** y hay **1**,
   `1 < 3`: **el veredicto sellado se sostiene**. Lo que queda es una diferencia de **prosa**,
   no de medición — la frase impresa dice «en toda la ventana vale» mientras su propio arreglo
   imprime 28 de 29.
3. **Convención del piso.** El campo `criterio.riz_T_C1 = 0.05` es el **techo del criterio Q1**
   (`riz <= 0.05`), no el riz medido: el riz medido en C1b es `0.005523` (§2), que es además el
   piso que usa el repositorio hermano `acustica-monotile`. El bin del punto 1 excede la banda
   **bajo cualquiera de los dos pisos** (con `0.005523` la banda sería `0.802647` y
   `|dif| = 1.5847017` sigue por encima), y la racha más larga sigue siendo 1 en ambos casos.
4. **Identidad del árbol.** El crudo declara `git_head = 7ffbcfa`; el `HEAD` de este repositorio
   ha avanzado desde entonces (`8138315` al cerrar este empaque). La identidad de los scripts sí
   se reproduce: los campos `malla_hat_sha16` y `fdtd2d_tm_sha16` del crudo coinciden con el
   `sha256[0:16]` medido hoy sobre los archivos (§5).

## 5. Instrumentos y huellas

`sha256[0:16]`, medido el `2026-10-13` (motor: `shasum -a 256` sobre el archivo en disco):

| Archivo | Líneas | `sha256[0:16]` |
|---|---|---|
| `malla_hat.py` | 910 | `914a8666f1d13ddf` |
| `fdtd2d_tm.py` | 501 | `56a8d5f896b9cfea` |
| `inyector_diel.py` | 767 | `01ffcf31c040368d` |

Los dos primeros **coinciden** con los campos que el propio crudo `acto3` declara
(`malla_hat_sha16`, `fdtd2d_tm_sha16`): dos mediciones independientes del mismo objeto.

## 6. Reproducción

Desde la raíz del repositorio:

```bash
# 1. Geometría (malla del monotilo Hat); no escribe nada sin --salida
python3 malla_hat.py --nivel 2 --n 384 --semilla 20261012 --json

# 2. Pieza C1 del instrumento (Q1..Q4); --json imprime el crudo
python3 fdtd2d_tm.py --json

# 3. Ingesta dieléctrica: sin banderas SOLO construye y verifica; no corre FDTD ni escribe
python3 inyector_diel.py

# 4. Acto 3.2 (contraste): corre U1 + las 5 semillas y escribe un crudo nuevo con hora UTC
python3 inyector_diel.py --correr
```

`inyector_diel.py --correr` **escribe un crudo nuevo** en `crudos/` con la hora UTC del momento;
el crudo sellado `crudo_acto3_20261012T205230Z.json` **no se sobrescribe**. Para correr sin
escribir, añadir `--sin-crudo`.

Los tres scripts se resuelven desde su propio directorio, así que las corridas anteriores
funcionan con cualquier directorio de trabajo. `--json` está disponible en los tres.

## 7. Estructura

```
fotonica-monotile/
  PREREGISTRO_INSTRUMENTO_FDTD_20261012.md   pre-registro sellado + veredicto del acto 1
  fdtd2d_tm.py                               solver Yee 2D TM
  malla_hat.py                               geometria del monotilo Hat
  inyector_diel.py                           ingesta dielectrica + acto 3.2
  crudos/                                    6 crudos de la pieza C1 + el crudo del acto 3.2
```

Todos los crudos se versionan en par `.json` + `.txt`: el `.json` es el objeto legible por
máquina, el `.txt` su impresión para lectura humana.

## 8. Cita y licencia

- `CITATION.cff` — el formulario canónico de cita (esquema CFF 1.2.0).
- `.zenodo.json` — los metadatos de depósito; el campo `license` viaja como **cadena escalar**.
- Licencia: **MIT** (ver `LICENSE`); la elección se declara en `CITATION.cff`.
- Repositorio hermano: [`acustica-monotile`](https://github.com/alexisromtz-lab/acustica-monotile)
  — el enlace que consume esta geometría.
- Geometría: commit `b498fcf` — el mismo que cita `acustica-monotile`; el `sha16` del generador `malla_hat.py` está sellado en la tabla de §6. <!-- HUELLA-GEOM-BILATERAL-FOTONICA-20261020 -->
