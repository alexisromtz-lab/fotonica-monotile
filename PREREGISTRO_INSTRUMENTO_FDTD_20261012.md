# PRE-REGISTRO DEL INSTRUMENTO FDTD 2D Y DEL PISO C1
## fotonica-monotile · acto 1 · sellado ANTES de la primera cifra

T_UTC de sellado: 2026-10-12T05:49Z
Casa de este archivo: `/Users/marlene/fotonica-monotile` (repo propio, `git init` 2026-10-12T05:49:17Z, rama `main`).
Casa de custodia: `/Users/marlene/platino-labs` (sin arbol git: es registro, no repo).

---

## §0 · Cómo se lee este archivo

Toda cifra de aquí en adelante viaja con **comando · unidad · hora UTC**. Un `0` de una
búsqueda o de una sonda es **el cero del metro**, no una ausencia del objeto: se declara
`0 (metro)` y nunca `no existe`.

Este archivo se escribe **antes** de correr nada. Su función es la de §3.c del pre-registro
de la casa de custodia (`PREREGISTRO_FOTONICA_MONOTILE_20261011.md`, sha16 `47c7c2ef5028ee40`):
fijar umbral, ventana, unidad y definición **antes** de ver la salida, para que elegirlas
después no sea posible.

---

## §1 · Las dos casas (declaradas)

| casa | ruta | qué guarda | git |
|---|---|---|---|
| custodia | `/Users/marlene/platino-labs` | actas, registro, punteros sha16 | **no** (RC 14 del gate A2.1, por diseño) |
| código | `/Users/marlene/fotonica-monotile` | solver, mallas, sondas, crudos | **sí** |

Regla heredada: **lo que no se admite es dejarlo sin declarar.** La casa de código nace
antes del primer script, no después.

---

## §2 · Universo del instrumento (fijado antes de medir)

UNIVERSO DECLARADO:
- **Unidades normalizadas:** constante de red `a = 1` · velocidad `c = 1` · `eps0 = mu0 = 1`.
  Frecuencia reducida `nu = omega a / (2 pi c)`; con `a = c = 1`, `nu` es directamente la
  frecuencia en unidades de `c/a`.
- **Malla Yee 2D, modo TM (`Ez`, `Hx`, `Hy`)**, `dx = dy = a/16`. Paso `S = c·dt/dx = 0.5`
  (límite de von Neumann en 2D = `1/sqrt(2) = 0.7071`; `0.5` queda por debajo).
- **Esponja absorbente:** 20 celdas en los cuatro bordes, perfil `borde^3`, `alfa = 0.30`
  por paso. Sin PML: es una esponja **declarada**, y su reflexión es lo que C1 mide.
- **Fuente:** pulso gaussiano `exp(-(t - t0)^2 / (2 sigma^2))`, `sigma = 0.5 a/c`,
  `t0 = 3 sigma`. Se inyecta como **corriente blanda** (se suma) en una sola columna
  `i = 2a`. Espectro analítico declarado `|S(omega)| = exp(-sigma^2 omega^2 / 2)`.
- **Sondas:** dos columnas a distancia fija de la fuente, `d1 = 4a` y `d2 = 8a`, en la fila
  central. **La posición es absoluta**, no relativa al borde: al crecer `L` el dominio sólo
  se alarga hacia la derecha.
- **Longitudes del sistema:** `L = 24a` (384 celdas) y `L = 36a` (576 celdas), dominio
  cuadrado `L x L`.
- **Ventana de medida:** `nu` en `[0.05, 0.55]`.
- **Resolución espectral:** `T_total = 72 a/c`, `Nt = 2304` pasos, `dnu = 1/T_total = 0.0139`.
  Ventana y duración del pulso son **la misma decisión** (trampa 5 del pre-registro): el
  pulso `sigma = 0.5` entrega `-10 dB` hasta `nu = 0.683`, es decir `0.683/0.55 = 1.24x`
  el borde alto de la ventana.
- **Umbral de planitud (declarado AHORA, antes de ver nada):** `0.05` en `|delta ln A|`.
- Motor: `numpy 2.4.6` · `python 3.14.7` · unidad de conteo = **celda**.

---

## §3 · Qué mide C1 y con qué criterio

C1 es el **piso de ruido del instrumento**: dominio **homogéneo `eps_r = 1`, misma malla,
mismo paso, misma esponja, misma fuente, mismas sondas**. No hay objeto en C1.

C1 entrega tres números, y los tres se declaran antes:

1. **Rizado dentro de una L:** `riz(L) = max_nu |ln A(L,nu) - media_nu ln A(L,·)|`, con
   `A(L,nu) = |FFT(sonda)| / |S(omega)|`. Es lo que un medio homogéneo **no** debería
   producir. Criterio: `riz <= 0.05`.
2. **Dependencia espuria con L:** `delta(L) = max_nu |ln A(36a,nu) - ln A(24a,nu)|`. Un medio
   homogéneo no cambia porque el borde esté más lejos, salvo por la reflexión de la esponja.
   **Este es el término que imita a la localización** y por eso es el número decisivo.
   Criterio: `delta <= 0.05`.
3. **Ley de reparto 2D (control positivo):** `A(d=8a)/A(d=4a)` debe valer `sqrt(4/8) = 0.7071`
   dentro de la ventana. Un control que **no puede diferir** no es evidencia; este sí puede
   fallar (si el solver no es Maxwell 2D, el reparto no sale `1/sqrt(r)`).

Y una cuarta medida, de cola: **cola de sonda** = `max|sonda|` para `t > 0.75·T_total`
dividido por el pico de la sonda. Es el nivel de reflexión que un decaimiento falso usaría
como disfraz.

**Regla dura (del pre-registro §4.1):** un decaimiento sólo cuenta **por encima de este
piso**. Si C1 no es plano, **se para y se repara el instrumento antes de medir el objeto** —
no se mide U1 con un instrumento sucio.

---

## §4 · Acompañante analítico (no es una corrida: es una fórmula)

La malla de Yee tiene una relación de dispersión propia:

```
(1/dt)^2 · sin^2(omega dt / 2) = (1/dx)^2 · [ sin^2(kx dx / 2) + sin^2(ky dy / 2) ]
```

Sobre el eje (`ky = 0`) y con `S = c dt/dx`, queda `sin(omega dt/2) = S · sin(k dx/2)`.
La velocidad de fase numérica `v_p = omega/k` se desvía de `c` y **esa desviación es el
piso de fase declarable**. Se imprime para `nu` en la ventana, y sirve de techo: si el
error de fase medido lo excede mucho, el rizado no es dispersión, es otra cosa.

---

## §5 · Las trampas, aplicadas a C1

1. **El piso numérico se disfraza de localización.** Por eso C1 existe y va **primero**.
2. **Error de escalera** en bordes dieléctricos: en C1 no hay dieléctrico; la trampa queda
   **abierta y declarada** para el acto siguiente (malla del *Hat*), donde se decide entre
   suavizado sub-píxel de `eps` o cruce contra FEM en unas `omega` declaradas.
3. **Un «hueco» sin definición previa es un rótulo.** Aquí no se mide hueco: se mide planitud.
   La definición de hueco (intervalo sin estados por encima de un umbral a una resolución)
   se sella cuando exista una DOS que mirar.
4. **Los aproximantes periódicos convergen o no:** `n`/`n+1` se publican en pares. No aplica a C1.
5. **Resolución espectral y duración del pulso son la misma decisión.** Tratada en §2.
6. **Doble precisión y CFL se declaran** (§2); **una frontera absorbente mal sintonizada
   fabrica decaimiento** ⇒ por eso C1 mide la cola de sonda.
7. **El instrumento se cuenta a sí mismo.** La ventana `[0.05, 0.55]` la elegimos nosotros:
   se declara, no se deduce del resultado.

---

## §6 · Predicción falsable de C1 (escrita antes de correr)

- **Q1.** `riz(L) <= 0.05` en la ventana para ambas `L`. Si falla, el instrumento tiene ruido
  propio por encima del umbral declarado.
- **Q2.** `delta <= 0.05`: el borde lejano **no** cambia la medida. Si `delta > 0.05`, entonces
  **cualquier** caída medida después con `L` es de la esponja, no del objeto ⇒ U1 y U2 quedan
  **vetados** hasta reparar.
- **Q3.** `A(8a)/A(4a) = 0.7071 +- 0.05` en la ventana. Si falla, el solver no reproduce el
  reparto cilíndrico 2D y no es Maxwell 2D válido en esta malla.
- **Q4.** El error de fase medido no excede `3x` el error de la fórmula de Yee de §4.

Un fallo de Q1..Q4 es un **resultado**, no un accidente: se publica como tal y detiene el
acto siguiente.

---

## §7 · Qué NO hace este acto (declarado, no insinuado)

- **No** construye la malla del *Hat* por inflación.
- **No** corre U1, U2, U3 ni C2. Ni un solo campo con dieléctrico.
- **No** toca `/Users/marlene/sfqed-paper` (ni copia: sólo se cita su núcleo como referencia
  de familia de instrumento; `sim/maxwell2d.py` sha16 `1db75ea7b0b8cba4`, 715 l, intacto).
- **No** afirma nada sobre el monotile. C1 no mide el objeto: mide **el metro**.

## §8 · Acto siguiente (uno, y ya está diseñado)

1. Malla del *Hat* por inflación, **con sus reglas de sustitución traídas de fuente y
   verificadas por matriz** (no de memoria), y control por construcción (conteo de fichas
   contra la matriz de sustitución).
2. U1 (Hat determinista) vs **U2 (Hat barajado: mismas fichas, misma geometría, misma
   fracción de llenado, orientaciones y posiciones barajadas)** — el nulo decisivo.
3. U3 (hexagonal periódica, mismas varillas, misma fracción de llenado).
4. C2 (cuasi­cristal conocido: Penrose) como control positivo del instrumento.
5. Sólo entonces: `T(omega)` vs `L` con ajuste exponencial contra ley de potencias en el
   mismo rango log, `IPR` vs `L`, estadística de espaciamiento, DOS y transporte anómalo.

**Condición de entrada al acto 2:** Q1..Q4 en verde.

---

## §2.b · ENMIENDA PRE-CORRIDA (T_UTC 2026-10-12T05:56Z · **cero figuras medidas todavía**)

La enmienda se escribe **antes de la primera corrida**. Su motivo es de diseño, no de resultado:
la §2 tal como quedó **no podía medir lo que dice medir**, y se declara por qué.

1. **Ventana de medida: `nu` en `[0.15, 0.55]`** (antes `[0.05, 0.55]`). El borde bajo no es
   libre: lo fija que la sonda más cercana quede a `>= 1` longitud de onda. Con `d = 8a` y
   `nu = 0.15`, `lambda = 6.67a` ⇒ `d = 1.20 lambda`. Por debajo de eso la sonda lee campo
   cercano y el «reparto cilíndrico» de Q3 no es una ley, es una vecindad.

2. **Sondas: posiciones absolutas `d = 8a` y `d = 16a`** (antes `4a` y `8a`). Razón idéntica:
   `k·d = 2 pi nu d / a` debe ser `>= 7.5` para que `H_0` esté en su régimen asintótico.

3. **Dos sub-corridas, no una.** El error de fase **no se puede medir con una fuente de columna
   finita**: en 2D la fase medida incluiría la corrección de campo cercano (`~1/(k r)`), que a
   `k r = 7.5` vale `~1–2 %`, es decir **~10x** el error de Yee que Q4 quiere acotar. Un criterio
   así no podría distinguir la dispersión numérica de la vecindad. Por eso:
   - **C1a · cuasi-1D**: `ny = 1`, esponja **sólo en x**, campo uniforme en `y` ⇒ onda plana.
     Aquí **no hay `H_0`, no hay campo cercano, no hay reparto**: la fase medida es
     directamente comparable con la fórmula de Yee de §4, y la amplitud no debe depender de `L`
     salvo por la reflexión de la esponja. Es el piso **más afilado**.
   - **C1b · 2D real**: `ny = nx`, esponja en los cuatro bordes, fuente de columna finita.
     Aquí viven el reparto cilíndrico (Q3) y la dependencia espuria con `L` en la geometría
     que el objeto va a usar.

4. **Convención de rama declarada** para `k_medido = (fase_desenvuelta + 2 pi m) / delta_d`: `m`
   se elige **una sola vez**, en el primer punto de la ventana, para que `k` quede junto a
   `2 pi nu / c` (vacío). Es **convención declarada**, no un ajuste contra la respuesta.

**§6.b · Q1..Q4 con esta geometría** (mismos umbrales; cero figuras vistas):

- **Q1** rizado dentro de una `L` `<= 0.05` en `C1a` y en `C1b`.
- **Q2** `delta = max_nu |ln A(36a) - ln A(24a)| <= 0.05` en `C1a` y en `C1b`.
- **Q3** `A(16a)/A(8a) = 0.7071 +- 0.05` en `C1b` (en `C1a` el valor teórico es `1.000`).
- **Q4** `max |error de fase medido| <= 3 · max |error de Yee|` en la ventana, medido en `C1a`.

---

## §2.c · ENMIENDA POST-CORRIDA (T_UTC 2026-10-12T06:37Z · **seis crudos medidos** · umbral `0.05`, ventana `nu` de §2.b.1, sondas de §2.b.2, Q1, Q2 y Q3 **NO se tocaron**)

UNIVERSO DECLARADO: este apéndice se mide sobre los **seis** crudos de `crudos/` (`.json` + `.txt`), el
árbol git de esta casa (**5** commits) y este mismo archivo; **unidad de conteo = celda**; hora de cada
cifra = la de su crudo. El acto 1 corrió **antes** de este apéndice: aquí no se predice, se declara lo
medido. §2.b sigue siendo verdadero **en su propia hora** («cero figuras medidas todavía») y **no se edita**.

### §2.c.1 · Veredicto del acto 1 (gen 6, medido)

`Q1 riz <= 0.05 en C1a y C1b: VERDE [C1a 0.005394 · C1b 0.005523]` ·
`Q2 delta = max_nu |ln A(36a) - ln A(24a)| <= 0.05 en C1a y C1b: VERDE [0.004548 · 0.006001]` ·
`Q3 reparto (2D 0.7071 · 1D 1.0000) +- 0.05: VERDE [C1a/24a 0.999907 · C1a/36a 1.000055 ·
C1b/24a 0.707449 · C1b/36a 0.707223]` ·
`Q4 fase medida 1.405e-03 <= 3 x Yee 1.463133e-03 = 4.389e-03: VERDE` (rama `m = 1`).
Acompañante analítico de §4: **1.463133e-03**. La condición de la l.160 de este archivo
(«**Condición de entrada al acto 2:** Q1..Q4 en verde») queda **cumplida**: el acto 2 queda **habilitado**,
no corrido. Cola de sonda del acto: `2.784e-03 · 2.636e-03 · 2.283e-02 · 2.283e-02`.

### §2.c.2 · Los SIETE defectos diagnosticados (mismo número que §2.b, otra causa)

Cada uno con su **testigo medido**; ninguno se repara tocando un umbral.

- **(a) Fila `j = 0` nunca escrita.** Con `ny = 1` el stencil de C1a escribía sólo `ez[1:, 1:]`. Testigo:
  gen 1, `cola de sonda (t > 0.75 T) = nan` en **las dos** corridas de C1a.
- **(b) Esponja multiplicativa de 1.25 a** = 0.19–0.69 longitudes de onda en la ventana: reflexión de orden
  `0.1–1`, es decir del tamaño del objeto. Reparado con **PML de campo partido (Berenger)**, 20 celdas por
  borde, orden 3, reflexión nominal `1.0e-06` ⇒ el absorbente pasa de objeto a frontera declarada.
- **(c) `FRACCION_COLA = 0.75`** era ciego al rebote que entra alrededor de `t = 28 a/c` ⇒ `T_COLA_MIN = 24.0 a/c`
  (1.5 x el cruce completo de la sonda lejana a 16a).
- **(d) Q1 tal como se selló era insatisfacible en C1b** (`A(8a)/A(4a)` con decaimiento cilíndrico `(k r)^(-1/2)`).
  Q1 se redefinió sobre `reparto`, **antes de correr** (§6.b).
- **(e) La rama de fase enmascaraba un cambio de signo.** La señal es **real** ⇒ la fase va como `-k d`; el
  código tomaba la diferencia y forzaba después `m >= 0`, con lo que `k` salía `+0.61` donde vale `+0.94`.
  Testigo: gen 2, `Q4 fase medida 7.495e+03`. Reparado con `fft[0] cdot conj(fft[1])` y la convención de
  §2.b.4 declarada. **Control**: gen 2 → gen 3, las **cuatro** colas de sonda quedaron **idénticas**
  (`2.784e-03 · 2.636e-03 · 2.283e-02 · 2.283e-02`) ⇒ la reparación tocó la fase, no el campo.
- **(f) La «cola de dominio» medía el PULSO PRIMARIO EN VUELO.** Con la ventana en `24 a/c` el máximo del
  dominio caía en `n = 768` ⇒ `t = 24.00 a/c`, **el propio primer paso de la ventana**. Testigo: gen 4,
  celdas `(375,0) n=768` · `(392,0) n=769` · `(352,18) n=768` · `(334,84) n=768`. Reparado con
  `T_COLA_DOM = 42.0 a/c` (= 1.17 x el cruce completo de `36 a/c`), declarado en docstring, constante,
  UNIVERSO e impresión. **Control de que se movió la métrica y no el objeto**: los **cuatro máximos
  finales** de gen 4 y gen 5 son **idénticos** — `2.841e-04 (74,0)` · `2.877e-04 (74,0)` · `3.721e-03 (9,352)` ·
  `5.265e-03 (6,17)` — mientras el residuo de ventana cayó `3.671e-02 → 2.780e-03` · `9.979e-01 → 2.727e-03` ·
  `8.501e-01 → 8.289e-03` · `8.533e-01 → 3.905e-02` (factores 13.2 · 366 · 102.6 · 21.9).
- **(g) La cola de dominio incluía las celdas de la PML**, donde el máximo es la onda incidente **entrando**
  al absorbente, no un residuo (el continente, no el objeto). Reparado con ventana **solo interior**:
  `vista = (slice(pml_x, nx - pml_x), slice(pml_y, ny - pml_y))`. Medido (gen 5 → gen 6):
  `C1b/24a 8.289e-03 → 7.728e-03` · `C1b/36a 3.905e-02 → 3.292e-02` · máximos finales
  `3.721e-03 → 2.522e-03` y `5.265e-03 → 2.899e-03` · **C1a sin cambio** (`2.780e-03` / `2.727e-03`;
  `2.841e-04` / `2.877e-04`): la corrección quitó el continente, no el objeto.

### §2.c.3 · Las seis generaciones de crudo (identidad medida, no atribuida)

El instrumento declarado dentro de cada crudo `sha16` **coincide** con `shasum -a 256` del archivo en su
commit (incluido HEAD en disco) ⇒ la identidad de cada corrida está **medida**, no inferida.

| gen | T_UTC | instrumento `sha16` | commit | ventana de dominio | veredicto |
|---|---|---|---|---|---|
| 1 | `055738Z` | `4251c1e0119afc8b` | `92834c4` | `t > 0.75 T` | **ROJO x4** |
| 2 | `061015Z` | `b19ea19ff70d080e` | (sin commit) | `24 a/c` | Q1–Q3 VERDE · **Q4 ROJO `7.495e+03`** |
| 3 | `061218Z` | `2b2056c6b83b1172` | `4e141af` | `24 a/c` | **VERDE x4** |
| 4 | `061815Z` | `c49c6a4de7224d4a` | `d64dc6d` | `24 a/c` (= `T_COLA_MIN`, sin constante propia) | VERDE x4 |
| 5 | `062245Z` | `406a51cfb85ccc10` | `d0bdc52` | `42 a/c` (`T_COLA_DOM`) | VERDE x4 |
| 6 | `062648Z` | `ef29596d00b873b5` | `0e93ff9` | `42 a/c`, **solo interior** | **VERDE x4** |

El gen 1 reporta `cola de dominio = nan` en C1a (el defecto (a)) y `1.097e-02` · `4.542e-02` en C1b; el gen 2
es la primera corrida de la v2 (reparación de fase aún sin aplicar) y no produjo commit propio.

### §2.c.4 · El control sin escritura en la casa (el residuo decae; no hay piso de PML)

Control en `/tmp/control_cola_decay.py`, **cero escrituras en esta casa** y **no está en custodia** (se cita
por su comando: `cd /tmp && python3 -B control_cola_decay.py`). Su cabecera declara
`UNIVERSO DECLARADO: 2 modos x L=24a x T_total en {72, 120, 180} a/c = 6 corridas`, y **no toca ningún
parámetro sellado**. Máximo final del **interior** con `T_total` = 72 / 120 / 180:

- `C1a: 2.841471e-04 → 2.259444e-04 → 1.874690e-04` (razones 0.7952 · 0.8297)
- `C1b: 2.521848e-03 → 1.445774e-03 → 1.102982e-03` (razones 0.5733 · 0.7629)

Decae **monótono** a lo largo de 108 `a/c` **sin piso**, con residuo de ventana y cola de sonda idénticos en
las tres `T` (`tau` de la razón medida = `48 a/c / |ln razon|` ⇒ `86–257 a/c`) ⇒ **no es un piso de reflexión
de la PML**. El primer punto del control (C1a `2.841471e-04` · C1b `2.521848e-03`) **reproduce** el crudo del
gen 6 ⇒ el control midió el **mismo objeto** que el acto.

### §2.c.5 · Lo que NO se tocó, y una hipótesis REFUTADA

- **No tocado:** el umbral `0.05` · la ventana `nu` de §2.b.1 · las sondas absolutas `8a` y `16a` de §2.b.2 ·
  las dos sub-corridas de §2.b.3 · Q1, Q2 y Q3. Las tres reparaciones ((e) fase, (f) ventana, (g) interior)
  **no cambiaron un solo veredicto**: la re-verificación en gen 5 y gen 6 volvió a dar Q1..Q4 **VERDE**.
- **Hipótesis refutada (por medición, no por argumento):** el componente de continua de la corriente blanda
  (`integral del pulso = sigma sqrt(2 pi) = 1.2533`; medido `≈ 1.25`) **no** explica el residuo. Media
  espacial de `Ez` en el interior / `pico2` (gen 6): `C1a -3.190e-07` · `+5.379e-07` · `C1b -8.242e-04` ·
  `-9.180e-04`, es decir `~1e-4` del residuo en C1a, y un campo estático **no** decaería como decae.
- **No adjudicado:** la **identidad** del residuo (contenido débil de `nu` baja absorbido despacio, esquina
  de la PML, o ambos). Lo medido es que **decae, no tiene piso y no crece con `L`**; adjudicarlo exigiría un
  acto propio, y no se hace aquí.
- **Custodia (medida, 06:36Z):** la casa de custodia `/Users/marlene/platino-labs` **no tiene árbol git**
  (gate A2.1 ⇒ RC 14, fail-closed por diseño, ya declarado en §1) y su carpeta
  `preregistro-fotonica-monotile-20261011/` contiene `PREREGISTRO_FOTONICA_MONOTILE_20261011.md`
  (`sha16 47c7c2ef5028ee40`, 17362 B) · 3 sondas · `crudos/`: **no contiene copia de este archivo** ⇒ no hay
  copia que parchear ni obsolescencia que declarar. El puntero a este apéndice desde la casa de custodia
  queda **PENDIENTE**: el gate A2.1 no puede crearlo (objetivo inexistente ⇒ RC 10; casa sin git ⇒ RC 14) y
  **no se escribe allí sin orden del dueño**.
