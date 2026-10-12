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
