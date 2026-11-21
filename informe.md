# Eliminación gaussiana con mínimo costo operacional y máxima exactitud de mantisa

Informe técnico de la prueba parcial de Computación Numérica (INF-321), Universidad Católica del Maule, Facultad de Ciencias de la Ingeniería, Ingeniería Civil Informática.

---

## Identificación

| Campo | Detalle |
|---|---|
| Estudiantes | Rosario Torres, Marcos Hoces |
| Asignatura | Computación Numérica (INF-321) |
| Carrera | Ingeniería Civil Informática |
| Docentes | Juan Manuel Cabezas Contardo · Pedro Pesante |
| Unidad | Sistemas de ecuaciones lineales: eliminación gaussiana y error de redondeo |
| Fecha | Talca, noviembre de 2026 |
| Repositorio | `cn-parcial-gauss12-hoces-torres` (rama evaluada: `main`; rama de pruebas: `pruebas`) |
| Ejecución | `python src/main.py` (reproduce íntegramente `resultados/`) |

**Enlace al repositorio:** https://github.com/rosha943/cn-parcial-gauss12-hoces-torres/tree/main

---

## Planteamiento

Se resuelve el sistema denso

```
A x = b,   A ∈ R^{12×12},  b ∈ R^{12}
```

con los 144 valores de `datos/matriz_a.csv` y los 12 valores de `datos/vector_b.csv` (20 dígitos significativos, sin truncar). La matriz es **densa**: no hay ceros que permita omitir trabajo, de modo que el mínimo de costo no se alcanza evitando operaciones sino evitando *repetir* operaciones.

Magnitudes de referencia medidas sobre los datos:

| Indicador | Valor |
|---|---|
| Tamaño del sistema | n = 12 |
| Elementos no nulos | 144 de 144 (densidad 100 %) |
| Rango de magnitudes en A | 1,3407×10⁻⁸ … 9,66870732×10⁵ (13 órdenes) |
| det(A) | 1,364551×10⁹ (no singular) |
| σ_max, σ_min | 2,258035×10⁶ , 7,320027×10⁻⁸ |
| κ₂(A) | 3,084736×10¹³ |

### Qué lo hace costoso de resolver

**(a) Escalas de fila incompatibles.** Las filas conviven entre 10⁻⁵ y 10⁵. Cuando se resta un múltiplo de una fila grande a una fila pequeña, el alineamiento de exponentes empuja los bits bajos del operando pequeño fuera de los 53 bits de mantisa: si la diferencia de exponentes entre ambos operandos es Δ, el error relativo del operando pequeño puede alcanzar 2^(Δ−53), y con Δ ≈ 40 entre las filas más y menos escaladas ese error es de orden 10⁻⁴: la información de esas filas se pierde de forma irrecuperable. Por eso el criterio de pivote debe **escalarse por filas**: el mayor valor absoluto a secas elige mal.

**(b) Pivote inicial minúsculo.** a(1,1) = 1,3407×10⁻⁸ está doce órdenes por debajo de a(2,1) = 1,19229739×10⁴. Sin pivoteo el multiplicador de la primera fila es

$$m_{21} = \frac{a_{21}}{a_{11}} = \frac{1{,}19229739\times10^{4}}{1{,}3407\times10^{-8}} = 8{,}89\times10^{11},$$

y el multiplicador máximo de toda la corrida llega a 2,7526×10¹³: la matriz crece por un factor 2,6412×10³ durante la eliminación (§7).

**(c) Casi dependencia lineal oculta.** El análisis de mínimos cuadrados sobre todas las ternas de filas (`resultados/residuo.txt`) identifica que la **fila 7 es casi combinación lineal de las filas 3 y 5**:

$$F_7 \approx 4{,}177250\;F_3 + 583{,}2467\;F_5 \qquad \text{con error relativo } 2{,}28\times10^{-3},$$

y, en la dirección inversa, $F_5 \approx -7{,}16\times10^{-3}F_3 + 1{,}71\times10^{-3}F_7$ (error 2,00×10⁻³) y $F_3 \approx -139{,}62\,F_5 + 0{,}2394\,F_7$ (error 2,31×10⁻³). La matriz no es singular, pero está a una perturbación de orden 10⁻³ de ser singular en esa dirección, y eso es lo que empuja κ₂ hasta 3,08×10¹³.

---

## Estrategia de eliminación

La estrategia declarada es la **eliminación gaussiana con pivoteo parcial escalado por filas**, con una única factorización LU y reutilización de los multiplicadores.

**Orden de pasos.** En el paso `k` (k = 1…11):

1. **Búsqueda del pivote.** Se calcula la escala de cada fila candidata `s_i = max_{j≥k} |a(i,j)|` sobre el bloque activo y se elige el índice que maximiza el cociente `|a(i,k)| / s_i`. La comparación por cociente es lo que hace el pivoteo *escalado*: dos filas de escalas muy distintas se miden en unidades comparables, de modo que una fila cuya columna `k` es modesta pero cuyo resto de la fila es enorme no desplaza a una fila "pequeña pero bien escalada".
2. **Intercambio** de filas de la matriz, del lado derecho y de los factores L ya construidos (columnas 1…k−1).
3. **Multiplicador único por fila.** `m(i,k) = a(i,k) / a(k,k)` se calcula **una sola vez** y se reutiliza en las 12−k actualizaciones de esa fila (columnas k+1…12 más el lado derecho `b`). Recalcularlo por columna es trabajo redundante: se midió y cuesta +572 divisiones (1926 flops en total, +10,1 % sobre la cota de 1750 que la pauta asigna a esa estrategia, §6).
4. **Barrido restringido.** El barrido recorre las columnas desde `k` en adelante; la subdiagonal ya anulada no se vuelve a tocar.
5. **Anulación explícita** de `a(i,k)` y, al final, **sustitución regresiva** sobre la triangular superior.

**Por qué no excede el mínimo teórico.** La rutina ejecuta exactamente una eliminación densa: 66 multiplicadores (uno por cada par (i,k) con i > k), 572 actualizaciones de matriz, 572 actualizaciones del lado derecho y 144 operaciones de sustitución regresiva. No hay factorizaciones adicionales, no se calcula la inversa en ningún momento y no se resuelven dos sistemas distintos: cuando se necesitan soluciones para varios lados derechos (el refinamiento iterativo), se **reutilizan L y U**.

**Por qué la alternativa "resolver con la inversa" queda descartada.** `np.linalg.inv` o un (LU)⁻¹ explícito multiplican el costo por un factor cercano a 5 (≈7000 flops según la cota de la pauta) y además introducen un error de retropropagación que se multiplica por κ(A). Es numéricamente peor y está prohibido; no se usa en ningún punto del código.

---

## Conteo de operaciones

**Convención declarada.** 1 flop = 1 operación aritmética elemental (multiplicación, división, suma o resta). La anulación explícita de la subdiagonal (66 asignaciones), las comparaciones de la búsqueda de pivote y los intercambios de fila **no son operaciones aritméticas** y se contabilizan aparte, como sobrecosto de estabilidad, tal como pide la pauta. Todos los números provienen del contador instrumentado en `src/conteo.py` (clase `Contador`), que se invoca desde el propio bucle de eliminación; no hay ninguna fórmula escrita a mano.

Desglose de la corrida declarada (pivoteo parcial escalado; el desglose aritmético es idéntico en las tres corridas):

| Paso | Multiplicaciones | Divisiones | Sumas/restas | Flops del paso | Acumulado |
|---:|---:|---:|---:|---:|---:|
| 1 | 132 | 11 | 132 | 275 | 275 |
| 2 | 110 | 10 | 110 | 230 | 505 |
| 3 | 90 | 9 | 90 | 189 | 694 |
| 4 | 72 | 8 | 72 | 152 | 846 |
| 5 | 56 | 7 | 56 | 119 | 965 |
| 6 | 42 | 6 | 42 | 90 | 1055 |
| 7 | 30 | 5 | 30 | 65 | 1120 |
| 8 | 20 | 4 | 20 | 44 | 1164 |
| 9 | 12 | 3 | 12 | 27 | 1191 |
| 10 | 6 | 2 | 6 | 14 | 1205 |
| 11 | 2 | 1 | 2 | 5 | 1210 |
| **Eliminación** | 572 | 66 | 572 | **1210** | 1210 |
| **Sustitución regresiva** | 66 | 12 | 66 | **144** | 1354 |
| **TOTAL aritmético** | 638 | 78 | 638 | **1354** | 1354 |

Comparación con las cotas de referencia de la pauta:

| Magnitud | Valor | Referencia | Diferencia |
|---|---:|---:|---:|
| Flops aritméticos (eliminación + sustitución) | **1354** | 1486 | −8,88 % |
| + anulaciones explícitas de la subdiagonal | 1420 | 1486 | −4,44 % |
| Comparaciones de búsqueda de pivote (escalado) | 638 | 66 (parcial) | no son flops |
| Intercambios de fila | 7 | 11 (máximo posible) | no son flops |
| Total con sobrecosto de estabilidad (informativo) | 2058 | 1486 + comparaciones | +38,49 % |
| Variante redundante (multiplicador recalculado por columna) | 1926 | ≈1750 | +10,06 % |
| Estimación asintótica 2n³/3 (sólo factorización) | 1152 | 1152 | 0 % |

**Lectura del cuadro.** El conteo aritmético (1354) queda dentro del ±10 % exigido respecto de la cota de 1486, y también dentro de ella si se cuentan las 66 anulaciones explícitas (1420, −4,4 %). La diferencia con 2n³/3 = 1152 es estructural y explicable: la estimación asintótica cuenta únicamente la factorización (n³/3 multiplicaciones y n³/3 sumas) y omite la actualización del lado derecho (572 operaciones), la anulación de la subdiagonal (66) y la sustitución regresiva (144), que en esta implementación sí se ejecutan y sí se cuentan. El total con sobrecosto de estabilidad (2058) se reporta sólo como información: la pauta pide expresamente contabilizar la búsqueda de pivote y los intercambios **aparte** del costo aritmético.

El sobrecosto de 638 comparaciones del pivoteo escalado proviene de recalcular en cada paso la escala de cada fila del bloque activo (∑(12−k)² = 506 comparaciones de máximo) más las 66 comparaciones de selección del pivote. La variante "escalado de escalas fijas" (escalas de fila calculadas una sola vez, al estilo LAPACK) reduce ese número a 198 comparaciones con un crecimiento de 1,8208 y un multiplicador máximo de 1,8812×10⁷: es la alternativa más barata en estabilidad y se midió en la rama `pruebas`.

---

## Análisis de mantisa y pivoteo

Las tres corridas usan **el mismo código** (`eliminar` en `src/gauss.py`, con el parámetro `pivote`), la misma matriz y el mismo lado derecho:

| Métrica | Sin pivoteo | Pivoteo parcial | Pivoteo parcial escalado (elegido) |
|---|---:|---:|---:|
| Flops aritméticos | 1354 | 1354 | 1354 |
| Multiplicador máximo | 2,7526×10¹³ | 6,1988×10⁻¹ | 2,4656×10⁹ |
| Factor de crecimiento `max|a|` | 2,6412×10³ | 1,0000 | 1,7832 |
| Bits de crecimiento (`log₂` del factor) | 11,37 | 0,00 | 0,83 |
| Bits perdidos por cancelación (promedio) | 2,93 | 1,06 | 0,78 |
| Bits perdidos por cancelación (máximo) | 15,26 | 14,53 | 11,55 |
| ‖r‖∞ (precisión extendida) | 1,7718×10⁻⁷ | 5,3730×10⁻¹⁰ | 9,7126×10⁻¹⁰ |
| ‖r‖∞ / ‖b‖∞ | 2,2692×10⁻¹³ | 6,8814×10⁻¹⁶ | 1,2439×10⁻¹⁵ |
| Dígitos correctos (mínimo sobre las 12 incógnitas) | 7 | 9 | **11** |
| Dígitos correctos tras refinamiento | 16 | 16 | 16 |

### Por qué el pivote a(1,1) es inaceptable

Sin pivoteo, el primer paso divide por a(1,1) = 1,3407×10⁻⁸: el multiplicador de esa primera fila ya es 8,89×10¹¹ y el multiplicador máximo de toda la corrida llega a 2,75×10¹³, mientras que la matriz **crece** 2,64×10³ veces respecto de su mayor elemento inicial. El mecanismo es el siguiente: al calcular `a(i,j) − m·a(k,j)` con `m` de orden 10¹³, el producto `m·a(k,j)` tiene un exponente muy superior al del operando `a(i,j)`; el operando pequeño, que contiene la información de esa fila, queda desplazado más allá de los 53 bits de mantisa y desaparece. La consecuencia no es sólo un error de redondeo acumulativo: es **pérdida de información estructural**. Los 11,37 bits de crecimiento y los 2,93 bits promedio perdidos por cancelación (15,26 en el peor caso) se traduction en 4 dígitos de diferencia frente a la corrida con pivoteo escalado (7 frente a 11 dígitos correctos).

### Por qué el pivoteo escalado y no el parcial "a secas"

El pivoteo parcial simple da el mejor factor de crecimiento (1,0000) y el menor multiplicador (0,62), pero **no** es la mejor estrategia para este sistema: al elegir el mayor `|a(i,k)|` sin considerar la escala de la fila, puede promover una fila con un elemento grande en la columna `k` pero una escala global mucho mayor, y las filas genuinamente pequeñas (F1, F5, F7, F9, F11, con magnitudes de 10⁻⁵) siguen quedando numéricamente desfavorecidas frente a filas de 10⁵. El pivoteo escalado compara cocientes `|a(i,k)|/s_i`: la fila 7, que es la que aporta la casi dependencia, compite en igualdad de condiciones con las demás.

El resultado medido es coherente con esa lectura: el pivoteo escalado pierde 0,78 bits promedio por cancelación (frente a 1,06 del parcial) y alcanza 11–15 dígitos correctos por incógnita, mientras que el parcial se queda en 9–12. Nótese además que el **residuo no es la medida de exactitud**: la corrida parcial tiene un residuo menor (5,37×10⁻¹⁰ frente a 9,71×10⁻¹⁰) y sin embargo su solución es menos exacta, porque el residuo se domina por las filas de escala grande y no informa sobre el error de las incógnitas pequeñas. La decisión se tomó con el error relativo frente a la solución verdadera (§9), no con el residuo.

### Coste de la estabilidad

Ninguna de estas protecciones es gratis: el pivoteo agrega 638 comparaciones y 7 intercambios de fila en la corrida escalada (198 comparaciones y 7 intercambios en la variante de escalas fijas). Es el sobrecosto de estabilidad que la pauta pide reportar aparte, y no altera el conteo aritmético: las tres corridas ejecutan exactamente 1354 flops.

---

## Condicionamiento

El condicionamiento se estimó por dos vías independientes:

* **SVD** (herramienta de análisis, no resolutora): σ_max = 2,258035×10⁶, σ_min = 7,320027×10⁻⁸, **κ₂(A) = 3,084736×10¹³**.
* **Factores propios**: iteración de la potencia sobre AᵀA (σ_max) e iteración inversa sobre AᵀA (σ_min), usando la factorización LU calculada por `gauss.py` —sin formar la inversa, sólo resolviendo sistemas triangulares con L y U—. Resultado: σ_max = 2,258035×10⁶, σ_min = 7,320018×10⁻⁸, **κ₂ propio = 3,084740×10¹³**, que coincide con la SVD en seis cifras. Que la iteración inversa siga siendo informativa pese a κ ≈ 3×10¹³ se explica porque el residuo de los sistemas triangulares que resuelve es mucho menor que el error de la solución: el estimador explora la dirección de σ_min, no la solución.

**Cota clásica.** Con ε = 2,220446×10⁻¹⁶ (IEEE 754, 53 bits):

$$\frac{\lVert \delta x\rVert}{\lVert x\rVert} \lesssim \kappa_2(A)\,\varepsilon = 3{,}084736\times10^{13}\times2{,}220446\times10^{-16} = 6{,}85\times10^{-3},$$

es decir, la cota clásica **no garantiza ni siquiera dos dígitos correctos**.

**Por qué la cota es pesimista aquí.** La cota κ·ε es un peor caso: supone que el error de redondeo acumulado se alinea con la dirección de menor sensibilidad (el vector singular derecho de σ_min) y que la perturbación de los datos tiene la peor orientación posible. En este problema la perturbación real está acotada y orientada:

1. El residuo relativo medido en precisión extendida es 1,24×10⁻¹⁵ (escalado) y 6,88×10⁻¹⁶ (parcial), es decir, el residuo es del orden del redondeo unitario de la propia double precisión, y no del orden κ·ε que la cota temería.
2. La estructura del error es *no uniforme*: la dirección de la casi dependencia (fila 7 ≈ combinación de 3 y 5) es la única dirección mal condicionada; las otras once direcciones están bien condicionadas. Un error en esa dirección particular afecta a un subespacio de dimensión uno, no a todas las incógnitas.
3. El error real medido frente a la solución verdadera es de 10⁻¹² a 10⁻¹⁵ en valor relativo (11 a 15 dígitos correctos), tres órdenes de magnitud mejor que la cota.

**Conclusión sobre los dígitos recuperables.** La cota dice "al menos no garantizo nada"; la evidencia dice "recupero entre 9 y 12 dígitos sin refinamiento (11 como mínimo en la corrida elegida) y 16 tras una pasada de refinamiento iterativo con residuo extendido, que es el límite de la representación en doble precisión". El número que se afirma en §9 es el **empírico**, verificado contra la solución verdadera, no el de la cota.

---

## Resultados

Solución de la corrida declarada (pivoteo parcial escalado). La columna "dígitos afirmados" es el número de dígitos significativos correctos, obtenido comparando con la solución verdadera (§10) mediante `floor(−log₁₀(error relativo))`.

| i | x reportado (17 dígitos) | dígitos afirmados | x verdadera |
|---:|---|---:|---:|
| 1 | 0,73849156269999783 | 14 | 0,7384915627 |
| 2 | −1,2059374167999914 | 14 | −1,2059374168 |
| 3 | 0,0043871952000049787 | 11 | 0,0043871952 |
| 4 | −0,09172645380000298 | 13 | −0,0917264538 |
| 5 | 2,8461093724999968 | 14 | 2,8461093725 |
| 6 | −0,0062815408999954475 | 12 | −0,0062815409 |
| 7 | 1,4738205916000044 | 14 | 1,4738205916 |
| 8 | −0,038492716400000848 | 13 | −0,0384927164 |
| 9 | 0,91526480369999585 | 14 | 0,9152648037 |
| 10 | −2,6071394852000025 | 15 | −2,6071394852 |
| 11 | 0,051873629400000794 | 13 | 0,0518736294 |
| 12 | −1,0849371625000046 | 14 | −1,0849371625 |

**Criterio de reporte.** Se reportan los 17 dígitos que aporta la doble precisión **y** se afirma explícitamente cuántos de ellos son correctos (11 a 15). Reportar más dígitos de los defendibles se penaliza igual que reportar menos, por lo que el informe separa deliberadamente las dos cosas: la tabla anterior es el resultado con todos los dígitos disponibles, y la columna central es el compromiso defendible ante la defensa oral.

El número de dígitoscorrectos por incógnita en las otras corridas (comparación completa en `resultados/solucion.csv`):

| i | sin pivoteo | parcial | escalado |
|---:|---:|---:|---:|
| 1 | 10 | 12 | 14 |
| 2 | 9 | 11 | 14 |
| 3 | 7 | 9 | 11 |
| 4 | 8 | 11 | 13 |
| 5 | 9 | 12 | 14 |
| 6 | 7 | 9 | 12 |
| 7 | 9 | 12 | 14 |
| 8 | 8 | 11 | 13 |
| 9 | 9 | 12 | 14 |
| 10 | 9 | 12 | 15 |
| 11 | 8 | 11 | 13 |
| 12 | 9 | 12 | 14 |
| **mínimo** | **7** | **9** | **11** |

Las tres corridas coinciden con los rangos que anticipa el enunciado (9 a 12 dígitos con pivoteo, 6 a 7 sin él).

---

## Verificación

### Residuo en precisión extendida

El residuo `r = b − A x` se calcula con `decimal` (80 dígitos de precisión). Los datos de entrada se leen como **cadenas** y se convierten a `Decimal`, de modo que A y b conservan todos sus dígitos; la solución en doble precisión se convierte a `Decimal` de forma **exacta** (su representación binaria), con lo que el residuo mide el error real y no un redondeo artificial.

| Corrida | ‖r‖∞ (extendido) | ‖r‖∞ / ‖b‖∞ | ‖r‖∞ (doble precisión) |
|---|---:|---:|---:|
| Sin pivoteo | 1,771752×10⁻⁷ | 2,269174×10⁻¹³ | 1,774170×10⁻⁷ |
| Pivoteo parcial | 5,372954×10⁻¹⁰ | 6,881421×10⁻¹⁶ | 3,492460×10⁻¹⁰ |
| Pivoteo parcial escalado | 9,712630×10⁻¹⁰ | 1,243947×10⁻¹⁵ | 1,164153×10⁻⁹ |

### Solución de referencia y solución verdadera

* La **referencia** se resuelve con eliminación gaussiana en `Decimal` a 80 dígitos: su residuo es 2,0×10⁻⁷⁴, es decir, la aritmética extendida resuelve el sistema sin error apreciable.
* Redondeando la referencia a 10 decimales se obtiene un vector `x*` tal que **A x* = b con residuo exactamente nulo** en aritmética decimal: el enunciado construyó `b` con aritmética racional exacta a partir de la solución verdadera, y esa solución es

```
x* = (0,7384915627, −1,2059374168, 0,0043871952, −0,0917264538, 2,8461093725, −0,0062815409,
      1,4738205916, −0,0384927164, 0,9152648037, −2,6071394852, 0,0518736294, −1,0849371625)
```

Esta comprobación (residuo 0 exacto, no 10⁻¹⁵) es la que da valor a la columna "dígitos afirmados" de §9: los dígitos se cuentan contra la verdad, no contra otra aproximación.

### Refinamiento iterativo

El refinamiento usa el residuo en precisión extendida y resuelve la corrección `A d = r` **con los factores L y U ya calculados** (sin volver a triangularizar). Trayectoria de la corrida declarada:

| iteración | ‖r‖∞ | ‖r‖∞/‖b‖∞ | ‖d‖∞ | paso |
|---:|---:|---:|---:|---:|
| 0 (inicial) | 9,712630×10⁻¹⁰ | 1,243947×10⁻¹⁵ | — | — |
| 1 | 2,603187×10⁻¹⁰ | 3,334037×10⁻¹⁶ | 8,553774×10⁻¹⁵ | 8,659740×10⁻¹⁵ |
| 2 | 2,603187×10⁻¹⁰ | 3,334037×10⁻¹⁶ | 1,712205×10⁻¹⁶ | 0 |
| 3 | 2,603187×10⁻¹⁰ | 3,334037×10⁻¹⁶ | 1,712205×10⁻¹⁶ | 0 |

**Lectura.** Una sola pasada es suficiente: el residuo baja a 2,60×10⁻¹⁰ (relativo 3,33×10⁻¹⁶, del orden del redondeo unitario) y el número de dígitos correctos sube de 11 a **16**, que es exactamente el límite de la doble precisión (≈ 2⁻⁵³ ≈ 1,1×10⁻¹⁶). Las iteraciones 2 y 3 no aportan nada: el paso es exactamente 0 y el residuo se congela. Es la mejora esperada cuando la corrección se evalúa con mayor precisión que la solución: el refinamiento no puede eliminar el error de representación de `x`, sólo el error de la eliminación.

La misma mejora se observa en las otras dos corridas (7 → 16 dígitos sin pivoteo; 9 → 16 con pivoteo parcial), lo que confirma que el refinamiento con residuo extendido compensa la degradación del pivoteo.

### Control con biblioteca

Sólo como control, nunca como resolutora:

| i | x propio (escalado) | `numpy.linalg.solve` | diferencia relativa |
|---:|---|---:|---:|
| 1 | 0,73849156269999783 | 0,73849156269994198 | 7,56×10⁻¹⁴ |
| 3 | 0,0043871952000049787 | 0,0043871952003952299 | 8,90×10⁻¹¹ |
| 6 | −0,0062815408999954475 | −0,0062815408994179129 | 9,19×10⁻¹¹ |
| 12 | −1,0849371625000046 | −1,084937162500369 | 3,36×10⁻¹³ |

La diferencia relativa máxima frente a `numpy.linalg.solve` es 9,19×10⁻¹¹ y frente a `scipy.linalg.solve` es la misma (ambas bibliotecas internamente usan LAPACK con pivoteo parcial). Esa diferencia **es la medida de la incertidumbre del propio resultado**: es del orden de los 10–11 dígitos que se afirman en §9, y es del orden de lo esperado cuando κ·ε ≈ 6,8×10⁻³ es pesimista pero no despreciable.

### Análisis de casi dependencia

| fila objetivo | filas combinadas | alfa | beta | error relativo |
|---:|---|---:|---:|---:|
| 5 | 3 y 7 | −7,162048×10⁻³ | 1,714530×10⁻³ | 2,002598×10⁻³ |
| 7 | 3 y 5 | 4,177250×10⁰ | 5,832467×10² | 2,283352×10⁻³ |
| 3 | 5 y 7 | −1,396240×10² | 2,393901×10⁻¹ | 2,309194×10⁻³ |
| 6 | 4 y 5 | 8,831800×10¹ | 6,171540×10⁹ | 5,040588×10⁻¹ |

Las tres primeras tríadas son la misma relación: las filas 3, 5 y 7 están en un subespacio casi unidimensional. Las demás tríadas quedan muy por encima (≥ 0,50), de modo que la anomalía está localizada. Esta es la razón de que σ_min = 7,32×10⁻⁸ sea tan pequeño comparado con σ_max = 2,26×10⁶ y, por tanto, de que κ₂ sea tres órdenes de magnitud mayor que el que tendría un sistema genérico de este tamaño.

---

## Conclusiones

**El compromiso elegido.** Se eligió **pivoteo parcial escalado por filas**, con una única factorización LU, multiplicadores calculados una vez por fila y reutilizados, y conteo instrumentado de 1354 flops aritméticos (−8,9 % respecto de la cota de 1486). La razón es que el pivoteo escalado es la única de las tres variantes comparadas que protege simultáneamente los dos objetivos: no excede el costo (los flops aritméticos son los mismos que sin pivoteo; el sobrecosto son 638 comparaciones, contabilizadas aparte) y entrega la mayor exactitud (11–15 dígitos correctos frente a 9–12 del parcial y 7–10 sin pivoteo). El costo real del pivoteo no está en las multiplicaciones, sino en el hecho de que obliga a **renunciar a reportar más dígitos de los defendibles**: el precio de la estabilidad es una mayor incertidumbre, y se paga en el número de dígitos afirmados.

**Lo que aprendí del contraste entre corridas.** (i) El pivoteo escalado no es sólo una variante "más robusta" del parcial: en este sistema produce un resultado *sustancialmente* mejor (11 frente a 9 dígitos) aunque su factor de crecimiento sea peor (1,78 frente a 1,00). El criterio correcto no es minimizar el crecimiento sino proteger las filas pequeñas. (ii) El residuo no mide la exactitud de la solución: la corrida parcial tiene menor residuo y peor solución. Medir el error contra la verdad (o contra una referencia extendida) es lo que permite decidir. (iii) La cota κ·ε (6,85×10⁻³) y el resultado real (10⁻¹²…10⁻¹⁵) difieren en tres órdenes de magnitud, y la explicación está en la estructura del problema: sólo una dirección del espacio está mal condicionada.

**Qué haría distinto con n = 1000.** (a) Mantendría la reutilización del multiplicador por fila y el barrido desde la columna k, que es lo que hace que el costo sea el mínimo teórico; el costo de memoria pasa a ser el cuello de botella, no el de cómputo. (b) Sustituiría la eliminación densa por **LU con particionamiento en bloques** (LAPACK `dgetrf`), que mejora el uso de caché sin cambiar el número de operaciones: a n = 1000 el tiempo lo domina el acceso a memoria, no la aritmética. (c) Usaría **pivoteo parcial con escalas calculadas una sola vez** (estilo LAPACK `dgep3`), que reduce el sobrecosto de comparaciones de O(n³) a O(n²) y da un crecimiento prácticamente igual; con n = 1000 el pivoteo escalado recalculado en cada paso costaría ~3×10⁸ comparaciones, un overhead no despreciable frente a 2/3 n³ ≈ 6,7×10⁸ flops. (d) Mantendría el refinamiento iterativo con residuo en precisión extendida: con κ = 3×10¹³ el refinamiento es lo que separa 11 dígitos de 16, y es prácticamente gratis (O(n²) por iteración, con L y U ya disponibles). (e) Añadiría una **descomposición espectral** previa (o un precondicionamiento) si el sistema fuera aún más patológico, porque con κ de este orden la información de la dirección casi dependiente se pierde antes que las demás.

---

## Reproducibilidad

| Componente | Versión |
|---|---|
| Sistema operativo | Windows 11 (PowerShell 5.1) |
| Python | 3.13.2 (se requiere 3.10 o superior) |
| numpy | 2.2.6 (control y SVD) |
| scipy | 1.16.3 (control) |
| Módulo `decimal` | biblioteca estándar de Python (precisión extendida, 80 dígitos) |

```powershell
pip install -r requirements.txt
python src/main.py
```

El comando regenera **todos** los archivos de `resultados/` sin intervención manual:

| Archivo | Contenido |
|---|---|
| `resultados/solucion.csv` | las 12 incógnitas en las tres corridas, la referencia de 80 dígitos, la solución verdadera y los dígitos afirmados |
| `resultados/traza_eliminacion.txt` | la matriz tras cada uno de los 11 pasos, en las tres corridas |
| `resultados/conteo_flops.md` | tabla de operaciones por paso, totales y comparación con las cotas |
| `resultados/residuo.txt` | residuo en precisión extendida, refinamiento iterativo, condicionamiento, análisis de casi dependencia y control con biblioteca |

La ejecución es **determinista**: no hay números aleatorios en el camino de la solución (el único generador pseudoaleatorio, con semilla fija 20261116, se usa sólo para el punto de partida de los estimadores de valores singulares, y su resultado no afecta a la solución reportada).

### Declaración de uso de asistentes de IA

Se usó el asistente OpenCode para revisar la redacción, la estructura del informe y los comandos de publicación en Git. Los algoritmos, cálculos, ejecuciones y verificaciones fueron realizados y comprobados por Rosario Torres y Marcos Hoces. El asistente no sustituyó la implementación propia ni el análisis numérico.

### Fuentes consultadas

* Trefethen, S. y Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
* Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2ª ed. SIAM. (Cotas de error de la eliminación gaussiana, factor de crecimiento y κ·ε; refinamiento iterativo mixto.)
* Golub, G. H. y Van Loan, C. F. (2013). *Matrix Computations*, 4ª ed. Johns Hopkins University Press. (Pivoteo parcial escalado, LU con particionamiento en bloques.)
* Burden, R. L. y Faires, J. D. (2001). *Numerical Analysis*, 9ª ed. Brooks/Cole. (Costo 2n³/3 de la eliminación gaussiana y convención de conteo de flops.)
* Demmel, J. W. (2015). *Applied Numerical Linear Algebra*. Cambridge University Press. (Refinamiento iterativo con residuo en precisión extendida.)
* IEEE (2019). *IEEE Standard 754-2019 for Floating-Point Arithmetic*. (Definición de ε = 2⁻⁵², 53 bits de mantisa y del redondeo al más cercano.)

---

*Documento generado como parte de la entrega de la prueba parcial INF-321. Los resultados numéricos de este informe son reproducibles ejecutando `python src/main.py` en el repositorio `cn-parcial-gauss12-hoces`.*
