# cn-parcial-gauss12-hoces-torres

Resuelve el sistema lineal denso de 12 incógnitas por eliminación gaussiana con pivoteo parcial escalado, minimizando el costo operacional y acotando el error de redondeo de la mantisa.

## Autoría y asignatura

**Rosario Torres, Marcos Hoces** · Computación Numérica (INF-321) · Ingeniería Civil Informática · Universidad Católica del Maule · noviembre de 2026.

Docentes: Juan Manuel Cabezas Contardo · Pedro Pesante.

## Requisitos

* Python 3.10 o superior (probado en 3.13.2).
* Dependencias: `numpy==2.2.6` y `scipy==1.16.3`. La precisión extendida usa `decimal` de la biblioteca estándar, sin dependencias adicionales.

```powershell
pip install -r requirements.txt
```

## Ejecución

```powershell
python src/main.py
```

El comando reproduce **todo** el contenido de `resultados/`, sin intervención manual:

| Archivo generado | Contenido |
|---|---|
| `resultados/solucion.csv` | las 12 incógnitas en las tres corridas, la referencia de 80 dígitos, la solución verdadera y los dígitos afirmados como correctos |
| `resultados/traza_eliminacion.txt` | la matriz tras cada uno de los 11 pasos de eliminación, en las tres corridas |
| `resultados/conteo_flops.md` | tabla de operaciones paso a paso, totales y comparación con las cotas de referencia |
| `resultados/residuo.txt` | residuo en precisión extendida, refinamiento iterativo, condicionamiento, análisis de casi dependencia y control con biblioteca |

## Estructura de carpetas

```
cn-parcial-gauss12-hoces-torres/
├─ readme.md            # esta tarjeta: qué es y cómo se ejecuta
├─ informe.md           # informe técnico con el análisis de error (10 secciones)
├─ requirements.txt     # dependencias con versión fija
├─ datos/               # A densa 12x12 (144 valores sin truncar) y b con 20 dígitos
├─ src/                 # gauss.py (eliminación), conteo.py (flops), verificacion.py, main.py
├─ resultados/          # salida generada por main.py
└─ .gitignore
```

* `datos/`: los datos de entrada en CSV (una fila por línea; en `matriz_a.csv`, una fila por fila de la matriz).
* `src/`: la implementación propia. La eliminación gaussiana y el contador de flops están escritos desde cero; numpy y scipy se usan sólo como control y para la SVD, nunca para resolver el sistema.
* `resultados/`: los archivos que produce la ejecución.

## Resultado resumido

Estrategia declarada: pivoteo parcial escalado, una sola factorización LU.

| i | x (10 dígitos) |
|---:|---:|
| 1 | 0,7384915627 |
| 2 | −1,205937417 |
| 3 | 0,0043871952 |
| 4 | −0,0917264538 |
| 5 | 2,846109372 |
| 6 | −0,0062815409 |
| 7 | 1,473820592 |
| 8 | −0,0384927164 |
| 9 | 0,9152648037 |
| 10 | −2,607139485 |
| 11 | 0,0518736294 |
| 12 | −1,084937163 |

| Magnitud | Valor |
|---|---:|
| Flops aritméticos (eliminación + sustitución regresiva) | 1354 |
| Comparaciones de búsqueda de pivote (escalado) | 638 |
| ‖r‖∞ del residuo en precisión extendida | 9,71×10⁻¹⁰ |
| ‖r‖∞ tras una pasada de refinamiento iterativo | 2,60×10⁻¹⁰ |
| κ₂(A) | 3,08×10¹³ |

## Informe

El análisis completo (planteamiento, estrategia, conteo de operaciones, mantisa y pivoteo, condicionamiento, resultados, verificación, conclusiones y reproducibilidad) está en **[informe.md](informe.md)**.

## Licencia y nota de uso académico

Trabajo académico de la prueba parcial de Computación Numérica (INF-321) de la Universidad Católica del Maule. Uso exclusivamente educativo y de evaluación.

Responsable de la documentación y la publicación: Rosario Torres.
