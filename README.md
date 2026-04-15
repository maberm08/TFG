# Valoración numérica de opciones bajo Black–Scholes mediante diferencias finitas

Este repositorio contiene una implementación en MATLAB de esquemas de diferencias finitas para la valoración de opciones bajo el modelo de Black–Scholes.

Se incluyen los tres esquemas clásicos estudiados en el trabajo:

- esquema explícito,
- esquema implícito,
- esquema de Crank–Nicolson.

El código está planteado de forma modular para poder reutilizar el mismo núcleo numérico en distintos contratos, modificando únicamente:

- el payoff terminal,
- las condiciones de contorno,
- el intervalo espacial.

## Estructura del repositorio

El repositorio se organiza en torno a dos casos:

- `BS/cte`: programas para coeficientes constantes.
- `BS/dept`: programas para coeficientes dependientes del tiempo.

### Carpeta `BS/cte`

Incluye las implementaciones y scripts de prueba para el caso en que \(r\), \(q\) y \(\sigma\) son constantes.

- `bs_explicito.m`: implementación del esquema explícito.
- `bs_implicito.m`: implementación del esquema implícito.
- `bs_cn.m`: implementación del esquema de Crank–Nicolson.
- `lu_tridiag.m`: factorización LU para matrices tridiagonales, reutilizada en los esquemas implícito y de Crank–Nicolson.
- `test_call_europea.m`: comparación de los tres esquemas para una call europea vanilla.
- `test_put_europea.m`: comparación de los tres esquemas para una put europea vanilla.

### Carpeta `BS/dept`

Incluye las implementaciones y scripts de prueba para el caso en que \(r\), \(q\) y \(\sigma\) dependen del tiempo.

- `bs_explicito_t.m`: implementación del esquema explícito con coeficientes dependientes de \(t\).
- `bs_implicito_t.m`: implementación del esquema implícito con coeficientes dependientes de \(t\).
- `bs_cn_t.m`: implementación del esquema de Crank–Nicolson con coeficientes dependientes de \(t\).
- `test_call_europea_t.m`: comparación de los tres esquemas para una call europea vanilla con coeficientes dependientes de \(t\).
- `test_put_europea_t.m`: comparación de los tres esquemas para una put europea vanilla con coeficientes dependientes de \(t\).

## Requisitos

Para ejecutar los ficheros base del método numérico basta con MATLAB base.

Los scripts de comparación que usan funciones cerradas de MATLAB requieren toolboxes adicionales:

- `blsprice` pertenece a **Financial Toolbox** y calcula precios de opciones europeas vanilla bajo Black–Scholes.
<!-- - `barrierbybls` pertenece a **Financial Instruments Toolbox** y calcula precios de opciones barrera europeas bajo Black–Scholes. -->

En el caso de coeficientes dependientes del tiempo, se ha implementado una función auxiliar análoga a `blsprice`, por lo que no es necesario un toolbox adicional para esa comparación.

## Uso básico

### 1. Call europea vanilla con coeficientes constantes
Ejecutar:

`BS/cte/test_call_europea.m`

Este script:

- define el payoff y las condiciones de contorno de una call europea,
- ejecuta los tres esquemas,
- compara el resultado con `blsprice`.

### 2. Put europea vanilla con coeficientes constantes
Ejecutar:

`BS/cte/test_put_europea.m`

Este script:

- define el payoff y las condiciones de contorno de una put europea,
- ejecuta los tres esquemas,
- compara el resultado con `blsprice`.

<!--### 3. Opción barrera knock-out con coeficientes constantes
Ejecutar:

`BS/cte/test_barrera_knockout.m`

Al comienzo del script se puede seleccionar:

- tipo de barrera: `down` o `up`,
- tipo de opción: `call` o `put`.

La comparación se realiza con `barrierbybls`.-->

### 3. Call europea vanilla con coeficientes dependientes del tiempo
Ejecutar:

`BS/dept/test_call_europea_t.m`

Este script:

- define funciones \(r(t)\), \(q(t)\) y \(\sigma(t)\),
- ejecuta los tres esquemas,
- compara el resultado con una función auxiliar análoga a `blsprice`.

### 4. Put europea vanilla con coeficientes dependientes del tiempo
Ejecutar:

`BS/dept/test_put_europea_t.m`

Este script:

- define funciones \(r(t)\), \(q(t)\) y \(\sigma(t)\),
- ejecuta los tres esquemas,
- compara el resultado con una función auxiliar análoga a `blsprice`.

## Notación

El código sigue la notación desarrollada en el texto teórico:

- intervalo espacial: \( S \in [a,b] \),
- paso espacial: \( h = \frac{b-a}{N} \),
- paso temporal: \( \tau = \frac{T}{M} \),
- nodos espaciales: \( S_i = a + ih \),
- nodos temporales: \( t_j = j\tau \),
- aproximación numérica: `F(i+1,j+1) ≈ f(S_i,t_j)`.

## Observaciones

- En el esquema explícito, la estabilidad depende de la elección de la malla.
- En los esquemas implícito y de Crank–Nicolson para coeficientes constantes se reutiliza una factorización LU tridiagonal.
- En el caso de coeficientes dependientes del tiempo, las matrices del sistema cambian con el nivel temporal.
- La comparación con funciones cerradas (blsprice y blsprice_tdep) permite validar la implementación en los casos europeos estudiados, pero no garantiza su correcta implementación.
