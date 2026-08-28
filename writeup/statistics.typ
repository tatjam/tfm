#import "@preview/physica:0.9.8": *
#import "utils.typ": *

== Variables aleatorias

El tratamiento matemático de los eventos aleatorios requiere asignar números a los posibles resultados de una situación aleatoria.

#definition[Sí denominamos $Omega$ al espacio de posibles resultados aleatorios, una *variable aleatoria $X = X(omega)$* es una función que asigna a cada $omega in Omega$ un valor numérico @xiuNumericalMethodsStochastic2010.]

Imaginemos el caso de tirar un dado, tenemos $Omega = {⚀, ⚁, ⚂, ⚃, ⚄, ⚅}$. Podriamos estudiar la variable aleatoria que asigna a cada resultado el valor de la cara del dado. Pero también podemos estudiar otra clase de eventos, como "el resultado es impar" o "el resultado está entre 2 y 5".

La generalización de estos se denomina una σ-álgebra que denominaremos con la letra $cal(F)$.

#definition[Una *σ-álgebra* es un conjunto de subconjuntos de $Omega$, tal que @xiuNumericalMethodsStochastic2010

  - $emptyset in cal(F)$ y $Omega in cal(F)$
  - Sí $A in cal(F)$ entonces $A^c in cal(F)$, siendo $A^c$ el complementario del evento $A$. Es decir, si la σ-álgebra contiene un evento, debe también contener su contrario.
  - Sí $A_1, A_2, ... in cal(F)$ , entonces $union.big_i A_i in cal(F)$ y $inter.big_i A_i in cal(F)$. Es decir, si la σ-álgebra contiene un cierto conjunto de eventos, deberá también contener su unión e intersección.
]


Con el ejemplo del dado, tenemos tres σ-álgebra elementales,

- $cal(F)_1 = {emptyset, Omega}$, que únicamente contiene el resultado "el dado tuvo un valor entre 1 y 6" o el elemento contrario.
- $cal(F)_2 = {emptyset, Omega, A_1, A_1^c}$, es decir, el mismo conjunto de antes pero aumentado con otro subconjunto de los resultados, por ejemplo, "el dado marca 3" o "el dado marca un numero impar", junto con su complementario.
- $cal(F)_3 = {A | A subset Omega}$, es decir, todos los posibles subconjuntos de de resultados. Este último conjunto se denomna el "conjunto potencia" de $Omega$, y se escribe $cal(F)_3 = 2^Omega$ @xiuNumericalMethodsStochastic2010.

Con un concepto básico de σ-álgebra podemos ya definir la probabilidad de un evento.

#definition[Un *espacio de probabilidad* es una tupla $(Omega, cal(F), P)$ dónde @xiuNumericalMethodsStochastic2010

  - $cal(F)$ es una σ-álgebra cualquiera de $Omega$
  - Para todo elemento $A in cal(F)$, existe $P(A): cal(F) -> [0, 1]$, que además cumple $P(Omega) = 1$ y $P(union.big_i A_i) = sum_i P(A_i)$
]

Sí además incluimos una variable aleatoria $X(omega)$

== Expansión en caos polinómico
