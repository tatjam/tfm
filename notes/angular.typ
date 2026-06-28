#import "@preview/bananote:0.1.2": *
#import "@preview/physica:0.9.8": *

#show: note.with(
  title: [Notas sobre incertidumbre en el círculo],
)

= Propagación asumiendo la variable angular como si fuese lineal

#box[
  == Variables del problema

  #abstract[Se introducen los elementos MEE y las ecuaciones que necesitamos de estos para propagar la dinámica del sistema (perturbado, o no) y convertir de vuelta a posiciones euclídeas.]
]

Consideremos el vector de estado para la propagación orbital en coordenadas MEE (modified equinoctial elements), que en relación a los elementos Keplerianos clásicos se escribe @walkerSetModifiedEquinoctial1985

$
  va(x) = vec(p, f, g, h, k, L) = vec(a (1 - e^2), e cos(omega + Omega), e sin(omega + Omega), tan(i / 2) cos(Omega), tan(i / 2) sin(Omega), Omega + omega + nu).
$

El vector se puede subdividir en tres grupos:

- Magnitud escalar, la variable $p$ es directamente proporcional al semi-eje mayor y por lo tanto actua como una variable de escala de la órbita
- Magnitudes "direccionales", las variables $f$, $g$, $h$ y $k$ esencialmente representan direcciones espaciales
- Magnitud angular, la variable $L$ es claramente un ángulo.


La ecuación general de propagación, incluyendo la fuerza del cuerpo central y una perturbación genérica en coordenadas CSN, es @walkerSetModifiedEquinoctial1985 #footnote[Notamos erratas en el manuscrito original, que fueron corregidas en una siguiente publicación de la revista, estas ecuaciones incluyen las correciones]:

$
  dv(p, t) &= (2 p C) / w sqrt(p / mu) \
  dv(f, t) &= sqrt(p / mu) (S sin(L) + (((w+1) cos(L) + f) C) / w - (g (h sin(L) - k cos(L)) N) / w) \
  dv(g, t) &= sqrt(p / mu) (-S cos(L) + (((w+1) sin(L) + g) C) / w - (f (h sin(L) - k cos(L)) N) / w) \
  dv(h, t) &= sqrt(p / mu) (s^2 N) / (2 w) sin(L) \
  dv(L, t) &= sqrt(mu p) (w / p)^2 + sqrt(p / mu) ((h sin(L) - k cos(L)) N) / w
$

dónde $w = 1 + f cos(L) + g sin(L)$ y $s^2 = 1 + h^2 + k^2$, y $C$, $S$ y $N$ son las acceleraciones perturbadoras.

En forma vectorial, escribimos las anteriores ecuaciones como

$
  dv(va(x), t) = F(va(x)) = F(p(t), f(t), g(t), h(t), L(t)).
$

Finalmente, las ecuaciones que nos permiten convertir $va(x)$ a coordenadas euclideanas son #footnote[En el anexo MEEToEuclid.nb se observa la derivación de estas expresiones] @jacobwilliamsDegenerateConicModified

$
  va(x)^e_p = p / (w(1 + h^2 + k^2)) mat(
    (1 + h^2 - k^2) cos(L) + 2 h k sin(L);
    (1 - h^2 + k^2) sin(L) + 2 h k cos(L);
    2(h sin(L) - k cos(L));
  )\
  va(x)^e_v = sqrt(mu / p) / (1 + h^2 + k^2) mat(
    (- & 1 - h^2 + k^2) sin(L) + 2 h k cos(L) + g k^2 + 2 f h k - g (1 + h^2);
    ( & 1 - h^2 + k^2) cos(L) - 2 h k sin(L) + f k^2 - 2 g h k + f (1 - h^2);
    & 2(h cos(L) + k sin(L) + f h + g k);
  )
$

separando en dos subvectores la parte de posición y la de velocidad, tal que $va(x)^e = mat(va(x)^e_p; va(x)^e_v)$.

Escribimos $va(x)^e = X (va(x))$ para referirnos a esta transformación.

#box[
  == Equivalencia física de las trayectorias ante desfases de $2 pi$ en $L$

  #abstract[Se demuestra que solucionar el sistema dinámico con un desfase $L' = L + 2 k pi$ para $k in ZZ$ resulta en trayectorias físicamente equivalentes en el espacio euclideo.]

]

Asumamos que el sistema de ecuaciones tiene solución única, y tomemos dos condiciones iniciales, $va(x_0)$ y $va(y_0)$ tal que


$
  va(y_0) = va(x_0) + mat(0, 0, 0, 0, 0, 2 k pi)^TT
$

para $k in ZZ$, es decir, una trayectoria desfasada un múltiplo de $2 pi$ en $L$. Notamos que $F(p(t), f(t), g(t), h(t), k(t), L(t)) = F(p(t), f(t), g(t), h(t), k(t), L(t) + 2 k pi)$ ya que $L$ solo aparece dentro de seno o coseno. Por lo tanto, concluimos que la solución del sistema cumple

$
  va(y)(t) = va(x)(t) + mat(0, 0, 0, 0, 0, 2 k pi).
$

De igual forma, ya que la transformación a coordenadas euclídeas solo depende del coseno o seno de $L$, el desfase en la solución va a representar el mismo estado físico, es decir

$
  va(y)^e (t) = va(x)^e (t).
$

#box[
  == Equivalencia numérica de las trayectorias ante desfases de $2 pi$ en $L$

  #abstract[Se deduce que, sí deseamos utilizar métodos numéricos convencionales (General Linear Methods), no debemos truncar los ángulos $L$, ya que en caso contrario sería posible obtener resultados incorrectos durante la integración. ]
]

Asumiremos que el sistema numérico respeta

$
  sin(x + 2 pi k) = sin(x) \
  cos(x + 2 pi k) = cos(x),
$

algo que no es cierto en aritmética de punto flotante excepto para $k$ pequeño, pero asumiremos suficientemente cierto para toda aplicación práctica.


=== Métodos GLM

Un GLM (General Lineal Method) es la representación más general de los métodos Runge-Kutta y los métodos multi-paso lineales @butcherGeneralLinearMethods1996. Estos, en su forma más general, utilizan dos conjuntos de vectores

- $r$ vectores de estado históricos del sistema ${va(y)_1^"(n)", va(y)_2^"(n)", ..., va(y)_r^"(n)"}$, escritos en forma de vector columna (es un vector de vectores) como

$
  y^"(n)" = mat(va(y)_1^"(n)", va(y)_2^"(n)", ..., va(y)_r^"(n)")^TT.
$

- $s$ puntos intermedios ${va(Y)_1, va(Y)_2, ..., va(Y)_s}$, escritos también en forma de vector columna (de nuevo un vector de vectores) como

$
  Y = mat(va(Y)_1, va(Y)_2, ..., va(Y)_s)^TT.
$

Los puntos intermedios se utilizan para la evaluación de $F$, formando otro vector de vectores


$
  F(Y) = mat(F(va(Y)_1), F(va(Y)_2), ..., F(va(T)_s))^TT.
$

Una vez escritas todas las magnitudes, la relación entre estas es

$
  va(Y)_i = h A_i F(Y) + U_i y^"(n)", \
  va(y)^"(n+1)"_i = h B_i F(Y) + V_i y^"(n)",
$

dónde $A_i$, $U_i$, $B_i$ y $V_i$ son matrices constantes (normalmente escritas como el producto de Kronecker de otra matriz y la matriz unitaria).

En todo caso, la estimación del siguiente estado depende de la suma de:

- Una combinación lineal de vectores estado anteriores.
- Una combinación lineal de evaluaciones de $F$ en puntos derivados de solucionar una ecuación formada por combinaciones lineales de los vectores estado anteriores y las propias evaluaciones de $F$ (si el método es implícito).

=== Wrap-around en las derivadas

Introducimos ahora el concepto de wrap-around. Consideremos la función $W mat(p, f, g, h, k, L) = mat(p, f, g, h, k, L mod 2 pi)$, afirmamos que

$
  F (sum_i a_i va(x)_i) = F (W(sum_i a_i va(x)_i)) = F (sum_i W(a_i va(x)_i)).
$

Sí $a_i in RR$. Para ello, notamos que

$
  sin(sum_i a_i x_i) = sin((sum_i a_i x_i) mod 2 pi) = sin (sum_i (a_i x_i mod 2 pi)),
$

y de igual forma para el coseno. Ya que la dependencia en $L$ de $F$ ocurre exclusivamente a través del seno y coseno, confirmamos la identidad.

Por otra parte, debemos ser cuidadosos, ya que

$
  sin(sum_i a_i (x_i mod 2 pi)) != sin(sum_i a_i x_i)
$

excepto sí $a_i in ZZ$. Ya que $a_i$ son reales, esta última expresión es incorrecta en su caso general.

Debido a lo anterior, el cálculo de las derivadas será equivalente siempre que se realicen las operaciones de truncado al círculo:
- o bien tras combinar linealmente los puntos
- o bien a cada punto tras ser escalado por su coeficiente

pero nunca realizando el truncado al círculo a los puntos para posteriormente ser escalados.

=== Wrap-around en la integración

Por lo anteriormente descrito, un método numérico GLM convencional no va a causar problemas en la evaluación de $F(Y)$, ya que $Y$ es una combinación lineal, sin truncamiento, de varios puntos.

El método GLM tampoco va a tener problemas a la hora de evaluar $va(y)^"(n+1)"$, ya que de nuevo, este se obtiene por combinación lineal de puntos sin ningún truncamiento.

Por otra parte, *resultaría incorrecto* realizar un truncamiento al círculo de cada vector $va(y)_i^"(n)"$, ya que posteriormente estos valores van a ser utilizados en una combinación lineal, y la suma ponderada de ángulos truncados al círculo no es equivalente al truncamiento al círculo de la suma ponderada de ángulos.

En conclusión, para asegurar que no hay problemas, un integrador convencional deberá operar con ángulos "desenrollados". Si deseamos utilizar ángulos truncados, *no siempre será apto utilizar un integrador convencional*, y sería necesario utilizar un integrador geométrico o alguna otra técnica.


#box[
  == Funciones de varias trayectorias con y sin desenrollamiento

  #abstract[
    Se concluye que integrar con el ángulo desenrollado no es un problema para
    - Los momentos estadísticos direccionales
    -
  ]
]

Consideraremos un sistema multi-trayectoria como un conjunto de soluciones de las ecuaciones dinámicas ${va(x)_i (t)}$, y asumiremos que estos se utilizan para deducir propiedades mediante ciertas funciones de las variables. Por ejemplo, la media se podría calcular para este tipo de sistemas como

$
  EE[va(x)(t)] approx 1 / n sum_(i = 1)^n va(x)_i (t),
$

pero en general se estudiará una propiedad

$
  phi(t) = phi(va(x)_1 (t), va(x)_2 (t), ..., va(x)_n (t)),
$

y asumiremos que los vectores de estado $va(x)_i (t)$ contienen el ángulo "desenrollado".

En su caso más general, $phi(t)$ puede depender de $L$ directamente, como en el ejemplo de la media anteriormente presentado. Pero, son más interesantes aquellas funciones que dependen de $L$ a través del seno y coseno. Por ejemplo, consideremos la media direccional

$
  EE_theta [L(t)] approx "Arg"(1 / n sum_(i = 1)^n e^(i L_i(t))) =
  "atan2"(1/n sum_(i=1)^n sin(L_i (t)), 1/n sum_(i=1)^n cos(L_i (t))),
$

o cualquier función que dependa de la posición euclídea del objeto (debido a que la transformación $X$ depende únicamente de seno y coseno de $L$).

Esto se generaliza a todos los momentos estadísticos, ya que $sin(n x)$ se puede escribir en función de $sin(x)$ y $cos(x)$ siempre que $n in NN$. Por otra parte, sí $n in RR$, esto no se cumple, y por lo tanto el desenrollamiento deberá ser tratado cuidadosamente.

#bibliography("../writeup/refs.bib")
