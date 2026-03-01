#import "@preview/physica:0.9.8": *
#import "utils.typ": *

#let cite_mech(x) = cite(
  <arnoldMathematicalMethodsClassical1989>,
  supplement: [pág. #x],
)

== Variedades diferenciables

Comenzamos la sección con un repaso breve de varios conceptos que serán necesarios en el resto del trabajo. Si bien las definiciones no tienen un rigor matemático completo, son suficientes para su uso en mecánica clásica. Todas ellas han sido tomadas y adaptadas de @arnoldMathematicalMethodsClassical1989.

#definition[Una *variedad topológica* se define como un conjunto de puntos, junto con una cantidad contable de cartas, tal que todo punto del conjunto esté representado por, al menos, una carta #cite_mech(77)] <def:manifold>

#definition[Una *carta* se puede definir como un subconjunto abierto $U$ de puntos de un espacio de coordenadas euclídeo, junto con una función biyectiva $f$ (uno a uno en ambas direcciones) que permite asignar a cada punto de $U$ un punto de la variedad #cite_mech(77). Representaremos la carta con la notación $(U, f)$.] <def:chart>

Equivalentemente, podemos hablar de una carta como el mapa de un cierto conjunto de puntos de la variedad al espacio euclídeo. Intuitivamente, cada carta nos permite asignar coordenadas a un cierto subconjunto de los puntos de la variedad. Esta última definición es más típica en los textos de matemáticas, mientras que la anterior es más apropiada en mecánica clásica.

#example[Carta trivial para el espacio euclídeo][
  El espacio euclídeo $RR^n = {(q_1, ..., q_n) | q_i in RR}$ es una variedad topológica trivialmente, tomando $U = RR^n$, y el mapa identidad, formando la carta $(RR^n, "Id")$.
]

#example[Carta trivial para la esfera][
  Un ejemplo menos trivial es la esfera $S^2 = {(x_1, x_2, x_3) in RR^3 | x_1^2 + x_2^2 + x_3^2 = 1}$, ya que $S^2$ no es un espacio euclídeo, no podemos utilizar la carta trivial que utilizamos en el ejemplo anterior.
]

#example[Las coordenadas esféricas NO son una carta para la esfera][
  Retomando el ejemplo anterior, consideremos $U = {(theta, phi) in RR^2 | theta in [0, pi), phi in [0, 2 pi)}$, y el mapa
  $
    x(theta, phi) = mat(x_1; x_2; x_3) = mat(
      sin theta cos phi;
      sin theta sin phi;
      cos theta
    ),
  $

  Si bien $x$ hace posible asignar a cada punto del rectángulo en $RR^2$ un punto de la esfera, el sistema no es una carta. Por una parte, incumplimos el requisito de que $U$ sea un subconjunto abierto de $RR^2$. Además, el mapa no es biyectivo. Por ejemplo, para $mat(x_1, x_2, x_3)^TT = (0, 0, 1)^TT$ tenemos $theta = 0$ y $phi in [0, 2 pi)$.
]

#definition[Una *variedad diferenciable* es una variedad topológica que puede ser totalmente cubierta por un atlas #cite_mech(77) #footnote[Convencionalmente, se define la variedad diferencial a partir del atlas, con cierta pérdida del valor intuitivo.]] <def:diff_manifold>

#definition[Un *atlas* es la unión de varias cartas compatibles #cite_mech(78).]
#definition[Dos cartas $(U, f)$ y $(V, g)$ de la misma variedad $M$ son *compatibles* sí, para los subconjuntos $X subset U$ y $Y subset V$, ambos con la misma imagen en $M$, la función que lleva desde cada punto de $X$ a $Y$ es diferenciable #cite_mech(78).]

En nuestro caso de estudio particular, podremos asumir siempre que $X$ y $Y$ son espacios de las mismas dimensiones, digamos $n$, y por lo tanto el mapa de $X$ a $Y$ serán $n$ funciones de $n$ variables. Por lo tanto, que el mapa sea diferenciable no es más que estas funciones sean independientemente diferenciables en cuanto a cada variable.

Además, podemos obtener el mapa notando que las cartas son biyectivas, por lo tanto, para un punto $x in X$, tenemos que $f(x) = m$, siendo $m$ un punto de la variedad, y $g^(-1)(m) = y$, siendo $y in Y$. Por lo tanto, el mapa de $X$ a $Y$ no es más que $g^(-1)(f(x))$, expresado de forma general como $g^(-1) comp f$.

#text(red)[Figura]

#example[Atlas para la esfera][
  El mapa estereográfico de la esfera respecto a un punto $p in S^2$ se genera proyectando cada punto de esta, excepto $p$, sobre el plano tangente a la esfera en el punto $q$ diametralmente opuesto a $p$:

  #text(red)[Figura]

  Reaizando el desarrollo con coordenadas, consideremos el sistema de coordenadas de la figura. Comenzamos con una carta $(U, f)$, formada por la proyección estereográfica de la esfera desde el punto $n = (0, 0, 1)^TT in S^2$ (el polo norte). $f^(-1): S^2 \\{n} -> U$ se puede obtener por trigonometría según

  $
    mat(u_1; u_2) = mat(
      x_1 / (1 - x_3);
      x_2 / (1 - x_3)
    ),
  $

  aplicable para toda $S^2$ menos $n$.

  De forma análoga, podemos obtener una segunda carta $(V, g)$ de igual forma, pero proyectando desde $s = (0, 0, -1)^TT in S^2$ (el polo sur), y obteniendo #box[$g^(-1): S^2 \\ {s} -> V$] según

  $
    mat(v_1; v_2) = mat(
      x_1 / (1 + x_3);
      x_2 / (1 + x_3)
    ).
  $

  Ahora, a excepción de los polos norte (punto $n$) y sur (punto $s$), todo el resto de puntos de $S^2$ son compartidos entre ambas cartas, y por lo tanto, para que estas formen un atlas, el mapa de una carta a otra debe ser diferenciable en todo punto compartido. Previamente hemos deducido que el mapa debe ser $g^(-1) comp f$. Tras algo de manipulación algebraica #footnote[Disponible en el anexo "Sphere.nb"] <annex:sphere> , este mapa queda como

  $
    mat(v_1; v_2) = mat(
      u_1 / (u_1^2 + u_2^2);
      u_2 / (u_1^2 + u_2^2);
    ),
  $

  cuyo Jacobiano @annex:sphere existe en todo punto menos $mat(0, 0)^T$.

  Concluimos por tanto que la esfera es una variedad diferenciable, ya que hemos encontrado un atlas que la cubre por completo.
] <ej:sphereAtlas>

== Fibrado tangente

La introducción del concepto de diferenciabilidad de una variedad es vital, ya que nos permite definir el espacio tangente a un punto de la variedad.

#definition[Consideremos una variedad diferenciable $M$, y dos curvas denominadas $phi: RR -> M$ y $psi: RR -> M$, así como una carta en la que ambas tomen coordenadas (en la carta) $vb(phi)(t)$ y $vb(psi)(t)$ respectivamente. Ambas se denominan *curvas equivalentes en un punto* sí $vb(phi)(0) = vb(psi)(0)$ y $lim_(t -> 0) (vb(phi)(t) - vb(psi)(t)) / t = 0$ #cite_mech(81).] <def:equivalentCurves>


#example[Dos curvas equivalentes en la esfera][
  Como ejemplo, consideremos las curvas cuya proyección en la carta estereográfica desde el polo norte $(U, f)$ del @ej:sphereAtlas, en particular

  $
    vb(phi)(t) = mat(phi_1(t); phi_2(t)) = mat(
      1 + t;
      1 + t^2;
    ) \
    vb(psi)(t) = mat(psi_1(t); psi_2(t)) = mat(
      1 + sin(t);
      2 - cos(t);
    ).
  $

  Notamos que $vb(phi)(0) = mat(1, 1)^TT$ y $vb(psi)(0) = mat(1, 1)^TT$, y a su vez

  $
    lim_(t -> 0) (vb(phi)(t) - vb(psi)(t)) / t =
    lim_(t -> 0) 1 / t mat(
      t - sin(t);
      t^2 + cos(t) - 1
    ) =
    lim_(t -> 0) mat(
      1 - sin(t)/t;
      t + (cos(t) - 1)/t
    ) =
    mat(0; 0).
  $

  Por lo tanto, ambas curvas son equivalentes para $t = 0$. En la @fig:sphere_equivalent_curves se pueden ver representadas ambas curvas en las dos cartas de la esfera y en la propia variedad @annex:sphere.

  #figure(
    image("img/motiv/sphere_equivalent_curves.svg", width: 80%),
    caption: [Dos curvas equivalentes en la carta estereográfica desde el polo norte (izquierda), desde el polo sur (centro) y en la esfera representada dentro de $RR^3$ (derecha).],
  ) <fig:sphere_equivalent_curves>
]

#definition[Un *vector tangente a una variedad diferenciable en un punto* es la clase de equivalencia de todas las curvas equivalentes en ese punto #cite_mech(81).] <def:tangentVector>

#note[Los vectores tangentes forman un espacio vectorial][
  Estas clases de equivalencia son, esencialmente, lo mismo que el vector derivada de las curvas en la carta en la que se proyectan. Por lo tanto,
  cada clase de equivalencia de curvas se asocia a un vector derivada en la carta, y se puede definir adición de estos vectores y multiplicación por escalares, formando un espacio vectorial.
]

#definition[El *espacio tangente a una variedad diferenciable $M$ en un punto $x$* se denomina $T M_x$, y es el espacio vectorial formado por todos los vectores tangentes a ese punto #cite_mech(81).] <def:tangentSpace>

#note[El espacio tangente es euclídeo][
  Debido a la construcción que hemos utilizado para formar el espacio tangente, utilizando las cartas que son espacios euclídeos, el espacio tangente a un punto siempre será un espacio euclídeo. Además, este espacio tendrá las mismas dimensiones que la variedad a la que pertenece.
]

#definition[El *fibrado tangente a una variedad diferenciable M* se denomina $T M$ y es la unión #footnote[La definición rigurosa de esta operación de unión se escapa del alcance del trabajo.] de todos los espacios tangentes de cada punto de esta, y es a su vez una variedad diferenciable con el doble de dimensiones que $M$ #cite_mech(81).]

#example[Fibrado tangente al círculo][
  Como ejemplo directo y fácil de visualizar, consideremos la construcción del fibrado tangente al círculo $S^1$. El círculo se puede representar por dos cartas, de forma análoga a la esfera, mediante proyección estereográfica (asumiendo radio unitario y centro en el origen). Denominemos $(U, f)$ y $(V, g)$ a ambas cartas, con

  $
    u = f^(-1)(x_1, x_2) = x_1 / (1 - x_2), \
    v = g^(-1)(x_1, x_2) = x_1 / (1 + x_2).
  $

  Notamos que estas dos cartas (la primera proyección estereográfica desde el "polo norte" del círculo, la segunda desde el "polo sur") se relacionan por

  $
    (g^(-1) comp f)(u) = tau (u) = v = 1 / u,
  $

  utilizando la identidad $x_1^2 + x_2^2 = 1$ y operando.

  Consideremos ahora dos curvas equivalentes, $phi: RR -> S^1$, $psi: RR -> S^1$, tal que $x = phi(0) = psi(0)$. Denominaremos a la proyección de las dos curvas en la primera carta #footnote[Sí $x$ no es representable en la primera carta, se podría usar la segunda equivalentemente.] $phi^U (t)$ y $psi^U (t)$. Recordando la definición de curva equivalente, afirmamos que ambas curvas son equivalentes sí $dot(phi)^U (0) = dot(psi)^U (0) = epsilon^u$, dónde el punto denota la derivada respecto al único argumento. Notamos que al ser funciones de una sola variable, las derivadas toman cualquier valor $epsilon^u in RR$, y por lo tanto obtenemos el espacio tangente a $x$ como $T_x S^1 tilde.equiv RR$.

  Definidos los espacios tangentes a cada punto, podemos empezar a construir el fibrado tangente. Primero, construyamos para cada carta el producto cartesiano de esta con su espacio tangente,

  $
    A = {(u, epsilon^u) | u in U, epsilon^u in RR}, \
    B = {(v, epsilon^v) | v in V, epsilon^v in RR}.
  $

  Por otra parte, podemos construir para cada conjunto de pares una función que nos lleva a puntos de $T S^1$. Definiremos, para $A$:

  - La *proyección natural* $p: A -> S^1$ tal que $p(u, epsilon^u) = f(u)$
  - El *diferencial* $dd(f): A -> T_u S^1$, tal que $dd(f)(u, epsilon^u)$ representa una familia de equivalencia de curvas que pasan por $u$, es decir $phi^U (0) = psi^U (0) = u$ con $dot(phi)^U (0) = dot(psi)^U (0) = epsilon^u$.
  - La función $F$ que combina las dos anteriores, es decir, asigna a cada $(u, epsilon^u)$ un punto de $M$ (mediante la proyección natural) y una familia de equivalencia de curvas que pasan por ese punto.
  - La carta $(A, F)$ que representa parte de $T S^1$

  Equivalentemente para $B$ podremos definir la función $G$, y formar la carta $(B, G)$. A continuación demostraremos que estas dos cartas forman un atlas para $T S^1$.


  Primero de todo, para que las cartas formen un atlas, deben ser compatibles. Es inmediato ver que la transformación de las proyecciones naturales de $A$ a $B$ existe, y es precisamente la función $tau$ previamente definida. Por otra parte, debemos demostrar la compatibilidad de sus diferenciales.

  De la definición de curva equivalente, y el hecho de estar trabajando en una variedad diferencial, podemos obtener que

  $
    dot(phi)^U (0) = dot(psi)^U (0) = epsilon^u & , quad "y" \
    dot(phi)^V (0) = dot(psi)^V (0) = epsilon^v & .
  $

  Las previas igualdades se pueden explotar para relacionar $epsilon^u$ y $epsilon^v$ si introducimos la relación entre ambas cartas $tau = (g^(-1) comp f)$,

  $
    phi^V (0) & = tau (phi^U (0)).
  $

  Tomando la derivada, obtenemos

  $
    dot(phi)^V (0) = dd(tau) (phi^U (0)) dot(phi)^U (0),
  $

  dónde identificamos cada término y obtenemos

  $
    epsilon^v = dd(tau) (u) thick epsilon^u.
  $

  Notamos el uso notacional, dónde $dd(tau)$ denota un mapa de un espacio tangente a otro, mientras que $dot(phi)$ es una derivada convencional respecto a una variable. En este ejemplo particular, ambas son derivadas de una sola variable. Para nuestras dos cartas,

  $
    epsilon^v = - epsilon^u / u^2,
  $

  concluyendo que ambas son compatibles gracias a las dos funciones $tau$ y $dd(tau)$. Resumiendo de forma gráfica:

  #text(red)[Figura]

] <ej:circleFibration>

#example[Fibrado tangente a la esfera][
  El mismo procedimiento realizado anteriormente es aplicable a la esfera $S^2$. Recordando el ejemplo @ej:sphereAtlas, la función que nos permite pasar de puntos de la carta $(U, f)$ a la carta $(V, g)$ es

  $
    mat(v_1; v_2) = (g^(-1) comp f)(u_1, u_2) = tau(u_1, u_2) = mat(
      u_1 / (u_1^2 + u_2^2);
      u_2 / (u_1^2 + u_2^2)
    ).
  $

  Por su parte, si consideramos dos curvas equivalentes que pasan por un cierto punto $x = phi(0) = psi(0)$, cuya representación en ambas cartas es denominado $u$ y $v$ respectivamente, tendremos

  $
    & dot(phi)^U (0) = dot(psi)^U (0) = mat(epsilon^u_1; epsilon^u_2) quad "y" \
    & dot(phi)^V (0) = dot(psi)^V (0) = mat(epsilon^v_1; epsilon^v_2).
  $

  La relación entre ambas proyecciones de ls mismas curvas es la función $tau$, y de forma análoga al ejemplo anterior, tenemos

  $
    mat(epsilon^v_1; epsilon^v_2) = dd(tau) (u_1, u_2) med mat(epsilon^u_1; epsilon^u_2),
  $

  dónde $dd(tau) (u_1, u_2)$ es el Jacobiano de $tau$.

  Al contrario de el caso del círculo, no es sencillo representar la variedad $T S^2$, ya que presenta cuatro dimensiones, pero si que podemos representar $tau$ y $dd(tau)$ para un solo punto:

  #text(red)[Figura]


] <ej:sphereFibration>

== Fibrado cotangente

Comencemos esta sección con un ejemplo introductorio, que sentará la motivación necesaria.

#example[Gradiente de una función escalar sobre una variedad][

  Consideremos ahora una función suficientemente suave definida sobre cierta variedad diferencial #box[$h: M -> RR$], y una curva $gamma: RR -> M$. Estudiemos la tasa de cambio de $h$ a medida que avanzamos en $gamma$, para un punto $t = 0$. Para indicar esta derivada utilizaremos la siguiente notación, dando énfasis en que se trata de una derivada convencional en cuanto a una única variable:

  $
    (h comp gamma)' (0) = lim_(t -> 0) ((h comp gamma)(t) - (h comp gamma)(0)) / t.
  $

  Tomemos dos cartas para la variedad diferencial, $(U, f)$ y $(V, g)$, compatibles y relacionadas por un mapa $tau: U -> V$, con su Jacobiano $dd(tau)$. La proyección de la curva $gamma$ en ambas cartas será entonces $gamma^U: RR -> U$ y $gamma^V: RR -> V$ respectivamente. La previamente introducida derivada se puede calcular de forma similar en ambas cartas según

  $
    (h comp gamma)' (0) = (h comp f comp gamma^U)' (0) = (h comp g comp gamma^V)' (0).
  $

  Podemos aplicar la ley de la cadena, por ejemplo, para la primera carta

  $
    (h comp f comp gamma^U)' (0) = dd((h comp f)) (u) med dot(gamma)^U (0),
  $

  dónde $u = f^(-1) (gamma (0))$. Análogamente, se puede obtener la misma relación para la segunda carta, quedando

  $
    dd((h comp f)) (u) med dot(gamma)^U (0) = dd((h comp g)) (v) med dot(gamma)^V (0)
  $ <eq:derivative_invariance>

  Podemos representar la anterior expresión gráficamente:

  #figure(
    image("img/motiv/1forms_motivation_1.svg", width: 80%),
    caption: [Igualdad presentada en la @eq:derivative_invariance de forma visual, para dos cartas de la esfera. Se aprecia que, sí bien los "vectores" que permiten obtener la derivada de la función (colores de la esfera y sus cartas, junto con las líneas de nivel) son muy diferentes, su "producto escalar" es idéntico.],
  ) <fig:1forms_motivation_1>

  Ahora, podemos relacionar las magnitudes a ambos lados de la @eq:derivative_invariance mediante la función $tau$. Recordando los anteriores ejemplos, tenemos que

  $
    dot(gamma)^U (0) = dd(tau) med dot(gamma)^V (0),
  $

  por lo que, sustituyendo en la @eq:derivative_invariance,

  $
    dd((h comp f))(u) med dd(tau) med dot(gamma)^V (0) = dd((h comp g))(v) med dot(gamma)^V (0).
  $

  Ya que la curva $gamma$ elegida es arbitraria, podemos concluir que

  $
    dd((h comp f))(u) med dd(tau) = dd((h comp g))(v), quad "es decir" \
    dd((h comp f))(u) = dd((h comp g))(v) med (dd(tau))^(-1)
  $

  Es aparente que, si bien los dos objetos que hemos dibujado como flechas en la @fig:1forms_motivation_1 son similares, son fundamentalmente diferentes, ya que la proyección del vector tangente $dot(gamma)^U$ se obtiene a partir de la otra proyección, $dot(gamma)^V$ a través de $dd(tau)$, pero $dd((h comp f))(u)$ se obtiene a partir de $dd((h comp g))(v)$ a partir de $dd(tau)^(-1)$.

  Estas nuevas entidades matemáticas "que hemos descubierto" no son vectores tangentes, ni pertenecen a $T M$, si bien están íntimamente relacionados con estos. En lo que sigue daremos nombre a esta entidad.
]

