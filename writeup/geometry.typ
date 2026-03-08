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

#definition[Consideremos una variedad diferenciable $M$, y dos curvas denominadas $phi: RR -> M$ y $psi: RR -> M$, así como una carta en la que ambas tomen coordenadas (en la carta) $phi^U (t)$ y $psi^U (t)$ respectivamente. Ambas se denominan *curvas equivalentes en un punto* sí $phi^U (0) = psi^U (0)$ y $lim_(t -> 0) (phi^U (t) - psi^U (t)) / t = 0$ #cite_mech(81).] <def:equivalentCurves>


#example[Dos curvas equivalentes en la esfera][
  Como ejemplo, consideremos las curvas cuya proyección en la carta estereográfica desde el polo norte $(U, f)$ del @ej:sphereAtlas, en particular

  $
    phi^U (t) = mat(phi^U_1(t); phi^U_2(t)) = mat(
      1 + t;
      1 + t^2;
    ) \
    psi^U (t) = mat(psi^U_1(t); psi^U_2(t)) = mat(
      1 + sin(t);
      2 - cos(t);
    ).
  $

  Notamos que $phi^U (0) = mat(1, 1)^TT$ y $psi^U (0) = mat(1, 1)^TT$, y a su vez

  $
    lim_(t -> 0) (phi^U (t) - psi^U (t)) / t =
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

  Consideremos ahora dos curvas equivalentes, $phi: RR -> S^1$, $psi: RR -> S^1$, tal que $x = phi(0) = psi(0)$. Denominaremos a la proyección de las dos curvas en la primera carta #footnote[Sí $x$ no es representable en la primera carta, se podría usar la segunda equivalentemente.] $phi^U (t)$ y $psi^U (t)$. Recordando la definición de curva equivalente, afirmamos que ambas curvas son equivalentes sí $(phi^U)' (0) = (psi^U) (0) = epsilon^u$, dónde la comilla denota la derivada convencional respecto a un único argumento. Notamos que al ser funciones de una sola variable, las derivadas toman cualquier valor $epsilon^u in RR$, y por lo tanto obtenemos el espacio tangente a $x$ como $T_x S^1 tilde.equiv RR$.

  Definidos los espacios tangentes a cada punto, podemos empezar a construir el fibrado tangente. Primero, construyamos para cada carta el producto cartesiano de esta con su espacio tangente. Al ser espacios planos, estos son trivialmente equivalentes a sus fibrados tangentes:

  $
    T U = {(u, epsilon^u) | u in U, epsilon^u in RR}, \
    T V = {(v, epsilon^v) | v in V, epsilon^v in RR}.
  $

  Por otra parte, podemos construir para cada conjunto de pares una función que nos lleva a puntos de $T S^1$. Definiremos, para $T U$:

  - La *proyección natural* $p: T U -> S^1$ tal que $p(u, epsilon^u) = f(u)$
  - El *diferencial* $dd(f)_u: T_u U -> T_f(u) S^1$, que asigna a cada coordenada del espacio tangente a $U$ un vector tangente de la variedad de forma lineal #footnote[La linealidad garantiza que exista un único $dd(f)_u$, hecho que no demostraremos.].
  - El *pushforward* $dd(f): T U -> T S^1$ que combina #footnote[Hemos introducido un ligero razonamiento circular al utilizar  $T S^1$, que es el objeto que estamos construyendo. Podríamos alternativamente definir $dd(f): T U -> union.sq.big_x T_x S^1$ y posteriormente equipar esta unión disjunta con la topología suave del fibrado tangente, evitando así el argumento circular.] las dos anteriores, es decir, #box[$dd(f)(u, epsilon^u) = dd(f)_u (epsilon^u)$].
  - La carta $(T U, dd(f))$, que representa parte de la variedad $T S^1$.

  Equivalentemente para $B$ podremos definir una proyección natural, un diferencial $dd(g)_v$, y con estas dos el pushforward $dd(g)$, y formar la carta $(T V, dd(g))$. A continuación demostraremos que estas dos cartas forman un atlas para $T S^1$ y por lo tanto que $T S^1$ es una variedad diferenciable.

  Es rápido ver que ambas cartas cubren por completo $T S^1$. Ahora, para que las cartas formen un atlas, deben ser compatibles. Es inmediato ver que la transformación de las proyecciones naturales de $T U$ a $T V$ existe, y es precisamente la función $tau$ previamente definida. Por otra parte, debemos demostrar la compatibilidad de sus diferenciales.

  De la definición de curva equivalente, y el hecho de estar trabajando en una variedad diferencial, podemos obtener que, para dos curvas equivalentes $phi$ y $psi$,

  $
    (phi^U)' (0) = (psi^U)' (0) = epsilon^u & , quad "y" \
    (phi^V)' (0) = (psi^V)' (0) = epsilon^v & .
  $

  Las previas igualdades se pueden explotar para relacionar $epsilon^u$ y $epsilon^v$ si introducimos la relación entre ambas cartas $tau = (g^(-1) comp f)$,

  $
    phi^V (0) & = tau (phi^U (0)).
  $

  Tomando la derivada convencional,

  $
    (phi^V)' (0) = tau' (phi^U (0)) (phi^U)' (0),
  $

  dónde identificamos cada término y obtenemos

  $
    epsilon^v = tau' (u) thick epsilon^u = - epsilon^u / u^2.
  $

  Demostrando por tanto que los diferenciales son compatibles, y que $(T U, dd(f))$ y $(T V, dd(g))$ forman un atlas de $T S^1$.

  Por último, consideremos la composición $dd(tau) = dd(g)^(-1) comp dd(f)$, tal que

  $
    dd(tau): T U -> T V.
  $

  Este objeto que nos lleva de las coordenadas de un vector tangente ($epsilon^u$) en cierta carta, a las coordenadas del mismo vector tangente ($epsilon^v$) en otra carta, se conoce también como pushforward. Es importante también notar que $dd(tau)$ actua sobre el par $(u, epsilon^u)$, no solo sobre $epsilon^u$. Para este ejemplo, podríamos escribir

  $
    dd(tau)(u, epsilon^u) = (1 / u, -epsilon^u / u^2).
  $

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
    & (phi^U)' (0) = (psi^U)' (0) = mat(epsilon^u_1; epsilon^u_2) quad "y" \
    & (phi^V)' (0) = (psi^V)' (0) = mat(epsilon^v_1; epsilon^v_2).
  $

  La relación entre ambas proyecciones de las mismas curvas es la función $tau$, y de forma análoga al ejemplo anterior, pero tomando explícitamente el Jacobiano (al tratarse de una función multivariable), tenemos

  $
    mat(epsilon^v_1; epsilon^v_2) = upright(J)_tau (u_1, u_2) med mat(epsilon^u_1; epsilon^u_2),
  $

  dónde $upright(J)_tau (u_1, u_2)$ es el Jacobiano de $tau$ evaluado en $(u_1, u_2)$.

  Equivalentemente, podríamos definir el pushforward y su actuación

  $
    dd(tau)(mat(u_1; u_2), mat(epsilon^u_1, epsilon^u_2)) = (mat(v_1; v_2), med upright(J)_tau (u_1, u_2) mat(epsilon^u_1; epsilon^u_2)).
  $

] <ej:sphereFibration>

#example[Notación de la velocidad de una curva][

  Consideraremos que la variedad $M$ tiene una carta, $(U, f)$, al igual que su fibrado tangente $T M$, con la carta $(T U, dd(f))$. La curva tendrá entonces su proyección $gamma^U: RR -> U$.

  Vamos a considerar ahora la siguiente definición

  $
    dot(gamma)^U: RR -> T U quad dot(gamma^U)(t) = (gamma^U (t), (gamma^U)' (t)),
  $

  junto con

  $
    dot(gamma): RR -> T M quad dot(gamma)(t) = (gamma(t), dd(f)_(gamma^U (t))((gamma^U)'(t))).
  $

  Si bien la notación de la definición es muy densa, encapsular el punto y su derivada en el mismo objeto nos permite usar la elegante relación

  $
    dot(gamma) = dd(f) comp dot(gamma)^U.
  $

  Es importante recordar que $dot(gamma)^U$ nos lleva a vectores con coordenadas claramente definidas (al ser parte de $U$ que es euclídeo), pero $dot(gamma)$ nos lleva a verdaderos vectores tangentes abstractos.

  Imaginemos ahora una segunda carta para $M$, $(V, g)$, y su respectiva carta para el fibrado tangente $(T V, dd(g))$, relacionados por $tau: U -> V$ y $dd(tau): T U -> T V$. Podemos entonces relacionar

  $
    dot(gamma)^V = dd(tau) comp dot(gamma)^U,
  $

  notación que justifica el nombre de pushforward, $dd(tau)$ "empuja" un vector tangente de una carta a otra. La misma notación se utiliza para transformaciones hacía la variedad diferencial, o incluso dentro de la misma variedad.
]


== Fibrado cotangente

Comencemos esta sección con un ejemplo introductorio, que sentará la motivación necesaria.

#example[Derivada de una función escalar a lo largo de una curva sobre una variedad diferenciable][

  Consideremos una función suficientemente suave definida sobre cierta variedad diferencial #box[$h: M -> RR$], y una curva $gamma: RR -> M$, junto con las previamente introducidas cartas $(U, f)$, $(V, g)$, $(T U, dd(f))$ y $(T V dd(g))$.

  La evaluación de $h$ a lo largo de la curva se puede escribir como $h comp gamma$, y por lo tanto debemos estudiar la derivada (convencional) de esta composición. Para ello, afirmamos que, ya que el valor de $h$ es invariante entre cartas, es decir,

  $
    (h comp gamma) = (h comp f comp gamma^U) = (h comp g comp gamma^V),
  $

  su derivada también lo debe ser,

  $
    (h comp gamma)' = (h comp f comp gamma^U)' = (h comp g comp gamma^V)'.
  $

  Estas igualdades de la derecha son funciones convencionales,

  $
    (h comp f comp gamma^U): RR ->^(gamma^U) underbrace(RR^n, U) ->^(h comp f) RR,
  $

  lo que motiva aplicar la regla de la cadena (si bien lo debemos hacer con coordenadas, utilizando matrices), para obtener

  $
    (h comp f comp gamma^U)'(t) = underbrace(dd((h comp f))(gamma^U (0)), mat(circle, circle, circle)) med underbrace((gamma^U)'(0), mat(circle; circle; circle)).
  $

  Con el fin de que la notación no obstaculize la intuición, vamos a escribir la anterior expresión como

  $
    (h comp f comp gamma^U)' (t) = omega^U med (gamma^U)' (0)
  $

  De forma análoga, se obtiene la misma construcción para la otra carta y podemos plantear (con la misma notación),

  $
    omega^U med (gamma^U)'(0) = omega^V med (gamma^V)' (0)
  $<eq:derivative_invariance>

  Gráficamente, si bien no del todo rigurosamente#footnote[La introducción del producto escalar es innecesaria, todo el desarrollo es posible sin métrica sobre las variedades], podemos entender ambas matrices como "vectores" en las cartas, y el producto como una proyección de uno en el otro.

  #figure(
    image("img/motiv/1forms_motivation_1.svg", width: 80%),
    caption: [Igualdad presentada en la @eq:derivative_invariance de forma visual, para dos cartas de la esfera. Se aprecia que, sí bien los "vectores" que permiten obtener la derivada de la función (colores de la esfera y sus cartas, junto con las líneas de nivel) son muy diferentes, su "producto escalar" es idéntico.],
  ) <fig:1forms_motivation_1>

  Consideremos ahora la relación $tau: U -> V$ entre las cartas (una función multivariable convencional), y su Jacobinao $upright(J) (tau)$, que nos permite escribir $(gamma^V)'(0) = upright(J) (tau) (gamma^U)'$. Tenemos entonces

  $
    omega^U med (gamma^U)'(0) = omega^V med upright(J) (tau) med (gamma^U)'(0),
  $

  ya que la curva $gamma$ elegida, y su derivada, son arbitrarias, podemos eliminarla de ambos lados y relacionar

  $
    omega^U = omega^V med upright(J)(tau),
  $

  es decir,

  $
    omega^V = omega^U med upright(J)(tau)^(-1) = (upright(J)(tau)^TT)^(-1) omega^U
  $

  Notamos que, por contraparte, los vectores tangentes se transforman según

  $
    (gamma^V)' = upright(J)(tau) (gamma^U)'.
  $

  Concluimos por tanto que $omega^U$ y $omega^V$ no son vectores en el sentido estricto, si no que son un nuevo objecto, que se transforma entre cartas de forma diferente, cuyo nombre introduciremos en breve.

]

#definition[Una *1-forma* o *covector* es una transformación lineal que lleva de un espacio vectorial a los reales, es decir, #box[$omega: RR^n -> RR$] #cite_mech(163).] <def:1form>

#note[Las 1-formas forman un espacio vectorial][
  Si definimos la suma de dos 1-formas $omega_1$ y $omega_2$, aplicadas a cierto vector $vb(v)$, como $(omega_1 + omega_2)(vb(v)) = omega_1(vb(v)) + omega_2(vb(v))$ y la multiplicación por un escalar #box[$(lambda omega)(vb(v)) = lambda omega(vb(v))$], podemos considerar que las 1-formas presentan estructura de espacio vectorial de tantas dimensiones como los vectores sobre los que se aplican #cite_mech(163).

  Este espacio vectorial se denomina como el *espacio dual* a $RR^n$, y generalmente se escribe como $(RR^n)^*$.
]

#example[Visualización de las 1-formas][
  Si podemos representar los vectores tangentes como flechas, las 1-formas tienen una representación clara: líneas de nivel @needhamVisualDifferentialGeometry2021a. A su vez, su expresión en coordenadas se puede enteder como una flecha perpendicular a estos, con su longitud proporcional al número de líneas de nivel por unidad de longitud.

  La actuación de la 1-forma en un vector finalmente es proporcional al número de líneas de nivel cortadas por el vector. Notamos que la actuación de la 1-forma sobre un vector es independiente de la carta utilizada.

  Para imaginar el efecto de las transformaciones entre cartas, imaginemos que los vectores y las líneas de nivel de las 1-formas se dibujan sobre una substancia elástica, y la transformación $tau$ deforma esta substancia.

  Representamos el efecto de una transformación tipo escala $tau = mat(2, 0; 0, 1)$:

  #figure(
    image("img/motiv/1form_scale.svg", width: 60%),
    caption: [Efecto de la transformación en un vector (azul) y una 1-forma (líneas de nivel grises). Apreciamos como las coordenadas del vector (en ejes anteriores a la transformación) cumplen $tau$, mientras que las coordenadas de la 1-forma cumplen $(tau^TT)^(-1) = mat(1/2, 0; 0, 1)$ para mantenerse perpendiculares a las líneas de nivel.],
  ) <fig:1forms-visualize-scale>

  Para apreciar el efecto de la transpuesta, consideremos la transformación tipo sesgo $tau = mat(1, -1; 0, 1)$:

  #figure(
    image("img/motiv/1form_sket.svg", width: 60%),
    caption: [En este caso,la transformación sesga el plano, teniendo $(tau^TT)^(-1) = mat(1, 0; 1, 1).$ Las coordenadas están escritas en el plano antes de la transformación.],
  ) <fig:1forms-visualize-skew>

  Es apreciable que tratar a las 1-formas como vectores lleva a confusión, ya que estos se transforman de formas geométricamente muy distintas. En cambio, la visualización como líneas de nivel es mucho más intuitiva.
]

Esta estructura de 1-forma es muy general, y nosotros nos fijamos en su aplicación sobre una variedad diferencial. Al igual que definimos para cada punto de la variedad el espacio tangente, es posible asignar a cada punto de una variedad un espacio de 1-formas denominado el espacio cotangente.

#definition[El *espacio cotangente* a un punto de una variedad diferencial $M$ es el espacio dual del espacio tangente de este y se denomina $T^*_x M$.]

Es decir, para cada punto $x in M$ tenemos un espacio tangente $T_x M$, y a su vez podemos definir el espacio de todas las 1-formas que nos llevan de este mismo espacio tangente a los reales, es decir $omega: T_x M -> RR$, denominado $T^*_x M$.

Equivalentemente, podríamos hablar de que el espacio cotangente es el espacio vectorial formado por todas las "funciones equivalentes" (funciones con mismo Jacobiano en una carta para ese punto) que se pueden definir para ese punto.

#example[Ejemplo de un covector del espacio cotangente][

  Retomemos el ejemplo anterior, notamos que, si bien hemos realizado el procedimiento con coordenadas (usando matrices), la actuación de una matriz fila en una matriz columna no es más que un mapa lineal. Identificamos entonces al objeto $omega^U$, con una 1-forma, es decir,

  $
    omega^U: T_u U -> RR,
  $

  dónde $u = gamma^U (0)$.

  Esta 1-forma pertenece entonces al espacio $T^*_u U$, y podríamos encontrar una infinidad de 1-formas similares.
]

De forma análoga al espacio tangente, podemos construir una variedad diferencial como unión de todos los espacios cotangentes.

#definition[El *fibrado cotangente a una variedad diferenciable M* se denomina $T^* M$ y es la unión de todos los espacios cotangentes de cada punto de esta, y es a su vez una variedad diferenciable con el doble de dimensiones que $M$.]
#example[Fibrado cotangente al círculo][
  Podemos realizar la construcción del fibrado cotangente de la esfera de forma análoga al fibrado tangente de la esfera, si bien es interesante realizar de nuevo los pasos para observar las diferencias entre ambos.

  Como en el @ej:circleFibration, tenemos dos cartas $(U, f)$ y $(V, g)$ para $S^1$ mediante la proyección estereográfica, relacionadas por $tau: U -> V$, junto con las cartas de $T S^1$, que hemos denominado $(T U, dd(f))$ y $(T V, dd(g))$, relacionadas por $dd(tau): T U -> T V$.

  Consideremos ahora una función $h: S^1 -> RR$, con su proyección en ambas cartas $h^U: U -> RR$ y $h^V: V -> RR$. Consideremos que existe un punto $x in S^1$ tal que $f(u) = x$, $g(v) = x$.

  Podemos entonces evaluar la derivada de $h$ en las cartas en este punto, que denominaremos

  $
    omega^u & = (h^U)'(u) in RR \
    omega^v & = (h^V)'(v) in RR.
  $

  Los diferentes valores que pueden tomar $omega^u$ y $omega^v$ son precisamente el espacio cotangente en ese punto.

  Ahora, podemos construir para cada carta de $S^1$ un conjunto formado por el producto de cada punto de la carta con su espacio cotangente,

  $
    T^* U & = {(u, omega^u) | u in U, omega^u in RR}, \
    T^* V & = {(v, omega^v) | v in V, omega^v in RR}
  $

  dónde en el caso del círculo, los covectores son también números reales. Ahora, de forma análoga, podemos introducir

  - La *proyección natural* $p: T^* U -> S^1$ tal que $p(u, omega^u) = f(u)$

]

