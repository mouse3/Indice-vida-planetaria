Este proyecto trata sobre construir un modelo matemático-físico que exprese la probabilidad de que exista o se cree vida en un planeta dados los datos radiométricos y aproximaciones de cuerpo negro de las emisiones fotónicas de la estrella que orbita.

La tierra absorbe fotones de alta energía y baja entropía (la energía es inversamente proporcional a la longitud de onda). Lo que hacen los organismos vivos es tomar esta energía en forma de fotones, aprovecharla para realizar procesos fisicoquímicos y, tras esto, liberar energía de manera inevitable en forma de calor, que serían fotones IR (infrarrojos) de baja energía y alta entropía. 

La vida se encuentra entre el límite inferior y el superior de una zona del espectro electromagnético (EM), el inferior lo delimita la energía necesaria para formar una reacción química y el superior lo limita la energía para romper un enlace dado.

Este rango se encuentra en la luz visible: los IR hacen vibrar a las moléculas pero son incapaces de formar reacciones termoquímicas y, los rayos X, tienen tanta energía que rompen los enlaces químicos.

Separamos a la vida en 2 grupos, la vida en la tierra se define dentro del grupo CHONPS y otra alternativa basada en Silicio, que es completamente hipotética.

Necesitamos obtener las energías de formación de los enlaces expresados en eV/enlace, en kJ/mol y en nm de fotón:

| Enlace | Contexto       | Rotura eV/en. | Rotura kJ/mol | Rotura en $\lambda$ (1240/eV) en nm |
| ------ | -------------- | ------------- | ------------- | ----------------------------------- |
| O-H    | Alcoholes, ac. | 4.79          | 463           | 258.9                               |
| O-H    | simple (agua)  | 4.81          | 464           | 257.8                               |
| C-C    | glúcidos       | 3.61          | 348           | 343.5                               |
| C-C    | simple         | 3.61          | 348           | 343.5                               |
| C-N    | aminas         | 3.16          | 305           | 392.4                               |
| C-N    | simple         | 3.16          | 306           | 392.4                               |
| P-O    | nucleótidos    | 3.64          | 351           | 340.7                               |
| P-O    | simple         | 3.64          | 351           | 340.7                               |
|        |                |               |               |                                     |
| Si-Si  | simple         | 2.34          | 226           | 529.9                               |
| Si-O   | simple         | 4.83          | 466           | 256.7                               |
| Si-C   | simple         | 3.18          | 307           | 389.9                               |


Tomando como "simple" el enlace dado con hidrógenos en los radicales sobrantes.

La energía que recibe un planeta es inversamente proporcional a la distancia con la que se encuentra a su estrella y proporcional al radio del planeta. $$A=4\pi r^2$$
Tenemos una función $$E_{espectral}(\lambda)=\frac{2hc^2}{\lambda^5}\frac{1}{e^{\frac{hc}{\lambda k_B T}}-1}$$
Donde la energía total emitida es una constante:
$$\int_0^\infty{E_{espectral}(\lambda)}d\lambda=\frac{2\pi^5k_B^4}{15c^2h^3}T^4$$


Un organismo puede vivir dentro de un umbral de temperaturas, ni tan caliente como para que su entropía interna sea muy grande ni tan pequeña como para que dejen de haber reacciones químicas.
$$k=AT^me^{-\frac{E_a}{RT}}$$
Donde $E_a$ es la energía de reacción mínima para que ocurra una reacción indicada.

El flujo de radiación se define como:
$$E/A=\frac{\int_a^bE_{espectral}(\lambda)d\lambda}{4\pi r^2}=\left[\frac{eV}{m^2}\right]$$
Donde $a$ y $b$ son las delimitaciones del espectro, desde que $\lambda_a$ a $\lambda_b$ se mide el flujo de radiación 


La probabilidad de que un fotón de energía E impulse una reacción fotoquímica sin destruir la estructura molecular se define como:
$$W(E, E_b)=\left(\frac{1}{1+e^{-\beta_1(E-E_{exc})}}\right)\cdot\left(\frac{1}{1+e^{\beta_2(E-E_{b})}}\right)$$
Si $\beta_1\neq\beta_2$ entonces $W(E, E_b)$ es una distribución con sesgo.
Donde $\beta_1$ y $\beta_2$ son parámetros de suavizado térmico/cuántico que determinan la pendiente de transición en los umbrales
1. $\frac{1}{1 + e^{-\beta_1 (E - E_{\text{exc}})}}$ representa la probabilidad de activación fotosintética o metabólica para $E \ge E_{\text{exc}}$.
2. $E_{exc}$ es energía mínima de excitación electrónica (HOMO-LUMO)
3. $\frac{1}{1 + e^{\beta_2 (E - E_b)}}$ representa el factor de supervivencia de la molécula frente a la fotodisolución para $E \le E_b=tabla$.

Se puede ver un ejemplo en este [link a desmos](https://www.desmos.com/calculator/7wxwj4henc?lang=es).

Ahora bien, como $W(E, E_b)$ depende de la energía de fotodisolución, concluimos que esta ventana fotoquímica óptima se define como la suma de las ventanas fotoquímicas de cada enlace de la tabla que expusimos:
$$\left[ \frac{Prob}{m^2}\right]=1-\prod_i{(1-(W_i(E, E_{b_i}))}=W_T(E/A)\in[0, 1]\subset\mathbb{R}$$
Donde es intuitivo que  $\beta_{1_i} = \beta_{2_i}$ por lo que mencionamos anteriormente, pero esta afirmación no es más que una hipótesis.

Si cada cuanto no transporta suficiente energía, un electrón no puede escapar del átomo por muchos cuantos que incidan sobre él; pero, si un electrón absorbe un cuanto suficientemente energético, escapa con una energía cinética que es, como máximo, igual a la diferencia entre la energía transportada por el cuando absorbido y la función del trabajo:
$$E_{c-max.}=h\nu-W_0$$

También tendríamos que definir la entropía de un flujo de fotones.
El flujo de energía de un cuerpo negro se  expresa por la ley de Stefan-Boltzmann $$J_E=\sigma T_{estrella}^4$$
El flujo de entropía $J_S$ se relaciona directamente con el flujo de energía tal que: $$J_S=\frac{3}{4}\frac{J_E}{T_{estrella}}$$
El balance Negentrópico (como bien bautizó Schrödinger) se define como $E_{in}=E_{out}$ donde las entropías de estos dos flujos son completamente distintas:
$$S_{in}=\frac{3}{4}\frac{E}{T_{estrella}}$$
$$S_{out}=\frac{3}{4}\frac{E}{T_{planeta}}$$
Como $T_{estrella}>T_{planeta}$ , entonces $S_{out}>S_{in}$. La diferencia exacta es la producción de entropía disponible se define como $$\Delta S=\frac{4E}{3}\left( \frac{1}{T_{planeta}}-\frac{1}{T_{estrella}}\right)$$

Entonces, con lo que tenemos ahora, la ventana fotoquímica $W_T$ expresa el modelo probabilístico de que los fotones puedan interactuar con la materia orgánica y este diferencial de entropía establece el límite energético de los procesos fisicoquímicos que la vida puede realizar en este planeta antes de irradiar calor al espacio.



Bien, ahora que tenemos el modelo base estructurado matemáticamente, queda visualizar otras variables independientes que influyen en la creación de vida en otros planetas:
1. El área que queda en permanente oscuridad y, por ende, los ángulos de rotación del planeta respecto a la órbita y su excentricidad .
2. La temperatura media (debido al efecto invernadero y albedo).

![[Imagen1.jpeg]]

Descubrimos la relevancia del primero con un experimento mental: dado que sabemos la existencia del prisma y el efecto refractario de la atmósfera, además del diagrama (no a escala) de la imagen anterior sabemos que, si $\theta=0 \vee 180º, \phi=0 \vee 180º$ entonces habrá una zona que estará en oscuridad persistente y, entre la zona donde es completamente de día y la zona que es completamente de noche, habrá una franja de habitabilidad. Pues en esta franja solo entrará las ondas más cortas de la luz de la estrella.

Entonces, dado esto, tenemos el siguiente problema:
"Halla el plano que contiene a la recta que contiene al vector
$$v=\left(\frac{1}{\tan \left(\phi \right)\sqrt{1+\frac{1}{\tan ^{2}\left(\phi \right)}+\tan ^{2}\left(\theta \right)}},\ \frac{1}{\sqrt{1+\frac{1}{\tan ^{2}\left(\phi \right)}+\tan ^{2}\left(\theta \right)}},\ \frac{\tan \left(\theta \right)}{\sqrt{1+\frac{1}{\tan ^{2}\left(\phi \right)}+\tan ^{2}\left(\theta \right)}}\right)$$
$\phi\to0 \therefore \tan(\phi),\tan^2(\phi)\to \infty$
y también contiene los puntos que cortan el plano ($x+z=z$) con la esfera (la figura geométrica) $x^{2}+y^{2}+z^{2}=R$ que, además, corta con el eje Y (la tangente a la trayectoria)"

Visualizándolo: Hay que hallar el plano $alpha$ (morado) 

![[Pasted image 20260901010016.png]]

1. Los puntos de corte del plano $x+z=z$ y la esfera se encuentran dentro del plano proyección, tal que $x, z=0$ y que, por tanto $y^2=R \therefore y_{0,1}=\pm \sqrt{R}$
2. Definimos los vectores punto: $(0, y_0, 0); (0, y_1, 0)$
3. Simplificando para el vector de rotación, nos queda que  $$\bar{v}=R(cot(\phi), 1, \tan(\theta))$$
4. El vector normal al plano es $v \times OP$, sabiendo que $z_{0,1}=0$, la componente y se anula: $$n=(-y_0\tan(\theta), 0, y_0\cot(\phi)$$
5. teniendo el origen en (0, 0, 0) y el nuevo vector normal, la ecuación general queda tal que $$-y_0\tan(\theta)x+y_0\cot(\phi)z=0$$
6. Simplificando un poco nos queda que $$-x\tan(\theta)+z\cot(\phi)=0$$
Seguimos explicando que cualquier punto inicial en la superficie del planeta expresado como $r_0=(x_0, y_0, z_0)$ cambia de posición con el tiempo gracias al guro sobre su eje. Normalizando el vector de rotación $\hat{v}=\frac{v}{|v|}$, sabemos que la posición del punto en función del tiempo de rotación $t_{rot}$ se calcula mediante la fórmula de rotación de Rodrigues alrededor del eje $\hat{v}$ con un ángulo barrido $\theta(t_{rot})=\omega_{rot}t_{rot}$, es decir $$r(t_{rot})=r_0\cos(\theta)+(\hat{v}\times r_0)\sin(\theta)+\hat{v}(\hat{v}\cdot r_0)(1-\cos(\theta))$$
que nos permite saber en qué coordenada x se encuentra un punto en cada instante para saber si cruza (o no) el plano del terminador x=0.
Como, por lo general, se cumple que $e>0$ (excentricidad), tenemos que incorporar la traslación orbital con la rotación; tomando en cuenta que el eje de rotación del planeta mantiene su orientación fija en el espacio a medida que recorre su órbita elíptica (es justo por esto por lo que existen las estaciones).
el vector unitario de rotación $\hat{v}$ permanece constante a medida que el planeta se desplaza por la elipse $R_{centro}(t_{orbit})$, donde $$R_{total}(t_{rot}, t_{orbit})=R_{centro}(t_{orbit})+r(t_{rot})$$
Para parametrizar $R_{centro}(t_{orbit})$ en función del tiempo orbital, situamos a la estrella en el origen del sistema inercial $(0, 0, 0)$. En una órbita elíptica con semieje mayor $a$ y excentricidad $e$, la posición del centro del planeta se expresa mediante la anomalía verdadera $f(t_{orbit})$
$$R_{centro}(t_{orbit})=(r(t)\cos(f), r(t)\sin(f), 0)$$ donde la distancia en cada instante se expresa como $r(t)=\frac{a(1-e^2)}{1+e\cos(f)}$ y la relación entre la posición y el periodo orbital se resuelve usando la ecuación de Kepler $$M=\frac{2\pi}{T_{orbit}}t_{orbit}=E-e\sin(E)$$
donde se relaciona con la anomalía a través de $$\tan{\frac{f}{2}}=\sqrt{\frac{1+e}{1-e}}\tan{\frac{E}{2}}$$
Entonces, la rotación (local) del punto sobre el eje fijo $\hat{v}$ avanza según el periodo rotacional $T_{\text{rot}}$, cuya velocidad angular es $\omega_{\text{rot}} = \frac{2\pi}{T_{\text{rot}}}$ y el ángulo barrido es $\theta(t_{\text{rot}}) = \omega_{\text{rot}} t_{\text{rot}}$ dentro de la fórmula de Rodrigues.

Al unificar ambos movimientos, cualquier punto de la superficie traza una trayectoria espiralada tridimensional alrededor del Sol. Para determinar si queda atrapado en la oscuridad perpetua, evaluamos si el vector de línea de visión desde dicho punto hacia el disco estelar (de radio $R_s$) intersecta permanentemente el volumen de sombra o umbra proyectado por el propio planeta (de radio $R_p$) a lo largo de todo el ciclo de rotación $T_{\text{rot}}$ para cada punto de la órbita $T_{\text{orbit}}$.

El área de la sombra es un casquete esférico cuya área se define como $$A=2\pi R^2(\cos(\rho))$$
El radio angular $\rho$ está determinado por la inclinación del vector rotación respecto a la dirección de los rayos de luz. A partir del coseno del ángulo calculado previamente, tenemos que $$\cos{\psi(t_{orbit})}=\frac{v(-R_{centro}(t_{orbit}))}{|v||R_{centro}(t_{orbit})}$$
y que $$\rho_{perm}=\frac{2}{\pi}-\max_{t_{orbit}}{\left|\frac{\pi}{2}-\phi(t_{orbit})\right|}$$

ergo, como $v=(\cot(\phi), 1, \tan{\theta})$ $$\rho_{perm}=\frac{2}{\pi}-\max_{t_{orbit}}{\left|\frac{\pi}{2}-\arccos\left({\frac{-\cot({\beta})R_x(t_{orbit})-R_y(t_{orbit})-\tan({\alpha})R_z(t_{orbit})}{\sqrt{\cot^2({\phi})+1+\tan^2(\theta)}|R_{centro}(t_{orbit})|}}|\right)\right|}$$
$$A=2\pi R^2\left( 1- \cos\left({\frac{2}{\pi}-\max_{t_{orbit}}{\left|\frac{\pi}{2}-\arccos\left({\frac{-\cot({\beta})R_x(t_{orbit})-R_y(t_{orbit})-\tan({\alpha})R_z(t_{orbit})}{\sqrt{\cot^2({\phi})+1+\tan^2(\theta)}|R_{centro}(t_{orbit})|}}\right)\right|}}\right)\right)$$
donde $R_{centro}(t_{orbit})=(R_x, R_y, R_z)$ describe la trayectoria elíptica en función del semieje mayor $a$, la excentricidad $e$ y la anomalía verdadera $f(t_{orbit})$



Una vez tenemos esto, podemos declarar que, solo en caso de que una estrella sobrepase la temperatura óptima, puede haber vida en las zonas circundantes a la zona de oscuridad permanente. Pero, en caso de que la estrella tenga una temperatura por debajo de los niveles óptimos, no podrá ser plausible la vida bajo ninguna circunstancia debido a la falla en la ventana fotoquímica.



Ahora bien, en cuanto a la temperatura, es muy conocida la ecuación de Johnson-Lewin , que define la tasa de actividad biológica $r(T)$ y la temperatura óptima resultante $$r(T) = \frac{A e^{\left(-\frac{E_a}{R T}\right)}}{1 + C e^{\left(-\frac{\Delta H_{\text{denat}}}{R T}\right)}}$$

Pero, para descripciones microbiológicas y/o astrobiológicas más completas, se suele utilizar la ecuación de Sharpe-Schoolfield:

$$r(T) = \frac{r_0 \cdot \frac{T}{T_{\text{ref}}} e^{\left[\frac{\Delta H_a^\ddagger}{R}\left(\frac{1}{T_{\text{ref}}} - \frac{1}{T}\right)\right]}}{1 + e^{\left[\frac{\Delta H_L}{R}\left(\frac{1}{T_L} - \frac{1}{T}\right)\right]} + e^{\left[\frac{\Delta H_H}{R}\left(\frac{1}{T_H} - \frac{1}{T}\right)\right]}}$$

Donde:
- $\Delta H^\ddagger_a$: entalpia de activación que modela el incremento (exponencial) de las colisiones efectivas en las reacciones químicas mientras la temperatura aumente
- $\Delta H_L$: entalpia asociada con el factor de atenuación en la inactivación de macromoléculas a bajas temperaturas
- $\Delta H_H$: entalpia de desnaturalización a altas temperaturas
- $R$: constante de los gases ideales
- $T_L$: temperatura de transición donde ocurre la inactivación de la enzima por bajas temperaturas
- $T_H$: temperatura de transición donde se da lugar a la inactivación térmica por altas temperaturas
- $T_{ref}$: temperatura de referencia para normalizar la tasa de actividad base

Entonces, de manera computacional, se han de insertar todas estas variables para cada tipo de reacción, por lo general: hidrolisis de ATP, la fijación del CO2 en rubisco, la nitrogenasa y la fermentación glucolítica.


De momento, esta será la base matemática para este índice