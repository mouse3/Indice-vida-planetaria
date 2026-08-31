# Índice de Vida Planetaria
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
$$W(E, E_b)=\left(\frac{1}{1+e^{-\beta_1(E-E_{exc})}}\right)\cdot\left(\frac{1}{1+e^{\beta_1(E-E_{b})}}\right)$$
Donde $\beta_1$ y $\beta_2$ son parámetros de suavizado térmico/cuántico que determinan la pendiente de transición en los umbrales:
- $\frac{1}{1 + e^{-\beta_1 (E - E_{\text{exc}})}}$ representa la probabilidad de activación fotosintética o metabólica para $E \ge E_{\text{exc}}$.
- $E_{exc}$ energía mínima de excitación electrónica (HOMO-LUMO)
- $\frac{1}{1 + e^{\beta_2 (E - E_b)}}$ representa el factor de supervivencia de la molécula frente a la fotodisolución para $E \le E_b=tabla$.

Ahora bien, como $W(E, E_b)$ depende de la energía de fotodisolución, concluimos que esta ventana fotoquímica óptima se define como la suma de las ventanas fotoquímicas de cada enlace de la tabla que expusimos:
$$\left[ \frac{Prob}{m^2}\right]=1-\prod_i{(1-(W_i(E, E_{b_i}))}=W_T(E/A)\in[0, 1]\subset\mathbb{R}$$
Donde $\beta_1$ y $\beta_2$, al ser parámetros de suavizado, se convierten en constantes universales para todo este sistema.
Cabe mencionar que, al tratarse de una suma, habrá regiones $E$ donde $\left[ \frac{Prob}{m^2}\right]\geq 1$ y, es justo por esto, por lo que deberíamos incluir una función de normalización porque una de tipo "corte" (que deja techos planos y pendientes muy inclinadas) es un tanto contraintuitivo.

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
