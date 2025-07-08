# Análisis Exploratorio de Datos de Transcriptómica Espacial para la Identificación de Poblaciones Celulares y Patrones Genéticos en el Microentorno Tumoral

Este repositorio contiene el código empleado para la validación de los algoritmos Louvain, BayesSpace, conGi, GraphST y stKeep; y la implementación de GraphST sobre cuatro muestras secuenciadas espacialmente usando la tecnología Visium 10x, de cáncer de ovario, incluidas en el proyecto [PITAGORAS ](https://cima.cun.es/investigacion/proyecto-pitagoras). Para ello se ha usado la métrica ARI, sobre el clustering en 12 cortes del dataset [DLPFC](https://github.com/LieberInstitute/HumanPilot).

Este Trabajo de Fin de Máster ha tenido los siguientes objetivos: 

- OC1. Realizar una revisión bibliográfica centrada en el microentorno tumoral y en los algoritmos disponibles para el análisis de datos de spaRNA-seq.
- OC2. Comparar el resultado de la identificación de dominios espaciales o clustering obtenido por diferentes algoritmos del estado del arte basados metodologías clásicas, en Machine Learning y Deep Learning, utilizando datos anotados de spaRNA-seq.
- OC3. Usar el mejor algoritmo de clustering resultante de la comparación anterior (OC.2) sobre datos de spaRNA-seq de tumores de ovario, integrando sobre los dominios o clústeres obtenidos la información genética y espacial de clones de receptores de las células T (TCR) en linfocitos infiltrados de tumor (TIL, del inglés, tumor infiltrating lymphocytes), con el fin de identificar patrones relevantes en la distribución celular y expresión génica.

<div align="center">
  <img src="https://github.com/user-attachments/assets/b7617747-0d33-4e52-834e-8de7af194fd1" alt="Ejemplo GraphST" width="500"/>
</div>

A continuación, se muestra un ejemplo de clustering para los cortes 151507, 151510, 151669 y 151676:

<div align="center">
  <img src="https://github.com/user-attachments/assets/24556198-432c-48f1-af48-64d143fd54e5" alt="Clustering por muestra" width="600"/>
</div>

Seguidamente, se aplicó el test de normalidad Shapiro-Wilk, ANOVA (dada su normalidad) y el test pareado de Tukey:

<div align="center">
  <img src="https://github.com/user-attachments/assets/8a330606-981a-450c-9588-42031fee358d" alt="Comparación estadística" width="600"/>
</div>
