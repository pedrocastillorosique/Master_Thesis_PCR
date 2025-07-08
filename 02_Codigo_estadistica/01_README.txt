Esta carpeta está compuesta de los siguientes scripts: 

- Cluster_validation.ipynb, en el que se calcula el ARI sobre los cortes para cada algoritmo, con el estadístico de noramlidad, ANOVA y pareado. Fnialmente se obtiene un diagrama
  de Box-Whisker. 
- TCR_spatial_single.ipynb, que se encarga de unificar los UMIs detectados en la librerá de TCR con aquellos de los que ha sido posible secuenciar la zona CDR3 VDJ del linfocito T.
  Además, sirve para realacionar los clonotipos de la librería de TCR de single cell RNA seq con la librería de TCR de Spatial Transcriptomics.