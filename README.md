# Estudio Bioinformático de la Variación de Isoformas Génicas en Cáncer de Mama (TCGA-BRCA)

Este repositorio contiene el pipeline computacional completo desarrollado como parte de mi **Trabajo de Fin de Grado (TFG) en Biología** para la **Universidad de Sevilla** (Departamento de Bioquímica Vegetal y Biología Molecular). 

* **Autor:** Abdelfatah Oumouss Aitmoussa
* **Tutor:** Dr. José Manuel García Heredia

---

## 🔬 Descripción del Proyecto
El cáncer de mama es una enfermedad con una alta heterogeneidad celular y transcriptómica. Tradicionalmente, los análisis genómicos evalúan cambios cuantitativos a nivel global de gen, ignorando alteraciones cualitativas cruciales como el *switching* (cambio de dominancia) de isoformas debido al splicing alternativo.

El objetivo general de este proyecto es implementar un flujo bioinformático robusto para caracterizar el panorama de splicing alternativo utilizando datos de secuenciación masiva (RNA-seq) del consorcio **The Cancer Genome Atlas (TCGA)**, aislando el ruido genético individual mediante un enfoque de pareado estricto paciente-control e identificando biomarcadores con potencial traslacional clínico (como el modelo genético del gen `DUT`).

---

## 🛠️ Estructura del Pipeline Bioinformático (`.Rmd`)
El script principal `TFG_TCGA_BRCA_Analysis.Rmd` está estructurado en bloques independientes y reproducibles que ejecutan las siguientes tareas moleculares y estadísticas:

1. **Carga y Preprocesamiento de Datos Brutos:** Carga de matrices de expresión analizadas por transcritos (procedentes de la plataforma Firebrowse del Broad Institute / formato RSEM normalizado).
2. **Clasificación y Filtro de Muestras:** Separación de muestras moleculares en tejido tumoral primario (`código 01`) y tejido normal adyacente (`código 11`) utilizando el barcode identificativo de TCGA.
3. **Pareado Estricto de Pacientes:** Identificación de pacientes con ambas condiciones para reducir la varianza interindividual. Transformación logarítmica $\log_2(x + 1)$ para la estabilización de la varianza.
4. **Análisis de Expresión Diferencial (DEA) a Nivel de Isoformas:** Identificación de isoformas sobreexpresadas y reprimidas aplicando umbrales asimétricos basados en un cambio biológico del $\pm40\%$ ($\log_2(1.4) = 0.485$ y $\log_2(0.6) = -0.737$) con control de falsos positivos mediante False Discovery Rate (FDR < 0.05).
5. **Estudio Dirigido y Análisis de Supervivencia (Modelo `DUT`):** Evaluación del impacto pronóstico del desbalance de isoformas en el gen `DUT` (comparando la dominancia de la isoforma `DUT-M` vs `DUT-N`) mediante curvas de supervivencia acumulativa de **Kaplan-Meier** y el test estadístico de **Log-Rank**.
6. **Análisis de Enriquecimiento Funcional (GSEA):** Mapeo de IDs biológicos a NCBI Entrez Gene e inferencia de implicaciones moleculares mediante Gene Ontology (GO) para las tres ontologías primarias (Procesos Biológicos, Componentes Celulares y Funciones Moleculares) utilizando `clusterProfiler`.

---

## 📊 Requisitos e Instalación

El código se ejecuta bajo el entorno de programación **R (v4.3.3 o superior)** en **RStudio**. El propio script está optimizado para verificar e instalar de forma automática los paquetes que falten de la red CRAN, no obstante, requiere que tengas configurado tu entorno correctamente.

### Librerías principales utilizadas:
* **Manipulación de datos:** `dplyr`, `readr`
* **Estadística molecular:** `diptest`, `survival`
* **Visualización de alta calidad:** `ggplot2`, `ggrepel`, `survminer`
* **Bioinformática de redes (Bioconductor):** `clusterProfiler`, `org.Hs.eg.db` *(Nota: asegúrate de tener instalado BiocManager)*

```R
# Para instalar los paquetes de Bioconductor de forma manual si es necesario:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
