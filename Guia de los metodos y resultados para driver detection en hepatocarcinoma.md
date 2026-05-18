# Guia de los metodos y resultados para driver detection en hepatocarcinoma

Se analizó un archivo `.maf` (Mutation Annotation Format), el cuál es un archivo que contiene anotaciones de 29,929 mutaciones con un total de 357 muestras únicas de donadores. 

Posterior a cargar el archivo, se removieron los hipermutantes, definidos por tener más de 500 mutaciones en una misma muestra, esto debido a que son factores de confusión para driver detection y pueden causar falsos positivos. 

En total se eliminaron dos muestras hipermutantes: `TCGA-4R-AA8I-01` y `TCGA-UB-A7MB-01` con 848 y 1169 mutaciones respectivamente 

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-16-14-50-39-image.png)

Figura 1: Numero de mutaciones en cada donante. La linea roja significa el nivel de significangia para definir un hipermutante (más de 500 mutaciones), en total se eliminaron 2 muestras

Posteriormente, se usó el programa `dndscv` para encontrar posibles gene divers dentro del hepatocarcinoma. Despues a correr el programa se detectaron 17 posibles gene drivers, esto debido a que tenían un p-valor ajustado global $<0.1$. 

Adicionalmente, no se sidntificó ningún gen bajo selección negativa. Se buscó información de los genes restantes además de consultar sus coeficientes de seleccion significativos a un nivel de significancia de 0.1, los resultados fueron estos:

### Well established oncogenes

Incluyen a *CTNNB1* y *NFE2L2* como oncogenes. Esto se refuerza teoricamente ya que tienen exlusivamente mutaciones missense signficiativas.

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-17-12-21-image.png)

Figura 2: Coeficiente de seleccion para los genes identificados como oncogenes para hepatocarcinoma. Se muestran los coeficientes de seleccion sigifnciativos, los cuales eran exlusivamente missense

Anotaciones:

- *CTNNB1*: Participa en la via de sañalizacion Wnt, en donde de forma downstream, se afctivan factores de transcripción, como el de la proliferacion celular con ayuda de *Myc*

- *NFE2L2* Participa en la via de señalización para la resistencia al estres oxidativo, regulado por *KEAP1*. *NFE2L2* libera señales de superivivencia cuando esta presente este tipo de estres

### Well-established tumor suppresors

Se identificaron 9 supresores de tumores, los cuales incluyen a: *TP53*, *AXIN1*, *ARID1A*, *KEAP1*, *TSC2*, *BAP1*, *RB1*, *ARID2* y *ACVR2A*. Al momento de revisar los coeficientes de selección siginificativos para estos genes, se encontró que tenían mayoritariamente mutaciónes de tipo Indeles, nonsense  y splicing, sin embargo, se identificaron dos coeficientes de seleccion significativos para las mutaciones missense.

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-17-21-56-image.png)

Figura 3: Coeficientes de seleccion signifcativos para supresores de tumores.

Anotaciones 

- *TP53*: El famoso guardian del genoma, cuando hay mucho daño celular mata a las celulas

- *AXIN1*: Degrada a la $\beta$catenina, apaga la via de señalizacion WnT, evitando que se transcriban señales de division

- *ARID1A*: Es un remodelador de la cromatina moviendo nuclesomas, perder este gen provoca una desregulacion de esto, fomentando la transcripcion

- *KEAP1*: Controla al gen *NFE2L2*, oncogen explicado previamente 

- *TSC2*: Apaga la via PI3K/AKT/mTOr, la via que le dice a las celulas vive 

- *BAP1*: Desubiquitina proteinas relacionadas con la reparacion del ADN, creguladoras del ciclo celular y reguacion trascricional, intercatua con BRCA1

- *RB1*: Primer supresor tumoral descubierto, estabiliza la cromatina y es un regulador negativo de la transcripcion

- *ACVR2A*: Participa en la regulacion del ciclo celular, su inactivacion lo desarrgula y aumenta la glucolisis

- *ARID2*: Tiene una funcion similar que ARID1A

- *IRX1*: Factor de trascricion involucrado en el desarrollo del sistema neviosos, puede ser supresor de tumor o un oncogen

### Putative oncogenic role status

Los genes restantes, corresponden a genes usados como biomarcadores, o con un rol no oncogenico ni supresor:

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-17-54-58-image.png)

- *ALB*: Albumina, regula la presion osmotica, transporta hormonas, lipidos, etc. Actualmente se considera un marcador de hepatocarcinoma, sin embargo estudios recientes lo asocian con agresividad

- *NLGN1*: Se observo que puede promover la invasion a traves de nervios en celulas canceringeas, gen que regula la angiogenesis. 

- *CRIP3*: Tiene una funcion dual, puede ser un supresor de tumores o un oncogen, se una como un biomarcador en algunos tipos de cancer

- *KRTAP3-3*: Se usa como un biomarcador para prognosticar la carcinogenesis (ayuda al trasporte celular y a la keratina del pelo) 

- *KRTAP22-1*: Gen de la queratina
  
  ### Filtrado de gene drivers
  
  Posterior a la anotación, se descartaron *IRX1* y *KRTAP22-1* esto debido a que tenian una notacion conflictiva o no asociada con el cancer, no tener ningun coeficiente signficiativo y debido a que sus q values globales estaban cercanos a  0.1
  
  <img src="file:///home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-17-54-22-image.png" title="" alt="" width="415">

Tabla 3. Lon coeficientes de seleccion estadisticamente signifcitavios se subrrayan con gris. Los oncogenes tienen exclusivamente coeficientes para mutaciones missense (verde), mientras que los todos supresores tumorales tiene coeficientes para truncar proteinas adicional a otros. Los genes en azul corresponden a genes con funcion oncogenica ni supresora de tumores

## Numero de mutaciones seleccionadas

Posterior a realizar el filtrado de los genes, se estimo el numero de mutaciones seleccionadas para cada gen donde ocurrio un coeficiente de seleccion significativo. Definimos una mutacion seleccionada como el numero de mutaciones que dan una ventaja genuina dentro de este gen. Para estimar el numero de mutaciones seleccionadas definimos:

$$
Muts_{\text{selected}} = Sel_\text{prop}*N_{obs}\\ \text{donde:}\\
\\
Sel_{\text{prop}}=\frac{w-1}{w}
$$

Esto se estimo para todos los genes que tenian coeficientes signifcicativos

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-19-23-57-image.png)

Opcion 2

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-19-54-11-image.png)

Es una grafica similar solo que se subrallan los genes que estan en la misma via: En verde es el complejo SWI/SNF de remodelacion de la cromatina.Morado es el pathway Wnt / $\beta$catetina de la transcripcion. En rosa el complejo KEAP1-NRF2 relacionado al estres oxidativo y en gris el ciclo celular

CTNNB1 y AXIN1: Son antagonistas, perder AXIN1 es lo mismo que ganar CTNNB1

KEAP y NF2EL2 igual son antagonistas

Perder TP3 es lo mismo que perder RB1





### Estimacion de hotspots

Se buscaron hostspots para los genes. Los mas importantes fueron:

#### CTNNB1

Los conteos fueron:

- proteina: S45P, candena: T133C en la posicion: 41266136 de tipo Missense: 8

- proteina:D32G, candena: A95G en la posicion: 41266098 de tipo Missense: 4

- proteina:H36P, candena: A107C en la posicion: 41266110 de tipo Missense: 4

- proteina: T41A, candena: A121G en la posicion: 41266124 de tipo Missense: 4

- proteina: K335I, cadena: A1004T en la posicion: 41268766 de tipo Missense: 3

- proteina: S33P, cadena: T97C en la posicion: 41266100 de tipo Missense: 3

- proteina: S45Y cadena: C134A en la posicion: 41266137 de tipo Missense: 3

Todas son cambios missense en la posicion para la cadena 41266XXX, rondando en los cambios como minimo de 32 a 45, esto corresponde al exon 3, el cual abarca delsde el aminoacido 26 al 80. Este eston es importante porque activa sitios de fosforiacion para la degradacion de la catetina. Mutaciones en este pueden hacer uq ete catenina no se degrade y sea una GOF





![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-20-19-59-image.png)



<img src="file:///home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-20-13-18-image.png" title="" alt="" width="530">



![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-20-20-14-image.png)



Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC5797067/#s3

#### TP53

- Proteina: R249S, cadena: G747T en la posicion: 7577534 de tipo Missense: 5

- Proteina: H193R, cadena: A578G en la posicion: 7578271 de tipo Missense: 4

R249S es particularmente interesante — es una **firma mutacional de exposición a Aflatoxina B1**, un carcinógeno hepático muy conocido.



<img src="file:///home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-20-39-46-image.png" title="" alt="" width="211"><img title="" src="file:///home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-20-40-13-image.png" alt="" width="160"><img title="" src="file:///home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-20-40-44-image.png" alt="" width="124">

R249S: Igual afecta que se una al DNA, es una de las mutaciones somaticas mas estudiadas chance por hepatitis B o por hongos D: [Aflatoxin levels and prevalence of TP53 aflatoxin-mutations in hepatocellular carcinomas in Mexico](https://ru.dgb.unam.mx/items/46f9c0b8-ab8c-4709-b276-8d6dd310120a) de la unam :D

H193R: Afecta al dominino de union al DNA: [OncoKB™ - MSK's Precision Oncology Knowledge Base](https://www.oncokb.org/gene/TP53/somatic/H193R)

### Usando sitednds para detectar hotspots

Mientras que `dndscv` detecta hotspots a nivel del gen, `sitednds` usa posiciondes del genoma. Al momento de correrlo, se confirmaron los mismos resultados que con dndsv. Los hostsopots seleccionados son reales 

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-20-51-42-image.png)

### Ahora con codondnds para detectar hotspots tmb

Ay no manchen me dio lo mismo :(

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-00-05-image.png)

Solo se detecto S37 con este programa, tambien esta en el exon 3

### Oncoplots

Se hiceron con maftools

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-04-01-image.png)

La titina es putative passenger. El concoplot confirma a los dirvers:

- TP53: Se destruye de cualquier forma, patron LOF

- CTNNB1: Se busca missense para una GOF, poco de los otros 

- ARID1A: Paton LOF

- AXIN: Patron LOG

- Probablemente sean passenger TTN, MUC16, PCLO, OBSCN, RYR2

- Raramente TP53 y CTNNB1 

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-08-18-image.png)



![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-17-59-image.png)





![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-18-17-image.png)

 Oncoplot de los genes significativos por dndsv

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-24-34-image.png)

CTNNB1 y AXIN1: Eran antagonistas y casi nunca coocurren. Un GOF de CTNNB1 es lo mismo que un LOF de AXIN1

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-30-24-image.png)

KEAP y NF2EL2 igual son antagonistas y nunca coocurren, un LOF de KEAP1 es lo mismo que un GOF de NFEL2

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-29-59-image.png)



Perder TP3 es lo mismo que perder RB1

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-28-46-image.png)

pasa algo similar con TP53 y CTNNB1

    ![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-28-12-image.png)

Oncoplots con notas

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-32-37-image.png)

Con pathways en lugar de genes

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-35-33-image.png)

Se afectan principalemtne dos vias, la de la intefgriad del genoma y la de wnt signaling

#### Summary

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-37-26-image.png)

Dominado por mutaciones missense por SNPs, No manchem aqui aparecen las firmas, en las boxplots dominan los missense, el gen mas mutado es la titna pero e spaasenger

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-40-11-image.png)

dominan C>T ~30% y C>A ~22% casi son iguales Ti Tv



#### Lolipop de TP53 y sus hotspots

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-45-31-image.png)

#### CTNNB1

![](/home/jcontreras/.var/app/com.github.marktext.marktext/config/marktext/images/2026-05-17-21-47-57-image.png)

La mayoria esta en una missense del N teerminal (donde se degradaban)

## Workflow

1. Remoción de hipermutantes (> 500 mutaciones)

2. Se corrió `dndscv` para identificar gene drivers

3. Identificacion de gene drivers (global q-value < 0.01)

4. Anotacion de los posibles gene drivers

5. Analisis de coeficientes de seleccion estadisticamente significativos para cada subtipo de gene drivers

6. Descarte de no drivers raros *IRX1* y *KRTAP22-1*

7. Estimacion de numero de mutaciones seleccionadas a partir de los coeficicentes de seleccion significativos

8. Estimacion de hotspots con dndsv, sitednds y codondnds

9. Creacion de oncoplots

10. Creacion de lolipoop plots para genes con hotspot





## Referencias:

1. Dantzer, C., Dif, L., Vaché, J. *et al.* Specific features of ß-catenin-mutated hepatocellular carcinomas. *Br J Cancer* **131**, 1871–1880 (2024). https://doi.org/10.1038/s41416-024-02849-7

2. Amaddeo G, Guichard C, Imbeaud S, Zucman-Rossi J. Next-generation sequencing identified new oncogenes and tumor suppressor genes in human hepatic tumors. Oncoimmunology. 2012 Dec 1;1(9):1612-1613. doi: 10.4161/onci.21480. PMID: 23264911; PMCID: PMC3525620.

3. Lauren V. Albrecht, Nydia Tejeda-Muñoz, Maggie H. Bui, Andrew C. Cicchetto, Daniele Di Biagio, Gabriele Colozza, Ernst Schmid, Stefano Piccolo, Heather R. Christofk, Edward M. De Robertis, GSK3 Inhibits Macropinocytosis and Lysosomal Activity through the Wnt Destruction Complex Machinery,Cell Reports,Volume 32, 
   https://doi.org/10.1016/j.celrep.2020.107973.

4. Wu RC, Wang TL, Shih IeM. The emerging roles of ARID1A in tumor suppression. Cancer Biol Ther. 2014 Jun 1;15(6):655-64. doi: 10.4161/cbt.28411. Epub 2014 Mar 11. PMID: 24618703; PMCID: PMC4049780.

5. Krishna, A., Meynert, A., Dolt, K.S. *et al.* Mutational scanning reveals oncogenic *CTNNB1* mutations have diverse effects on signaling. *Nat Genet* **58**, 366–375 (2026). https://doi.org/10.1038/s41588-025-02496-5

6. Moon EJ, Giaccia A. Dual roles of NRF2 in tumor prevention and progression: possible implications in cancer treatment. Free Radic Biol Med. 2015 Feb;79:292-9. doi: 10.1016/j.freeradbiomed.2014.11.009. Epub 2014 Nov 29. PMID: 25458917; PMCID: PMC4339613.

7. Salahshor S, Woodgett JR. The links between axin and carcinogenesis. J Clin Pathol. 2005 Mar;58(3):225-36. doi: 10.1136/jcp.2003.009506. Erratum in: J Clin Pathol. 2005 Dec;58(12):1344. PMID: 15735151; PMCID: PMC1770611.

8. Cui G, Zhou Y, Liao W, Cossu E, Evert M, Zhang S, Wang J, Deng S, Yonemura A, David L, Xu M, Doumergue JM, Lu X, Yang L, Li J, Liang B, Wang H, Xu H, Zhong S, Deng Y, Calvisi DF, Zhao J, Chen X, Wang X. Inhibition of GSK3 and TSC2 Mediates the Oncogenic Activity of AKT in Hepatocellular Carcinoma. Cancer Res. 2025 Dec 15;85(24):5049-5065. doi: 10.1158/0008-5472.CAN-25-1615. PMID: 40939191; PMCID: PMC12500318.

9. Han A, Purwin TJ, Aplin AE. Roles of the BAP1 Tumor Suppressor in Cell Metabolism. Cancer Res. 2021 Jun 1;81(11):2807-2814. doi: 10.1158/0008-5472.CAN-20-3430. Epub 2021 Jan 14. PMID: 33446574; PMCID: PMC8178170.

10. National Center for Biotechnology Information. (2026). RB1 RB transcriptional corepressor 1, Homo sapiens (human). NCBI Gene. https://www.ncbi.nlm.nih.gov/gene/5925

11. Yasukawa K, Shimada S, Akiyama Y, Taniai T, Igarashi Y, Tsukihara S, Tanji Y, Umemura K, Kamachi A, Nara A, Yamane M, Akahoshi K, Shimizu A, Soejima Y, Tanabe M, Tanaka S. ACVR2A attenuation impacts lactate production and hyperglycolytic conditions attracting regulatory T cells in hepatocellular carcinoma. Cell Rep Med. 2025 Apr 15;6(4):102038. doi: 10.1016/j.xcrm.2025.102038. Epub 2025 Mar 25. PMID: 40139191; PMCID: PMC12047472.

12. Wodziński D, Wosiak A, Pietrzak J, Świechowski R, Jeleń A, Balcerczak E. Does the expression of the ACVR2A gene affect the development of colorectal cancer? Genet Mol Biol. 2019 Jan-Mar;42(1):32-39. doi: 10.1590/1678-4685-GMB-2017-0332. Epub 2019 Mar 11. PMID: 30856244; PMCID: PMC6428132.

13. Shukla, P., Dange, P., Mohanty, B.S. *et al.* ARID2 suppression promotes tumor progression and upregulates cytokeratin 8, 18 and β-4 integrin expression in *TP53*-mutated tobacco-related oral cancer and has prognostic implications. *Cancer Gene Ther* **29**, 1908–1917 (2022). https://doi.org/10.1038/s41417-022-00505-x

14. Fu X, Yang Y, Zhang D. Molecular mechanism of albumin in suppressing invasion and metastasis of hepatocellular carcinoma. Liver Int. 2022 Mar;42(3):696-709. doi: 10.1111/liv.15115. Epub 2021 Dec 7. PMID: 34854209; PMCID: PMC9299813.

15. Bizzozero L, Pergolizzi M, Pascal D, Maldi E, Villari G, Erriquez J, Volante M, Serini G, Marchiò C, Bussolino F, Arese M. Tumoral Neuroligin 1 Promotes Cancer-Nerve Interactions and Synergizes with the Glial Cell Line-Derived Neurotrophic Factor. Cells. 2022 Jan 14;11(2):280. doi: 10.3390/cells11020280. PMID: 35053395; PMCID: PMC8774081.

16. Judaki AA, Shirinpoor M, Farahani M, Aldaghi T, Arefi-Oskouie A, Nazari E. Bioinformatics and machine learning reveal novel prognostic biomarkers in head and neck squamous cell carcinoma. J Appl Genet. 2025 Oct 1. doi: 10.1007/s13353-025-01018-7. Epub ahead of print. PMID: 41028529.
