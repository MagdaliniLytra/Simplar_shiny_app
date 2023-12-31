# Simplar_shiny_app

The central focus of this project entailed the conceptualization
and development of an application geared towards easing the
workflow for biologists possessing a scientific foundation in
RNA-Seq analyses but exhibiting limited or no proficiency in
programming.

Simplar is a flexible user-friendly tool for the analysis and visualization
of gene expression data. Within the graphical user interface (GUI), users are enabled
to upload pertinent files, subsequently triggering a series of analytical procedures
tailored for the examination of data derived from RNA-Seq
experiments. These sequential analytical steps encompass the
Trimming of raw sequences originating from RNA-Seq
experiments, Alignment to the reference genome,
Quantification of sequences with subsequent conversion to
gene expression levels, Differential Expression Analysis (DGEA), and
Gene Set Enrichment Analysis (GSEA).

The outcomes of these meticulous analyses are formatted to
meet publication standards and are conveniently accessible to
users through the graphical user interface (GUI). Additionally,
users are afforded the option to download these results in
diverse file formats, including .csv, .pdf, .png, among others.
It is imperative to note that all the aforementioned processes
and analyses are meticulously executed within the R
programming language, utilizing established packages widely
recognized and endorsed by the research community.



Running the app

To be able to use the app, users should have pre-installed R, R studio and R tools to their computers.
Users have to download from github the "WWW" folder, the "main_script.r", the "global2v.r", 
the "Shiny_App_Manual.RMD" and save them in the same directory.
The users have to open the R script “main_script.r” within R studio and press the “Run app” button.
The app will open either on a new R preview windows or on a new browser window, depending on users’ settings.


Packages

Trimming: Rbowtie2(AdapterRemoval)

Mapping: Rsubread

Quantification: Rsubread(Feature Counts)

Differential Expression: DESeq2

Gene Set Enrichment: AnnotationHub, clusterProfiler
