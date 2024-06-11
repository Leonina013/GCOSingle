# Sys.setenv("plotly_username"=" your_plotly_username")
# Sys.setenv("plotly_api_key"="your_api_key")
## test repo

# wd <- dirname(rstudioapi::getActiveDocumentContext()$path) # set wd as the current folder
# print(wd == getwd())
# print(wd)
# print(getwd())
# if (!wd == getwd()) {
#   setwd(wd)
# }

print("start loading")
start.load <- Sys.time() ### time

if (length(find.package(package = "shiny", quiet = T)) > 0) {
  library(shiny)
} else {
  print("Package shiny not installed")
  install.packages("shiny")
  print("Package shiny installed")
  library(shiny)
}

if (length(find.package(package = "cyjShiny", quiet = T)) > 0) {
  library(cyjShiny)
} else {
  remotes::install_github("cytoscape/cyjShiny")
  library(cyjShiny)
}

if (length(find.package(package = "shinythemes", quiet = T)) > 0) {
  library(shinythemes)
} else {
  print("Package shinythemes not installed")
  install.packages("shinythemes")
  print("Package shinythemes installed")
  library(shinythemes)
}

if (length(find.package(package = "rstudioapi", quiet = T)) > 0) {
  library(rstudioapi)
} else {
  install.packages("rstudioapi")
  library(rstudioapi)
}

if (!length(find.package(package = "rlang", quiet = T)) > 0) {
  install.packages("rlang")
}

#################################
if (length(find.package(package = "RColorBrewer", quiet = T)) > 0) {
  library(RColorBrewer)
} else {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if (length(find.package(package = "kohonen", quiet = T)) > 0) {
  library(kohonen)
} else {
  install.packages("kohonen")
  library(kohonen)
}

if (length(find.package(package = "sparsepca", quiet = T)) > 0) {
  library(sparsepca)
} else {
  install.packages("sparsepca")
  library(sparsepca)
}

if (length(find.package(package = "randomForest", quiet = T)) > 0) {
  library(randomForest)
} else {
  install.packages("randomForest")
  library(randomForest)
}

if (length(find.package(package = "cluster", quiet = T)) > 0) {
  library(cluster)
} else {
  install.packages("cluster")
  library(cluster)
}

if (length(find.package(package = "ComplexHeatmap", quiet = T)) > 0) {
  library(ComplexHeatmap)
} else {
  install.packages("ComplexHeatmap")
  library(ComplexHeatmap)
}

## for t-sne
# if (length(find.package(package = "reticulate", quiet = T)) > 0) {
#   library(reticulate)
# } else {
#   install.packages("reticulate")
#   library(reticulate)
# }

if (length(find.package(package = "Rtsne", quiet = T)) > 0) {
  library(Rtsne)
} else {
  install.packages("Rtsne")
  library(Rtsne)
}

######################################For Report Generation#################
if (length(find.package(package = "webshot2", quiet = T)) > 0) {
  library(webshot2)
} else {
  remotes::install_github("rstudio/webshot2")
  library(webshot2)
}
library(gridExtra)
library(plotly)
if (!require("processx")) install.packages("processx")
library(png)

if (length(find.package(package = "capture", quiet = T)) > 0) {
  library(capture)
} else {
  remotes::install_github("dreamRs/capture")
  library(capture)
}

####################### Dependencies For RAFSIL ###################################
# if (length(find.package(package = "RAFSIL", quiet = T)) > 0) {
#   library(RAFSIL)
# }

if (length(find.package(package = "gridGraphics", quiet = T)) > 0) {
  library(gridGraphics)
} else {
  install.packages("gridGraphics")
  library(gridGraphics)
}

if (length(find.package(package = "gridExtra", quiet = T)) > 0) {
  library(gridExtra)
} else {
  install.packages("gridExtra")
  library(gridExtra)
}

if (length(find.package(package = "tidyverse", quiet = T)) > 0) {
  library(tidyverse)
} else {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

if (length(find.package(package = "ggpubr", quiet = T)) > 0) {
  library(ggpubr)
} else {
  install.packages("ggpubr")
  library(ggpubr)
}

if (length(find.package(package = "networkD3", quiet = T)) > 0) {
  library(networkD3)
} else {
  install.packages("https://cran.r-project.org/src/contrib/Archive/networkD3/networkD3_0.2.10.tar.gz", repo=NULL, type="source")
  library(networkD3)
}

if (length(find.package(package = "data.tree", quiet = T)) > 0) {
  library(data.tree)
} else {
  install.packages("data.tree")
  library(data.tree)
}


if (length(find.package(package = "bubbles", quiet = T)) > 0) {
  library(bubbles)
} else {
  devtools::install_github("jcheng5/bubbles", upgrade = FALSE)
  library(bubbles)
}

###################################################################################
####################### Dependencies For Uniprot ###################################

if (length(find.package(package = "UniprotR", quiet = T)) > 0) {
  library(UniprotR)
} else {
  install.packages("UniprotR")
  library(UniprotR)
}

if (length(find.package(package = "scales", quiet = T)) > 0) {
  library(scales)
} else {
  install.packages("scales")
  library(scales)
}

###################################################################################

####################### Dependencies For Pathway Enrichment ###################################

if (length(find.package(package = "gprofiler2", quiet = T)) > 0) {
  library(gprofiler2)
} else {
  install.packages("gprofiler2")
  library(gprofiler2)
}

###################################################################################

####################### Dependencies For Protein Interactions ###################################

if (length(find.package(package = "httr", quiet = T)) > 0) {
  library(httr)
} else {
  install.packages("httr")
  library(httr)
}

if (length(find.package(package = "curl", quiet = T)) > 0) {
  library(curl)
} else {
  install.packages("curl")
  library(curl)
}

if (length(find.package(package = "later", quiet = T)) > 0) {
  library(later)
} else {
  install.packages("later")
  library(later)
}

if (length(find.package(package = "qdapTools", quiet = T)) > 0) {
  library(qdapTools)
} else {
  install.packages("qdapTools")
  library(qdapTools)
}

if (length(find.package(package = "alakazam", quiet = T)) > 0) {
  library(alakazam)
} else {
  install.packages("https://cran.r-project.org/src/contrib/Archive/alakazam/alakazam_1.0.0.tar.gz", repo=NULL, type="source")
  library(alakazam)
}

if (length(find.package(package = "msa", quiet = T)) > 0) {
  library(msa)
} else {
  BiocManager::install("msa", update = FALSE)
  library(msa)
}

if (length(find.package(package = "ape", quiet = T)) > 0) {
  library(ape)
} else {
  install.packages("ape")
  library(ape)
}

if (length(find.package(package = "seqinr", quiet = T)) > 0) {
  library(seqinr)
} else {
  install.packages("seqinr")
  library(seqinr)
}

if (length(find.package(package = "qdapRegex", quiet = T)) > 0) {
  library(qdapRegex)
} else {
  install.packages("qdapRegex")
  library(qdapRegex)
}

###################################################################################

####################### Dependencies For Co-expression ###################################

if (length(find.package(package = "shinyjs", quiet = T)) > 0) {
  library(shinyjs)
} else {
  install.packages("shinyjs")
  library(shinyjs)
}

###################################################################################

####################### Dependencies For Microarray ###################################

if (length(find.package(package = "devtools", quiet = T)) > 0) {
  library(devtools)
} else {
  install.packages("devtools")
  library(devtools)
}

if (length(find.package(package = "remotes", quiet = T)) > 0) {
  library(remotes)
} else {
  devtools::install_github("r-lib/remotes")
  library(remotes)
}

if (length(find.package(package = "maEndToEnd", quiet = T)) > 0) {
  suppressPackageStartupMessages({library("maEndToEnd")})
} 
# else {
#   remotes::install_github("b-klaus/maEndToEnd", ref="master")
#   suppressPackageStartupMessages({library("maEndToEnd")})
# }

if (length(find.package(package = "oligoClasses", quiet = T)) > 0) {
  library(oligoClasses)
} else {
  print("Package oligoClasses not installed")
  BiocManager::install("oligoClasses")
  print("Package oligoClasses installed")
  library(oligoClasses)
}

if (length(find.package(package = "ArrayExpress", quiet = T)) > 0) {
  library(ArrayExpress)
} else {
  print("Package ArrayExpress not installed")
  BiocManager::install("ArrayExpress")
  print("Package ArrayExpress installed")
  library(ArrayExpress)
}

if (length(find.package(package = "pd.hugene.1.0.st.v1", quiet = T)) > 0) {
  library(pd.hugene.1.0.st.v1)
} else {
  print("Package pd.hugene.1.0.st.v1 not installed")
  BiocManager::install("pd.hugene.1.0.st.v1")
  print("Package pd.hugene.1.0.st.v1 installed")
  library(pd.hugene.1.0.st.v1)
}

if (length(find.package(package = "hugene10sttranscriptcluster.db", quiet = T)) > 0) {
  library(hugene10sttranscriptcluster.db)
} else {
  print("Package hugene10sttranscriptcluster.db not installed")
  BiocManager::install("hugene10sttranscriptcluster.db")
  print("Package hugene10sttranscriptcluster.db installed")
  library(hugene10sttranscriptcluster.db)
}

if (length(find.package(package = "arrayQualityMetrics", quiet = T)) > 0) {
  library(arrayQualityMetrics)
} else {
  print("Package arrayQualityMetrics not installed")
  BiocManager::install("arrayQualityMetrics")
  print("Package arrayQualityMetrics installed")
  library(arrayQualityMetrics)
}

if (length(find.package(package = "limma", quiet = T)) > 0) {
  library(limma)
} else {
  print("Package limma not installed")
  BiocManager::install("limma")
  print("Package limma installed")
  library(limma)
}

if (length(find.package(package = "topGO", quiet = T)) > 0) {
  library(topGO)
} else {
  print("Package topGO not installed")
  BiocManager::install("topGO")
  print("Package topGO installed")
  library(topGO)
}

if (length(find.package(package = "ReactomePA", quiet = T)) > 0) {
  library(ReactomePA)
} else {
  print("Package ReactomePA not installed")
  BiocManager::install("ReactomePA")
  print("Package ReactomePA installed")
  library(ReactomePA)
}

###################################For GEO import###############################
library(xml2)
if (length(find.package(package = "GEOquery", quiet = T)) > 0) {
  library(GEOquery)
} else {
  BiocManager::install("GEOquery")
  library(GEOquery)
}
if (length(find.package(package = "umap", quiet = T)) > 0) {
  library(umap)
} else {
  install.packages("umap")
  library(umap)
}
if (length(find.package(package = "maptools", quiet = T)) > 0) {
  library(maptools)
} else {
  install.packages("maptools", repos="http://R-Forge.R-project.org")
  library(maptools)
}

#### Other libaries ####
library(BiocParallel)
library(SeuratObject)
library(Seurat)
library(patchwork)
library(hdf5r)
library(SingleR)
library(celldex)
library(pheatmap)
library(harmony)
library(monocle3)
library(SeuratWrappers)
library(singleCellHaystack)
library(markdown)

###################################################################################
########################### Style files for Cytoscape.js ################

styles <- c(
  "generic style"="./www/style/basicStyle.js",
  "red-yellow"="./www/style/red-yellow.js",
  "red-pink" = "./www/style/red-pink.js",
  "green-blue"="./www/style/green-blue.js",
  "green-blue(ppi)"="./www/style/green-blue(ppi).js")

##########################################################################

#
# ## sourcing util files
source("./www/utils.R")
source("./www/PhyscochemicalSep.R")

#
loadPkg()

id_to_name <- read.csv(paste0("./www/TransTable_Human.csv"))


#################### Complex Enrichment ##########################

complexes <- load(paste0("./www/allComplexes.RData")) #allComplexes is masked under complexes
up_corum_mapping <- read.csv(paste0("./www/UniProt_CORUM_Mapping.csv"))

##################################################################


end.load <- Sys.time()
print("loading time")
print(end.load - start.load)

##### UI from here ###########
ui <- tagList(
  shinyjs::useShinyjs(),
  navbarPage(
    id = "navbar",
    theme = shinytheme("flatly"),
    title = "",
    tabPanel(
      "GeneCloudOmics",
      br(),
      sidebarLayout(
        sidebarPanel(
          img(
            src = "GeneCloudOmics-logo.png",
            width = "100%", height = "100%"
          )
        ),
        mainPanel(
          h2("Welcome to GeneCloudOmics", align = "center",style = "color:#73C6B6;font-weight: bold;"),
          p("The Biostatistical Tool for Gene Expression Data Analysis", align = "center"),
          h4(span("GeneCloudOmics", style = "color:#73C6B6;font-weight: bold;")," an easy-to-use web server for high-throughput gene expression analysis a comprehensive range of data analytics tools in one package that no other current standalone software or web-based tool can do currently"),
          h4(span("GeneCloudOmics", style = "color:#73C6B6;font-weight: bold;"),"  provides the user access to 23 different data analytical and bioinformatics tasks including reads normalization, scatter plots, linear/non-linear correlations, PCA, clustering (hierarchical, k-means, t-SNE, SOM), differential expression analyses, pathway enrichments, evolutionary analyses, pathological analyses, and protein-protein interaction (PPI) identifications."),
          h4(span("Supported Transcriptome data:", style = "color:#73C6B6;font-weight: bold;"), " RNA-Seq and Microarray "),
          h4(span("Data Preprocessing:", style = "color:#73C6B6;font-weight: bold;"), " GeneCloudOmics performs raw data normalization using four normalization methods RPKM, 
      FPKM, TPM and RUV. The raw vs. normalized data are visualized as boxplots and violin plots."),
          h4(span("Differential Gene Expression (DGE) Analysis:", style = "color:#73C6B6;font-weight: bold;"), " GDE using five methods EdgeR, DESeq2 and NOISeq."),
          h4(span("Bio-statistical Analysis:", style = "color:#73C6B6;font-weight: bold;"), " GeneCloudOmics provides the user with the following bio-statistical analyses: 
      Pearson and Spearman rank correlations, PCA, k-means and hierarchical clustering, 
      Shannon entropy and noise (square of the coefficient of variation), t-SNE, random forest and SOM analyses. All analyses include proper high-resolution visualization."),
          h4(span("Bioinformatics Analysis of Gene and Protein sets:", style = "color:#73C6B6;font-weight: bold;"), " For the differential expressed genes (DEG), GeneCloudOmics provides 
      the users with multiple bioinformatics tools to investigate their 
      gene/protein list including gene ontology (GO), pathway enrichment analysis, PPI, co-expression, gene/protein function, subcellular localization, complex enrichment, protein domains, tissue expression, sequence properties (acidity, hydrophobicity and charge),
       evolutionary analysis (gene tree, phylogenetic tree and species/chromosome location)  and pathological analysis (diseases that these genes/proteins are involved in). The analyses include proper high-resolution visualization, when applicable."),
          h4( #span("For questions and bug reporting, please write to", a("Mohamed Helmy", href="mailto:mohamed_helmy@bii.a-star.edu.sg", target="_blank"), "or ", a("Kumar Selvarajoo", href="mailto:kumar_selvarajoo@bii.a-star.edu.sg", target="_blank")),
            withTags({
              div(class="header", checked=NA,
                  h4("User manual can be found ", 
                     a("here", href= "https://github.com/buithuytien/GeneCloudOmics/blob/master/GeneCloudOmics_Help_1.0.pdf", target="_blank")), # "GeneCloudOmics_Help_1.0.pdf"
                  h4(span("For questions and bug reporting, please write to", 
                          a("Mohamed Helmy", href="mailto:mohamed_helmy@bii.a-star.edu.sg", target="_blank"), "or ", 
                          a("Kumar Selvarajoo", href="mailto:kumar_selvarajoo@bii.a-star.edu.sg", target="_blank")),
                  )
                  
              )
            })
          )
        )
      )
    ),
    navbarMenu('Preprocessing',
               tabPanel(
                 "RNASeq Data",
                 value = "active_tab_rnaseq",
                 tabsetPanel(
                   id = "Rnaseq_pre",
                   tabPanel(
                     "Upload data",
                     sidebarPanel(
                       radioButtons(
                         "file_type", "Choose File Type",
                         c("Raw file (read count)" = "raw", "Normalised file" = "norm")
                       ),
                       conditionalPanel(
                         condition = "input.file_type=='raw'", # raw
                         withTags({
                           div(class="header", checked=NA,
                               p("Example ", 
                                 a(target="_blank", href="https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Yeast%20Biofilm%202%20-%205%20genotypes/Yeast-biofilm2-raw.csv", "csv"),
                                 a(target="_blank", href="https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_raw.png", "image"))
                           )
                         }),
                         fileInput("file1", "Choose Raw Counts (required)"),
                         
                         withTags({
                           div(class="header", checked=NA,
                               p("Example ", 
                                 a("csv", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Yeast%20Biofilm%202%20-%205%20genotypes/Yeast-biofilm2-length.csv"),
                                 a("image", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_gene_length.png")), # ADD EXAMPLE
                           )
                         }),
                         fileInput("length1", "Choose Gene Length (optional)"), # gene id + length
                         
                         withTags({
                           div(class="header", checked=NA,
                               p("Example ", 
                                 a("csv", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Zebra%20fish%20microarray/zfGenes_neg_control.csv", target="_blank"),
                                 a("image", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_negative_control_genes.png", target="_blank")), # ADD EXAMPLE
                           )
                         }),
                         fileInput("spikes1", "Choose negative controls (eg. ERCC Spike-in) (optional)")
                       ),
                       conditionalPanel(
                         condition = "input.file_type=='norm'", # normalized
                         withTags({
                           div(class = "header",
                               p("Example ", 
                                 a("csv", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Yeast%20Biofilm%202%20-%205%20genotypes/Yeast-biofilm2-normalized-tpm.csv", target="_blank"),
                                 a("image", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_normalised.png", target="_blank")), # ADD EXAMPLE
                           )
                         }),
                         fileInput("file2", "Choose Normalized Expression (required)")
                         # helpText("* Format requirement: CSV file. Gene names in rows and genotypes in columns, following the usual format of files deposited in the GEO database.")
                       ),
                       
                       withTags({
                         div(class = "header",
                             p("Example ", 
                               a("csv",   href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Yeast%20Biofilm%202%20-%205%20genotypes/Yeast-Biofilm2-meta.csv", target="_blank"),
                               a("image", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_metadata.png", target="_blank")), # ADD EXAMPLE
                         )
                       }),
                       fileInput("metafile1", "Choose Meta Data File (required)"),
                       actionButton("submit_input", "Submit")
                     ),
                     
                     
                     
                     mainPanel(
                       # h3("Welcome to GeneCloudOmics --"),
                       # h3("A Biostatistical tool for Transcriptomics Analysis"),
                       # img(
                       #   src = "GeneCloudOmics-logo.png",
                       #   width = 570, height = 370
                       # )
                     )
                   ),
                   tabPanel(
                     "Preprocessing",
                     value="active_tab_preprocessing_rnaseq",
                     sidebarPanel(
                       h4("Filtering"),
                       splitLayout(
                         numericInput("min_val", "Min. value", min = 0.1, step = 0.1, value = 1.0),
                         numericInput("min_col", "Min. columns", min = 1, value = 2)
                       ),
                       conditionalPanel(
                         condition = "input.file_type=='raw'",
                         radioButtons(
                           "norm_method", "Normalisation method",
                           c(
                             "None (Black)" = "None",
                             "RPKM (Gene length input required)" = "RPKM", 
                             "FPKM (Gene length input required)" = "FPKM",
                             "TPM (Gene length input required)" = "TPM",
                             "RUV (Negative control genes input required)" = "RUV"
                           )
                         )
                       ),
                       actionButton("submit_preprocessing", "Submit"),
                       conditionalPanel(
                         condition = "input.preprocessing_tabs == 'Data table' ",
                         br(),
                         br(),
                         downloadButton("download_norm_data", "Download table (csv)")
                       )
                     ),
                     mainPanel(
                       h3("Preprocessing Rnaseq data"),
                       tabsetPanel(
                         type = "tabs", id = "preprocessing_tabs",
                         tabPanel(
                           "RLE plot",
                           conditionalPanel(
                             condition = "$('html').hasClass('shiny-busy')",
                             div(img(src = "load.gif", width = 240, height = 180),
                                 h4("Processing ... Please wait"),
                                 style = "text-align: center;"
                             )
                           ),
                           conditionalPanel(
                             condition = "!$('html').hasClass('shiny-busy')",
                             plotOutput("RLE.plot2")
                           ),
                           
                           conditionalPanel(
                             condition = "input.file_type=='raw'",
                             conditionalPanel(
                               condition = "$('html').hasClass('shiny-busy')",
                               div(img(src = "load.gif", width = 240, height = 180),
                                   h4("Processing ... Please wait"),
                                   style = "text-align: center;"
                               )
                             ),
                             conditionalPanel(
                               condition = "!$('html').hasClass('shiny-busy')",
                               plotOutput("RLE.plot")
                             )
                           )
                         ),
                         tabPanel(
                           "Violin Plot",
                           conditionalPanel(
                             condition = "$('html').hasClass('shiny-busy')",
                             div(img(src = "load.gif", width = 240, height = 180),
                                 h4("Processing ... Please wait"),
                                 style = "text-align: center;"
                             )
                           ),
                           conditionalPanel(
                             condition = "!$('html').hasClass('shiny-busy')",
                             plotlyOutput("violin_plot2")
                           ),
                           conditionalPanel(
                             condition = "!$('html').hasClass('shiny-busy')",
                             plotlyOutput("violin_plot")
                           )
                         ),
                         tabPanel(
                           "Data table",
                           h3("Normalized data"),
                           DT::dataTableOutput("norm_table")
                         ),
                         tabPanel(
                           "Description table",
                           h3("Data description"),
                           DT::dataTableOutput("meta_table")
                         )
                       )
                     )
                   )
                 )
               ),
               tabPanel(
                 "Microarray Data",
                 value = "active_tab_micro",
                 sidebarPanel(
                   withTags({
                     div(class = "header",
                         p("Example data ", a("here", href = "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2967/", target="_blank")), # ADD EXAMPLE ( have to change )
                     )
                   }),
                   fileInput("file_micro", "Choose Microarray Data"),
                   downloadButton("downloadMicroRaw", "Download Raw Data as CSV"),
                   br(), br(),
                   downloadButton("downloadMicroMeta", "Download Meta Data as CSV")
                 ),
                 
                 mainPanel(
                   h3("Preprocessing Microarray Data"),
                   conditionalPanel(
                     condition = "$('html').hasClass('shiny-busy')",
                     div(img(src = "load.gif", width = 240, height = 180),
                         h4("Processing ... Please wait"),
                         style = "text-align: center;"
                     )
                   )
                 )
               ),
               tabPanel(
                 "GEO Data Import",
                 value = "active_tab_geo",
                 sidebarPanel(
                   tabsetPanel(type="tabs", id = "geo_tab",
                               tabPanel("GEO DATA",
                                        div(class= "header","For Example: GSE153941"),
                                        textInput("geo_acc_no", "Enter Accession Number", value = "", width = NULL, placeholder = NULL),
                                        radioButtons("file_type_button","FILE TYPE",
                                                     c("RnaSeq","Microarray","Auto")),
                                        actionButton("submit_geo_acc_no", "Submit")),
                               tabPanel("PREPROCESSING",value="geo_pre",radioButtons("file_name_button","SELECT FILE",
                                                                                     c("a")),
                                        
                                        withTags({
                                          div(class="header", checked=NA,
                                              p("Example ", 
                                                a("csv", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Yeast%20Biofilm%202%20-%205%20genotypes/Yeast-biofilm2-length.csv"),
                                                a("image", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_gene_length.png")), # ADD EXAMPLE
                                          )
                                        }),
                                        fileInput("length2", "Choose Gene Length (optional)"), # gene id + length
                                        
                                        withTags({
                                          div(class="header", checked=NA,
                                              p("Example ", 
                                                a("csv", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Zebra%20fish%20microarray/zfGenes_neg_control.csv", target="_blank"),
                                                a("image", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_negative_control_genes.png", target="_blank")), # ADD EXAMPLE
                                          )
                                        }),
                                        fileInput("spikes2", "Choose negative controls (eg. ERCC Spike-in) (optional)"),
                                        actionButton("submit_geo_preprocessing", "Submit"),
                               )
                               
                   )),
                 
                 
                 mainPanel(
                   h3("GEO DATA IMPORT"),
                   
                   tabsetPanel(type ="tabs", 
                               tabPanel(
                                 uiOutput("help_text_geo")
                               )),
                   h4("Overview of GEO Imported Data"),
                   tabsetPanel(type ="tabs",  tabPanel("Box Plot",plotOutput("geo_box_plot")),
                               tabPanel("Expression Density", plotOutput("geo_expr_plot")),
                               tabPanel("Mean variance", plotOutput("geo_mean_plot")),
                               tabPanel("UMAP", plotOutput("geo_umap_plot")))
                 )
               )
               
    ),
    navbarMenu('Transcriptome Analysis',
               tabPanel(
                 "    Scatter    ",
                 sidebarPanel(
                   selectInput(inputId = "scatter.x", label = "X-axis", choices = ""),
                   selectInput(inputId = "scatter.y", label = "Y-axis", choices = ""),
                   radioButtons(
                     "trans", "Transformation:",
                     c("None", "Natural log", "log2", "log10")
                   ),
                   checkboxInput("regline", "Display regression line", value = FALSE),
                   actionButton("submit_scatter", "Plot"),
                   br(),
                   br(),
                   downloadButton("downloadscatter", "Download as PDF"),
                   #h6("Download all pairs of samples in one PDF (this may take some time to run) :"),
                   br(),
                   br(),
                   actionButton("add_scatter","Add to report"),
                   downloadButton("downloadscatter_collage", "Download collage")
                 ),
                 mainPanel(
                   h3("Heatscatter"),
                   uiOutput("help_text_scatter"),
                   plotlyOutput("scatter.plot"),
                   #      tags$script('
                   # document.getElementById("downloadscatter").onclick = function() {
                   # var gd = document.getElementById("scatter.plot");
                   # Plotly.Snapshot.toImage(gd, {format: "png"}).once("success", function(url) {
                   #   var a = window.document.createElement("a");
                   #   a.href = url; 
                   #   a.type = "image/png";
                   #   a.download = "plot.png";
                   #   document.body.appendChild(a);
                   #   a.click();
                   #   document.body.removeChild(a);                      
                   # });
                   # }
                   # ')
                 )
               ),
               tabPanel(
                 "Distribution Fit",
                 sidebarPanel(
                   conditionalPanel(
                     condition = "input.dist_tabs=='Distribution Fit'",
                     selectInput(inputId = "dist.var", label = "Choose a column", choices = colnames("dataset")),
                     checkboxGroupInput("distributions", "Distributions:",
                                        choices = c("Log-normal", "Log-logistic", "Pareto", "Burr", "Weibull", "Gamma"), selected = c("Log-normal", "Pareto")
                     ),
                     radioButtons("dist_zoom", "Zoom to see fit", c("slider", "text input")),
                     conditionalPanel(
                       condition = "input.dist_zoom=='slider'",
                       sliderInput("dist_range", "Range:",
                                   min = 0.1, max = 1000, step = 1,
                                   value = c(0.1, 1000)
                       )
                     ),
                     actionButton("submit_distfit","Plot"),
                     conditionalPanel(
                       condition = "input.dist_zoom=='text input'",
                       textOutput("dist_range_allowed"),
                       numericInput("dist_range_min", "min", value = 0.1, min = 0.1, max = 1000),
                       numericInput("dist_range_max", "max", value = 1000, min = 0.1, max = 1000)
                     ),
                     downloadButton("downloaddist", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.dist_tabs=='AIC table'",
                     downloadButton("downloaddistaic", "Download as CSV")
                   ),
                   checkboxInput("checkbox_distfit", label = "Distribution Fit", value = FALSE),
                   checkboxInput("checkbox_aic", label = "AIC table", value = FALSE),
                   actionButton("add_distfit","Add to report")
                 ),
                 mainPanel(
                   h3("Distribution Fit"),
                   tabsetPanel(
                     type = "tabs", id = "dist_tabs",
                     tabPanel(
                       "Distribution Fit",
                       uiOutput("help_text_dis_fit"),
                       plotOutput("dist.plot")),
                     tabPanel(
                       "AIC table",
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         div(tableOutput("dist.aic"), style = "font-size:80%")
                       )
                     )
                   )
                 )
               ),
               tabPanel(
                 "  Correlation  ",
                 sidebarPanel(
                   radioButtons(
                     "cor_method", "Method:",
                     c("Pearson correlation", "Spearman correlation")
                   ),
                   actionButton("submit_corr","Plot"),
                   br(),
                   br(),
                   conditionalPanel(
                     condition = "input.cor_tabs == 'Correlation heatmap'",
                     downloadButton("downloadcorrplot", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.cor_tabs == 'Correlation plot'",
                     downloadButton("downloadcorrplot2", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.cor_tabs == 'Correlation matrix'",
                     downloadButton("downloadcorrmat", "Download as CSV")
                   ),
                   checkboxInput("checkbox_corrheatmap", label = "Correlation Heatmap", value = FALSE),
                   checkboxInput("checkbox_corrplot", label = "Correlation Plot", value = FALSE),
                   checkboxInput("checkbox_corrmat", label = "Correlation Matrix", value = FALSE),
                   actionButton("add_correlation","Add to report"),
                   
                 ),
                 mainPanel(
                   conditionalPanel(
                     condition = "input.cor_method=='Pearson correlation'",
                     h3("Pearson correlation")
                   ),
                   conditionalPanel(
                     condition = "input.cor_method=='Spearman correlation'",
                     h3("Spearman correlation")
                   ),
                   tabsetPanel(
                     type = "tabs", id = "cor_tabs",
                     tabPanel(
                       "Correlation heatmap",
                       uiOutput("help_text_correlation"),
                       plotOutput("corr.plot")),
                     tabPanel("Correlation plot", plotOutput("corr.plot2")),
                     tabPanel("Correlation matrix", div(tableOutput("corr.matrix"), style = "font-size:80%"))
                   )
                 )
               ),
               tabPanel(
                 "PCA",
                 sidebarPanel(
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-2D plot'",
                     selectInput(inputId = "pca.x", label = "X-axis", choices = ""),
                     selectInput(inputId = "pca.y", label = "Y-axis", choices = "")
                   ),
                   selectInput(inputId = "gene_size", label = "Gene sample size", choices = ""),
                   radioButtons(
                     "gene_order", "Gene sample order (wrt column 1)",
                     c("Descending (highest to lowest)" = "Descending", "Ascending (lowest to highest)" = "Ascending", "Random")
                   ),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-2D plot' || input.pca_tabs == 'PCA-3D plot'",
                     checkboxInput("pca_cluster", strong("Kmeans clustering on columns"), FALSE),
                     conditionalPanel(
                       condition = "input.pca_cluster == true",
                       sliderInput("pca_cluster_num", "Number of clusters:", value = 1, min = 1, max = 1, step = 1),
                       checkboxInput("pca_text", strong("Display sample name"), FALSE)
                     )
                   ),
                   ######################################
                   radioButtons(
                     "pca_type", "Type of PCA",
                     c("PCA" = "PCA", "Sparse PCA" = "SPCA")
                   ),
                   ######################################
                   conditionalPanel(
                     condition = "input.gene_order=='Random'",
                     helpText("* Click multiple times to resample"),
                     actionButton("pca_refresh", "Resample", style = "background-color: #337ab7;border-color:#337ab7"),
                     br(), br()
                   ),
                   actionButton("submit_pca","Plot"),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA variance'",
                     downloadButton("downloadpcavar", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-2D plot'",
                     downloadButton("downloadpca2d", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-3D plot'",
                     downloadButton("downloadpca3d", "Download as PDF")
                   ),
                   checkboxInput("checkbox_pcavar", label = "PCA Variance", value = FALSE),
                   checkboxInput("checkbox_pca2d", label = "PCA-2D Plot", value = FALSE),
                   checkboxInput("checkbox_pca3d", label = "PCA3D Plot", value = FALSE),
                   actionButton("add_pca","Add to report")
                 ),
                 mainPanel(
                   h3("PCA"),
                   tabsetPanel(
                     type = "tabs", id = "pca_tabs",
                     tabPanel("PCA variance", 
                              uiOutput("help_text_PCA"),
                              plotlyOutput("pcavar.plot")),
                     tabPanel("PCA-2D plot", plotlyOutput("pca2d.plot")),
                     tabPanel("PCA-3D plot", plotlyOutput("pca3d.plot"))
                   )
                 )
               ),
               tabPanel(
                 "DE Analysis",
                 sidebarPanel(
                   radioButtons("n_rep", "Replicates?", choices = c("Multiple" = 1, "Single" = 0)),
                   conditionalPanel(
                     condition = "input.n_rep=='1'",
                     radioButtons("de_method1", "DE Method", choices = c("EdgeR", "DESeq2", "NOISeq"))
                   ),
                   conditionalPanel(
                     condition = "input.n_rep=='0'",
                     radioButtons("de_method0", "DE Method", choices = c("NOISeq"))
                   ),
                   h5("Choose 2 experiment conditions for DE analysis"),
                   selectInput("f1", "Condition 1", choices = ""),
                   selectInput("f2", "Condition 2", choices = ""),
                   
                   h5("DE criteria"),
                   splitLayout(
                     numericInput("p_val", "FDR", min = 0.01, max = 1, value = 0.05, step = 0.01),
                     numericInput("fc", "Fold Change", min = 1, value = 2, step = 0.1)
                   ),
                   fluidRow(
                     column(
                       4,
                       actionButton("submit_DE", "Plot")
                     ),
                     column(
                       6,
                       conditionalPanel(
                         condition = "input.DE_tabs=='DE genes' ",
                         downloadButton("download_de_table", "Download table (csv)")
                       ),
                       conditionalPanel(
                         condition = "input.DE_tabs=='Volcano plot' ",
                         downloadButton("download_volcano", "Download plot (PDF)")
                       ),
                       conditionalPanel(
                         condition = "input.DE_tabs=='Dispersion plot' ",
                         downloadButton("download_dispersion", "Download plot (PDF)")
                       ),
                       br(),
                       br(),
                       checkboxInput("checkbox_degenes", label = "DE genes", value = FALSE),
                       checkboxInput("checkbox_volcano", label = "Volcano plot", value = FALSE),
                       checkboxInput("checkbox_dispersion", label = "Dispersion plot", value = FALSE),
                       actionButton("add_de_analysis","Add to report")
                       
                     )
                   )
                 ),
                 mainPanel(
                   h3("DE Analysis"),
                   tabsetPanel(
                     type = "tabs", id = "DE_tabs",
                     tabPanel(
                       "DE genes",
                       uiOutput("help_text_DE_anal"),
                       # h3("Differential Expression Analysis"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         DT::dataTableOutput("DE_table")
                       )
                     ),
                     tabPanel(
                       "Volcano plot", # for DESeq and edgeR
                       h6("Volcano plot is only available for edgeR and DESeq2 methods"),
                       conditionalPanel(
                         condition = "input.n_rep=='1' && input.method1!='NOISeq'",
                         conditionalPanel(
                           condition = "$('html').hasClass('shiny-busy')",
                           div(img(src = "load.gif", width = 240, height = 180),
                               h4("Processing ... Please wait"),
                               style = "text-align: center;"
                           )
                         ),
                         conditionalPanel(
                           condition = "!$('html').hasClass('shiny-busy')",
                           plotOutput("volcano_plot")
                         )
                       ),
                       conditionalPanel(
                         condition = "input.method0=='NOISeq' || input.method1=='NOISeq'",
                         h6("Volcano Plot is only applicable to DESeq2 and edgeR")
                       )
                     ),
                     tabPanel(
                       "Dispersion plot", # for edgeR
                       h6("Dispersion plot is only available for edgeR and DESeq2 methods"),
                       conditionalPanel(
                         condition = "input.n_rep=='1' && input.method1!='NOISeq'",
                         # h3("Dispersion plot"),
                         conditionalPanel(
                           condition = "$('html').hasClass('shiny-busy')",
                           div(img(src = "load.gif", width = 240, height = 180),
                               h4("Processing ... Please wait"),
                               style = "text-align: center;"
                           )
                         ),
                         conditionalPanel(
                           condition = "!$('html').hasClass('shiny-busy')",
                           plotOutput("dispersion_plot")
                         )
                       ),
                       conditionalPanel(
                         condition = "input.method0=='NOISeq' || input.method1=='NOISeq'",
                         h6("Dispersion Plot is only applicable to DESeq2 and edgeR")
                       )
                     )
                   )
                 )
               ),
               tabPanel(
                 "Heatmap",
                 sidebarPanel(
                   conditionalPanel(
                     condition = "input.heatmap_tabs=='Heatmap'",
                     
                     radioButtons("heatmap_de_ind", label = "Choose data", choices = c("Indenpendent" = "ind", "DE result" = "de")),
                     numericInput("numOfCluster", "Number of clusters on rows", value = 2, min = 2, max = 30, step = 1),
                     conditionalPanel(
                       condition = "input.heatmap_de_ind == 'ind' ",
                       # selectInput('numOfGeno',"Number of genotypes (mutants)",choices=c(1)),
                       splitLayout(
                         numericInput("fold", "Fold change", value = 2, min = 1, step = 1),
                         numericInput("fold_ncol", "min. column", value = 2, min = 1, step = 1)
                       )
                       
                       # uiOutput("refGeno"),
                       # radioButtons('heatmap_value',"Values",
                       #              c('Fold change','Log fold change'))
                     ),
                     actionButton("heatmap_plot", "Plot"),
                     
                     downloadButton("downloadheatmap", "Download as PDF")
                     
                     # conditionalPanel(
                     #   condition = "input.heatmap_de_ind == 'ind' ",
                     #   h5('Specify names of the genotypes'),
                     #   uiOutput("expand_genonames")
                     # )
                   ),
                   
                   conditionalPanel(
                     condition = "input.heatmap_tabs=='Gene clusters'",
                     uiOutput("heatmap_display"),
                     conditionalPanel(
                       condition = "input.display_cluster=='ALL'",
                       downloadButton("downloadclusters", "Download as CSV")
                     )
                   ),
                   checkboxInput("checkbox_heatmap", label = "Heatmap", value = FALSE),
                   checkboxInput("checkbox_gene_cluster", label = "Gene Clusters", value = FALSE),
                   actionButton("add_heatmap","Add to report")
                   
                 ),
                 mainPanel(
                   h3("Heatmap"),
                   tabsetPanel(
                     type = "tabs", id = "heatmap_tabs",
                     tabPanel(
                       "Heatmap",
                       uiOutput("help_text_heatmap"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         plotOutput("heatmap.plot")
                       )
                     ),
                     tabPanel("Gene clusters", dataTableOutput("cluster.info"))
                   )
                 )
               ),
               
               ######## NOISE ######
               #############################################
               tabPanel(
                 "Noise",
                 sidebarPanel(
                   radioButtons("noise_situation", "Select desired noise plot between", choices = c("replicates" = "a", "genotypes (average of replicates)" = "b", "genotypes (no replicate)" = "c")),
                   conditionalPanel(
                     condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
                     textInput("noise_numOfRep", "Number of replicates", value = 1),
                     helpText("* Please order the sample columns in input file properly. Replicates of the same genotype should be put in adjacent columns.")
                   ),
                   conditionalPanel(
                     condition = "input.noise_situation=='b'",
                     uiOutput("noise_anchor_choices")
                   ),
                   conditionalPanel(
                     condition = "input.noise_situation=='c'",
                     selectInput("noise_anchor_c", "Anchor genotype", choices = "")
                   ),
                   radioButtons(
                     "noise_graph_type", "Graph type:",
                     c("Bar chart", "Line chart")
                   ),
                   actionButton("noise_plot", "Plot"),
                   downloadButton("downloadnoise", "Download as PDF"),
                   
                   conditionalPanel(
                     condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
                     h5("Specify names of the genotypes"),
                     uiOutput("expand_genonames_noise")
                   ),
                   
                   actionButton("add_noise","Add to report")
                 ),
                 mainPanel(
                   h3("Noise"),
                   uiOutput("help_text_Noise"),
                   conditionalPanel(
                     condition = "$('html').hasClass('shiny-busy')",
                     div(img(src = "load.gif", width = 240, height = 180), h4("Processing ... Please wait"), style = "text-align: center;")
                   ),
                   conditionalPanel(
                     condition = "!$('html').hasClass('shiny-busy')",
                     plotlyOutput("noise.plot")
                   )
                 )
               ),
               
               
               ###### ENTROPY #############
               #########################################
               tabPanel(
                 "Entropy",
                 sidebarPanel(
                   checkboxInput("tsflag", strong("Time series data"), FALSE),
                   conditionalPanel(
                     condition = "input.tsflag==true",
                     textInput("entropy_timepoints", "Number of time points"),
                     helpText("* Please order the sample columns in input file properly. Time series data of the same genotype should be put in adjacent columns.")
                   ),
                   radioButtons(
                     "entropy_graph_type", "Graph type:",
                     c("Bar chart", "Line chart")
                   ),
                   actionButton("submit_entropy","Plot"),
                   downloadButton("downloadentropy", "Download as PDF"),
                   conditionalPanel(
                     condition = "input.tsflag==true",
                     h5("Specify names of the genotypes"),
                     uiOutput("expand_genonames_entropy")
                   ),
                   
                   actionButton("add_entropy","Add to report")
                 ),
                 mainPanel(
                   h3("Shannon entropy"),
                   uiOutput("help_text_Entropy"),
                   plotlyOutput("entropy.plot")
                 )
               ),
               
               ################### Support Vector Machine ##################### 
               
               # tabPanel('SVM',
               #          sidebarPanel(
               #            actionButton("submit_svm","Submit")),
               #          mainPanel(
               #            h3('SVM Plot'),
               #            tabsetPanel(
               #              type = "tabs", id = "SVM_tabs",
               #              tabPanel("SVM Classification",
               #                       uiOutput("help_text_SVM"),
               #                       plotOutput('svm_plot')),
               #              tabPanel("Raw Plot", plotOutput("svm_df_plot"))
               #            )
               #          )),
               
               #####################################################################
               
               tabPanel(
                 't-SNE',
                 sidebarPanel(
                   splitLayout(
                     numericInput("perplexity_value","Perplexity value", min=1, value=2),
                     numericInput("no_of_pca","No. of PCs", min=1, value=30)
                     # numericInput("no_of_clusters","No. of clusters", min=2, value=2)
                   ),
                   radioButtons('tsne2_trans',"Transformation:",
                                c('None', 'log10')),
                   checkboxInput("tsne_cluster", strong("Kmeans clustering on columns"), FALSE),
                   conditionalPanel(
                     condition = "input.tsne_cluster == true",
                     numericInput("tsne_cluster_num", "Number of clusters:", min = 1, value = 2),
                     
                   ),
                   checkboxInput("tsne_text", strong("Display sample name"), FALSE),
                   actionButton("submit_tsne2","Plot"),
                   conditionalPanel(
                     condition = "input.tsne_tabs=='t-SNE plot",
                     downloadButton("download_tsne2", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.tsne_tabs=='t-SNE table'",
                     downloadButton("download_tsne", "Download as CSV")
                   ),
                   checkboxInput("checkbox_tsne_plot", label = "t-SNE Plot", value = FALSE),
                   checkboxInput("checkbox_tsne_table", label = "t-SNE table", value = FALSE),
                   actionButton("add_tsne","Add to report")
                   
                 ),
                 mainPanel(
                   h3('t-SNE Plot'),
                   
                   tabsetPanel(
                     type = "tabs", id = "tsne_tabs",
                     tabPanel("t-SNE plot", 
                              uiOutput("help_text_tsne"),
                              plotlyOutput('tsne2.plot')),
                     tabPanel("t-SNE table", 
                              DT::dataTableOutput("tsne_table") )
                   )
                 )),
               
               
               # Random forest tab is temporarily hidden.
               # to activate the RF panel: uncomment the below lines
               # tabPanel(
               #   "Random Forest",
               #   sidebarPanel(
               #     radioButtons(
               #       "analysis_type", "Choose Analysis Type",
               #       c("RF clustering" = "rf", "RAFSIL" = "rafsil")
               #     ),
               #     conditionalPanel(
               #       condition = "input.analysis_type=='rf'",  #rf
               #       splitLayout(
               #         numericInput("num_trees", "No. of trees", min = 1, value = 25),
               #         numericInput("num_clusters", "No. of clusters", min = 1, value = 2)
               #       ),
               #       radioButtons(
               #         "rf_trans", "Transformation:",
               #         c("None", "log10")
               #       ),
               #       actionButton("submit_rf", "Submit")
               #     ),
               #     conditionalPanel(
               #       condition = "input.analysis_type=='rafsil'",  #rafsil
               #       actionButton("submit_rafsil", "Submit")
               #     )
               #     # conditionalPanel(
               #     #          condition = "input.rf_tabs == 'RF plot'",
               #     #          downloadButton("downloadrfplot", "Download as PDF")
               #     #        ),
               #     #        conditionalPanel(
               #     #          condition = "input.rf_tabs == 'RF matrix'",
               #     #          downloadButton("downloadrfmatrix", "Download as PDF")
               #     #        )
               #   ),
               #   mainPanel(
               #     h3("Clustering With Random Forest"),
               #     tabsetPanel(type = "tabs",id="rf_tabs",
               #                 tabPanel("RF plot",
               #                          uiOutput("help_text_rf"),
               #                          plotlyOutput("rf.plot")),
               #                 tabPanel("RAFSIL plot", plotOutput("RAFSIL.plot")),
               #                 tabPanel("RF matrix", div(tableOutput('rf.matrix'), style = "font-size:80%"))
               #     )
               #   )
               # ),
               
               
               tabPanel(
                 "SOM",
                 sidebarPanel(
                   selectInput(inputId = "som_samples", label = "Samples used", choices = ""),
                   splitLayout(
                     numericInput("som_grid_h", "No. of horizontal grids", min = 1, value = 2),
                     numericInput("som_grid_v", "No. of vertical grids", min = 1, value = 2)
                   ),
                   numericInput("som_cluster_size", "No. of clusters (for cluster plot)", min = 2, value = 2),
                   radioButtons(
                     "som_trans", "Transformation:",
                     c("None", "log10")
                   ),
                   actionButton("submit_som", "Plot"),
                   br(),
                   br(),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Property plot'",
                     downloadButton("downloadProperty", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Count plot'",
                     downloadButton("downloadCount", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Codes plot'",
                     downloadButton("downloadCodes","Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Distance plot'",
                     downloadButton("downloadDistance", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Cluster plot'",
                     downloadButton("downloadCluster","Download as PDF")
                   ),
                   checkboxInput("checkbox_property", label = "Property plot", value = FALSE),
                   checkboxInput("checkbox_count", label = "Count plot", value = FALSE),
                   checkboxInput("checkbox_codes", label = "Codes plot", value = FALSE),
                   checkboxInput("checkbox_distance", label = "Distance plot", value = FALSE),
                   checkboxInput("checkbox_cluster", label = "Cluster plot", value = FALSE),
                   actionButton("add_som","Add to report"),
                   
                 ),
                 mainPanel(
                   h3("SOM Analysis"),
                   tabsetPanel(
                     type = "tabs", id = "som_tabs",
                     tabPanel("Property plot", 
                              uiOutput("help_text_SOM"),
                              plotOutput("som_property.plot")),
                     tabPanel("Count plot", plotOutput("som_count.plot")),
                     tabPanel("Codes plot", plotOutput("som_codes.plot")),
                     tabPanel("Distance plot", plotOutput("som_dist.plot")),
                     tabPanel("Cluster plot", plotOutput("som_cluster.plot"))
                   )
                 ),
                 tags$style(type = 'text/css', 
                            '.navbar { font-size: 17px;}'
                 )
               )
               
    ),
    ###############################################
    ###############################################
    ###############################################
    navbarMenu(
      "Gene Set Analysis",
      
      
      ########## Pathway Enrichment ##############
      #########################################
      tabPanel(
        "Gene Pathways Enrichment",
        tags$head(tags$style("#path_enri_visu{height:95vh !important;}")),
        sidebarLayout(
          sidebarPanel(
            withTags({
              div(class = "header",
                  p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gPro_gene_names.csv", target="_blank")),
              )
            }),
            fileInput("file_path_enri_gene", "Upload genes CSV file"),
            textInput("text_path_enri_gene", "Enter gene id"),
            actionButton("submit_path_enri_gene", "Submit"),br(),br(),
            selectInput("loadStyleFile_path_gene", "Select Style: ", choices=styles),
            # selectInput(inputId = "overlap_min", label = "Minimum Overlap", choices = ""),
            hidden(sliderInput("overlap_min_gene", "Minimum Overlap",
                               min = 0, max = 100,
                               value = 1)),
            hidden(selectInput("doLayout_path_gene", "Select Layout:",
                               choices=c("",
                                         "cose",
                                         "cola",
                                         "circle",
                                         "concentric",
                                         "breadthfirst",
                                         "grid",
                                         "random",
                                         "dagre",
                                         "cose-bilkent"))),
            hidden(actionButton("sfn_path_gene", "Select First Neighbor")),
            br(),br(),
            
            #             actionButton("fit_path_gene", "Fit Graph"),br(),br(),
            #             actionButton("fitSelected_path_gene", "Fit Selected"),br(),br(),
            #             actionButton("clearSelection_path_gene", "Clear Selection"), br(),br(),
            #             actionButton("removeGraphButton_path_gene", "Remove Graph"), br(),br(),
            #             actionButton("addRandomGraphFromDataFramesButton_path_gene", "Add Random Graph"),br(),br(),
            #             actionButton("getSelectedNodes_path_gene", "Get Selected Nodes"), br(),br(),
            #             htmlOutput("selectedNodesDisplay_path_gene"),
            #             checkboxInput("checkbox_plot", label = "Plot", value = FALSE),
            #             checkboxInput("checkbox_visualization", label = "Visualizationt", value = FALSE),
            #             actionButton("add_gene_path_enrich","Add to report"),
            
            
            hidden(actionButton("fit_path_gene", "Fit Graph")),br(),br(),
            hidden(actionButton("fitSelected_path_gene", "Fit Selected")),br(),br(),
            hidden(actionButton("clearSelection_path_gene", "Clear Selection")), br(),br(),
            hidden(actionButton("removeGraphButton_path_gene", "Remove Graph")), br(),br(),
            hidden(actionButton("addRandomGraphFromDataFramesButton_path_gene", "Add Random Graph")),br(),br(),
            hidden(actionButton("getSelectedNodes_path_gene", "Get Selected Nodes")), br(),br(),
            hidden(htmlOutput("selectedNodesDisplay_path_gene")),
            width=2
            
          ),
          mainPanel(
            h3("Pathways Enrichment"),
            tabsetPanel(
              type = "tabs", id = "path_enri_tab_gene",
              tabPanel("Plot",
                       uiOutput("help_text_path_enri_gene"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         plotlyOutput("path_enri.plot_gene", width = "100%", height = "100%")
                       ), 
              ),
              tabPanel(
                "Visualization",
                conditionalPanel(
                  condition = "$('html').hasClass('shiny-busy')",
                  div(img(src = "load.gif", width = 240, height = 180),
                      h4("Processing ... Please wait"),
                      style = "text-align: center;"
                  )
                ),
                conditionalPanel(
                  condition = "!$('html').hasClass('shiny-busy')",
                  cyjShinyOutput('path_enri_visu_gene', height=350)
                ),
              )
            )
          )
        )),
      
      
      ###### Tissue Expression #############
      #########################################
      tabPanel(
        "Tissue Expression",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv", target="_blank")),
            )
          }),
          fileInput("file_prot_expr", "Upload UniProt accession CSV file"),
          textInput("text_prot_expr","Enter Uniprot accession numbers"),
          actionButton("submit_prot_expr", "Submit"),br(),br(),
          downloadButton("prot_expr_download", "Download as CSV"),
          
          actionButton("add_tissue_exp","Add to report")
        ),
        mainPanel(
          h3("Tissue Expression"),
          uiOutput("help_text_prot_exp"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_expr_table")
          ), 
        )),
      
      
      ########## Co-expression #############
      #########################################
      tabPanel(
        "Co-expression",
        sidebarPanel(
          selectInput("organismID", "Choose organism:", 
                      choices= list(
                        "Arabidopsis_thaliana",
                        "Caenorhabditis_elegans",
                        "Danio_rerio",
                        "Drosophila_melanogaster",
                        "Escherichia_coli",
                        "Homo_sapiens",
                        "Mus_musculus",
                        "Rattus_norvegicus",
                        "Saccharomyces_cerevisiae"
                      )),
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_names.csv", target="_blank")),
            )
          }),
          fileInput("file_gene", "Upload genes CSV file"),
          textInput("text_gene","Enter gene ids"),
          actionButton("genemania_submit", "Submit")
        ),
        mainPanel(
          h3("Co-expression"),
          uiOutput("help_text_gene_mania"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            div(id = "hide_link", 
                p("Please click", htmlOutput("linkCo")))
            %>% shinyjs::hidden()
          ),
        ))
      #########################################
      
      #########################################
      
    ),
    navbarMenu(
      "Protein Set Analysis",
      tabPanel("Upload a Protein Set",
               sidebarPanel(fileInput("file_protein_set", "Upload UniProt accession CSV file"),
                            textInput("text_protein_set", "Enter UniProt accession csv file"),
                            actionButton("submit_protein_set", "Submit")),
               mainPanel(
                 h3("Upload a Protein Set"),
                 uiOutput("help_text_protein_set")
               )),
      
      ########## Gene Ontology #############
      #########################################
      tabPanel(
        "Gene ontology",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv", target="_blank")),
            )
          }),
          fileInput("file_uniprot", "Upload UniProt accession CSV file"),
          textInput("text_uniprot", "Enter UniProt accession csv file"),
          actionButton("submit_uniprot", "Submit"),br(),br(),
          conditionalPanel(
            condition = "input.uniprot_tabs == 'Biological process'",
            downloadButton("download_bio_plot", "Download Plot"),br(),br(),
            downloadButton("download_bio_pro", "Download Table")
          ),
          conditionalPanel(
            condition = "input.uniprot_tabs == 'Molecular function'",
            downloadButton("download_mole_plot", "Download Plot"),br(),br(),
            downloadButton("download_mole_func", "Download Table")
          ),
          conditionalPanel(
            condition = "input.uniprot_tabs == 'Cellular component'",
            downloadButton("download_cell_plot", "Download Plot"),br(),br(),
            downloadButton("download_cell_comp", "Download Table")
          ),
          checkboxInput("checkbox_bio_proc_plot", label = "Biological Process Plot", value = FALSE),
          checkboxInput("checkbox_bio_proc_table", label = "Biological Process Table", value = FALSE),
          checkboxInput("checkbox_mol_func_plot", label = "Molecular Function Plot", value = FALSE),
          checkboxInput("checkbox_mol_func_table", label = "Molecular Function Table", value = FALSE),
          checkboxInput("checkbox_cell_comp_plot", label = "Cellular Component Plot", value = FALSE),
          checkboxInput("checkbox_cell_comp_table", label = "Cellular Component Plot", value = FALSE),
          actionButton("add_gene_onto","Add to report"),
        ),
        mainPanel(
          h3("Gene ontology"),
          tabsetPanel(
            type = "tabs", id = "uniprot_tabs",
            tabPanel("Biological process",
                     uiOutput("help_text_bio_pr"),
                     conditionalPanel(
                       condition = "$('html').hasClass('shiny-busy')",
                       div(img(src = "load.gif", width = 240, height = 180),
                           h4("Processing ... Please wait"),
                           style = "text-align: center;"
                       )
                     ),
                     conditionalPanel(
                       condition = "!$('html').hasClass('shiny-busy')",
                       plotOutput("uniprotbioplot"),
                       shiny::dataTableOutput("uniprot_biotable")),
            ),
            tabPanel("Molecular function",
                     conditionalPanel(
                       condition = "$('html').hasClass('shiny-busy')",
                       div(img(src = "load.gif", width = 240, height = 180),
                           h4("Processing ... Please wait"),
                           style = "text-align: center;"
                       )
                     ),
                     conditionalPanel(
                       condition = "!$('html').hasClass('shiny-busy')",
                       plotOutput("uniprot_molcplot"),
                       shiny::dataTableOutput("uniprot_molctable")
                     ),
            ),
            tabPanel("Cellular component",
                     conditionalPanel(
                       condition = "$('html').hasClass('shiny-busy')",
                       div(img(src = "load.gif", width = 240, height = 180),
                           h4("Processing ... Please wait"),
                           style = "text-align: center;"
                       )
                     ),
                     conditionalPanel(
                       condition = "!$('html').hasClass('shiny-busy')",
                       plotOutput("uniprot_celplot"),
                       shiny::dataTableOutput("uniprot_celtable")
                     ),
            )
          )
        )), ##### Gene ontlogy closing 
      
      
      ###### Protein Interaction #############
      #########################################
      tabPanel(
        "Protein Interactions",
        tags$head(tags$style("#cyjShiny{height:95vh !important;}")),
        sidebarLayout(
          sidebarPanel(
            withTags({
              div(class = "header",
                  p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv", target="_blank")),
              )
            }),
            fileInput("file_prot_Int", "Upload UniProt accession CSV file"),
            textInput("text_prot_Int","Enter UniProt accession numbers"),
            actionButton("submit_prot_Int", "Submit"),br(),br(),
            selectInput("loadStyleFile", "Select Style: ", choices=styles),
            selectInput("doLayout", "Select Layout:",
                        choices=c("",
                                  "cose",
                                  "cola",
                                  "circle",
                                  "concentric",
                                  "breadthfirst",
                                  "grid",
                                  "random",
                                  "dagre",
                                  "cose-bilkent")),
            # selectInput("showCondition", "Select Condition:", choices=rownames(output$tbl.lfc)),
            # selectInput("selectName", "Select Node by ID:", choices = c("", sort(tbl.nodes$id))),
            actionButton("sfn", "Select First Neighbor"),
            br(),br(),
            actionButton("fit", "Fit Graph"),br(),br(),
            actionButton("fitSelected", "Fit Selected"),br(),br(),
            actionButton("clearSelection", "Clear Selection"), br(),br(),
            actionButton("removeGraphButton", "Remove Graph"), br(),br(),
            actionButton("addRandomGraphFromDataFramesButton", "Add Random Graph"),br(),br(),
            actionButton("getSelectedNodes", "Get Selected Nodes"), br(),br(),
            htmlOutput("selectedNodesDisplay"),
            checkboxInput("checkbox_pp_visu", label = "Visualization", value = FALSE),
            checkboxInput("checkbox_pp_interact", label = "Protein Interactions", value = FALSE),
            checkboxInput("checkbox_prot_name", label = "Protein Names", value = FALSE),
            actionButton("add_pp_inter","Add to report"),
          ),
          mainPanel(
            h3("Protein-Protein Interactions"),
            tabsetPanel(
              type = "tabs", id = "prot_inte_tab",
              tabPanel("Visualization",
                       uiOutput("help_text_p_inte"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         cyjShinyOutput('cyjShiny', height=350)
                       )),
              tabPanel("Protein Interaction",
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         DT::dataTableOutput("prot_int_table")
                       )),
              tabPanel("Protein Name",
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         DT::dataTableOutput("prot_name_table")
                       ))
            ),
            
          )
          #  mainPanel(cyjShinyOutput('cyjShiny', height=400),
          #           width=10,
          #           tabPanel("Protein Name",
          #           DT::dataTableOutput("prot_name_table"))
          #  )
        )),   #### PPI Closing
      
      
      ###### Protein Function #############
      #########################################
      tabPanel(
        "Protein Function",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv", target="_blank")),
            )
          }),
          fileInput("file_prot_func", "Upload UniProt accession CSV file"),
          textInput("text_prot_func","Enter UniProt accession numbers"),
          actionButton("submit_prot_func", "Submit"),br(),br(),
          downloadButton("prot_func_download", "Download as CSV"),
          actionButton("add_prot_func","Add to report"),
        ),
        mainPanel(
          h3("Protein Function"),
          uiOutput("help_text_prot_fn"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_func_table")
          ),
        )),  ### Protein function closing 
      
      
      ###### Subcellular Localization #############
      #########################################
      tabPanel(
        "Subcellular Localization",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv", target="_blank")),
            )
          }),
          fileInput("file_prot_local", "Upload UniProt accession CSV file"),
          textInput("text_prot_local","Enter UniProt accession numbers"),
          actionButton("submit_prot_local", "Submit"),br(),br(),
          downloadButton("prot_local_download", "Download as CSV"),
          actionButton("add_subcell_loc","Add to report"),
          
        ),
        mainPanel(
          h3("Subcellular Localization"),
          uiOutput("help_text_sub_loc"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_local_table")
          ), 
        )),  ### Subcellular closing 
      
      ########## Protein Domains ##############
      #########################################
      tabPanel(
        "Protein Domains",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv", target="_blank")),
            )
          }),
          fileInput("file_prot_domain", "Upload UniProt accession CSV file"),
          textInput("text_prot_domain","Enter UniProt accession numbers"),
          actionButton("submit_prot_domain", "Submit"),br(),br(),
          downloadButton("prot_domain_download", "Download as CSV"),
          actionButton("add_prot_dom","Add to report")
        ),
        mainPanel(
          h3("Protein Domains"),
          uiOutput("help_text_pro_dom"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_domain_table")
          ), 
        )),   ### Protein domain closing 
      
      
      ######### Protein Sequences #############
      #########################################
      tabPanel(
        "Protein properties",
        sidebarPanel(
          fileInput("file_prot_seq", "Upload UniProt accession CSV file"),
          textInput("text_prot_seq","Enter UniProt accession numbers"),
          actionButton("submit_prot_Seq", "Submit"),br(),br(),
          shinyjs::hidden(downloadButton('downloadData', 'Download Sequence FASTA')),
          checkboxInput("checkbox_seq_charge", label = "Sequence Charge", value = FALSE),
          checkboxInput("checkbox_seq_acid", label = "Sequence Acidity", value = FALSE),
          checkboxInput("checkbox_seq_grav_ind", label = "Sequence Gravy Index", value = FALSE),
          checkboxInput("checkbox_physio_prop", label = "All Physiochecmical Properties", value = FALSE),
          actionButton("add_prot_seq","Add to report"),
          
        ),
        mainPanel(
          h3("Protein Sequences"),
          tabsetPanel(type = "tabs",
                      tabPanel("Sequence information",uiOutput("help_text_prot_seq")),
                      tabPanel(
                        "Sequence charge", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("ChargePlot") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Sequence acidity", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("AcidityPlot") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Sequence gravy index",
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("GravyPlot") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "All physicochemical properties", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("SequencePlot" , width="900px",height="750px")
                        ),
                      )
                      #################################################################################
                      
          ),
          
        )), # For protein properties tab panel 
      
      #Here we start evolution tab panel 
      tabPanel(
        "Evolutionary Analysis",
        sidebarPanel(
          fileInput("file_prot_seq_evol", "Upload UniProt accession CSV file"),
          textInput("text_prot_seq_evol","Enter UniProt accession numbers"),
          actionButton("submit_prot_seq_evol","Submit"),
          checkboxInput("checkbox_prot_gene", label = "Protein's genes tree", value = FALSE),
          checkboxInput("checkbox_prot_chrom", label = "Protein's chromosomal location", value = FALSE),
          checkboxInput("checkbox_prot_evol", label = "Evolutionary analysis", value = FALSE),
          actionButton("add_prot_evol_analysis","Add to report"),
          
          
        ),
        mainPanel(
          h3("Protein Evolutionary analysis"),
          tabsetPanel(type = "tabs",
                      tabPanel("Evolutionary information",uiOutput("help_text_prot_seq_evol")),
                      tabPanel(
                        "Protein's genes tree", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          radialNetworkOutput("GenePlot", width="900px",height="750px") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Protein's chromosomal location", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("Chromo", width="1200px",height="1500px") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Evolutionary analysis", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("Phylogenetic", width="900px",height="750px") # , width="900px",height="750px"
                        ),
                      )
                      
          ),
          
        )), #For protein evolution tab panel 
      
      tabPanel(
        "Pathological Analysis",
        sidebarPanel(
          fileInput("file_prot_seq_Patho", "Upload UniProt accession CSV file"),
          textInput("text_prot_seq_Patho","Enter UniProt accession numbers"),
          actionButton("submit_prot_seq_Patho","Submit"),
          checkboxInput("checkbox_dise_role", label = "Protein's disease role", value = FALSE),
          checkboxInput("checkbox_dise_dist", label = "Protein's disease distribution", value = FALSE),
          actionButton("add_prot_path_analysis","Add to report"),
        ),
        mainPanel(
          h3("Protein pathological analysis"),
          tabsetPanel(type = "tabs",
                      tabPanel("Pathological information",uiOutput("help_text_prot_seq_Patho")),
                      tabPanel(
                        "Protein's disease role", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          dataTableOutput("DisaeseTable") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Protein's disease distribution", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          bubblesOutput("DiseasePlot") #  , width="750px",height="750px"
                        ),
                      )
                      
          ),
          
        )), #For protein pathology tab panel 
      
      
      #########################################
      tabPanel(
        "Protein Complex Enrichment",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv", target="_blank")),
            )
          }),
          fileInput("file_complex_prot", "Upload UniProt accession CSV file"),
          textInput("text_complex_prot","Enter UniProt accession numbers"),
          actionButton("submit_complex_prot", "Submit"),br(),br(),
          downloadButton("complex_download_prot", "Download as CSV"),
          actionButton("add_prot_comp_enrich","Add to report"),
        ),
        mainPanel(
          h3("Complex Enrichment"),
          uiOutput("help_text_complex_en_prot"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            shiny::dataTableOutput("complex_table_prot")
          ),
        )), ###Complex enrichment
      
      
      ########## Pathway Enrichment ##############
      #########################################
      tabPanel(
        "Pathways Enrichment",
        tags$head(tags$style("#path_enri_visu{height:95vh !important;}")),
        sidebarLayout(
          sidebarPanel(
            withTags({
              div(class = "header",
                  p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gPro_gene_names.csv", target="_blank")),
              )
            }),
            fileInput("file_path_enri_prot", "Upload genes CSV file"),
            textInput("text_path_enri_prot", "Enter gene id"),
            actionButton("submit_path_enri_prot", "Submit"),br(),br(),
            selectInput("loadStyleFile_path_prot", "Select Style: ", choices=styles),
            # selectInput(inputId = "overlap_min", label = "Minimum Overlap", choices = ""),
            hidden(sliderInput("overlap_min_prot", "Minimum Overlap",
                               min = 0, max = 100,
                               value = 1)),
            hidden(selectInput("doLayout_path_prot", "Select Layout:",
                               choices=c("",
                                         "cose",
                                         "cola",
                                         "circle",
                                         "concentric",
                                         "breadthfirst",
                                         "grid",
                                         "random",
                                         "dagre",
                                         "cose-bilkent"))),
            hidden(actionButton("sfn_path_prot", "Select First Neighbor")),
            br(),br(),
            
            #             actionButton("fit_path_prot", "Fit Graph"),br(),br(),
            #             actionButton("fitSelected_path_prot", "Fit Selected"),br(),br(),
            #             actionButton("clearSelection_path_prot", "Clear Selection"), br(),br(),
            #             actionButton("removeGraphButton_path_prot", "Remove Graph"), br(),br(),
            #             actionButton("addRandomGraphFromDataFramesButton_path_prot", "Add Random Graph"),br(),br(),
            #             actionButton("getSelectedNodes_path_prot", "Get Selected Nodes"), br(),br(),
            #             htmlOutput("selectedNodesDisplay_path_prot"),
            #             checkboxInput("checkbox_plot_prot", label = "Plot", value = FALSE),
            #             checkboxInput("checkbox_visualization_prot", label = "Visualizationt", value = FALSE),
            #             actionButton("add_prot_path_enrich","Add to report"),
            
            hidden(actionButton("fit_path_prot", "Fit Graph")),br(),br(),
            hidden(actionButton("fitSelected_path_prot", "Fit Selected")),br(),br(),
            hidden(actionButton("clearSelection_path_prot", "Clear Selection")), br(),br(),
            hidden(actionButton("removeGraphButton_path_prot", "Remove Graph")), br(),br(),
            hidden(actionButton("addRandomGraphFromDataFramesButton_path_prot", "Add Random Graph")),br(),br(),
            hidden(actionButton("getSelectedNodes_path_prot", "Get Selected Nodes")), br(),br(),
            hidden(htmlOutput("selectedNodesDisplay_path_prot")),
            width=2
            
          ),
          mainPanel(
            h3("Pathways Enrichment"),
            tabsetPanel(
              type = "tabs", id = "path_enri_tab_prot",
              tabPanel("Plot",
                       uiOutput("help_text_path_enri_prot"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         plotlyOutput("path_enri.plot_prot")
                       ), 
              ),
              tabPanel(
                "Visualization",
                conditionalPanel(
                  condition = "$('html').hasClass('shiny-busy')",
                  div(img(src = "load.gif", width = 240, height = 180),
                      h4("Processing ... Please wait"),
                      style = "text-align: center;"
                  )
                ),
                conditionalPanel(
                  condition = "!$('html').hasClass('shiny-busy')",
                  cyjShinyOutput('path_enri_visu_prot', height=350)
                ),
              )
            )
          )
        )) #Pathway enrichemnt analysis 
      
    ), ####Nav bar closing 
    navbarMenu('RDS Data Analysis',
               tabPanel("Documentation", value=1,
                        uiOutput('markdown')
               ),
               tabPanel(
                 "Upload RDS",
                 value = "rds_tab",
                 tabsetPanel(
                   id = "upload_rds_tab",
                   tabPanel(
                     "Upload Data",
                     sidebarPanel(
                       fileInput("rds_file", "Upload RDS File", accept = ".rds")
                     ),
                     mainPanel()
                   )
                 )
               ),
               tabPanel("Double Marker", value=2,
                        br(),
                        div(style="display: inline-block;vertical-align:top; width: 19%;",
                            selectInput("dataset", "Numeric Analysis Type:",
                                        c('Numeric Metadata', 'Genes','PCs'))),
                        div(style="display: inline-block;vertical-align:top; width: 19%;",
                            selectInput("reduction_double", "Reduction:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 19%;",
                            selectInput("categorical", "Identity:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 19%;",
                            selectInput("numeric", "Primary Numeric:", "")),
                        div(style="display: inline-block;vertical-align:top; width: 19%;",
                            selectInput('numeric2', 'Secondary Numeric', "")),
                        actionButton("add_doublemarker", "Add to report"),
                        mainPanel(width = 12,
                                  br(),
                                  br(),
                                  #h3(textOutput("caption")),
                                  plotOutput("MarkerGenePlot"),
                                  plotOutput("ViolinPlot"),
                                  plotOutput("CategoricalPlot")
                        )
               ),
               tabPanel("Single Marker", value=3,
                        br(),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("dataset_single", "Numeric Analysis Type:",
                                        c('Numeric Metadata', 'Genes','PCs'))),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("reduction_single", "Reduction:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("categorical_single", "Identity:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("numeric_single", "Primary Numeric:", "")),
                        actionButton("add_singlemarker", "Add to report"),
                        mainPanel(width = 12,
                                  br(),
                                  br(),
                                  #h3(textOutput("caption")),
                                  plotOutput("MarkerGenePlotSingle"),
                                  plotOutput("ViolinPlotSingle"),
                                  plotOutput("CategoricalPlotSingle")
                        )
               ),
               tabPanel("Multiple Feature Plot", value=5,
                        br(),
                        selectInput("multiple_feature_categorical_plot", "Identity:",
                                    NULL),
                        selectInput("multiple_feature_reduction", "Reduction:",
                                    NULL),
                        selectizeInput("multiple_feature_list", "Primary Numeric: \n 
                                                  - Csv format works best here if pasted in from premade lists. \n
                                                  - Optimal for >5 and <16 input. \n
                                                  - To be most efficient when removing entries hold SHIFT and click all, then delete.", "", 
                                       options = list(
                                         maxItems=16,
                                         delimiter = ',',
                                         create = I("function(input, callback){
                                              return {
                                                value: input,
                                                text: input
                                               };
                                            }")),
                                       selected = NULL, multiple = TRUE), ## and switch multiple to True,
                        actionButton("add_multiplefeature", "Add to report"),
                        mainPanel(width = 12,
                                  br(),
                                  br(),
                                  plotOutput("MultipleFeatureCategoricalPlot"),
                                  plotOutput("MultipleFeaturePlot",  height = "1000px")
                        )
               ),
               tabPanel("Cluster Tree", value=6,
                        br(),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("identity_tree", "Identity:",
                                        NULL)),
                        actionButton("add_clustertree", "Add to report"),
                        mainPanel(width = 12,
                                  br(),
                                  br(),
                                  #h3(textOutput("caption")),
                                  plotOutput("ClusterTree"),
                        )
               ),
               tabPanel("Seperated Feature", value=7,
                        br(),
                        div(style="display: inline-block;vertical-align:top; width: 20%;",
                            selectInput("dataset_seperated", "Numeric Analysis Type:",
                                        c('Genes', 'Numeric Metadata','PCs'))),
                        div(style="display: inline-block;vertical-align:top; width: 20%;",
                            selectInput("reduction_seperated", "Reduction:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 20%;",
                            selectInput("identity_seperated", "Cell Type/Cluster:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 20%;",
                            selectInput("identity_seperated2", "Identity:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 20%;",
                            selectInput("numeric_seperated", "Primary Numeric:", "")),
                        actionButton("add_separatedfeature", "Add to report"),
                        mainPanel(width = 12,
                                  br(),
                                  br(),
                                  #h3(textOutput("caption")),
                                  plotOutput("SeperatedFeature", height = "500px"),
                                  plotOutput("SeperatedDim"),
                                  plotOutput("SeperatedViolin", width="2000px"),
                                  tableOutput("SeperatedCounts")
                                  
                        )
               ),
               tabPanel("Seperated Categorical", value=8,
                        br(),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("reduction_seperated_categorical", "Reduction:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("identity_seperated_categorical", "Identity:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("identity2_seperated_categorical", "Secondary Identity:", "")),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("split_seperated_categorical", "Split By:",
                                        NULL)),
                        actionButton("add_separatecategorical", "Add to report"),
                        mainPanel(width = 12,
                                  br(),
                                  br(),
                                  #h3(textOutput("caption")),
                                  plotOutput("SeperatedIdentityCategorical", height = "500px"),
                                  plotOutput("SeperatedIdentity2Categorical"),
                                  plotOutput("SeperatedCountsCategorical")
                                  
                        )
               ),
               tabPanel("Marker Table", value=9,
                        br(),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("identity_table", "Identity:",
                                        NULL)),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("markers_table", "Get markers for:", "", multiple = TRUE)),
                        div(style="display: inline-block;vertical-align:top; width: 24%;",
                            selectInput("compare_table", "Compare to (blank is all other groups):", "", multiple = TRUE)),
                        mainPanel(width = 12,
                                  br(),
                                  br(),
                                  #h3(textOutput("caption")),
                                  tableOutput("markers")
                        )
               )
    ),
    navbarMenu('Single Cell Analysis',
               tabPanel(
                 "Cleaner",
                 value = "cleaner_tab",
                 tabsetPanel(
                   id = "analysis_pipeline",
                   tabPanel(
                     "Clean data",
                     sidebarPanel(
                       fileInput("file1", "Upload CSV File", accept = ".csv"),
                       downloadButton("downloadData", "Download Clean Data")
                     ),
                     mainPanel()
                   )
                 )
               ), 
               tabPanel(
                 "Raw Count/Expression Matrix", 
                 value = "expression_tab", 
                 tabsetPanel(
                   id = "expression_mat", 
                   tabPanel(
                     "Read Raw Count/Expression Matrix",
                     sidebarPanel(
                       fileInput("file2", "Upload CSV File", accept = ".csv")
                       
                     ),
                     mainPanel()
                   )
                 )
               ), 
               tabPanel(
                 "CCA Integration", 
                 value = "CCA_tab", 
                 tabsetPanel(
                   id = "CCA_integration", 
                   tabPanel(
                     "Upload Data",
                     sidebarPanel(
                       fileInput("file3", "Upload CSV File", accept = ".csv", multiple = TRUE), 
                       uiOutput("vlnplot3_selectize"), 
                       uiOutput("vlnplot4_selectize"), 
                       uiOutput("vlnplot5_selectize"), 
                       uiOutput("featureplot1_selectize"), 
                       uiOutput("featureplot2_selectize"), 
                       actionButton("add_cca", "Add to report")
                     ),
                     mainPanel(
                       plotOutput("vlnplot"),
                       plotOutput("combineplot1"), 
                       plotOutput("elbowplot2"), 
                       plotOutput("combineplot2"), 
                       plotOutput("combineplot3"), 
                       plotOutput("combineplot4"), 
                       plotOutput("combineplot5"), 
                       plotOutput("combineplot6"), 
                       plotOutput("combineplot7"),
                       plotOutput("dimplot1"), 
                       plotOutput("vlnplot3"), 
                       plotOutput("jackstrawplot"), 
                       plotOutput("vlnplot4"), 
                       plotOutput("vlnplot5"), 
                       plotOutput("featureplot1"), 
                       plotOutput("dotplot1"), 
                       plotOutput("dimplot2"), 
                       plotOutput("featureplot2")
                     )
                   )
                 )
               ), 
               tabPanel(
                 "Harmony Integration", 
                 value = "harmony_tab", 
                 tabsetPanel(
                   id = "harmony_integration", 
                   actionButton("add_harmony", "Add to report"),
                   mainPanel(
                     plotOutput("elbowplot1"), 
                     plotOutput("before"), 
                     plotOutput("after")
                   )
                 )
               ), 
               tabPanel(
                 "Non Integration", 
                 value = "nonintegration_tab", 
                 tabsetPanel(
                   id = "non_integration", 
                   tabPanel(
                     "Filter Genes",
                     sidebarPanel(
                       uiOutput("vlnplot16_selectize"), 
                       uiOutput("vlnplot17_selectize"), 
                       uiOutput("featureplot9_selectize"), 
                       actionButton("add_nonintegration", "Add to report")
                     ),
                     mainPanel(
                       plotOutput("vlnplot6"), 
                       plotOutput("combineplot8"), 
                       plotOutput("combineplot9"), 
                       plotOutput("vizdimloadings1"), 
                       plotOutput("dimplot3"), 
                       plotOutput("dimheatmap1"), 
                       plotOutput("dimheatmap2"), 
                       plotOutput("jackstrawplot2"), 
                       plotOutput("elbowplot3"), 
                       plotOutput("dimplot_umap"), 
                       plotOutput("vlnplot16"), 
                       plotOutput("vlnplot17"), 
                       plotOutput("featureplot9")
                     )
                   )
                 )
               ), 
               tabPanel(
                 "Auto Annotation", 
                 value = "autoannotation_tab", 
                 tabsetPanel(
                   id = "auto_annotation", 
                   tabPanel(
                     "Filter Genes",
                     sidebarPanel(
                       uiOutput("vlnplot7_selectize"), 
                       uiOutput("vlnplot8_selectize"), 
                       uiOutput("vlnplot9_selectize"), 
                       uiOutput("vlnplot10_selectize"), 
                       uiOutput("featureplot5_selectize"), 
                       uiOutput("ridgeplot2_selectize"), 
                       actionButton("add_auto", "Add to report")
                     ),
                     mainPanel(
                       plotOutput("combineplot11"), 
                       plotOutput("combineplot12"), 
                       plotOutput("combineplot13"), 
                       plotOutput("scoreheatmap1"), 
                       plotOutput("deltadistribution"), 
                       plotOutput("pheatmap1"), 
                       plotOutput("dimplot4"), 
                       plotOutput("scoreheatmap2"), 
                       plotOutput("dimplot5"), 
                       plotOutput("dimplot6"), 
                       plotOutput("vlnplot7"), 
                       plotOutput("vlnplot8"), 
                       plotOutput("vlnplot9"), 
                       plotOutput("vlnplot10"), 
                       plotOutput("featureplot5"), 
                       plotOutput("doheatmap1"), 
                       plotOutput("ridgeplot2"), 
                       plotOutput("doheatmap2")
                     )
                   )
                 )
               ), 
               tabPanel(
                 "Differential Expression", 
                 value = "diff_expression_tab", 
                 tabsetPanel(
                   id = "differential_expression", 
                   tabPanel(
                     "Filter Genes",
                     sidebarPanel(
                       uiOutput("de_selectize"), 
                       uiOutput("idents1_selectize"), 
                       actionButton("add_diffexp", "Add to report")
                     ),
                     mainPanel(
                       plotOutput("v1")
                     )
                   )
                 )
               ), 
               tabPanel(
                 "Pseudobulk", 
                 value = "pseudobulk_tab",
                 tabsetPanel(
                   id = "pseudobulk",
                   tabPanel(
                     "Filter Genes",
                     sidebarPanel(
                       uiOutput("pseudobulk_selectize"), 
                       uiOutput("idents2_selectize"), 
                       actionButton("add_pseudobulk", "Add to report")
                     ),
                     mainPanel(
                       plotOutput("vb1")
                     )
                   )
                 )
               ), 
               tabPanel(
                 "Compare between SCDE and PSCDE", 
                 value = "compare_tab", 
                 tabsetPanel(
                   id = "compare",
                   tabPanel(
                     "Filter Genes",
                     sidebarPanel(
                       uiOutput("vlnplot14_selectize"), 
                       uiOutput("vlnplot15_selectize"), 
                       actionButton("add_comparison", "Add to report")
                     ),
                     mainPanel(
                       plotOutput("vlnplot11"),
                       plotOutput("vlnplot12"),
                       plotOutput("vlnplot13"),
                       plotOutput("vlnplot14"),
                       plotOutput("vlnplot15")
                     )
                   )
                 )
               ), 
               tabPanel(
                 "Monocle3 Workflow", 
                 value = "monocle3_tab", 
                 tabsetPanel(
                   id = "monocle3", 
                   tabPanel(
                     "Select Genes",
                     sidebarPanel(
                       uiOutput("monocle_selectize"), 
                       actionButton("add_monocle3", "Add to report")
                     ),
                     mainPanel(
                       plotOutput("a1_monocle"),
                       plotOutput("a2_monocle"),
                       plotOutput("clusters"),
                       plotOutput("plot_cell"),
                       plotOutput("plot_cell2"),
                       plotOutput("pseudotime"),
                       plotOutput("featureplot7"),
                       plotOutput("featureplot8")
                     )
                   )
                 )
               ), 
               tabPanel(
                 "Save", 
                 value = "save_tab", 
                 tabsetPanel(
                   id = "save_report", 
                   tabPanel(
                     "Download Final Data",
                     sidebarPanel(
                       downloadButton("downloadReport", "Download Data")
                     ),
                     mainPanel()
                   )
                 )
               )
    ),
    tabPanel('Generate Report',fluidPage(
      
      uiOutput("error_text_report"),
      tags$div( h2("GeneCloudOmics Report",align = "center"), id = 'placeholder'),
      capture_pdf(
        selector = "#placeholder",
        filename = "results",
        icon("camera"), "Download as PDF"
      )
    ))
    
  )
  
)

####################################################

server <- function(input, output, session) {
  # <<<<<<< develop
  
  value_var<- reactiveValues()
  value_var$geo_file_type<-"none"
  
  
  # >>>>>>> master
  gene_mania_link <- reactiveVal("https://genemania.org")
  count_fasta <- reactiveVal(0)
  count_id <- reactiveVal(0)
  ####################download report#############
  output$error_text_report <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          You have not added anything to the report. Please perform the analysis first.
          </b>
        </p>
      </center>
    ")
  })
  
  observeEvent(input$add_scatter, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("Scatter Plot",align = "center"), plotOutput("scatter", height = 500),br(),br()))
      ),
    )
    output$scatter <- renderPlot({
      scatterplot()
    })
  })
  observeEvent(input$add_distfit,{
    hide("error_text_report")
    if(input$checkbox_distfit == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Distribution Fit",align = "center"), plotOutput("dis_fit", height = 500),br(),br()))
        )
      )
      output$dis_fit <- renderPlot({
        distplot()
      })
      
    }
    if(input$checkbox_aic == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("AIC Table",align = "center"), div(tableOutput("dist_aic"), style = "font-size:80%"),br(),br()))
        )
      )
      output$dist_aic <- renderTable({
        distaic()
      })
      
    }
    
  })
  observeEvent(input$add_correlation, {
    hide("error_text_report")
    if(input$checkbox_corrheatmap == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Correlation Heatmap",align = "center"), plotOutput("corr_hm", height = 500),br(),br()))
        )
      )
      output$corr_hm <- renderPlot({
        corrplot1()
      })
      
    }
    if(input$checkbox_corrplot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Correlation Plot",align = "center"), plotlyOutput("corr_m", height = 500),br(),br()))
        )
      )
      output$corr_m <- renderPlot({
        corrplot2()
      })
      
      
    }
    if(input$checkbox_corrmat == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Correlation Matrix",align = "center"), div(tableOutput("corr_mat"), style = "font-size:80%"),br(),br()))
        )
      )
      output$corr_mat <- renderTable({
        cor_df()
      })
      
    }
  })
  observeEvent(input$add_pca, {
    hide("error_text_report")
    if(input$checkbox_pcavar == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          fluidRow(
            column(width=2),
            column(width= 8,h3("PCA Variance",align = "center"), plotlyOutput("pca_var", height = 500),br(),br()))
        )
      )
      output$pca_var <- renderPlotly({
        pcavarplot()
      })
      
    }
    if(input$checkbox_pca2d == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("PCA 2D Plot",align = "center"), plotlyOutput("pca_2d", height = 500),br(),br()))
        )
      )
      output$pca_2d <- renderPlotly({
        pca2dplot()
      })
      
      
    }
    if(input$checkbox_pca3d == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          h3("PCA 3D Plot"),
          fluidRow(plotlyOutput("pca_3d"),br(),br())
        )
      )
      output$pca_3d <- renderPlotly({
        pca3dplot()
      })
      
      
    }
  })
  observeEvent(input$add_de_analysis, {
    hide("error_text_report")
    if(input$checkbox_volcano == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Volcano Plot",align = "center"), plotOutput("vol_plot", height = 500),br(),br()))
        )
      )
      output$vol_plot <- renderPlot({
        volcano_plot()
      })
      
    }
    if(input$checkbox_dispersion == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Dispersion Plot",align = "center"), plotlyOutput("dis_plot", height = 500),br(),br()))
        )
      )
      output$dis_plot <- renderPlot({
        dispersion_plot()
      })
      
      
    }
    if(input$checkbox_corrmat == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("DE genes",align = "center"), div(DT::dataTableOutput("DE_gene"), style = "font-size:80%"),br(),br()))
        )
      )
      output$DE_gene <- DT::renderDataTable({
        res.df <- de_no_filt()
        p_val <- input$p_val
        fc <- input$fc
        rep_number <- input$n_rep
        if (input$submit_DE > 0) {
          res.df.filt <- de_filt(res.df, p_val, fc, rep_number)
          res.df.filt
        }
      })
      
    }
  })
  observeEvent(input$add_heatmap,{
    hide("error_text_report")
    if(input$checkbox_heatmap == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Heatmap",align = "center"), plotOutput("heat_plot", height = 500),br(),br()))
        )
      )
      output$heat_plot <- renderPlot({
        mapplot()
      })
      
    }
    if(input$checkbox_gene_cluster == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Gene Clusters",align = "center"), div(dataTableOutput("cluster_gene"), style = "font-size:80%"),br(),br()))
        )
      )
      output$cluster_gene <- DT::renderDataTable({
        clusternum <- input$display_cluster
        gl <- plotHeatmap()[[3]] # getCluster()
        if (!is.null(gl)) {
          if (clusternum == "ALL") {
            gl
          } else {
            clusternum <- as.numeric(clusternum)
            dplyr::filter(gl, cluster == clusternum)
          }
        }
        
      })
    }
  })
  observeEvent(input$add_noise, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("Noise",align = "center"), plotlyOutput("noise_p", height = 500),br(),br()))
      )
    )
    output$noise_p <- renderPlotly({
      noisePlot()
    })
  })
  observeEvent(input$add_entropy, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("Shannon Entropy",align = "center"), plotlyOutput("entropy_plot", height = 500),br(),br()))
      )
    )
    output$entropy_plot <- renderPlotly({
      entropyPlot()
    })
  })
  observeEvent(input$add_tsne,{
    hide("error_text_report")
    if(input$checkbox_tsne_plot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("t-SNE Plot",align = "center"), plotlyOutput("tsne_p", height = 500),br(),br()))
        )
      )
      output$tsne_p <- renderPlotly({
        li <- tsne2plot()
        p <- li[[1]]
        tsne_table <- li[[2]]
        p
      })
      
    }
    if(input$checkbox_tsne_table == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("t-SNE Table",align = "center"), div(tableOutput("tsne_data"), style = "font-size:80%"),br(),br()))
        )
      )
      output$tsne_data <- renderTable({
        tsne_table <- tsne2plot()[[2]] # get table
        
        tsne_table
      })
      
    }
    
  })
  observeEvent(input$add_rf, {
    hide("error_text_report")
    if(input$checkbox_rf_plot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("RF Plot",align = "center"), plotlyOutput("rf", height = 500),br(),br()))
        )
      )
      output$rf <- renderPlotly({
        rfplot()
      })
      
    }
    if(input$checkbox_rafsil_plot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("RAFSIL Plot",align = "center"), plotlyOutput("raf_sil", height = 500),br(),br()))
        )
      )
      output$raf_sil <- renderPlotly({
        rafsilplot()
      })
      
      
    }
    if(input$checkbox_rf_matrix == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("RF Matrix",align = "center"), div(tableOutput("rf_mat"), style = "font-size:80%"),br(),br()))
        )
      )
      output$rf_mat <- renderTable({
        rf_matrix()
      })
      
    }
  })
  observeEvent(input$add_som, {
    hide("error_text_report")
    if(input$checkbox_property == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("SOM Analysis - ( Property Plot )",align = "center"), plotOutput("som_prop", height = 500),br(),br()))
        )
      )
      output$som_prop <- renderPlot({
        sompropertyplot()
      })
      
    }
    if(input$checkbox_count == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("SOM Analysis - ( Count Plot )",align = "center"), plotOutput("som_co", height = 500),br(),br()))
        )
      )
      output$som_co <- renderPlot({
        somcountplot()
      })
      
    }
    if(input$checkbox_codes == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("SOM Analysis - ( Codes Plot )",align = "center"), plotOutput("som_cod", height = 500),br(),br()))
        )
      )
      output$som_cod <- renderPlot({
        somcodesplot()
      })
      
    }
    if(input$checkbox_distance == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("SOM Analysis - ( Distance Plot )",align = "center"), plotOutput("som_dis", height = 500),br(),br()))
        )
      )
      output$som_dis <- renderPlot({
        somdistplot()
      })
      
    }
    if(input$checkbox_cluster == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("SOM Analysis - ( Cluster Plot )",align = "center"), plotOutput("som_clus", height = 500),br(),br()))
        )
      )
      output$som_clus <- renderPlot({
        somclusterplot()
      })
      
    }
    
  })
  observeEvent(input$add_gene_path_enrich,{
    hide("error_text_report")
    if(input$checkbox_plot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Pathways Enrichment (Plot)",align = "center"), plotlyOutput("path_enr_plot_gene", height=500),br(),br()))
        )
      )
      output$path_enr_plot_gene <- renderPlotly({
        df_path_enri_id_gene()
        gene_name <- as.data.frame(df_path_enri_id_gene())
        gene_name[,1] <- as.character(gene_name[,1])
        
        ggplotly(Pathway.Enr(gene_name[,1]), tooltip = c("text"))
      })
      
    }
    if(input$checkbox_visualization == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Pathways Enrichment (Visualization)",align = "center"), cyjShinyOutput('path_enri_vis', height=350),br(),br()))
        )
      )
      output$path_enri_vis <- renderCyjShiny({
        
        print("visualization")
        df_path_enri_id_gene()
        Enrich <- gost(df_path_enri_id_gene(),evcodes = T, sources = c('KEGG', 'REAC'))
        Pathway <- Construct.COPathway(Enrich, input$overlap_min_gene)
        nodes_tot <- c(unique(Pathway[,1],unique(Pathway[,2])))
        
        
        path_enri.nodes <- data.frame(id=nodes_tot,
                                      type=nodes_tot,
                                      stringsAsFactors=FALSE)
        
        path_enri.edges <- data.frame(source=Pathway[,1],
                                      target=Pathway[,2],
                                      interaction=Pathway[,1],
                                      stringsAsFactors=FALSE)
        
        graph.json <- dataFramesToJSON(path_enri.edges, path_enri.nodes)
        cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
        
      })
      
      
    }
    
  })
  observeEvent(input$add_tissue_exp, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("Scatter Plot",align = "center"),div(tableOutput("tissue_expr_table"), style = "font-size:80%"),br(),br()))
      )
    )
    output$tissue_expr_table <- renderTable({
      df_expr_table()
    })
  })
  observeEvent(input$add_gene_onto, {
    hide("error_text_report")
    if(input$checkbox_bio_proc_plot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Gene Ontlogy (Biological Process Plot)",align = "center"), plotOutput("gene_bio_proc", height = 500),br(),br()))
        )
      )
      output$gene_bio_proc <- renderPlot({
        GO_df <- plotUniprot()
        PlotGOBiological(GO_df,20)
      })
      
    }
    if(input$checkbox_bio_proc_table == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Gene Ontology (Biological Process Table)",align = "center"), div(tableOutput("bio_proc_table"), style = "font-size:80%"),br(),br()))
        )
      )
      output$bio_proc_table <- renderTable({
        GO_df <- plotUniprot()
        BiologicalDF <- Goparse(GO_df, 3)
        BiologicalDF <- na.omit(BiologicalDF)
        download_bio_table <- BiologicalDF
        BiologicalDF
      })
      
    }
    if(input$checkbox_mol_func_plot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Gene Ontlogy (Molecular Function Plot)",align = "center"), plotOutput("gene_mol_fun", height = 500),br(),br()))
        )
      )
      output$gene_mol_fun <- renderPlot({
        GO_df <- plotUniprot()
        Plot.GOMolecular(GO_df,20)
      })
      
    }
    if(input$checkbox_mol_func_table == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Gene Ontology (Molecular Function Table)",align = "center"), div(tableOutput("mol_func_table"), style = "font-size:80%"),br(),br()))
        )
      )
      output$mol_func_table <- renderTable({
        GO_df <- plotUniprot()
        MolecularDF <- Goparse(GO_df, 4)
        MolecularDF <- na.omit(MolecularDF)
        download_mol_table <- MolecularDF
        MolecularDF
      })
      
    }
    if(input$checkbox_cell_comp_plot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Gene Ontlogy (Cellular Component Plot)",align = "center"), plotOutput("gene_cell_comp", height = 500),br(),br()))
        )
      )
      output$gene_cell_comp <- renderPlot({
        GO_df <- plotUniprot()
        Plot.GOSubCellular(GO_df,20)
      })
      
    }
    if(input$checkbox_cell_comp_table == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h4("Gene Ontology (Cellular Component Table)",align = "center"), div(tableOutput("cell_comp_table"), style = "font-size:80%")))
        )
      )
      output$cell_comp_table <- renderTable({
        GO_df <- plotUniprot()
        CellularDF <- Goparse(GO_df, 5)
        CellularDF <- na.omit(CellularDF)
        download_cel_table <- CellularDF
        CellularDF
      })
      
    }
  })
  observeEvent(input$add_pp_inter,{
    hide("error_text_report")
    if(input$checkbox_pp_visu == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein-Protein Interactions (Visualization)",align = "center"), cyjShinyOutput("pp_visu", height=350),br(),br()))
        )
      )
      output$pp_visu <- renderCyjShiny({
        
        print(" renderCyjShiny invoked")
        print("graph.json:")
        
        
        print("running...")
        
        
        # tryCatch({
        
        Accessions <- df_prot_int_id()
        print("Please Wait... Fetching interaction data. It may take a while")
        protein_interaction_df <- getInteraction(Accessions)
        df_interaction(protein_interaction_df)
        print("Fetched...")
        
        #migrating rowId to first colunm 
        # protein_interaction_df <- cbind(ID = rownames(protein_interaction_df),protein_interaction_df)
        # rownames(protein_interaction_df) <- 1:nrow(protein_interaction_df)
        
        #making nodes
        nodes <- as.character(protein_interaction_df[,1])
        for (i in 1:nrow(protein_interaction_df))
        {
          if(!(is.na(protein_interaction_df[i,2])))
          {
            data_df <- strsplit(as.character(protein_interaction_df[i,2]),"; ")
            for(j in data_df)
            {
              nodes <- c(nodes,j)
            }
          }
        }
        
        print(nodes)
        
        print("Please Wait... Fetching Gene Names. It may take a while")
        protein_gene_name <- getGeneNames(nodes)
        df_names(protein_gene_name)
        print("........................")
        print(as.character(protein_gene_name[,1]))
        print("Fetched...")
        edge_source <- character()
        edge_target <- character()
        
        for (i in 1:nrow(protein_interaction_df))
        {
          if(!(is.na(protein_interaction_df[i,2])))
          {
            data_df <- strsplit(as.character(protein_interaction_df[i,2]),"; ")
            for(j in data_df)
            {
              edge_source <- c(edge_source,rep(as.character(protein_gene_name[as.character(protein_interaction_df[i,1]),1]),length(j)))
              print(as.character(protein_gene_name[j,1]))
              edge_target <- c(edge_target,as.character(protein_gene_name[j,1]))
            }
          }
        }
        
        tbl.nodes <- data.frame(id=as.character(protein_gene_name[,1]),
                                type=as.character(protein_gene_name[,1]),
                                stringsAsFactors=FALSE)
        
        
        tbl.edges <- data.frame(source=edge_source,
                                target=edge_target,
                                interaction=edge_target,
                                stringsAsFactors=FALSE)
        
        # }, error = function(error_condition) {
        #   print("using defauslt value")
        # })
        
        graph.json <- dataFramesToJSON(tbl.edges, tbl.nodes)
        
        print(fromJSON(graph.json))
        cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
        
        
      })
      
      
    }
    if(input$checkbox_pp_interact == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein-Protein Interactions(Protein Interactions) ",align = "center"), div(tableOutput("pp_inter_table"), style = "font-size:80%"),br(),br()))
        )
      )
      output$pp_inter_table <- renderTable({
        protein_interaction_df <- df_interaction()
        protein_gene_name <- df_names()
        print(protein_interaction_df)
        print("here")
        print(class(protein_interaction_df))
        if(df_names() == 0)
        {
          
          p_int_formatted <- data.frame()
          
        } else {
          
          protein_interaction_df[,1] <- as.character(protein_interaction_df[,1])
          
          p_int_formatted <- data.frame()
          count = 0
          n = 1
          for ( id in protein_interaction_df[,1])
          {
            count = count + 1
            if(!is.null(protein_interaction_df[,2]))
            {
              a = strsplit(as.character(protein_interaction_df[,2]),"; ")
              
              for(int_with in a[[count]])
              {
                p_int_row <- data.frame(id = as.character(paste0(as.character(lookup(id, as.data.frame(id_to_name), missing="Not found"))," ( ", id," )")),
                                        Interacts_With = as.character(paste0(as.character(lookup(int_with, as.data.frame(id_to_name), missing="Not found"))," ( ", int_with," )")),
                                        row.names = n)
                p_int_formatted <- rbind(p_int_formatted,p_int_row)
                n = n + 1
              }
            }
          }
          
          # for(i in 1:nrow(protein_interaction_df))
          # {
          #     protein_interaction_df[i,1] <- paste0(protein_interaction_df[i,1],
          #                               ' (',
          #                               protein_gene_name[protein_interaction_df[i,1],1],
          #                               ')')
          # }
          # print(protein_interaction_df)
          # colnames(protein_interaction_df)[2] <- "Interacts With"
          
        }
        
        p_int_formatted
      })
      
    }
    if(input$checkbox_prot_name == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein-Protein Interactions(Protein Names)",align = "center"), div(tableOutput("pp_name_table"), style = "font-size:80%"),br(),br()))
        )
      )
      output$pp_name_table <- renderTable({
        protein_gene_name <- df_names()
        if(protein_gene_name == 0)
        {
          protein_gene_name <- data.frame()
        } else {
          
          protein_gene_name <- cbind(ID = rownames(protein_gene_name),protein_gene_name)
          rownames(protein_gene_name) <- 1:nrow(protein_gene_name)
          colnames(protein_gene_name)[2] <- "Names"
          
        } 
        protein_gene_name
      })
      
    }
    
  })
  observeEvent(input$add_prot_func, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("Protein Function",align = "center"),div(tableOutput("prot_func_table"), style = "font-size:80%"),br(),br()))
      )
    )
    output$prot_func_table <- renderTable({
      df_func_table()
    })
  })
  observeEvent(input$add_subcell_loc, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("SubCellular Localization",align = "center"),div(tableOutput("sub_cell_table"), style = "font-size:80%"),br(),br()))
      )
    )
    output$sub_cell_table <- renderTable({
      df_local_table()
    })
  })
  observeEvent(input$add_prot_dom, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("Protein Domain",align = "center"),div(tableOutput("prot_dom_table"), style = "font-size:80%"),br(),br()))
      )
    )
    output$prot_dom_table <- renderTable({
      df_domain_table()
    })
  })
  observeEvent(input$add_prot_seq, {
    hide("error_text_report")
    if(input$checkbox_seq_charge == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein Sequences (Sequence Charge)",align = "center"), plotOutput("plot_seq_char", height = 500),br(),br()))
        )
      )
      output$plot_seq_char <- renderPlot({
        if (!is.null(df_prot_seq()))
        {
          hide("help_text_prot_seq")
          if(is.null(Seqdata))
          {
            Proteins <- df_prot_seq()
            Seqdata <<- GetSequences(Proteins)
          }
          PlotCharge(Seqdata)
        }
      })
      
    }
    if(input$checkbox_seq_acid == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein Sequences (Sequence Acidity)",align = "center"), plotOutput("plot_seq_acid", height = 500),br(),br()))
        )
      )
      output$plot_seq_acid <- renderPlot({
        if (!is.null(df_prot_seq()))
        {
          hide("help_text_prot_seq")
          if(is.null(Seqdata))
          {
            Proteins <- df_prot_seq()
            Seqdata <<- GetSequences(Proteins)
          }
          PlotAcidity(Seqdata)
        }
      })
      
    }
    if(input$checkbox_seq_grav_ind == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein Sequences (Sequence Gravy Index)",align = "center"), plotOutput("plot_seq_grav", height = 500),br(),br()))
        )
      )
      output$plot_seq_grav <- renderPlot({
        if (!is.null(df_prot_seq()))
        {
          hide("help_text_prot_seq")
          if(is.null(Seqdata))
          {
            Proteins <- df_prot_seq()
            Seqdata <<- GetSequences(Proteins)
          }
          PlotGravy(Seqdata)
        }
      })
      
    }
    if(input$checkbox_physio_prop == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein Sequences (All Physiochemical Properties)",align = "center"), plotOutput("plot_physio_prop", height = 500),br(),br()))
        )
      )
      output$plot_physio_prop <- renderPlot({
        if (!is.null(df_prot_seq()))
        {
          hide("help_text_prot_seq")
          if(is.null(Seqdata))
          {
            Proteins <- df_prot_seq()
            Seqdata <<- GetSequences(Proteins)
          }
          PlotPhysicochemical(Seqdata)
        }
      })
      
    }
  })
  observeEvent(input$add_prot_evol_analysis, {
    hide("error_text_report")
    if(input$checkbox_prot_gene == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein Evolutionary analysis (Protein's Gene Trees)",align = "center"), radialNetworkOutput("prot_gene_plot", width="500px",height="500px"),br(),br()))
        )
      )
      output$prot_gene_plot <- renderRadialNetwork(
        {
          if (!is.null(df_prot_seq_evol()))
          {
            if (is.null(GenesObj))
            {
              Proteins <- df_prot_seq_evol()
              GenesObj <- GetNamesTaxa(Proteins)
            }
            ConstructGenes(GenesObj)
          }
        }
      )
      
      
    }
    if(input$checkbox_prot_chrom == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein Evolutionary analysis (Protein's chromosomal location)",align = "center"), plotOutput("plot_prot_chrom", height = 500),br(),br()))
        )
      )
      output$plot_prot_chrom <- renderPlot(
        if (!is.null(df_prot_seq_evol()))
        {
          if(is.null(GenesObj))
          {
            Proteins <- df_prot_seq_evol()
            GenesObj <- GetNamesTaxa(Proteins)
          }
          PlotChromosomeInfo(GenesObj)
        }
      )
      
    }
    if(input$checkbox_prot_evol == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein Evolutionary analysis (Evolutionary analysis)",align = "center"), plotOutput("plot_evol_ana", height = 500),br(),br()))
        )
      )
      output$plot_evol_ana <- renderPlot(
        {
          if(!is.null(df_prot_seq_evol()))
            if(is.null(Seqdata))
            {
              Proteins <- df_prot_seq_evol()
              Seqdata <<- GetSequences(Proteins)
            }
          ConstructPhylogeny(Seqdata)
        }
      )
      
      
    }
    
  })
  observeEvent(input$add_prot_path_analysis, {
    hide("error_text_report")
    
    if(input$checkbox_dise_role == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein pathological analysis (Protein's disease role)",align = "center"), div(tableOutput("prot_dis_table"), style = "font-size:80%"),br(),br()))
        )
      )
      output$prot_dis_table <- renderTable({
        if(!is.null(df_prot_seq_Patho()))
        {
          Proteins <- df_prot_seq_Patho()
          Pathodata <- GetPathology_Biotech(Proteins)
          DiseaseTable <- Get.diseases(Pathodata) 
        }
      }, escape = F)
      
    }
    if(input$checkbox_dise_dist == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Protein pathological analysis (Protein's disease distribution)",align = "center"), bubblesOutput("prot_dis_plot"),br(),br()))
        )
      )
      output$prot_dis_plot <- renderBubbles({
        if(!is.null(df_prot_seq_Patho()))
        {
          if(!is.null(DiseaseTable))
          {
            Plot.NDiseases(DiseaseTable)
          }
          else {
            Proteins <- df_prot_seq_Patho()
            Pathodata <- GetPathology_Biotech(Proteins)
            DiseaseTable <- Get.diseases(Pathodata)
            Plot.NDiseases(DiseaseTable)
          }
        }
      })
      
    }
    
  })
  observeEvent(input$add_prot_comp_enrich, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        
        fluidRow(
          column(width=2),
          column(width= 8,h3("Complex Enrichment",align = "center"),div(tableOutput("comp_enrich_table"), style = "font-size:80%"),br(),br()))
      )
    )
    output$comp_enrich_table <- renderTable({
      df_com_table()
    })
  })
  observeEvent(input$add_prot_path_enrich,{
    hide("error_text_report")
    if(input$checkbox_plot_prot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Pathways Enrichment (Plot)",align = "center"), plotlyOutput("path_enr_plot_prot", height=500),br(),br()))
        )
      )
      output$path_enr_plot_prot <- renderPlotly({
        df_path_enri_id_prot()
        gene_name <- as.data.frame(df_path_enri_id_prot())
        gene_name[,1] <- as.character(gene_name[,1])
        
        ggplotly(Pathway.Enr(gene_name[,1]), tooltip = c("text"))
      })
      
    }
    if(input$checkbox_visualization_prot == TRUE){
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          
          fluidRow(
            column(width=2),
            column(width= 8,h3("Pathways Enrichment (Visualization)",align = "center"), cyjShinyOutput('path_enri_vis_prot', height=350),br(),br()))
        )
      )
      
      output$path_enri_vis_prot <- renderCyjShiny({
        
        print("visualization")
        df_path_enri_id_prot()
        
        print(df_path_enri_id_prot()[,1])
        
        Enrich <- gost(df_path_enri_id_prot()[,1],evcodes = T, sources = c('KEGG', 'REAC'))
        
        Pathway <- Construct.COPathway(Enrich, input$overlap_min_prot)
        
        nodes_tot <- c(unique(Pathway[,1],unique(Pathway[,2])))
        
        
        path_enri.nodes <- data.frame(id=nodes_tot,
                                      type=nodes_tot,
                                      stringsAsFactors=FALSE)
        
        path_enri.edges <- data.frame(source=Pathway[,1],
                                      target=Pathway[,2],
                                      interaction=Pathway[,1],
                                      stringsAsFactors=FALSE)
        
        graph.json <- dataFramesToJSON(path_enri.edges, path_enri.nodes)
        cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
        
      })
      
    }
    
  })
  
  ########################################
  ##### Increases the Upload Limit #######
  ########################################
  
  options(shiny.maxRequestSize=25000*1024^2)
  
  ########################################
  ##### get variable names for input #####
  ########################################
  
  observe({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_norm()
    } else if (type == "raw") {
      DS <- df_raw()
    }
    
    nms <- colnames(DS)
    updateSelectInput(session, "scatter.x", choices = nms, selected = nms[1])
    updateSelectInput(session, "scatter.y", choices = nms, selected = nms[2])
    updateSelectInput(session, "dist.var", choices = nms)
    col_num <- ncol(DS)
    updateSliderInput(session, "pca_cluster_num", max = col_num - 1)
    genotype_num <- NULL
    if (is.null(DS) == FALSE) {
      for (i in 2:col_num) {
        if (col_num %% i == 0) {
          genotype_num <- c(genotype_num, i)
        }
      }
    }
    updateSelectInput(session, "numOfGeno", choices = genotype_num)
    updateSelectInput(session, "noise_anchor_c", choices = nms)
    
    ##############################################
    ##############################################
    updateSelectInput(session, "som_samples", choices = c("All", nms), selected = "All")
    updateSelectInput(session, "overlay.x1", choices = nms, selected = nms[1])
    updateSelectInput(session, "overlay.y1", choices = nms, selected = nms[2])
    updateSelectInput(session, "overlay.x2", choices = nms, selected = nms[3])
    updateSelectInput(session, "overlay.y2", choices = nms, selected = nms[4])
    ##############################################
    ##############################################
    
    observeEvent(input$start_rnaseq, {
      updateNavbarPage(session, inputId = "navbar", selected = "active_tab_rnaseq")
    })
    
    observeEvent(input$start_micro, {
      updateNavbarPage(session, inputId = "navbar", selected = "active_tab_micro")
    })
    
    ### preprocessing tab
    f <- group_names()
    f <- unique(as.character(f))
    if (is.null(f)) {
      hideTab(inputId = "preprocessing_tabs", target = "Description table")
      # hideTab(inputId="preprocessing_tabs", target="Description table")
    } else {
      showTab(inputId = "preprocessing_tabs", target = "Description table")
      updateSelectInput(session, "f1", choices = f, selected = f[1])
      updateSelectInput(session, "f2", choices = f, selected = f[2])
    }
    
    ### gene expression range for distribution fit ###
    if (is.null(DS) == FALSE) {
      DS_dist <- distfit_df()
      print(DS_dist)
      range_min <- min(DS_dist)
      range_max <- max(DS_dist)
      updateSliderInput(session, "dist_range", max = round(range_max), value = c(0.1, range_max))
      updateNumericInput(session, "dist_range_min", min = 0.000001, max = round(range_max), value = 0.1)
      updateNumericInput(session, "dist_range_max", min = 0.000001, max = round(range_max), value = round(range_max))
    }
    
    ### gene sample size choices for PCA ###
    # print("line 647 check input$submit_preprocessing")
    # v=input$submit_preprocessing
    if (input$submit_preprocessing > 0) {
      if (type == "norm") {
        DS_filt <- df_shiny()
      } else if (type == "raw") {
        DS_filt <- df_raw_shiny()
      }
    } else {
      DS_filt <- DS
    }
    
    i <- 1
    min_size <- 25
    samplesize <- NULL
    while (i * min_size < length(DS_filt[, 1])) {
      samplesize <- c(samplesize, i * min_size)
      i <- i * 2
    }
    if (is.null(samplesize)) {
      samplesize <- c(samplesize, length(DS_filt[, 1]))
    } else if (samplesize[length(samplesize)] != length(DS_filt[, 1])) {
      samplesize <- c(samplesize, length(DS_filt[, 1]))
    }
    updateSelectInput(session, "gene_size", choices = samplesize, selected = samplesize[length(samplesize)])
    
    ### pca choices for PCA-2D ###
    pcchoices <- NULL
    if (is.null(DS) == FALSE) {
      for (i in 1:ncol(DS)) {
        pcchoices <- c(pcchoices, paste("PC", i, sep = ""))
      }
    }
    updateSelectInput(session, "pca.x", choices = pcchoices, selected = pcchoices[1])
    updateSelectInput(session, "pca.y", choices = pcchoices, selected = pcchoices[2])
    
  })
  
  
  observeEvent(input$submit_input, {
    type <- input$file_type
    if (type == "norm") {
      DS <- df_norm()
      lengths <- 0
    } else if (type == "raw") {
      DS <- df_raw()
      lengths <- gene_length()
      # if( length(intersect(rownames(lengths), rownames(DS))) < 1000 )
      #   length <- NULL
    }
    f <- group_names()
    spikes <- neg_control()
    
    # if any NULL value, throw error. TO CHANGE TO BE MORE SPECIFIC
    input.list <- list(DS, f)
    input.null <- sapply(input.list, is.null)
    names(input.null) <- c("Expression/Counts", "Meta Data")
    
    if (any(input.null)) {
      index.null <- which(input.null)
      errors <- paste(names(input.null)[index.null], collapse = ", ")
      # print(errors)
      showModal(modalDialog(
        type = "Error",
        paste("Please check these input:", errors, "and try again!")
      ))
    } else {
      updateTabsetPanel(session, inputId = "Rnaseq_pre", selected = "Preprocessing")
    }
    
    # update input
    updateNumericInput(session, "min_col", max = ncol(DS)) # update max column nunmber in filtering
    if (is.null(spikes)) {
      updateRadioButtons(session, "norm_method", choices = c(
        "None (Black)" = "None",
        "RPKM (Blue)" = "RPKM", "FPKM (Dark cyan)" = "FPKM",
        "TPM (Dark green)" = "TPM",
        "Upper Quartile (Brown)" = "RUV"
      ))
      # c("None",'RPKM','FPKM','TPM',"Upper Quartile"="RUV")
    } else {
      updateRadioButtons(session, "norm_method", choices = c(
        "None (Black)" = "None",
        "RPKM (Blue)" = "RPKM", "FPKM (Dark cyan)" = "FPKM",
        "TPM (Dark green)" = "TPM",
        "RUV (Brown)" = "RUV"
      ))
    }
    if (is.null(lengths) & !(is.null(spikes))) {
      updateRadioButtons(session, "norm_method", choices = c("None (Black)" = "None", "RUV (Brown)" = "RUV"))
    } else if (is.null(lengths) & (is.null(spikes))) {
      updateRadioButtons(session, "norm_method", choices = c("None (Black)" = "None", "Upper Quartile (Brown)" = "RUV"))
    }
    
    if (is.null(f)) {
      hideTab(inputId = "navbar", target = "DE Analysis")
    } else {
      showTab(inputId = "navbar", target = "DE Analysis")
    }
    # if(is.null(f)){
    #   hideTab(inputId="preprocessing_tabs", target="Description table")
    # } else {
    #   showTab(inputId="preprocessing_tabs", target="Description table")
    # }
    
    hide("help_text_scatter")
    hide("help_text_dis_fit")
    hide("help_text_correlation")
    hide("help_text_PCA")
    hide("help_text_DE_anal")
    hide("help_text_heatmap")
    hide("help_text_Noise")
    hide("help_text_Entropy")
    # hide("help_text_SVM")
    hide("help_text_tsne")
    hide("help_text_rf")
    hide("help_text_SOM")
    
  })
  
  
  ############################ Micro array ###############################
  
  
  output$downloadMicroRaw <- downloadHandler(
    filename = function() {
      paste("Microarray_Raw", ".csv", sep = "")
    },
    content = function(file) {
      raw_data <- df_micro()
      
      micro_raw_data <- as.data.frame(raw_data@assayData[["exprs"]][1:2000,])
      micro_raw_data <- tibble::rowid_to_column(micro_raw_data, "Geneid")
      
      n <- ncol(micro_raw_data)
      for(i in 2:n)
      {
        colnames(micro_raw_data)[i] <- paste0('new_',colnames(micro_raw_data)[i])
      }
      
      write.csv(micro_raw_data, file, row.names = FALSE)
    }
  )
  
  output$downloadMicroMeta <- downloadHandler(
    filename = function() {
      paste("Microarray_Meta", ".csv", sep = "")
    },
    content = function(file) {
      raw_data <- df_micro()
      
      micro_raw_data <- as.data.frame(raw_data@assayData[["exprs"]][1:2000,])
      micro_meta_data <- as.data.frame(raw_data@phenoData@data[["Factor.Value.disease."]])
      micro_meta_data <- add_column(micro_meta_data, as.data.frame(colnames(micro_raw_data)), .after = 0)
      
      names(micro_meta_data)[1] <- "Id"
      names(micro_meta_data)[2] <- "Types"
      
      n <- nrow(micro_meta_data)
      micro_meta_data[,1] <- as.character(micro_meta_data[,1])
      for(i in 1:n)
      {
        micro_meta_data[i,1] <- paste0("new_",micro_meta_data[i,1])
      }
      
      write.csv(micro_meta_data, file, row.names = FALSE)
    }
  )
  
  
  ####################################################################
  
  # observeEvent(input$submit_preprocessing, {
  #   type <- input$file_type
  #   if(type=='norm'){
  #     DS <- df_shiny()
  #   }else if(type=='raw'){
  #     DS <- df_raw_shiny()
  #   }
  #   ### gene sample size choices for PCA ###
  #   i <- 1
  #   min_size <- 25
  #   samplesize <- NULL
  #   while(i*min_size<length(DS[,1])){
  #     samplesize <- c(samplesize,i*min_size)
  #     i <- i*2
  #   }
  #   if(is.null(samplesize)){
  #     samplesize <- c(samplesize,length(DS[,1]))
  #   }else if(samplesize[length(samplesize)]!=length(DS[,1])){
  #     samplesize <- c(samplesize,length(DS[,1]))
  #   }
  #   updateSelectInput(session,"gene_size", choices = samplesize,selected = samplesize[length(samplesize)])
  # })
  
  ######################################
  ######### read in / get data #########
  ######################################
  
  #####################
  ## get data #########
  #####################
  
  # get normalized counts
  df_norm <- reactive({ # get normalized counts
    if (is.null(input$file2)) {
      return(NULL)
    }
    parts <- strsplit(input$file2$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    type <- tolower(type)
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    ds <- read.csv(input$file2$datapath)
    ds <- na.omit(ds)
    ds <- ds[!duplicated(ds[, 1]), ] # remove duplicated gene names
    
    row_names <- ds[, 1]
    DS <- data.frame(ds)
    if (ncol(DS) <= 1) {
      showModal(modalDialog(
        title = "Error",
        "Please check normalised data file format (Eg_normalised.png) and try again!"
      ))
      return(NULL)
    }
    DS <- DS[, -1]
    row.names(DS) <- row_names
    for (i in 1:ncol(DS)) {
      if (class(DS[, i]) != "numeric" & class(DS[, i]) != "integer") {
        showModal(modalDialog(
          title = "Error",
          "Please check normalised data file format (Eg_normalised.png) and try again!"
        ))
        return(NULL)
      }
    }
    return(DS)
  })
  
  # get raw counts
  df_raw <- reactive({
    print("in df_raw")
    print(value_var$geo_file_type)
    if (is.null(input$file1) && value_var$geo_file_type != "rnaseq") {
      return(NULL)
    }
    else if (!is.null(input$file1)){
      parts <- strsplit(input$file1$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      raw_ds <- read.csv(input$file1$datapath)
      
      
    }
    else if(value_var$geo_file_type == "rnaseq"){
      
      if(file.exists(file.path(getwd(),input$file_name_button))){
        print("Reading geo_file")
        
        raw_ds <- read.table(file.path(getwd(),input$file_name_button) ,header=TRUE, stringsAsFactors = FALSE)
        print(raw_ds)
        #value_var$geo_file_type<-"none"
        
        
      }else{
        return(NULL)
      }
      
      
    } # remove duplicated gene names
    else{
      return(NULL)
    }
    raw_ds <- na.omit(raw_ds)
    raw_ds <- raw_ds[!duplicated(raw_ds[, 1]), ]  
    # raw_ds <- as.data.frame(raw_ds)
    if (ncol(raw_ds) <= 1) {
      showModal(modalDialog(
        title = "Error",
        "Data file must contain at least 2 columns. Please check raw data format and try again!"
      ))
      return(NULL)
    }
    
    row_names <- raw_ds[, 1]
    rownames(raw_ds) <- row_names
    raw_DS <- raw_ds[, -1] # remove the first column, which is gene Id
    for (i in 1:ncol(raw_DS)) {
      
      if (class(raw_DS[, i]) != "numeric" & class(raw_DS[, i]) != "integer") {
        showModal(modalDialog(
          title = "Error",
          "Raw counts must be integer. Please check raw data formate and try again!"
        ))
        return(NULL)
      }
    }
    print(raw_DS)
    return(raw_DS)
  })
  
  ######################## Micro array ##########################
  
  
  df_micro <- reactive({
    print("running")
    print("In Micro array")
    print(value_var$geo_file_type)
    if (is.null(input$file_micro)&& value_var$geo_file_type != "microarray") {
      return(NULL)
    }
    else if(!is.null(input$file_micro)){
      parts <- strsplit(input$file_micro$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      
      if (type != "zip") {
        showModal(modalDialog(
          title = "Error",
          "Please input a zip file!"
        ))
        return(NULL)
      }
      
      unzip(input$file_micro$datapath,exdir = parts[[1]][1])
    }else if( value_var$geo_file_type == "microarray"){
      parts <- strsplit(input$file_name_button, ".", fixed = TRUE)
      untar(file.path(getwd(),input$file_name_button),exdir = parts[[1]][1])
      
    }else {
      return(NULL)
    }
    fol_name <- print(list.files(parts[[1]][1]))
    micro_data_dir <- paste0(parts[[1]][1],"/",fol_name)
    print(micro_data_dir)
    
    sdrf_location <- file.path(micro_data_dir, "E-MTAB-2967.sdrf.txt")
    SDRF <- read.delim(sdrf_location)
    
    rownames(SDRF) <- SDRF$Array.Data.File
    SDRF <- AnnotatedDataFrame(SDRF)
    
    raw_data <- oligo::read.celfiles(filenames = file.path(micro_data_dir, 
                                                           SDRF$Array.Data.File),
                                     verbose = FALSE, phenoData = SDRF)
    
    print(raw_data)
    
    return(raw_data)
    
  })
  
  
  
  ############################## GEO IMPORT###################
  
  output$help_text_geo <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Please enter the GEO accession number to begin analysis.
          </b>
        </p>
      </center>
    ")
  })
  getDirListing <- function(url) {
    # Takes a URL and returns a character vector of filenames
    a <- xml2::read_html(url)
    fnames = grep('^G',xml_text(xml_find_all(a,'//a/@href')),value=TRUE)
    return(fnames)
  }
  
  getFileUrl <- function(GEO,filetype){
    geotype <- toupper(substr(GEO,1,3))
    fileinfo <- list()
    stub = gsub('\\d{1,3}$','nnn',GEO,perl=TRUE)
    if(geotype=='GSM') {
      url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/%s/",stub,GEO,filetype)
    }
    if(geotype=='GSE') {
      url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/%s/",stub,GEO,filetype)
    }
    if(geotype=='GPL') {
      url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/%s/",stub,GEO,filetype)
    }
    return(url)
  }
  
  getFiles <- function(url) {
    
    fnames <- try(getDirListing(url),silent=TRUE)
    
    if(inherits(fnames,'try-error')) {
      message('No supplemental files found.')
      message('Check URL manually if in doubt')
      message(url)
      return(NULL)
    }
    return(fnames)
    
    
  }
  
  downloadFile <- function(url,fname){
    storedir <- getwd()
    #suppressWarnings(dir.create(storedir <- file.path(getwd(),GEO)))
    if(!file.exists(file.path(getwd(),input$file_name_button))){
      download.file(paste(file.path(url,fname),'tool=geoquery',sep="?"),
                    destfile=file.path(storedir,fname),
                    mode='wb',
                    method=getOption('download.file.method.GEOquery'))
    }
    #acc_data <- read.table(file.path(storedir,fname))
    #return(acc_data)
    
  }
  
  check_file_type <- function(){
    parts <- strsplit(input$file_name_button, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    print("before")
    print(value_var$geo_file_type)
    if (type == "gz"){
      
      value_var$geo_file_type<-"rnaseq"
      output$help_text_geo <- renderUI({
        HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          GEO Data Import Complete. Please go to the preprocessing tab of RnaSeq and proceed with the analysis.
          </b>
        </p>
      </center>
    ")
      })
    }
    else if(type == "tar"){
      file_list <- untar(file.path(getwd(),input$file_name_button),list=TRUE)
      print(file_list[1])
      parts <- strsplit(file_list[1], ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])-1]
      print(type)
      if(type == "cel" || type == "CEL"){
        value_var$geo_file_type<-"microarray"
        output$help_text_geo <- renderUI({
          HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Its Microarray file. Please go to the preprocessing tab of Microarray and proceed with the analysis.
          </b>
        </p>
      </center>
    ")
        })
      }
      #df_micro()
    }
    
    
  }
  
  observeEvent(input$submit_geo_acc_no, {
    
    url<-getFileUrl(input$geo_acc_no,"suppl")
    fname <- getFiles(url)
    updateRadioButtons(session, "file_name_button",
                       choices = fname,
                       selected = fname[1]
    )
    
    updateTabsetPanel(session, "geo_tab",
                      selected = "geo_pre")
    
  })
  
  observeEvent(input$submit_geo_preprocessing,{
    value_var$geo_file_type<-"none"
    url<-getFileUrl(input$geo_acc_no,"suppl")
    print(input$file_name_button)
    downloadFile(url,input$file_name_button)
    check_file_type()
    df_geo <- read.table(file.path(getwd(),input$file_name_button) ,header=TRUE, stringsAsFactors = TRUE)
    output$geo_box_plot <- renderPlot({# load series and platform data from GEO
      
      # box-and-whisker plot
      par(mar=c(7,4,2,1))
      title <- paste (input$geo_acc_no, "/",input$file_name_button , sep ="")
      boxplot(df_geo, boxwex=0.7, notch=FALSE, main=title, outline=FALSE, las=2)})
    
    output$geo_expr_plot <- renderPlot({
      
      par(mar=c(4,4,2,1))
      title <- paste (input$geo_acc_no, "/", input$file_name_button, " value distribution", sep ="")
      plotDensities(df_geo, main=title, legend=F)
      
    })
    output$geo_mean_plot <- renderPlot({
      
      # mean-variance trend
      ex <- na.omit(df_geo) # eliminate rows with NAs
      plotSA(lmFit(ex), main="Mean variance trend")
      
    })
    
    output$geo_umap_plot <- renderPlot({
      
      ex <- df_geo[!duplicated(df_geo), ]  # remove duplicates
      ump <- umap(t(ex), n_neighbors = 5, random_state = 123)
      plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", pch=20, cex=1.5)
      # point labels without overlaps
      pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
      
    })
  })
  
  
  ###############################################################
  
  # get gene length
  gene_length <- reactive({
    if (is.null(input$length1)&&is.null(input$length2)) {
      return(NULL)
    }
    else if (!is.null(input$length1)){
      lengths_df <- read.csv(input$length1$datapath)
    }else{
      lengths_df <- read.csv(input$length2$datapath)
    }
    lengths_df2 <- data.frame("len" = lengths_df[, 2])
    rownames(lengths_df2) <- as.character(lengths_df[, 1])
    return(lengths_df2)
  })
  
  # get spikes / negative control genes
  neg_control <- reactive({
    if (is.null(input$spikes1)&&is.null(input$spikes2)) {
      return(NULL)
    }else if (!is.null(input$spikes1)){
      spikes <- read.csv(input$spikes1$datapath, header = F)
    }else{
      spikes <- read.csv(input$spikes2$datapath, header = F)
    }
    spikes <- as.character(spikes[, 1])
    # print(spikes[1:10])
    return(spikes)
  })
  
  # get meta data table
  group_names <- reactive({
    # if no data
    if (is.null(input$metafile1)) {
      return(NULL)
    }
    
    # read in group names (metadata)
    groups <- read.csv(input$metafile1$datapath)
    group_colnames <- as.character(groups[, 1])
    
    type <- input$file_type
    if (type == "norm") {
      DS <- df_norm()
    } else if (type == "raw") {
      DS <- df_raw()
    }
    col_names <- colnames(DS) # columm names of DS in order
    
    # check if groups and column names are similar
    if (!all(col_names %in% group_colnames) || ncol(groups) < 2) {
      showNotification(type = "error", "group names and DS column names not similar")
      return(NULL)
    }
    
    if (ncol(groups) == 2) {
      f <- groups[match(col_names, groups[, 1]), ] [, 2] # arrange f in the same order as col_names
    } else {
      f <- groups[match(col_names, groups[, 1]), ] [, 2]
      for (i in 3:ncol(groups)) {
        f <- paste0(f, "_", groups[, i])
      }
    }
    f <- as.factor(make.names(f))
    # return(as.factor(f))
    return(f)
  })
  
  ### Gene ontology
  
  gene_list <- reactive({
    if (is.null(input$filego)) {
      return(NULL)
    }
    parts <- strsplit(input$filego$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    type <- tolower(type)
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    ds <- read.csv(input$filego$datapath, header = FALSE)
    if (ncol(ds) >= 2) {
      col1 <- ds[-1, 1]
    } else if (ncol(ds) == 1) {
      col1 <- ds[, 1]
    } else {
      showModal(modalDialog(
        title = "Error",
        "No data found! Please check required data format and try again!"
      ))
      return(NULL)
    }
    gene_list <- as.character(col1)
    print("gene list from gene_list")
    print(head(gene_list))
    return(gene_list)
  })
  
  bg_list <- reactive({
    if (is.null(input$filebg)) {
      return(NULL)
    }
    parts <- strsplit(input$filebg$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    type <- tolower(type)
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    ds <- read.csv(input$filebg$datapath, header = FALSE)
    if (ncol(ds) > 1) {
      col1 <- ds[-1, 1]
    } else if (ncol(ds) == 1) {
      col1 <- ds[, 1]
    } else {
      showModal(modalDialog(
        title = "Error",
        "No data found! Please check required data format and try again!"
      ))
      return(NULL)
    }
    bg_list <- as.character(col1)
    return(bg_list)
  })
  
  ####################################
  ########## PREPROCESSING ###########
  ####################################
  
  # filter normalized counts
  df_shiny <- eventReactive(input$submit_preprocessing, {
    DS_norm <- df_norm()
    min_val <- input$min_val
    min_col <- input$min_col
    keep <- rowSums(DS_norm >= min_val) >= min_col
    DS <- DS_norm[keep, ]
    # DS <- apply(DS_norm, 1, function(x) length(x[x>min_val])>=min_col)
    return(DS)
  })
  
  # filter raw counts
  df_raw_filt <- eventReactive(input$submit_preprocessing, {
    DS_raw <- df_raw()
    min_val <- input$min_val
    min_col <- input$min_col
    keep <- rowSums(DS_raw >= min_val) >= min_col
    DS_filt <- DS_raw[keep, ]
    # DS_filt <- apply(DS_raw, 1, function(x) length(x[x>min_val])>=min_col)
    return(DS_filt)
  })
  
  # normalizing raw counts
  df_raw_shiny <- reactive({
    raw_DS <- df_raw_filt() # get filtered raw counts
    method <- input$norm_method
    
    if (method %in% c("TPM", "RPKM", "FPKM")) {
      if (is.null(input$length1)&&is.null(input$length2)) {
        showModal(modalDialog(
          title = "Error",
          "Please Enter a Gene Length File first!"
        ))
      }
      
      lengths_df <- gene_length()
      merge_DS <- merge(raw_DS, lengths_df, by = "row.names")
      rownames(merge_DS) <- merge_DS[, 1]
      merge_DS <- merge_DS[, -1]
      raw_DS <- merge_DS[, -ncol(merge_DS)]
      lengths <- merge_DS[, ncol(merge_DS)]
      # print("length")
      # print(head(merge_DS))
    }
    # print("from line 981 df_raw_shiny")
    # print(method)
    # print("raw_DS")
    # print(head(raw_DS[,1:4]))
    # print("dimension of raw_DS")
    # print(dim(raw_DS))
    
    if (method == "TPM") {
      tpm.matrix <- apply(raw_DS, 2, function(x) tpm(x, lengths))
      tpm.df <- data.frame(tpm.matrix)
      return(tpm.df)
    } else if (method == "RPKM") {
      rpkm.matrix <- edgeR::rpkm(raw_DS, lengths)
      rpkm.df <- data.frame(rpkm.matrix)
      return(rpkm.df)
    } else if (method == "FPKM") {
      fpkm.matrix <- apply(raw_DS, 2, function(x) fpkm(x, lengths))
      fpkm.df <- data.frame(fpkm.matrix)
      return(fpkm.df)
    } else if (method == "None") {
      return(raw_DS)
    } else if (method == "RUV") {
      if (is.null(input$spikes1)&&is.null(input$spikes2)) {
        showModal(modalDialog(
          title = "Error",
          "Please Enter a Negative Control File first!"
        ))
      }
      spikes <- neg_control()
      if (!is.null(spikes)) {
        spikes <- intersect(spikes, rownames(raw_DS))
      }
      # f <- group_names()
      # if( is.null(spikes) )
      #   spikes <- getEmpirical(rawDS,f)
      set1 <- RUVg.apply(raw_DS, spikes)
      RUV.df <- as.data.frame(normCounts(set1))
      return(RUV.df)
    }
  })
  
  ### for distribution fitting
  distfit_df <- reactive({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    for (i in 1:ncol(DS)) {
      DS <- DS[which(DS[, i] > 0), ]
      DS <- na.omit(DS)
    }
    return(DS)
  })
  
  
  ######### ANALYSIS FROM HERE ############
  ######## RLEplot and Preprocessing ###########
  #############################################
  
  
  RLE.plot <- reactive({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    set1 <- newSeqExpressionSet(as.matrix(DS))
    norm_method_name <- input$norm_method
    colors <- c("RPKM" = "blue", "FPKM" = "darkcyan", "TPM" = "darkgreen", "RUV" = "Brown", "Upper Quartile" = "Brown")
    if (norm_method_name != "None" & input$submit_preprocessing != 0) {
      spikes <- neg_control()
      if (norm_method_name == "RUV" & is.null(spikes)) {
        norm_method_name <- "Upper Quartile"
      }
      par(mar = c(7, 4, 4, 4) + 1.2)
      plotRLE(set1,
              ylim = c(-1.5, 1.5), outline = FALSE, col = colors[norm_method_name],
              las = 2,
              hjust = 1,
              main = paste(norm_method_name, "Normalized")
      )
    }
  })
  
  
  violin_plot <- reactive({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    
    norm_method_name <- input$norm_method
    
    if (norm_method_name != "None" & input$submit_preprocessing != 0) {
      
      df <- as.data.frame(DS)
      df <- setNames(stack(df),c("norm_type","Genotype"))
      df$norm_type <- log(df$norm_type+1)
      
      p <- ggplot(df, aes(x=Genotype, y=norm_type)) + 
        geom_violin(trim=FALSE) + 
        scale_color_brewer(palette="Dark2") +
        labs(title=paste(norm_method_name,"Normalized",sep = ' '), y = paste("log(",norm_method_name,"+1)",sep = '') )+
        stat_summary(fun.data=mean_sdl, mult=1, 
                     geom="pointrange", color="red") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      p
      
      ggplotly(p, tooltip = c("text"))
      
    }
  })
  
  
  output$RLE.plot <- renderPlot({
    RLE.plot()
  })
  
  output$violin_plot <- renderPlotly({
    violin_plot()
  })
  
  output$RLE.plot2 <- renderPlot({ # for raw data
    start.rle <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      raw_DS <- df_shiny()
      main_title <- "Input data"
    } else if (type == "raw") {
      raw_DS <- df_raw()
      main_title <- "Raw data"
    }
    
    set1 <- newSeqExpressionSet(as.matrix(raw_DS))
    if (input$submit_preprocessing != 0) {
      par(mar = c(7, 4, 4, 4) + 1.2)
      plotRLE(set1, ylim = c(-1.5, 1.5), outline = FALSE, main = main_title, las = 2)
    }
    end.rle <- Sys.time()
    print("time for RLE plot and preprocessing")
    print(end.rle - start.rle)
  })
  
  violin_plot2 <- reactive({
    type <- input$file_type
    if (type == "norm") {
      raw_DS <- df_shiny()
      main_title <- "Input data"
    } else if (type == "raw") {
      raw_DS <- df_raw()
      main_title <- "Raw data"
    }
    
    df <- as.data.frame(raw_DS)
    df <- setNames(stack(df),c("norm_type","Genotype"))
    df$norm_type <- log(df$norm_type+1)
    
    p <- ggplot(df, aes(x=Genotype, y=norm_type)) + 
      geom_violin(trim=FALSE) + 
      scale_color_brewer(palette="Dark2") +
      labs(title="Raw Data", y = paste("log(Raw Data+1)",sep = '') ) +
      stat_summary(fun.data=mean_sdl, mult=1, 
                   geom="pointrange", color="red") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p
    
    ggplotly(p, tooltip = c("text"))
  })
  
  output$violin_plot2 <- renderPlotly({
    violin_plot2()
  })
  
  output$norm_table <- DT::renderDataTable({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    # if(input$submit_preprocessing != 0)
    DS # with filtering and normalization
  })
  
  output$meta_table <- DT::renderDataTable({
    f <- group_names()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (!is.null(f)) {
      meta_df <- data.frame("Column names" = colnames(DS), "Description" = f)
      meta_df
    }
  })
  
  output$download_norm_data <- downloadHandler(
    filename = function() {
      method <- input$norm_method
      paste(method, "normalized.csv")
    },
    content = function(file) {
      type <- input$file_type
      if (type == "norm") {
        DS <- df_shiny()
      } else if (type == "raw") {
        DS <- df_raw_shiny()
      }
      write.csv(DS, file, row.names = F)
    }
  )
  
  ############################
  ######## scatter ###########
  ############################
  
  # input scatter data
  plotScatter <- reactive({
    scatter.start <- Sys.time()
    trans <- input$trans
    x <- input$scatter.x
    y <- input$scatter.y
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    
    if (trans == "None") {
      scatter.data <- DS
    } else if (trans == "Natural log") {
      scatter.data <- log1p(DS)
    } else if (trans == "log2") {
      scatter.data <- log2(DS + 1)
    } else if (trans == "log10") {
      scatter.data <- log10(DS + 1)
    }
    scatter.end <- Sys.time()
    print("Scatter plot time")
    print(scatter.end - scatter.start)
    return(list(x, y, scatter.data))
  })
  
  
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  scatterplot <- function() {
    # get values from list
    li <- plotScatter()
    xval <- li[[1]]
    yval <- li[[2]]
    scatter.data <- li[[3]]
    
    # data frame needed for ggplot
    df <- data.frame(t1 = scatter.data[, xval], t2 = scatter.data[, yval])
    
    # get 2d kernel density estimate
    df$density <- get_density(df$t1, df$t2, n = 100)
    
    # plot heat scatter w/ ggplot
    p <- ggplot(df, aes(x = t1, y = t2, color = density, text = paste(xval, ": ", round(t1, 4), "\n", yval, ": ", round(t2, 4), sep = ""), group = 1)) +
      geom_point(shape = 19, size =    0.25) +
      scale_color_viridis()
    
    # modify label and fill defaults
    p <- p + xlab(xval) + ylab(yval) + labs(color = "KDE", title = paste("Scatter Plot, R=", round(cor(scatter.data[, xval], scatter.data[, yval]), 3)))
    
    # if checkbox is ticked, display regression line
    if (input$regline == TRUE) {
      p <- p + geom_smooth(method = lm, se = FALSE, size = 0.5, color = "blue")
    }
    
    
    hide("help_text_scatter")
    # add interactivity w/ plotly
    p
    ggplotly(p, tooltip = c("text"))
    
    
  }
  
  scatterplot_collage <- function() {
    li <- plotScatter()
    scatter.data <- li[[3]]
    par(mfrow = c(3, 3))
    for (i in 1:ncol(scatter.data)) {
      for (j in i:ncol(scatter.data)) {
        d <- kde2d(scatter.data[, i], scatter.data[, j])
        ColorLevels <- round(seq(min(d$z), max(d$z), length = 5), 4)
        heatscatter(x = scatter.data[, i], y = scatter.data[, j], xlab = colnames(scatter.data)[i], ylab = colnames(scatter.data)[j], main = "")
        legend("topleft", paste("R=", round(cor(scatter.data[, i], scatter.data[, j]), 3)), bty = "n")
        legend("bottomright", title = "KDE", legend = ColorLevels, pch = 19, col = LSD::colorpalette("heat"))
        if (i != j) {
          lines(lowess(scatter.data[, i], scatter.data[, j]), col = "black")
        }
      }
    }
  }
  
  observeEvent(input$submit_scatter, {
    output$scatter.plot <- renderPlotly({
      scatterplot()
    })
    
  })
  
  output$downloadscatter_collage <- downloadHandler(
    filename = function() {
      paste("heatscatter_collage", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      scatterplot_collage()
      dev.off()
    }
  )
  
  output$downloadscatter <- downloadHandler(
    filename = function() {
      paste("heatscatter", ".pdf", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(widget = scatterplot(), file = "scatterplot.html")
      webshot(url = "scatterplot.html", file = file)
    }
  )
  
  output$help_text_scatter <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Scatter plot compares global gene expression
          between 2 conditions.<br>Color density is calculated based on 2D Gaussian kernel density.
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ######## distfit ###########
  ############################
  
  
  output$downloaddist <- downloadHandler(
    filename = function() {
      paste("distribution_fit", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      distplot()
      dev.off()
    }
  )
  
  output$dist_range_allowed <- renderText({
    DS <- distfit_df()
    paste("Suggested range: ( 0", " ~ ", round(max(DS)), " ]", sep = "")
  })
  
  plotDist <- reactive({
    dist.start <- Sys.time()
    dis <- input$distributions
    var <- input$dist.var
    DS <- distfit_df()
    fits <- list()
    distrs <- NULL
    numcol <- c(0, 0, 0, 0, 0, 0)
    dist_zoom <- input$dist_zoom
    if (dist_zoom == "slider") {
      fit_range <- input$dist_range
    } else if (dist_zoom == "text input") {
      fit_range <- c(input$dist_range_min, input$dist_range_max)
    }
    if ("Log-normal" %in% dis) {
      fit_ln <- fitdist(DS[, var], "lnorm")
      fits <- c(fits, list(fit_ln))
      distrs <- c(distrs, "Log-normal")
      numcol[1] <- 1
    }
    if ("Log-logistic" %in% dis) {
      fit_ll <- fitdist(DS[, var], "llogis", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fits <- c(fits, list(fit_ll))
      distrs <- c(distrs, "Log-logistic")
      numcol[2] <- 1
    }
    if ("Pareto" %in% dis) {
      fit_P <- fitdist(DS[, var], "pareto", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fits <- c(fits, list(fit_P))
      distrs <- c(distrs, "Pareto")
      numcol[3] <- 1
    }
    if ("Burr" %in% dis) {
      fit_B <- fitdist(DS[, var], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1), lower = c(0, 0, 0))
      fits <- c(fits, list(fit_B))
      distrs <- c(distrs, "Burr")
      numcol[4] <- 1
    }
    if ("Weibull" %in% dis) {
      fit_W <- fitdist(DS[, var], "weibull", lower = c(0, 0))
      fits <- c(fits, list(fit_W))
      distrs <- c(distrs, "Weibull")
      numcol[5] <- 1
    }
    if ("Gamma" %in% dis) {
      fit_G <- fitdist(DS[, var], "gamma", lower = c(0, 0), start = list(scale = 1, shape = 1))
      fits <- c(fits, list(fit_G))
      distrs <- c(distrs, "Gamma")
      numcol[6] <- 1
    }
    dist.end <- Sys.time()
    print("Distribution fitting time")
    print(dist.end - dist.start)
    return(list(fits, distrs, numcol, var, fit_range))
  })
  
  distaic <- reactive({
    dist.start <- Sys.time()
    DS <- distfit_df()
    AIC.df <- as.data.frame(matrix(nrow = ncol(DS), ncol = 6))
    rownames(AIC.df) <- colnames(DS)
    colnames(AIC.df) <- c("Log-normal", "Log-logistic", "Pareto", "Burr", "Weibull", "Gamma")
    for (i in 1:nrow(AIC.df)) {
      fit_ln <- fitdist(DS[, i], "lnorm")
      fit_ll <- fitdist(DS[, i], "llogis", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fit_P <- fitdist(DS[, i], "pareto", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fit_B <- fitdist(DS[, i], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1), lower = c(0, 0, 0))
      fit_W <- fitdist(DS[, i], "weibull", lower = c(0, 0))
      fit_G <- fitdist(DS[, i], "gamma", lower = c(0, 0), start = list(scale = 1, shape = 1))
      fits <- list(fit_ln, fit_ll, fit_P, fit_B, fit_W, fit_G)
      AIC.df[i, ] <- gofstat(fits)$aic
    }
    for (i in 1:nrow(AIC.df)) {
      AIC.df$min.AIC[i] <- colnames(AIC.df)[which.min(AIC.df[i, 1:6])]
    }
    dist.end <- Sys.time()
    print("distribution fitting time")
    print(dist.end - dist.start)
    
    return(AIC.df)
  })
  
  distplot <- function() {
    li <- plotDist()
    fits <- li[[1]]
    distrs <- li[[2]]
    numcol <- li[[3]]
    var <- li[[4]]
    fit_range <- li[[5]]
    line_types <- c(1, 2, 3, 4, 5, 6) # par lty
    if (length(fits) != 0) {
      cdfcomp(fits,
              xlogscale = TRUE, ylogscale = TRUE,
              ylab = "CDF", xlab = "Expression levels (log)", xlim = c(fit_range[1], fit_range[2]), ylim = c(10^-3, 1),
              legendtext = distrs, cex = 0.5, main = var, fitcol = rainbow(6)[which(numcol == 1)], fitlty = line_types[which(numcol == 1)]
      )
    }
    
  }
  
  output$downloaddistaic <- downloadHandler(
    filename = function() {
      paste("aic", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(distaic(), file, row.names = TRUE)
    }
  )
  
  observeEvent(input$submit_distfit, {
    output$dist.plot <- renderPlot({
      distplot()
    })
    output$dist.aic <- renderTable({
      distaic()
    },
    rownames = TRUE
    )
    
  })
  
  
  output$help_text_dis_fit <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Fitting the selected probability distribution(s) to transcriptome-wide data of the selected sample.<br>
           Cumulative distribution function will be showned, with black lines being the empirical (transcriptomic) data, 
           and colored lines being the probability distribution(s) with best fitted parameters. <br>
           The AIC (Akaike information criterion) table provides a comparison in goodness-of-fit to transcriptomic data
           among selected probability distribution(s).
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ####### correlation ########
  ############################
  
  COR <- function(d, i, myMethod) {
    Result2 <- cor(x = d[, i], y = d[, i], method = myMethod)
    return(format(round(Result2, 5), nsmall = 5))
  }
  
  cor_df <- reactive({
    cor.start <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    
    
    method <- input$cor_method
    if (method == "Pearson correlation") {
      Cor2 <- data.frame(COR((DS), 1:length(DS), "pearson"))
      
    } else if (method == "Spearman correlation") {
      Cor2 <- data.frame(COR((DS), 1:length(DS), "spearman"))
    }
    Cor2 <- na.omit(Cor2)
    cor.end <- Sys.time()
    print("correlation time")
    print(cor.end - cor.start)
    return(Cor2)
  })
  
  observeEvent(input$submit_corr, {
    output$corr.plot <- renderPlot({
      corrplot1()
    })
    
    output$corr.plot2 <- renderPlot({
      corrplot2()
    })
    
    output$corr.matrix <- renderTable(
      {
        cor_df()
      },
      rownames = TRUE
    )
    #shinyjs::show("downloadcorrAll")
  })
  corrplot1 <- function() {
    corr <- as.matrix(cor_df())
    corr <- apply(corr, 2, as.numeric)
    rownames(corr) <- rownames(cor_df())
    if (ncol(corr) <= 20) {
      fontsize <- 1
    } else {
      fontsize <- 20 / ncol(corr)
    }
    corrplot(corr, method = "shade", shade.col = NA, tl.col = "black", cl.lim = c(min(corr), 1), is.corr = FALSE,tl.cex = fontsize)
    hide("help_text_correlation")
  }
  
  corrplot2 <- function() {
    corr <- as.matrix(cor_df())
    corr <- apply(corr, 2, as.numeric)
    rownames(corr) <- rownames(cor_df())
    if (ncol(corr) <= 20) {
      fontsize <- 1
    } else {
      fontsize <- 20 / ncol(corr)
    }
    corrplot(corr, type = "upper", tl.col = "black", cl.lim = c(min(corr), 1), is.corr = FALSE, tl.cex = fontsize)
    
  }
  
  output$downloadcorrplot <- downloadHandler(
    filename = function() {
      paste("corrheatmap", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      corrplot1()
      dev.off()
    }
  )
  
  output$downloadcorrplot2 <- downloadHandler(
    filename = function() {
      paste("corrplot", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      corrplot2()
      dev.off()
    }
  )
  
  output$downloadcorrmat <- downloadHandler(
    filename = function() {
      paste("correlation", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(cor_df(), file, row.names = TRUE)
    }
  )
  
  
  output$help_text_correlation <- renderUI({
    HTML("
    <br>
    <br>
      <center>
      <ol>
        <li>
        <p>
          <b>
          <h4>Pearson correlation</h4><br>
          Pearson correlation measures linear relationship between two vectors, where r = 1 if the two vectors are identical, 
          and r = 0 if there are no linear relationships between the vectors.<br>
          The correlation coefficient r between two vectors (e.g. transcriptome in two different samples), 
          containing n observations (e.g. gene expression values), is defined by (for large n):<br>
          <img src='https://www.statisticssolutions.com/wp-content/uploads/2019/09/ehtsht.png' alt='Pearson-correlation' border='0'><br>
          where x<sub>i</sub> and y<sub>i</sub> are the ith observation in the vectors X and Y, respectively.
          </b>
        </p>
        </li>
        <br>
        <li>
        <p>
          <b>
          <h4>Spearman correlation</h4><br>
          Spearman correlation is a non-parametric test that measrues the degree of association 
          between two vectors  (e.g. transcriptome in two different samples).  The Spearman rank correlation 
          test does not carry any assumptions about the distribution of the data and is the appropriate correlation 
          analysis when the variables are measured on a scale that is at least ordinal.
          The following formula is used to calculate the Spearman rank correlation:<br>
          <img src='https://www.dataanalytics.org.uk/wp-content/uploads/2019/06/Figure-8.5.png' alt='Spearman-correlation' border='0'><br>
          ρ = Spearman rank correlation<br>
          r<sub>x,i</sub>, r<sub>y,i</sub> = ranks of corresponding variables (or gene)<br>
          n = number of observations
          </b>
        </p>
        </li>
        </ol>
      </center>
    ")
  })
  
  ############################
  #######     PCA     ########
  ############################
  
  refreshDS1 <- eventReactive(input$pca_refresh, {
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    DS1 <- DS[sample(nrow(DS), nrow(DS), replace = FALSE), ]
    return(DS1)
  })
  
  plotPCA <- reactive({ # process and return data
    pca.start <- Sys.time()
    type <- input$file_type
    ############################
    pca_type <- input$pca_type
    ############################
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    order <- input$gene_order
    size <- input$gene_size
    x <- input$pca.x
    y <- input$pca.y
    cluster_flag <- input$pca_cluster
    rindex <- as.numeric(substring(x, 3))
    cindex <- as.numeric(substring(y, 3))
    if (order == "Ascending") {
      DS1 <- DS[order(DS[, 1]), ]
    } else if (order == "Descending") {
      DS1 <- DS[rev(order(DS[, 1])), ]
    } else if (order == "Random") {
      DS1 <- refreshDS1()
    }
    
    DSample <- head(DS1, n = size)
    
    ##### PCA & Sparse PCA #####
    if (pca_type == "PCA") {
      PR <- prcomp(t(DSample), center = TRUE)
      print("Normal PCA selected")
    }
    else if (pca_type == "SPCA") {
      PR <- spca(t(DSample), scale = FALSE, center = TRUE, max_iter = 10)
      PR$x <- PR$scores
      print("Sparse PCA selected")
    }
    
    col_val_x <- as.numeric(gsub("[^[:digit:]]", "", x))
    col_val_y <- as.numeric(gsub("[^[:digit:]]", "", y))
    #####################
    
    
    PCA.var <- PR$sdev^2
    PCA.var.per <- round(PCA.var / sum(PCA.var) * 100, 1)
    xlabel <- paste(colnames(PR$x)[rindex], " - ", PCA.var.per[rindex], "%", sep = "")
    ylabel <- paste(colnames(PR$x)[cindex], " - ", PCA.var.per[cindex], "%", sep = "")
    if (cluster_flag == TRUE) {
      num <- as.numeric(input$pca_cluster_num)
      ####################################################################################
      kmeans.data <- data.frame(x = PR$x[, col_val_x], y = PR$x[, col_val_y])
      print(kmeans.data)
      ####################################################################################
      set.seed(1)
      kmeans.result <- kmeans(kmeans.data, num)
      return(list(PR, PCA.var, PCA.var.per, rindex, cindex, xlabel, ylabel, cluster_flag, kmeans.result))
    }
    pca.end <- Sys.time()
    print("pca time")
    print(pca.end - pca.start)
    return(list(PR, PCA.var, PCA.var.per, rindex, cindex, xlabel, ylabel, cluster_flag))
  })
  
  pcavarplot <- function() {
    li <- plotPCA()
    PCA.var.per <- li[[3]] / 100
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    pcchoices <- NULL
    for (i in 1:length(PCA.var.per)) {
      pcchoices <- c(pcchoices, paste("PC", i, sep = ""))
    }
    xform <- list(
      categoryorder = "array",
      categoryarray = pcchoices
    )
    p <- plot_ly(
      x = pcchoices,
      y = PCA.var.per,
      name = "PCA variance",
      type = "bar"
    ) %>% layout(xaxis = xform)
    hide("help_text_PCA")
    return(p)
  }
  
  pca2dplot <- function() {
    li <- plotPCA()
    PR <- li[[1]]
    rindex <- li[[4]]
    cindex <- li[[5]]
    xlabel <- li[[6]]
    ylabel <- li[[7]]
    cluster_flag <- li[[8]]
    if (cluster_flag == FALSE) {
      p <- plot_ly(
        x = PR$x[, rindex],
        y = PR$x[, cindex],
        type = "scatter",
        mode = "markers"
      ) %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel))
    } 
    else if (cluster_flag == TRUE) {
      kmeans.result <- li[[9]]
      text_flag <- input$pca_text
      if (text_flag == TRUE) {
        p <- plot_ly(
          x = PR$x[, rindex],
          y = PR$x[, cindex],
          type = "scatter",
          color = as.character(kmeans.result$cluster),
          mode = "markers",
          colors = "Set1"
        ) %>%
          hide_colorbar() %>%
          add_trace(
            x = PR$x[, rindex],
            y = PR$x[, cindex],
            type = "scatter",
            mode = "text",
            text = names(kmeans.result$cluster),
            textposition = "top right"
          ) %>%
          layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel), showlegend = FALSE)
      } else if (text_flag == FALSE) {
        p <- plot_ly(
          x = PR$x[, rindex],
          y = PR$x[, cindex],
          type = "scatter",
          color = as.character(kmeans.result$cluster),
          mode = "markers",
          text = names(kmeans.result$cluster),
          colors = "Set1"
        ) %>%
          hide_colorbar() %>%
          layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel), showlegend = FALSE)
      }
    }
  }
  
  pca3dplot <- function() {
    li <- plotPCA()
    PR <- li[[1]]
    xlabel <- "PC1"
    ylabel <- "PC2"
    zlabel <- "PC3"
    cluster_flag <- li[[8]]
    if (cluster_flag == FALSE) {
      p <- plot_ly(
        x = PR$x[, 1],
        y = PR$x[, 2],
        z = PR$x[, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 5)
      ) %>% layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel), zaxis = list(title = zlabel)))
    } else if (cluster_flag == TRUE) {
      kmeans.result <- li[[9]]
      text_flag <- input$pca_text
      if (text_flag == TRUE) {
        p <- plot_ly(
          x = PR$x[, 1],
          y = PR$x[, 2],
          z = PR$x[, 3],
          type = "scatter3d",
          mode = "text",
          text = names(kmeans.result$cluster),
          color = as.character(kmeans.result$cluster),
          textfont = list(size = 10),
          textposition = "top right"
        ) %>% layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel), zaxis = list(title = zlabel)), showlegend = FALSE)
      } else if (text_flag == FALSE) {
        p <- plot_ly(
          x = PR$x[, 1],
          y = PR$x[, 2],
          z = PR$x[, 3],
          type = "scatter3d",
          color = as.character(kmeans.result$cluster),
          mode = "markers",
          marker = list(size = 5),
          text = names(kmeans.result$cluster),
          colors = "Set1"
        ) %>%
          hide_colorbar() %>%
          layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel), zaxis = list(title = zlabel)), showlegend = FALSE)
      }
    }
  }
  
  
  observeEvent(input$submit_pca, {
    output$pcavar.plot <- renderPlotly({
      pcavarplot()
    })
    
    output$pca2d.plot <- renderPlotly({
      pca2dplot()
    })
    
    output$pca3d.plot <- renderPlotly({
      pca3dplot()
    })
    
    
  })
  
  output$downloadpcavar <- downloadHandler(
    filename = function() {
      paste("pca_variance", ".pdf", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(widget = pcavarplot(), file = "pcavariance.html")
      webshot(url = "pcavariance.html", file = file)
    }
  )
  
  output$downloadpca2d <- downloadHandler(
    filename = function() {
      paste("pca2d", ".pdf", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(widget = pca2dplot(), file = "pca2d.html")
      webshot(url = "pca2d.html", file = file)
    }
  )
  
  output$downloadpca3d <- downloadHandler(
    filename = function() {
      paste("pca3d", ".pdf", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(widget = pca3dplot(), file = "pca3d.html")
      webshot(url = "pca3d.html", file = file)
    }
  )
  
  output$help_text_PCA <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Principal Components Analysis (PCA) is an multivariate statistical technique for simplifying high-dimentional 
          data sets (Basilevsky 1994). Given m observations on n variables, the goal of PCA is to reduce the dimensionality 
          of the data matrix by finding r new variables, where r is less than n. Termed principal components, these r new 
          variables together account for as much of the variance in the original n variables as possible while remaining 
          mutually uncorrelated and orthogonal. Each principal component is a linear combination of the original variables, 
          and so it is often possible to ascribe meaning to what the components represent.
          A PCA analysis of transcriptomic data consider the genes as variables, creating a set of “principal gene components” that indicate the features of genes that best explain the experimental responses they produce.
          </b>
        </p>
        <p>
          <b>
          A PCA analysis of transcriptomic data consider the genes as variables, creating a set of 
          “principal gene components” that indicate the features of genes that best explain the experimental responses they produce.
          </b>
        </p>
        <p>
          <b>
          To compute the principal components, the n eigenvalues and their corresponding eigenvectors are calculated 
          from the n×n covariance matrix of conditions. Each eigenvector defines a principal component. A component can be 
          viewed as a weighted sum of the conditions, where the coefficients of the eigenvectors are the weights. 
          The projection of gene i along the axis defined by the i<sup>th</sup> principal component is:<br>
          <img src='https://i.ibb.co/n6M8n79/Screenshot-from-2020-09-05-15-33-31.png' alt='Screenshot-from-2020-09-05-15-33-31' border='0'><br>
          Where <var>v<sub>tj</sub></var> is the t<sup>th</sup> coefficient for the i<sup>th</sup> principal component; <var>a<sub>it</sub></var> is the 
          expression measurement for gene i under the t<sup>th</sup> condition. 
          A′ is the data in terms of principal components. Since V is an orthonormal matrix, A′ is a rotation of the data from the original space of 
          observations to a new space with principal component axes.
          </b>
        </p>
        <p>
          <b>
          The variance accounted for by each of the components is its associated eigenvalue; 
          it is the variance of a component over all genes. Consequently, the eigenvectors with large eigenvalues 
          are the ones that contain most of the information; eigenvectors with small eigenvalues are uninformative.
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ######## DE analysis #######
  ############################
  group_names_de <- reactive({
    f <- group_names()
    f1 <- input$f1
    f2 <- input$f2
    if (is.null(f) || is.null(f1) || is.null(f2)) {
      return(NULL)
    }
    
    type <- input$file_type
    if (type == "norm") {
      raw_DS <- df_shiny()
    } else if (type == "raw") {
      raw_DS <- df_raw_filt()
    }
    f.df <- data.frame("f" = f)
    rownames(f.df) <- colnames(raw_DS)
    # f.df.slect <- subset(f.df,f %in% c(f1,f2) )
    # f.df.slect <- rbind( subset(f.df,f %in% f1),  subset(f.df,f %in% f2) )
    # f.df.slect2 <- f.df.slect;
    # f_new <- f.df.slect[,1]
    # f.df.slect2[,1] <- droplevels(f_new, except = levels(f_new)%in%f_new)
    return(f.df)
  })
  
  df_raw_de <- reactive({
    f_de <- group_names_de()
    if (is.null(f_de)) {
      return(NULL)
    }
    type <- input$file_type
    rep_number <- input$n_rep
    if (rep_number == 1) {
      de_type <- input$de_method1
    } else {
      de_type <- input$de_method0
    }
    
    if (type == "norm") {
      raw_DS <- df_shiny() # filtered and normalized
    } else if (type == "raw") {
      if (de_type == "NOISeq") {
        raw_DS <- df_raw_shiny() # filtered and normalized
      } else {
        raw_DS <- df_raw_filt() # filtered and UN-NORMALIZED
      }
    }
    # raw_DS_de <- raw_DS[,rownames(f_de)]
    return(raw_DS)
  })
  
  # all helper functions are in utils.R file
  de_no_filt <- eventReactive(input$submit_DE, { # return as table object
    start.de.table <- Sys.time()
    f_de <- group_names_de() # for edgeR >> f=f_de[,1]; for the rest factors=f_de
    if (is.null(f_de)) {
      return(NULL)
    }
    DS_de <- df_raw_de() # with only 2 conditions for DE analysis
    p_val <- 1 # input$p_val
    fc <- 1 # input$fc
    f1 <- input$f1
    f2 <- input$f2
    rep_number <- input$n_rep # either 0 = no replicates or 1 = have replicates
    
    
    spikes <- neg_control()
    norm_method <- input$norm_method
    if (!is.null(spikes) & norm_method == "RUV") {
      set1 <- RUVg.apply(DS_de, spikes)
      W_1 <- pData(set1)$W_1
    } else {
      W_1 <- NULL
    }
    # print("from de_no_filt")
    # print(pData(set1))
    # print(W_1)
    
    if (rep_number == 1) { # have replicates
      de_type <- input$de_method1
      if (de_type == "EdgeR") {
        res <- edgerApply(DS = DS_de, f = f_de[, 1], W_1 = W_1, f1 = f1, f2 = f2) # edgeR, return edgeR object
        res.df <- edgerFilter(res, FC = fc, p_val = p_val) # fitler edgeR object result
      } else if (de_type == "DESeq2") {
        res <- deseqApply(DS = DS_de, f.df = f_de, W_1 = W_1, f1 = f1, f2 = f2) # DESeq, return DESeq object
        res.df <- deseqFilter(res, FC = fc, p_val = p_val) # fitler DESeq object result
      } else if (de_type == "NOISeq") {
        res <- noiseqbioApply(DS = DS_de, f.df = f_de, f1 = f1, f2 = f2) # NOISeqbio, return NOIseq object
        res.df <- noiseqbioFilter(res, FC = fc, p_val = p_val) # filter return NOIseq object
      }
    } else { # no replicates
      de_type <- input$de_method0 # NOISeq
      res <- noiseqsimApply(DS = DS_de, f.df = f_de, f1 = f1, f2 = f2) # NOISeqbio, return NOIseq object
      res.df <- noiseqsimFilter(res, FC = fc)
    }
    end.de.table <- Sys.time()
    print("de table time")
    print(end.de.table - start.de.table)
    return(res.df)
  })
  
  
  de_filt <- function(res.df, p_val, fc, rep_number) {
    # res.df <- de_no_filt()
    if (is.null(res.df)) {
      return(NULL)
    }
    # p_val <- input$p_val
    # fc <- input$fc
    # rep_number <- input$n_rep
    
    if (rep_number == 1) {
      res.df.filt <- filter(res.df, FDR <= p_val, log2FCabs >= log2(fc))
    } else {
      res.df.filt <- filter(res.df, FDR <= p_val, log2FCabs >= log2(fc))
    }
    return(res.df.filt)
  }
  
  output$DE_table <- DT::renderDataTable({
    res.df <- de_no_filt()
    p_val <- input$p_val
    fc <- input$fc
    rep_number <- input$n_rep
    if (input$submit_DE > 0) {
      res.df.filt <- de_filt(res.df, p_val, fc, rep_number)
      res.df.filt
    }
    
    hide("help_text_DE_anal")
  })
  
  ##### volcano plot ######
  volcano_plot <- eventReactive(input$submit_DE, {
    volcano.start.time <- Sys.time()
    rep_number <- input$n_rep
    if (rep_number == 0) {
      return(NULL)
    }
    de_type <- input$de_method1
    if (de_type == "NOISeq") {
      return(NULL)
    }
    
    res <- de_no_filt() # de result, no filter
    if (is.null(res)) {
      return(NULL)
    }
    
    p_val <- input$p_val
    fc <- input$fc
    res$Gene <- rownames(res)
    res <- na.omit(res)
    # plot
    ymax <- quantile(-log10(res$PValue), c(0.98))
    xmax <- quantile(abs(res$log2FC), c(0.98))
    # if (ymax > 5) ymax <- 5
    # print("from volcano plot - range(res$PValue)")
    # print(range(res$PValue))
    with(res, plot(log2FC, -log10(PValue), pch = 20, main = "Volcano plot", 
                   xlim = c(-xmax, xmax), ylim = c(0, ymax))) # xlim=c(-5,5),ylim=c(0,ymax)
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(res, FDR < p_val), points(log2FC, -log10(PValue), pch = 20, col = "red"))
    with(subset(res, abs(log2FC) > log2(fc)), points(log2FC, -log10(PValue), pch = 20, col = "orange"))
    with(subset(res, FDR < p_val & abs(log2FC) > log2(fc)), points(log2FC, -log10(PValue), pch = 20, col = "green"))
    legend("topleft",
           bty = "n", col = c("red", "orange", "green", "black"), pch = 19,
           legend = c("FDR < FDR limit", "FC > FC limit", "Both", "Other")
    )
    # library(calibrate)
    # with(subset(res, FDR<.05 & abs(log2FC)>1), textxy(log2FC, -log10(PValue), labs=Gene, cex=.8))
    volcano.end.time <- Sys.time()
    print("volcano time")
    print(volcano.end.time - volcano.start.time)
    
  })
  
  output$volcano_plot <- renderPlot({
    volcano_plot()
  })
  
  ##### dispersion plot ######
  dispersion_plot <- eventReactive(input$submit_DE, {
    dispersion.start.time <- Sys.time()
    f_de <- group_names_de() # for edgeR >> f=f_de[,1]; for the rest factors=f_de
    if (is.null(f_de)) {
      return(NULL)
    }
    rep_number <- input$n_rep # either 0 or 1
    if (rep_number == 0) {
      return(NULL)
    }
    DS_de <- df_raw_de() # with only 2 conditions for DE analysis
    p_val <- 1 # input$p_val
    fc <- 1 # input$fc
    de_type <- input$de_method1
    if (de_type == "EdgeR") {
      edgerDisp(DS_de, f_de[, 1])
    } else if (de_type == "DESeq2") {
      deseqDisp(DS_de, f_de)
    }
    dispersion.end.time <- Sys.time()
    print("dispersion time")
    print(dispersion.end.time - dispersion.start.time)
    
    
  })
  
  output$dispersion_plot <- renderPlot({
    dispersion_plot()
  })
  
  ########## download buttons DE analysis ###########
  output$download_de_table <- downloadHandler(
    filename = function() {
      paste0("DE analysis", ".csv")
    },
    content = function(file) {
      res.df <- de_no_filt()
      p_val <- input$p_val
      fc <- input$fc
      rep_number <- input$n_rep
      res.df.filt <- de_filt(res.df, p_val, fc, rep_number)
      write.csv(res.df.filt, file, row.names = F)
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("Volcano", ".pdf")
    },
    content = function(file) {
      pdf(file)
      volcano_plot()
      dev.off()
    }
  )
  
  output$download_dispersion <- downloadHandler(
    filename = function() {
      paste0("Dispersion plot", ".pdf")
    },
    content = function(file) {
      pdf(file)
      dispersion_plot()
      dev.off()
    }
  )
  
  output$help_text_DE_anal <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          DE analysis identifies the genes that are statistically different in expression levels between the 2 selected conditions. Two important threshold are:
          </b>
        </p>
        <p>
          <b>
            1. The lower bound of expression fold change  between the 2 selected conditions<br>
            2. The upper bound of hypothesis test p-value
          </b>
        </p>
        <p>
          <b>
            GeneCloudOmics implements 3 popular methods to identify DE genes:<br>
              1. <a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8'>DESeq2</a><br>
              2. <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/'>EdgeR</a><br>
              3. <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666377/'>NOISeq</a>
          </b>
        </p>
      </center>
    ")
  })
  
  # output$download_heatmap <- downloadHandler(
  #   filename = function(){
  #     paste0("Heatmap",".pdf")
  #   },
  #   content = function(file) {
  #     pdf(file)
  #     print(heatmap_plot())
  #     dev.off()
  #   }
  # )
  
  ############################
  ######### heatmap ##########
  ############################
  
  ####### heatmap renderUI commented #########
  # output$expand_genonames <- renderUI({
  #   type <- input$file_type
  #   if(type=='norm'){
  #     DS <- df_shiny()
  #   }else if(type=='raw'){
  #     DS <- df_raw_shiny()
  #   }
  #   if(ncol(DS)==input$numOfGeno){
  #     lapply(1:input$numOfGeno, function(i) {
  #       textInput(paste('type',i,sep=""), paste('Type',i,sep=" "),value = colnames(DS)[i])
  #     })
  #   }else{
  #     lapply(1:input$numOfGeno, function(i) {
  #       textInput(paste('type',i,sep=""), paste('Type',i,sep=" "))
  #     })
  #   }
  # })
  #
  # output$refGeno <- renderUI({
  #   selectInput('heatmap_anchor',"Reference genotype",choices=c(1:input$numOfGeno))
  # })
  #
  output$heatmap_display <- renderUI({
    display <- "ALL"
    for (i in 1:input$numOfCluster) {
      display <- c(display, i)
    }
    selectInput("display_cluster", "Display cluster", choices = display)
  })
  ################
  
  setOneWithinFold <- function(arr) { # logFC
    fold <- as.numeric(input$fold)
    for (i in 1:length(arr)) {
      if ((arr[i] <= (fold)) & (arr[i] >= (1 / fold))) {
        arr[i] <- 1
      }
    }
    return(arr)
  }
  
  plotHeatmap <- eventReactive(input$heatmap_plot, { # process and return data
    heatmap.start.time <- Sys.time()
    
    type <- input$file_type
    value <- input$heatmap_value
    de_type <- input$heatmap_de_ind
    
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    names <- NULL
    # for(i in 1:input$numOfGeno){
    #   id <- paste('type',i,sep="")
    #   names <- c(names,input[[id]])
    # }
    
    # numOfGeno <- input$numOfGeno
    # ref <- as.numeric(input$heatmap_anchor)
    clusterNum <- input$numOfCluster
    
    if (de_type == "ind") {
      fold <- as.numeric(input$fold)
      fold_ncol <- input$fold_ncol
      DS2 <- deWithoutStats(DS, FC = fold, n_col = fold_ncol)
      de_genes <- rownames(DS2)
      # print("from heatmap YT version")
      # print("de genes")
      # print(head(de_genes))
      # print(paste(length(de_genes),"genes"))
    } else if (de_type == "de") {
      res.df <- de_no_filt()
      if (is.null(res.df)) {
        return(NULL)
      }
      p_val <- input$p_val
      fc <- input$fc
      rep_number <- input$n_rep
      res.df.filt <- de_filt(res.df, p_val, fc, rep_number)
      de_genes <- res.df.filt$Gene
      # print("from line heatmap de result from DE analysis")
      # print("res.df.filt")
      # print(head(res.df.filt))
    }
    
    de_genes_exp <- DS[rownames(DS) %in% de_genes, ]
    DS3 <- t(scale(t(de_genes_exp)))
    DS3 <- na.omit(DS3)
    # print("from line 1894 - heatmap de type")
    # print("DS3")
    # print(head(DS3))
    
    set.seed(110)
    set.seed(1)
    a <- ComplexHeatmap::Heatmap(DS3,
                                 name = "Normalized expression",
                                 col = colorRamp2(c(min(DS3), 0, max(DS3)), c("red", "black", "green")),
                                 row_names_gp = gpar(fontsize = 1),
                                 row_dend_gp = gpar(fontsize = 1),
                                 row_title_gp = gpar(fontsize = 10),
                                 cluster_columns = FALSE,
                                 row_dend_width = unit(3, "cm"),
                                 split = clusterNum, clustering_distance_rows = "pearson",
                                 show_heatmap_legend = TRUE,
                                 show_row_names = FALSE, show_column_names = T,
                                 heatmap_legend_param = list(title = "Normalized expression")
    )
    set.seed(110)
    rcl.list <- ComplexHeatmap::row_order(a)
    DS3.1 <- as.matrix(rownames(DS3))
    
    # Cluster <- NULL
    # for(i in 1:length(rcl.list)){
    #   for(j in 1:length(rcl.list[[i]])){
    #     pair <- c(i,DS3.1[rcl.list[[i]][j]])
    #     Cluster <- rbind(Cluster,pair)
    #   }
    # }
    # Cluster <- data.frame(Cluster,row.names = NULL)
    # colnames(Cluster) <- c("cluster","GeneID")
    
    rcl.list2 <- rcl.list
    for (i in 1:length(rcl.list)) {
      rcl.list2[[i]] <- rownames(DS3)[rcl.list[[i]]]
    }
    
    for (i in 1:length(rcl.list2)) {
      genes <- rcl.list2[[i]]
      group_name <- rep(i, length(genes))
      Cluster_i <- data.frame("GeneID" = genes, "Cluster" = i)
      if (i == 1) {
        Cluster <- Cluster_i
      } else {
        Cluster <- rbind(Cluster, Cluster_i)
      }
    }
    
    print("line 2089, Cluster")
    print(head(Cluster))
    print(paste0("de_type = ", de_type))
    print("length of input DS")
    print(dim(DS3))
    
    # end heat map analysis
    heatmap.end.time <- Sys.time()
    print("heat map time")
    print(heatmap.end.time - heatmap.start.time)
    return(list(a, DS3, Cluster))
  })
  
  getCluster <- eventReactive(input$heatmap_plot, {
    set.seed(110)
    ll <- plotHeatmap()
    a <- ll[[1]]
    DS3 <- ll[[2]]
    rcl.list <- ComplexHeatmap::row_order(a)
    DS3.1 <- as.matrix(rownames(DS3))
    
    Cluster <- NULL
    for (i in 1:length(rcl.list)) {
      for (j in 1:length(rcl.list[[i]])) {
        pair <- c(i, DS3.1[rcl.list[[i]][j]])
        Cluster <- rbind(Cluster, pair)
      }
    }
    Cluster <- data.frame(Cluster, row.names = NULL)
    colnames(Cluster) <- c("cluster", "gene.id")
    
    return(Cluster[, c("gene.id", "cluster")])
  })
  
  mapPlot <- function() {
    myHeatmap <- plotHeatmap()[[1]]
    myHeatmap <- draw(myHeatmap)
    
    hide("help_text_heatmap")
  }
  
  output$heatmap.plot <- renderPlot({
    mapPlot()
  })
  
  output$downloadheatmap <- downloadHandler(
    filename = function() {
      paste("heatmap", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      p <- mapPlot()
      dev.off()
    }
  )
  
  
  output$cluster.info <- DT::renderDataTable({
    clusternum <- input$display_cluster
    gl <- plotHeatmap()[[3]] # getCluster()
    if (!is.null(gl)) {
      if (clusternum == "ALL") {
        gl
      } else {
        clusternum <- as.numeric(clusternum)
        dplyr::filter(gl, cluster == clusternum)
      }
    }
    
  })
  
  output$downloadclusters <- downloadHandler(
    filename = function() {
      paste("genelist", ".csv", sep = "")
    },
    content = function(file) {
      gl <- plotHeatmap()[[3]]
      write.csv(gl, file, row.names = FALSE)
    }
  )
  
  output$help_text_heatmap <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Hierarchical clustering is used to find the groups of co-expressed genes. 
          The clustering is performed on normalized expressions of differentially expressed genes using Ward clustering method. 
          Normalized expression of the jth gene at time ti is defined as<br>
          <img src='https://i.ibb.co/tJgCMVD/Screenshot-from-2020-09-05-16-37-39.png' alt='Screenshot-from-2020-09-05-16-37-39' border='0'><br>
          where <var>x<sub>j</sub>(t<sub>i</sub>)</var> is the expression of the j<sup>th</sup> gene at time t<sub>i</sub>, 
          <var>x<sub>j</sub>&#772</var> is the mean expression across all time points, and <var>&#963<sub>j</sub></var> is the standard deviation.
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ########## noise ###########
  ############################
  
  SQCO <- function(MT) {
    if (ncol(MT) == 1) {
      res <- matrix(0)
      return(res)
    }
    temp <- NULL
    for (i in 1:nrow(MT)) {
      m <- sum(MT[i, ]) / length(MT[i, ])
      v <- stats::var(MT[i, ])
      if (m != 0) {
        temp <- c(temp, v / (m * m))
      }
    }
    res <- matrix(sum(temp) / length(temp))
    return(res)
  }
  
  output$expand_genonames_noise <- renderUI({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS) / numOfRep
    
    lapply(1:numOfGeno, function(i) {
      textInput(paste("noisetype", i, sep = ""), paste("Type", i, sep = " "), value = colnames(DS)[(i - 1) * numOfRep + 1])
    })
  })
  
  output$noise_anchor_choices <- renderUI({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS) / numOfRep
    names <- NULL
    for (i in 1:numOfGeno) {
      id <- paste("noisetype", i, sep = "")
      names <- c(names, input[[id]])
    }
    selectInput("noise_anchor_b", "Anchor genotype", choices = names)
  })
  noisePlot <- function(){
    
    noise.start.time <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS) / numOfRep
    graph <- input$noise_graph_type
    names <- NULL
    for (i in 1:numOfGeno) {
      id <- paste("noisetype", i, sep = "")
      names <- c(names, input[[id]])
    }
    
    situation <- input$noise_situation
    if (situation == "a") {
      DS1 <- list()
      for (j in 1:numOfGeno) {
        DS1[[j]] <- as.matrix(DS[, ((j - 1) * numOfRep + 1):(j * numOfRep)])
      }
      Noise <- NULL
      for (y in 1:numOfGeno) {
        Noise <- c(Noise, SQCO(DS1[[y]]))
      }
      xform <- list(
        categoryorder = "array",
        categoryarray = names
      )
      if (graph == "Bar chart") {
        p <- plot_ly(
          x = names,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      } else if (graph == "Line chart") {
        p <- plot_ly(
          x = names,
          y = Noise,
          type = "scatter",
          mode = "lines+markers"
        ) %>% layout(xaxis = xform, yaxis = list(range = c(0, max(Noise) + 0.001)))
      }
    } else if (situation == "b") {
      DS_ave <- NULL
      for (j in 1:numOfGeno) {
        part_DS <- as.matrix(DS[, ((j - 1) * numOfRep + 1):(j * numOfRep)])
        DS_ave <- cbind(DS, data.frame(matrixStats::rowMeans2(part_DS)))
      }
      anchor <- input$noise_anchor_b
      names <- NULL
      for (i in 1:numOfGeno) {
        id <- paste("noisetype", i, sep = "")
        names <- c(names, input[[id]])
      }
      anchor_index <- match(anchor, names)
      Noise <- NULL
      for (i in 1:numOfGeno) {
        if (i != anchor_index) {
          Noise <- c(Noise, SQCO(cbind(DS_ave[, anchor_index], DS_ave[, i])))
        }
      }
      names_wo_anchor <- NULL
      for (i in 1:numOfGeno) {
        if (i != anchor_index) {
          id <- paste("noisetype", i, sep = "")
          names_wo_anchor <- c(names_wo_anchor, input[[id]])
        }
      }
      xform <- list(
        categoryorder = "array",
        categoryarray = names_wo_anchor
      )
      if (graph == "Bar chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      } else if (graph == "Line chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "scatter",
          mode = "lines+markers"
        ) %>% layout(xaxis = xform, yaxis = list(range = c(0, max(Noise) + 0.001)))
      }
    } else if (situation == "c") {
      anchor <- input$noise_anchor_c
      names <- colnames(DS)
      anchor_index <- match(anchor, names)
      Noise <- NULL
      for (i in 1:ncol(DS)) {
        if (i != anchor_index) {
          Noise <- c(Noise, SQCO(cbind(DS[, anchor_index], DS[, i])))
        }
      }
      names_wo_anchor <- names[-anchor_index]
      xform <- list(
        categoryorder = "array",
        categoryarray = names_wo_anchor
      )
      if (graph == "Bar chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      } else if (graph == "Line chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "scatter",
          mode = "lines+markers"
        ) %>% layout(xaxis = xform, yaxis = list(range = c(0, max(Noise) + 0.001)))
      }
    }
    noise.end.time <- Sys.time()
    print("noise time")
    print(noise.end.time - noise.start.time)
    
    hide("help_text_Noise")
    return(p)
    
  }
  observeEvent(input$noise_plot, {
    output$noise.plot <- renderPlotly({
      noisePlot()
    })
    
  })
  
  
  
  output$downloadnoise <- downloadHandler(
    filename = function() {
      paste("noise", ".pdf", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(widget = noisePlot(), file = "noise.html")
      webshot(url = "noise.html", file = file)
    }
  )
  
  output$help_text_Noise <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          To quantify between gene expressions scatter of all replicates in one experimental condition, we 
          computed transcriptome-wide average noise for each cell type, defined as<br>
          <img src='https://i.ibb.co/kcxrCzv/Screenshot-from-2020-09-05-16-14-27.png' alt='Screenshot-from-2020-09-05-16-14-27' border='0'><br>
          where <var>n</var> is the number of genes and <var>n<sub>i</sub><sup>2</sup></var> is the pairwise noise of the i<sup>th</sup>
          gene (variability between any two replicates), defined as<br>
          <img src='https://i.ibb.co/dp21hsK/Screenshot-from-2020-09-05-16-14-57.png' alt='Screenshot-from-2020-09-05-16-14-57' border='0'><br>
          where <var>m</var> is the number of replicates in each condition and <var>n<sub>ijk</sub><sup>2</sup></var> is the expression noise of the i<sup>th</sup> gene, 
          defined by the variance divided by the squared mean expression in the pair of replicates (j,k).<br>
          Citation: <a href='https://www.nature.com/articles/srep07137'>https://www.nature.com/articles/srep07137</a> (Kumar’s embryonic development paper)
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ######### entropy ##########
  ############################
  
  computeBin <- function(arr) { # Doane's rule
    n <- length(arr)
    gx <- moments::skewness(arr)
    sigmag <- sqrt(6 * (n - 2) / ((n + 1) * n + 3))
    bin <- 1 + log2(n) + log2(1 + abs(gx) / sigmag)
    return(bin)
  }
  
  getBinCounts <- function(arr) {
    vec <- entropy::discretize(arr, computeBin(arr), r = range(arr))
    return(vec)
  }
  
  output$expand_genonames_entropy <- renderUI({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    tp <- as.numeric(input$entropy_timepoints)
    numOfGeno <- ncol(DS) / tp
    
    lapply(1:numOfGeno, function(i) {
      textInput(paste("entropytype", i, sep = ""), paste("Type", i, sep = " "), value = colnames(DS)[(i - 1) * tp + 1])
    })
  })
  
  entropyPlot <- reactive({
    entropy.start.time <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (is.null(DS) == FALSE) {
      tsflag <- input$tsflag
      graph <- input$entropy_graph_type
      names <- colnames(DS)
      xform <- list(
        categoryorder = "array",
        categoryarray = names
      )
      entropy.vector <- NULL # entropy of each column
      for (i in 1:length(DS)) {
        binCount <- getBinCounts(DS[, i])
        entropy <- entropy.empirical(binCount, unit = "log2")
        entropy.vector <- c(entropy.vector, entropy)
      }
      if (tsflag == FALSE) {
        if (graph == "Bar chart") {
          p <- plot_ly(
            x = names,
            y = entropy.vector,
            type = "bar"
          ) %>% layout(xaxis = xform)
        } else if (graph == "Line chart") {
          p <- plot_ly(
            x = names,
            y = entropy.vector,
            type = "scatter",
            mode = "lines+markers"
          ) %>% layout(xaxis = xform)
          # yaxis=list(range = c(0, max(ent)+0.002))
        }
      } else if (tsflag == TRUE) {
        tp <- as.numeric(input$entropy_timepoints)
        numOfGeno <- ncol(DS) / tp
        names <- NULL
        for (i in 1:numOfGeno) {
          id <- paste("entropytype", i, sep = "")
          names <- c(names, input[[id]])
        }
        time_index <- c(1:tp)
        ent <- data.frame(time_index)
        for (j in 1:numOfGeno) {
          part_ent <- entropy.vector[(tp * j - (tp - 1)):(tp * j)]
          ent <- cbind(ent, part_ent)
        }
        if (graph == "Bar chart") {
          p <- plot_ly(x = ent[, 1], y = ent[, 2], name = names[1], type = "bar")
          for (i in 1:(numOfGeno - 1)) {
            p <- add_trace(p, y = ent[, i + 2], name = names[i + 1], type = "bar")
          }
          p <- layout(p, xaxis = list(title = "Time"), yaxis = list(title = "Entropy"))
        } else if (graph == "Line chart") {
          p <- plot_ly(x = ent[, 1], y = ent[, 2], name = names[1], type = "scatter", mode = "lines+markers")
          for (i in 1:(numOfGeno - 1)) {
            p <- add_trace(p, y = ent[, i + 2], name = names[i + 1], type = "scatter", mode = "lines+markers")
          }
          p <- layout(p, xaxis = list(title = "Time"), yaxis = list(title = "Entropy"))
        }
      }
      entropy.end.time <- Sys.time()
      print("entropy time")
      print(entropy.end.time - entropy.start.time)
      
      
      hide("help_text_Entropy")
      return(p)
    }
  })
  
  observeEvent(input$submit_entropy, {
    output$entropy.plot <- renderPlotly({
      entropyPlot()
    })
    
  })
  
  
  output$downloadentropy <- downloadHandler(
    filename = function() {
      paste("entropy", ".pdf", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(widget = entropyPlot(), file = "entropy.html")
      webshot(url = "entropy.html", file = file)
    }
  )
  
  output$help_text_Entropy <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Shannon entropy (Shannon, 1948) measures the disorder of a high-dimensional
          system, where higher values indicate increasing disorder.  Entropy of each transcriptome, X, is defined as<br>
          <img src='https://i.ibb.co/5W0KwMP/Screenshot-from-2020-09-05-16-31-03.png' alt='Screenshot-from-2020-09-05-16-31-03' border='0'><br>
          where p(xi) is the probability of gene expression value x=xi.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ########## New-Features ###########
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ############   SVM     ############
  ###################################
  ###################################
  
  # data_svm <- eventReactive(input$submit_svm, {
  #   print("Running...")
  #   svm.start <- Sys.time()
  #   type <- input$file_type
  #   
  #   if (type == "norm") {
  #     DS <- df_shiny()
  #   }
  #   else if (type == "raw") {
  #     DS <- df_raw_shiny()
  #   }
  #   
  #   print("1")
  #   
  #   x1.1 <- as.data.frame(DS[,1:3])
  #   x1.11 <- as.data.frame(DS[1:100,1:3])
  #   x2.2 <- as.data.frame(DS[,4:6])
  #   x2.22 <- as.data.frame(DS[1:100,4:6])
  #   
  #   print("2")
  #   
  #   x1.1 <-  setNames(stack(x1.1),c("x1.1","colName"))
  #   x1.11 <-  setNames(stack(x1.11),c("x1.1","colName"))
  #   x2.2 <-  setNames(stack(x2.2),c("x2.2","colName"))
  #   x2.22 <-  setNames(stack(x2.22),c("x2.2","colName"))
  #   
  #   mut_type <- as.data.frame(x1.1[,2])
  #   mut_type1 <- as.data.frame(x1.11[,2])
  #   
  #   dat <- data.frame(matrix(ncol = 3, nrow = nrow(x1.1)))
  #   dat1 <- data.frame(matrix(ncol = 3, nrow = nrow(x1.11)))
  #   x <- c("x1.1", "x2.2", "y")
  #   colnames(dat) <- x
  #   colnames(dat1) <- x
  #   
  #   print("3")
  #   
  #   dat[,1] <- x1.1[,1]
  #   dat[,2] <- x2.2[,1]
  #   dat[,3] <- mut_type[,1]
  #   
  #   print("4")
  #   
  #   dat1[,1] <- x1.11[,1]
  #   dat1[,2] <- x2.22[,1]
  #   dat1[,3] <- mut_type1[,1]
  #   
  #   dat[,3] <- as.factor(as.numeric(dat[,3]))
  #   dat1[,3] <- as.factor(as.numeric(dat1[,3]))
  #   
  #   return(dat1)
  #   
  # })
  # 
  # plotSVM <- function() {
  #   
  #   dat <- data_svm()
  #   
  #   print("training")
  #   svmfit <- svm(y~., data = dat, kernel = "radial", cost = 10, gamma = 1)
  #   print("done")
  #   plot(svmfit , dat )
  #   
  # }
  # 
  # plotSVM_df <- function() {
  #   
  #   dat <- data_svm()
  #   
  #   ggplot(data = dat, aes(x = x2.2, y = x1.1, color = y, shape = y)) + 
  #     geom_point(size = 2) +
  #     scale_color_manual(values=c("#000000","#FF0000","#00BA00")) +
  #     theme(legend.position = "none")
  #   
  # }
  # 
  # output$svm_plot <- renderPlot({
  #   plotSVM()
  # })
  # 
  # output$svm_df_plot <- renderPlot({
  #   plotSVM_df()
  # })
  # 
  # output$help_text_SVM <- renderUI({
  #   HTML("<h3><b>To be implemented</b></h3>")
  # })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  ############   t-SNE    ###########
  ############# Python ##############
  ###################################
  # data for t-sne
  plotTSNE2 <- function() {
    
    tsne2_trans <- input$tsne2_trans
    type <- input$file_type
    perplexity_value <- input$perplexity_value
    no_of_pca <- input$no_of_pca
    # no_of_clusters <- input$no_of_clusters
    
    
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){ 
      DS <- df_raw_shiny()
    }
    if(tsne2_trans=='None'){
      tsne2.data <- t(apply(DS, MARGIN = 1, scale)); 
      colnames(tsne2.data) <- colnames(DS)
    }
    else if(tsne2_trans=='log10'){
      tsne2.data <- log10(DS+1)
    }
    
    tsne_cluster_flag <- input$tsne_cluster # ture or false
    
    return (list(tsne2.data, perplexity_value, no_of_pca, tsne_cluster_flag)) #, no_of_clusters
  }
  
  tsne2plot <-eventReactive(input$submit_tsne2,{
    set.seed(13)
    tsne2.start <- Sys.time()
    # get data 
    li <- plotTSNE2()
    tsne2.data <- li[[1]]
    perplexity_value <- li[[2]]
    no_of_pca <- li[[3]]
    tsne_cluster_flag <- li[[4]]
    tsne_text_flag <- input$tsne_text # display sample name or not
    # no_of_clusters <- li[[4]] 
    
    # get tsne value
    tsne_val <- Rtsne(t(tsne2.data), 
                      dims = 2,
                      initial_dims = no_of_pca,
                      perplexity = perplexity_value,
                      theta = 0.0)
    tsne_df <- data.frame(
      TSNE1 = tsne_val$Y[, 1],
      TSNE2 = tsne_val$Y[, 2],
      Sample = colnames(tsne2.data)
    )
    
    if(!tsne_cluster_flag){
      # plotting
      p <- plot_ly(data = tsne_df, x = ~TSNE1, y = ~TSNE2, text = ~Sample) %>% 
        add_trace(type = "scatter", mode = 'markers', opacity = 0.5)
      
    } else { # tsne_cluster_flag == TRUE
      tsne_cluster_num <- as.numeric(input$tsne_cluster_num)
      set.seed(1)
      tsne_kmeans_result <- kmeans(tsne_df[,1:2], tsne_cluster_num)
      tsne_df$cluster <- factor(tsne_kmeans_result$cluster, levels = 1:max(tsne_kmeans_result$cluster) )
      
      # plotting
      p <- plot_ly(data = tsne_df, x = ~TSNE1, y = ~TSNE2, text = ~Sample, color = ~cluster ) %>%
        add_trace(type = "scatter", mode = 'markers', opacity = 0.5)
    }
    
    if(tsne_text_flag){
      p <- p %>% hide_colorbar() %>%
        add_trace(type = "scatter", mode = 'text', textposition = "top right", showlegend = FALSE)
    }
    
    tsne2.end <- Sys.time()
    print("t-SNE plot time")
    print(tsne2.end - tsne2.start)
    hide("help_text_tsne")
    return(list(p, tsne_df))
  })
  
  output$tsne2.plot <- renderPlotly({
    li <- tsne2plot()
    p <- li[[1]]
    tsne_table <- li[[2]]
    p
  })
  
  output$tsne_table <- DT::renderDataTable({
    tsne_table <- tsne2plot()[[2]] # get table
    
    tsne_table
  })
  
  output$download_tsne <- downloadHandler(
    filename = function() {
      paste("tsne_list", ".csv", sep = "")
    },
    content = function(file) {
      gl <- tsne2plot()[[2]]
      write.csv(gl, file, row.names = FALSE)
    }
  )
  tsneplothtml<- function(){
    tsne2plot()[[1]]
  }
  
  output$download_tsne2 <- downloadHandler(
    filename = function() {
      paste("tsne_plot", ".pdf", sep = "")
    },
    content = function(file) {
      htmlwidgets::saveWidget(widget = tsneplothtml(), file = "tsneplot.html")
      webshot(url = "tsneplot.html", file = file)
    }
  )
  
  output$help_text_tsne <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          t-SNE (t-distributed stochastic neighbor embedding) (Hinton, 2008) uses the local relationships between points 
          to create a low-dimensional mapping. This allows it to capture non-linear structure. </b>
          t-SNE creates a probability distribution using the Gaussian distribution that defines 
          the relationships between the points in high-dimensional space. t-SNE uses the Student 
          t-distribution to recreate the probability distribution in low-dimensional space. This 
          prevents the crowding problem, where points tend to get crowded in low-dimensional 
          space due to the curse of dimensionality. </b>
          t-SNE optimizes the embeddings directly using gradient descent. The cost function is 
          non-convex, thus there is the risk of getting stuck in local minima. t-SNE uses multiple 
          tricks to try to avoid this problem.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  #######   Random Forest    ########
  ###################################
  ###################################
  # data for random forest
  plotRF <- eventReactive(input$submit_rf, {
    rf.start <- Sys.time()
    rf_trans <- input$rf_trans
    type <- input$file_type
    num_trees <- input$num_trees
    num_clusters <- input$num_clusters
    
    if (type == "norm") {
      DS <- df_shiny()
    }
    else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (rf_trans == "None") {
      rf.data <- DS
    }
    else if (rf_trans == "log10") {
      rf.data <- log10(DS + 1)
    }
    rf.end <- Sys.time()
    print("Random forest plot time")
    print(rf.end - rf.start)
    hide("help_text_rf")
    return(list(rf.data, num_trees, num_clusters))
  })
  
  
  plotRAFSIL <- eventReactive(input$submit_rafsil, {
    rf.start <- Sys.time()
    rf_trans <- input$rf_trans
    type <- input$file_type
    
    if (type == "norm") {
      DS <- df_shiny()
    }
    else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (rf_trans == "None") {
      rf.data <- DS
    }
    else if (rf_trans == "log10") {
      rf.data <- log10(DS + 1)
    }
    f <- group_names()
    if (!is.null(f)) {
      meta_df <- data.frame("Column names" = colnames(DS), "Description" = f)
      meta_df
    }
    
    meta_df <- meta_df %>% remove_rownames %>% column_to_rownames(var="Column.names")
    meta_df$Description <- as.numeric(as.factor(meta_df$Description))
    
    meta_df <- as.matrix(meta_df)
    DS <- as.matrix(DS)
    
    rf.end <- Sys.time()
    print("RFSIL plot time")
    print(rf.end - rf.start)
    hide("help_text_rf")
    return(list(meta_df,DS))
  })
  
  rafsilplot <- function() {
    
    tryCatch({
      # get data
      t_list <- plotRAFSIL()
      ord = order(t_list[[1]]) ; t_list[[2]]=t_list[[2]][,ord] ; t_list[[1]] = t_list[[1]][ord] ; rm(ord)
      
      #- run RAFSIL1 with 50 forests
      res.r1 = RAFSIL(t(t_list[[2]]),nrep = 50, method="RAFSIL1")
      res.r2 = RAFSIL(t(t_list[[2]]),           method="RAFSIL2")
      
      #- retriev the dissimilarities
      dis.r1  = res.r1$D
      dis.r2  = res.r2$D
      dis.cor = sqrt((1 - cor(t_list[[2]],method="spearman"))/2)
      
      par(mfrow=c(1,2))
      par(mai=c(.1,.1,.5,.1))
      plotTSNE(dis.r1,labels=t_list[[1]],is_distance=FALSE,verbose=TRUE,perplexity=5)
      mtext("rafsil-1 / embedding", line=1)
      plotTSNE(dis.r2,labels=t_list[[1]],is_distance=FALSE,verbose=TRUE,perplexity=5)
      mtext("rafsil-2 / embedding", line=1)
    }, error = function(error_condition) {
      plot_exception("RAFSIL cannot be applied on this dataset.\nPlease use random forest clustering instead")
    }) 
    
  }
  
  ####################################################################################
  
  plot_exception <-function(
    ...,
    sep=" ",
    type=c("message","warning","cat","print"),
    color="auto",
    console=TRUE,
    size = 6){      
    type=match.arg(type)
    txt = paste(...,collapse=sep)
    if(console){
      if(type == "message") message(txt)
      if(type == "warning") warning(txt)
      if(type == "cat") cat(txt)
      if(type == "print") print(txt)
    }
    if(color =="auto") color <- if(type == "cat") "black" else "red"
    if(txt == "warning") txt <- paste("warning:",txt)
    print(ggplot2::ggplot() +
            ggplot2::geom_text(ggplot2::aes(x=0,y=0,label=txt),color=color,size=size) + 
            ggplot2::theme_void())
    invisible(NULL)
  }
  
  ####################################################################################
  
  
  rfplot <- function() {
    # get data
    li <- plotRF()
    rf.data <- li[[1]]
    num_trees <- li[[2]]
    num_clusters <- li[[3]]
    
    # unsupervised random forest on data
    print("Running random forest...")
    rf.data <- t(rf.data)
    rf_out <- randomForest(rf.data, type = unsupervised, ntree = num_trees, proximity = TRUE)
    print("Done!")
    
    mds_out <- cmdscale(1 - rf_out$proximity, eig = TRUE, k = 2)
    clusters_pam <- pam(1 - rf_out$proximity, k = num_clusters, diss = TRUE)
    shape_lvl <- c(1:num_clusters)
    shape_legend <- factor(clusters_pam$clustering, levels = shape_lvl)
    
    # print the proximity matrix
    print(table(clusters_pam$clustering, shape_legend))
    print(str(mds_out))
    print(mds_out$points)
    print(rownames(mds_out$points))
    
    # plot the graph
    df <- data.frame(x = mds_out$points[, 1], y = mds_out$points[, 2], color = shape_legend, shape = shape_legend)
    p <- ggplot(data = df, aes(x = x, y = y, color = color, shape = shape, text = paste("x: ", round(x, 4), "\n", "y: ", round(y, 4), "\n", "Name: ", rownames(mds_out$points), "\n", "Cluster: ", shape, sep = ""), group = 1)) +
      geom_point(size = 1.40) + theme_bw()
    p <- p + theme(legend.position = "none")
    p
    
    # add interactivity w/ plotly
    ggplotly(p, tooltip = c("text"))
  }
  
  
  rf_matrix <- reactive({
    li <- plotRF()
    rf.data <- li[[1]]
    num_trees <- li[[2]]
    num_clusters <- li[[3]]
    
    # unsupervised random forest on data
    print("Running random forest...")
    rf.data <- t(rf.data)
    rf_out <- randomForest(rf.data, type = unsupervised, ntree = num_trees, proximity = TRUE)
    print("Done!")
    
    mds_out <- cmdscale(1 - rf_out$proximity, eig = TRUE, k = 2)
    
    return(mds_out$points)
  })
  
  
  output$rf.plot <- renderPlotly({
    rfplot()
  })
  
  output$RAFSIL.plot <- renderPlot({
    rafsilplot()
  })
  
  output$rf.matrix <- renderTable({
    rf_matrix()
  },rownames=TRUE)
  
  output$help_text_rf <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Clustering belongs to unsupervised learning, in which each sample is clustered 
          into different classes, based on their similarity (usually based on Euclidean distance).
          Random forest algorithm is used to generate a proximity matrix - a rough estimate of the 
          distance between samples based on the proportion of times the samples end up in the same 
          leaf node of the decision tree. The proximity matrix is converted to a dist matrix which 
          is then input to the hierarchical clustering algorithm.
          </b>
        </p>
        <p>
          <b>
          Implementation adapted from - <a href ='https://nishanthu.github.io/articles/ClusteringUsingRandomForest.html'>
          https://nishanthu.github.io/articles/ClusteringUsingRandomForest.html</a>
          </b>
        </p>
      </center>
    ")
  })
  
  # output$downloadrfplot <- downloadHandler(
  #   filename = function(){
  #     paste("randomforestplot",".pdf",sep="")
  #   },
  #   content = function(file){
  #     pdf(file) 
  #     rfplot()
  #     dev.off()
  #   }
  # )
  
  # output$downloadrfmatrix <- downloadHandler(
  #   filename = function(){
  #     paste("randomforestmatrix",".pdf",sep="")
  #   },
  #   content = function(file){
  #     pdf(file) 
  #     rf_matrix()
  #     dev.off()
  #   }
  # )
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ############   SOM    #############
  ###################################
  ###################################
  # data for SOM
  plotSOM <- eventReactive(input$submit_som, {
    som.start <- Sys.time()
    som_trans <- input$som_trans
    sample_choice <- input$som_samples
    grid_h <- input$som_grid_h
    grid_v <- input$som_grid_v
    plot_type <- input$som_plot_type
    cluster_size <- input$som_cluster_size
    type <- input$file_type
    
    if (type == "norm") {
      DS <- df_shiny()
    }
    else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (som_trans == "None") {
      som.data <- DS
    }
    else if (som_trans == "log10") {
      som.data <- log10(DS + 1)
    }
    
    # Use all samples or individual
    if (sample_choice == "All") {
      som.data <- som.data
    }
    else {
      som.data <- som.data[, sample_choice]
    }
    
    # some parameters
    som.data <- as.matrix(som.data)
    set.seed(1)
    som_grid <- somgrid(xdim = grid_h, ydim = grid_v, topo = "hexagonal")
    som_model <- som(som.data, grid = som_grid)
    
    som.end <- Sys.time()
    print("SOM plot time")
    print(som.end - som.start)
    return(list(som_model, cluster_size))
  })
  
  
  sompropertyplot <- function() {
    # get data
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: property
    colors <- function(n, alpha = 'Set1') {
      rev(brewer.pal(n, alpha))
    }
    # use codes vectors (weight) for property plot
    #shinyjs::show("downloadProperty")
    hide("help_text_SOM")
    plot(som_model, type = "property", property = getCodes(som_model), main = "Property", palette.name = colors)
    
  }
  
  somcountplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: count
    colors <- function(n, alpha = 'Set2') {
      rev(brewer.pal(n, alpha))
    }
    #shinyjs::show("downloadCount")
    # show how many genes are mapped to each node
    plot(som_model, type = "count", main = "Count", palette.name = colors)
  }
  
  somcodesplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: codes
    # shows codebook vectors of genes
    #shinyjs::show("downloadCodes")
    plot(som_model, type = "codes", main = "Codes")
  }
  
  somdistplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: distance
    colors <- function(n, alpha = 'Set3') {
      rev(brewer.pal(n, alpha))
    }
    #shinyjs::show("downloadDistance")
    # show how close genes are from each other when they are mapped
    plot(som_model, type = "dist.neighbours", main = "Distance", palette.name = colors)
  }
  
  somclusterplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    cluster_size <- li[[2]]
    
    # plot type: cluster
    colors <- function(n, alpha = 1) {
      rev(heat.colors(n, alpha))
    }
    
    # define colors from RColorBrewer
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    # use hierarchical clustering to cluster the SOM
    som.hc <- cutree(hclust(object.distances(som_model, "codes")), cluster_size)
    #shinyjs::show("downloadCluster")
    plot(som_model, type = "mapping", bgcol = col_vector[som.hc], main = "Clusters")
    add.cluster.boundaries(som_model, som.hc)
  }
  
  
  output$som_property.plot <- renderPlot({
    sompropertyplot()
  })
  output$som_count.plot <- renderPlot({
    somcountplot()
  })
  output$som_codes.plot <- renderPlot({
    somcodesplot()
  })
  output$som_dist.plot <- renderPlot({
    somdistplot()
  })
  output$som_cluster.plot <- renderPlot({
    somclusterplot()
  })
  
  output$downloadProperty <- downloadHandler(
    filename = function(){
      paste("SOMProperty",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      sompropertyplot()
      dev.off()
    }
  )
  
  output$downloadCount <- downloadHandler(
    filename = function(){
      paste("SOMCount",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somcountplot()
      dev.off()
    }
  )
  
  output$downloadCodes <- downloadHandler(
    filename = function(){
      paste("SOMCodes",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somcodesplot()
      dev.off()
    }
  )
  
  output$downloadDistance <- downloadHandler(
    filename = function(){
      paste("SOMDistance",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somdistplot()
      dev.off()
    }
  )
  
  output$downloadCluster <- downloadHandler(
    filename = function(){
      paste("SOMCluster",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somclusterplot()
      dev.off()
    }
  )
  
  output$help_text_SOM <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          A self-organizing map (SOM) produces a two-dimensional, discretized representation 
          of the high-dimensional gene expression matrix, and is therefore a dimensionality 
          reduction technique. Self-organizing maps apply uses a neighborhood function to 
          preserve the topological properties of the input gene expression matrix.
          </b>
        </p>
        <p>
          <b>
          Each data point (1 sample) in the input gene expression matrix recognizes 
          themselves by competeting for representation. SOM mapping steps starts 
          from initializing the weight vectors.From there a sample vector is 
          selected randomly and the map of weight vectors is searched to find 
          which weight best represents that sample. Each weight vector has neighboring 
          weights that are close to it. The weight that is chosen is rewarded by being able 
          to become more like that randomly selected sample vector. The neighbors of that 
          weight are also rewarded by being able to become more like the chosen sample vector. 
          This allows the map to grow and form different shapes. Most generally, they form square/rectangular/hexagonal/L shapes in 2D feature space.
          </b>
        </p>
        <p>
          <b>
          Citation: <a href ='https://doi.org/10.1016/S0925-2312(98)00037-X'>https://doi.org/10.1016/S0925-2312(98)00037-X</a>
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ######## Gene-Set Analysis ########
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  ###### Complex Enrichment ########
  ###################################
  ###################################
  
  download_com_table <- reactiveVal(0)
  
  
  df_complex <- reactive({
    print("running...")
    if (is.null(input$file_complex_prot)&& is.null(input$text_complex_prot)) {
      return(NULL)
    }
    else if(!is.null(input$file_complex_prot)){
      parts <- strsplit(input$file_complex_prot$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_complex_prot$datapath)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      
      Acessions<-strsplit(input$text_complex_prot,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_complex_prot," ")
        
      }
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    
    return(Accessions)
    
  })
  
  
  df_com_table <- function(){
    
    gene_id <- df_com_id()
    
    corum_table <- data.frame()
    n <- 1
    
    for (id in gene_id) {
      check <- lookup(id, as.data.frame(up_corum_mapping))
      print(class(check))
      if(!is.na(check))
      {
        for (c_id in as.matrix(check)) {
          row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
          c_row <- data.frame(Uniprot_id = sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name),
                              Corum_id = c_id,
                              Complex_Name = as.character(allComplexes[paste0(c_id),"Complex_Name"]),
                              Complex_comment = allComplexes[paste0(c_id),"Complex_comment"],
                              row.names = n)
          corum_table <- rbind(corum_table, c_row)
          n = n + 1
        }
        
      } else {
        c_row <- data.frame(Uniprot_id = paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )"),
                            Corum_id = "No Match",
                            Complex_Name = "No Match",
                            Complex_comment = "No Match",
                            row.names = n)
        corum_table <- rbind(corum_table, c_row)
        n = n + 1
      }
    }
    download_com_table(corum_table)
    return(corum_table)
    
  }
  
  output$complex_table_prot <- shiny::renderDataTable({
    df_com_table()
  }, escape = FALSE)
  
  
  df_com_id <- eventReactive(input$submit_complex_prot, {
    hide("help_text_complex_en")
    df <- df_complex()
    return(df)
  })
  
  
  output$complex_download_prot <- downloadHandler(
    filename = function() {
      paste("complex", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_com_table(), file, row.names = FALSE)
    }
  )
  
  output$help_text_complex_en_prot <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page performs a Complex enrichment through the <a href ='http://mips.helmholtz-muenchen.de/corum/'>CORUM database</a> 
          of a given set of UniProt accessions and links the results to <a href ='http://UniProt.org'>UniProt.org</a>.
          </b>
        </p>
      </center>
    ")
  })
  
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ######## Protein Function #########
  ###################################
  ###################################
  
  
  download_prot_func <- reactiveVal(0)
  
  df_prot_func <- reactive({
    print("running...")
    if (is.null(input$file_prot_func) && is.null(input$text_prot_func)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_func)){
      parts <- strsplit(input$file_prot_func$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_prot_func$datapath)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      
      Acessions<-strsplit(input$text_prot_func,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_prot_func," ")
        
      }
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  df_func_table <- function() {
    
    Accessions <- df_func_id()
    
    print("fetching...")
    df <- GetProteinFunction(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Function" = df[,"Function..CC."]
    )
    download_prot_func(output_table)
    return(output_table)
    
  }
  
  output$prot_func_table <- DT::renderDataTable({
    df_func_table()
  }, escape = FALSE)
  
  df_func_id <- eventReactive(input$submit_prot_func, {
    hide("help_text_prot_fn")
    df <- df_prot_func()
    return(df)
  })
  
  output$prot_func_download <- downloadHandler(
    filename = function() {
      paste("Protein-Function", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_func(), file, row.names = FALSE)
    }
  )
  
  output$help_text_prot_fn <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Protein function information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  ###################################
  ###################################
  ####### Tissue Expression ########
  ###################################
  ###################################
  
  download_prot_expr <- reactiveVal(0)
  
  df_prot_expr <- reactive({
    print("running...")
    if (is.null(input$file_prot_expr)&& is.null(input$text_prot_expr)) {
      return(NULL)
    }else if(!is.null(input$file_prot_expr)){
      parts <- strsplit(input$file_prot_expr$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_prot_expr$datapath)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }else{
      Acessions<-strsplit(input$text_prot_expr," ")
      print(Acessions)
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        print(x)
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    
    return(Accessions)
    
  })
  
  df_expr_table <- function() {
    
    Accessions <- df_expr_id()
    print("fetching...")
    df <- GetExpression(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Tissue Specificity" = df[,"Tissue.specificity"]
    )
    
    download_prot_expr(output_table)                            
    return(output_table)
    
  }
  
  output$prot_expr_table <- DT::renderDataTable({
    df_expr_table()
  }, escape = FALSE)
  
  df_expr_id <- eventReactive(input$submit_prot_expr, {
    hide("help_text_prot_exp")
    df <- df_prot_expr()
    return(df)
  })
  
  output$prot_expr_download <- downloadHandler(
    filename = function() {
      paste("Protein-Expression", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_expr(), file, row.names = FALSE)
    }
  )
  
  output$help_text_prot_exp <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Tissue Expression information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  ###################################
  ###################################
  #### Subcellular Localization #####
  ###################################
  ###################################
  
  
  download_prot_local <- reactiveVal(0)
  
  df_prot_local <- reactive({
    print("running...")
    if (is.null(input$file_prot_local) && is.null(input$text_prot_local)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_local)){
      parts <- strsplit(input$file_prot_local$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_prot_local$datapath)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      Acessions<-strsplit(input$text_prot_local,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_prot_local," ")
        
      }
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  df_local_table <- function() {
    
    Accessions <- df_local_id()
    print("fetching...")
    df <- GetSubcellular_location(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Subcellular Location" = df[,"Subcellular.location..CC."]
    )
    download_prot_local(output_table)                            
    return(output_table)
    
  }
  
  output$prot_local_table <- DT::renderDataTable({
    df_local_table()
  }, escape = FALSE)
  
  df_local_id <- eventReactive(input$submit_prot_local, {
    hide("help_text_sub_loc")
    df <- df_prot_local()
    return(df)
  })
  
  output$prot_local_download <- downloadHandler(
    filename = function() {
      paste("Subcellular-Localization", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_local(), file, row.names = FALSE)
    }
  )
  
  output$help_text_sub_loc <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Subcellular LocalizationSubcellular Localization information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  #############Enrichment area################
  showGraphNodes.Enrch <- function(gene = "_gene")
  {
    shinyjs::show(id = paste0("loadStyleFile_path", gene))
    shinyjs::show(id = paste0("overlap_min", gene))
    shinyjs::show(id = paste0("doLayout_path", gene))
    shinyjs::show(id = paste0("sfn_path", gene))
    shinyjs::show(id = paste0("fit_path", gene))
    shinyjs::show(id = paste0("fitSelected_path", gene))
    shinyjs::show(id = paste0("clearSelection_path", gene))
    shinyjs::show(id = paste0("removeGraphButton_path", gene))
    shinyjs::show(id = paste0("addRandomGraphFromDataFramesButton_path", gene))
    shinyjs::show(id = paste0("getSelectedNodes_path", gene))
    shinyjs::show(id = paste0("selectedNodesDisplay_path", gene))
  }
  
  
  observe(
    if (input$path_enri_tab_prot == "Visualization")
      showGraphNodes.Enrch(gene = "_prot")
  )
  observe(
    if (input$path_enri_tab_gene == "Visualization")
      showGraphNodes.Enrch(gene = "_gene")
  )
  ############################################
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  ###################################
  ###################################
  ######## Protein Domains ##########
  ###################################
  ###################################
  
  
  download_prot_domain <- reactiveVal(0)
  
  df_prot_domain <- reactive({
    print("running...")
    if (is.null(input$file_prot_domain)&& is.null(input$text_prot_domain)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_domain)){
      parts <- strsplit(input$file_prot_domain$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_prot_domain$datapath)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      
      Acessions<-strsplit(input$text_prot_domain,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_prot_domain," ")
        
      }
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  df_domain_table <- function() {
    
    
    Accessions <- df_domain_id()
    print("fetching...")
    df <- GetFamily_Domains(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Protein Families" = df[,"Protein.families"],
      "Protein Domain" = df[,"Domain..FT."]
    )
    download_prot_domain(output_table)                           
    return(output_table)
    
  }
  
  output$prot_domain_table <- DT::renderDataTable({
    df_domain_table()
  }, escape = FALSE)
  
  df_domain_id <- eventReactive(input$submit_prot_domain, {
    hide("help_text_pro_dom")
    df <- df_prot_domain()
    return(df)
  })
  
  output$prot_domain_download <- downloadHandler(
    filename = function() {
      paste("Protein-Domains", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_domain(), file, row.names = FALSE)
    }
  )
  
  output$help_text_pro_dom <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Protein Domains information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ###### Pathways Enrichment ########
  ###################################
  ###################################
  
  
  pathway_enri_df <- reactiveVal(0)
  pathway_enri_nodes <- reactiveVal(0)
  
  df_path_enri_gene <- reactive({
    print("running pathway gene...")
    if (is.null(input$file_path_enri_gene) && is.null(input$text_path_enri_gene) ) {
      return(NULL)
    }else if(!is.null(input$file_path_enri_gene)){
      parts <- strsplit(input$file_path_enri_gene$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_path_enri_gene$datapath)
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      
      Acessions<-strsplit(input$text_path_enri_gene," ")
      print(Acessions)
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        print(x)
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    
    return(Accessions)
  })
  
  
  df_path_enri_prot <- reactive({
    print("running...")
    if (is.null(input$file_path_enri_prot)&& is.null(input$text_path_enri_prot)) {
      return(NULL)
    }
    else if(!is.null(input$file_path_enri_prot)){
      parts <- strsplit(input$file_path_enri_prot$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_path_enri_prot$datapath)
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      Acessions<-strsplit(input$text_path_enri_prot,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_path_enri_prot," ")
        
      }
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        print(x)
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  
  df_path_enri_id_gene <- eventReactive(input$submit_path_enri_gene,{
    print("running")
    
    hide("help_text_path_enri")
    df <- df_path_enri_gene()
    
    return(df)
  })
  
  df_path_enri_id_prot <- eventReactive(input$submit_path_enri_prot,{
    print("running")
    
    hide("help_text_path_enri")
    df <- df_path_enri_prot()
    return(df)
  })
  
  
  output$help_text_path_enri_gene <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page performs Pathways Enrichment from for a given set of genes using
          <a href ='https://biit.cs.ut.ee/gprofiler/gost'>g:Profiler</a>.
          </b>
        </p>
      </center>
    ")
  })
  
  output$help_text_path_enri_prot <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page performs Pathways Enrichment from for a given set of genes using
          <a href ='https://biit.cs.ut.ee/gprofiler/gost'>g:Profiler</a>.
          </b>
        </p>
      </center>
    ")
  })
  
  
  plot_path_enri_gene <- function() {
    df_path_enri_id_gene()
    gene_name <- as.data.frame(df_path_enri_id_gene())
    gene_name[,1] <- as.character(gene_name[,1])
    
    path_list <- gost(gene_name[,1],exclude_iea = TRUE,evcodes = TRUE ,sources = "GO:BP")
    path_df <- path_list[[1]]
    pathway_enri_df(path_df)
    prot_num <- data.frame()
    for(i in 1:nrow(path_df))
    {
      prot_num <- rbind(prot_num,nrow(as.data.frame(strsplit(path_df[i,"intersection"],","))))
    }
    
    path_enrich_df <- data.frame(
      "term_name" = path_df[,"term_name"],
      "intersection" = prot_num[,1]
    )
    
    pathway_enri_nodes(path_enrich_df)
    
    path_enrich_df <- path_enrich_df[order(path_enrich_df$intersection),]
    
    
    bar_plot <- ggplot(data=path_enrich_df, aes(x=reorder(path_enrich_df$term_name , path_enrich_df$intersection), y=path_enrich_df$intersection)) +
      geom_bar(stat="identity", fill="steelblue" , alpha = 0.7) + xlab("Molecular function") + ylab("Number of Genes") +
      geom_text(aes(label = path_enrich_df$intersection), vjust = -0.03) + theme(axis.text.x = element_text(angle = 90 , hjust = 1 , vjust = 0.2))+
      theme_minimal() +coord_flip() + theme_bw()+theme(text = element_text(size=12, face="bold", colour="black"),axis.text.x = element_text(vjust=2))
    
    
    return(bar_plot)
  }
  
  plot_path_enri_prot <- function() {
    df_path_enri_id_prot()
    gene_name <- as.data.frame(df_path_enri_id_prot())
    gene_name[,1] <- as.character(gene_name[,1])
    
    path_list <- gost(gene_name[,1],exclude_iea = TRUE,evcodes = TRUE ,sources = "GO:BP")
    path_df <- path_list[[1]]
    pathway_enri_df(path_df)
    prot_num <- data.frame()
    for(i in 1:nrow(path_df))
    {
      prot_num <- rbind(prot_num,nrow(as.data.frame(strsplit(path_df[i,"intersection"],","))))
    }
    
    path_enrich_df <- data.frame(
      "term_name" = path_df[,"term_name"],
      "intersection" = prot_num[,1]
    )
    
    pathway_enri_nodes(path_enrich_df)
    
    path_enrich_df <- path_enrich_df[order(path_enrich_df$intersection),]
    
    
    bar_plot <- ggplot(data=path_enrich_df, aes(x=reorder(path_enrich_df$term_name , path_enrich_df$intersection), y=path_enrich_df$intersection)) +
      geom_bar(stat="identity", fill="steelblue" , alpha = 0.7) + xlab("Molecular function") + ylab("Number of Genes") +
      geom_text(aes(label = path_enrich_df$intersection), vjust = -0.03) + theme(axis.text.x = element_text(angle = 90 , hjust = 1 , vjust = 0.2))+
      theme_minimal() +coord_flip() + theme_bw()+theme(text = element_text(size=12, face="bold", colour="black"),axis.text.x = element_text(vjust=2))
    
    
    return(bar_plot)
  }
  
  
  output$path_enri.plot_gene <- renderPlotly({
    df_path_enri_id_gene()
    gene_name <- as.data.frame(df_path_enri_id_gene())
    gene_name[,1] <- as.character(gene_name[,1])
    
    ggplotly(Pathway.Enr(gene_name[,1]), tooltip = c("text"))
  })
  
  output$path_enri.plot_prot <- renderPlotly({
    df_path_enri_id_prot()
    gene_name <- as.data.frame(df_path_enri_id_prot())
    gene_name[,1] <- as.character(gene_name[,1])
    
    ggplotly(Pathway.Enr(gene_name[,1]), tooltip = c("text"))
  })
  
  #visualization
  
  observeEvent(input$fit_path_gene, ignoreInit=TRUE, {
    fit(session, 80)
  })
  
  observeEvent(input$fit_path_prot, ignoreInit=TRUE, {
    fit(session, 80)
  })
  
  
  observeEvent(input$showCondition_gene, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_gene)
    values <- as.numeric(pathway_enri_nodes()[,2])
    node.names <- pathway_enri_nodes()[,1]
    print(values)
    setNodeAttributes(session, attributeName="lfc", nodes=node.names, values)
  })
  
  
  observeEvent(input$showCondition_prot, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_prot)
    values <- as.numeric(pathway_enri_nodes()[,2])
    node.names <- pathway_enri_nodes()[,1]
    print(values)
    setNodeAttributes(session, attributeName="lfc", nodes=node.names, values)
  })
  
  
  observeEvent(input$loadStyleFile_path_gene,  ignoreInit=TRUE, {
    if(input$loadStyleFile_path != ""){
      tryCatch({
        loadStyleFile(input$loadStyleFile_path_gene)
      }, error=function(e) {
        msg <- sprintf("ERROR in stylesheet file '%s': %s", input$loadStyleFile_path_gene, e$message)
        showNotification(msg, duration=NULL, type="error")
      })
      later(function() {updateSelectInput(session, "loadStyleFile", selected=character(0))}, 0.5)
    }
  })
  
  
  
  observeEvent(input$loadStyleFile_path_prot,  ignoreInit=TRUE, {
    if(input$loadStyleFile_path_prot != ""){
      tryCatch({
        loadStyleFile(input$loadStyleFile_path_prot)
      }, error=function(e) {
        msg <- sprintf("ERROR in stylesheet file '%s': %s", input$loadStyleFile_path_prot, e$message)
        showNotification(msg, duration=NULL, type="error")
      })
      later(function() {updateSelectInput(session, "loadStyleFile", selected=character(0))}, 0.5)
    }
  })
  
  
  
  observeEvent(input$doLayout_path_gene,  ignoreInit=TRUE,{
    if(input$doLayout_path_gene != ""){
      strategy <- input$doLayout_path_gene
      doLayout(session, strategy)
      later(function() {updateSelectInput(session, "doLayout", selected=character(0))}, 1)
    }
  })
  
  observeEvent(input$doLayout_path_prot,  ignoreInit=TRUE,{
    if(input$doLayout_path_prot != ""){
      strategy <- input$doLayout_path_prot
      doLayout(session, strategy)
      later(function() {updateSelectInput(session, "doLayout", selected=character(0))}, 1)
    }
  })
  
  
  
  observeEvent(input$sfn_path_gene,  ignoreInit=TRUE,{
    selectFirstNeighbors(session)
  })
  
  observeEvent(input$sfn_path_prot,  ignoreInit=TRUE,{
    selectFirstNeighbors(session)
  })
  
  
  
  observeEvent(input$fitSelected_path_gene,  ignoreInit=TRUE,{
    fitSelected(session, 100)
  })
  
  observeEvent(input$fitSelected_path_prot,  ignoreInit=TRUE,{
    fitSelected(session, 100)
  })
  
  
  
  observeEvent(input$getSelectedNodes_path_gene, ignoreInit=TRUE, {
    output$selectedNodesDisplay_path_gene <- renderText({" "})
    getSelectedNodes(session)
  })
  
  observeEvent(input$getSelectedNodes_path_prot, ignoreInit=TRUE, {
    output$selectedNodesDisplay_path_prot <- renderText({" "})
    getSelectedNodes(session)
  })
  
  
  
  observeEvent(input$clearSelection_path_gene,  ignoreInit=TRUE, {
    clearSelection(session)
  })  
  
  observeEvent(input$clearSelection_path_prot,  ignoreInit=TRUE, {
    clearSelection(session)
  })  
  
  
  
  observeEvent(input$removeGraphButton_path_gene, ignoreInit=TRUE, {
    removeGraph(session)
  })
  
  observeEvent(input$removeGraphButton_path_prot, ignoreInit=TRUE, {
    removeGraph(session)
  })
  
  
  observeEvent(input$addRandomGraphFromDataFramesButton_path_gene, ignoreInit=TRUE, {
    source.nodes <-  LETTERS[sample(1:5, 5)]
    target.nodes <-  LETTERS[sample(1:5, 5)]
    tbl.edges <- data.frame(source=source.nodes,
                            target=target.nodes,
                            interaction=rep("generic", length(source.nodes)),
                            stringsAsFactors=FALSE)
    all.nodes <- sort(unique(c(source.nodes, target.nodes, "orphan")))
    tbl.nodes <- data.frame(id=all.nodes,
                            type=rep("unspecified", length(all.nodes)),
                            stringsAsFactors=FALSE)
    addGraphFromDataFrame(session, tbl.edges, tbl.nodes)
  })
  
  observeEvent(input$addRandomGraphFromDataFramesButton_path_prot, ignoreInit=TRUE, {
    source.nodes <-  LETTERS[sample(1:5, 5)]
    target.nodes <-  LETTERS[sample(1:5, 5)]
    tbl.edges <- data.frame(source=source.nodes,
                            target=target.nodes,
                            interaction=rep("generic", length(source.nodes)),
                            stringsAsFactors=FALSE)
    all.nodes <- sort(unique(c(source.nodes, target.nodes, "orphan")))
    tbl.nodes <- data.frame(id=all.nodes,
                            type=rep("unspecified", length(all.nodes)),
                            stringsAsFactors=FALSE)
    addGraphFromDataFrame(session, tbl.edges, tbl.nodes)
  })
  
  
  
  # observeEvent(input$selectedNodes, {
  #       newNodes <- input$selectedNodes;
  #       output$selectedNodesDisplay <- renderText({
  #          paste(newNodes)
  #          })
  #       })
  
  pathway_overlap <- reactiveVal(0)
  
  new_source_var <- reactiveVal(0)
  new_target_var <- reactiveVal(0)
  overlap_wt <- reactiveVal(0)
  
  output$path_enri_visu_gene <- renderCyjShiny({
    
    print("visualization")
    df_path_enri_id_gene()
    
    Enrich <- gost(df_path_enri_id_gene()[,1],evcodes = T, sources = c('KEGG', 'REAC'))
    
    Pathway <- Construct.COPathway(Enrich, input$overlap_min_gene)
    
    
    nodes_tot <- c(unique(Pathway[,1],unique(Pathway[,2])))
    
    
    path_enri.nodes <- data.frame(id=nodes_tot,
                                  type=nodes_tot,
                                  stringsAsFactors=FALSE)
    
    path_enri.edges <- data.frame(source=Pathway[,1],
                                  target=Pathway[,2],
                                  interaction=Pathway[,1],
                                  stringsAsFactors=FALSE)
    
    graph.json <- dataFramesToJSON(path_enri.edges, path_enri.nodes)
    cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
    
  })
  
  Construct.COPathway <- function(EnrichmentObject, threshold = 1)
  {
    
    PathwayNetwork <- data.frame()
    PathwayDF <- EnrichmentObject[["result"]] 
    for (i in 1:nrow(PathwayDF))
    {
      if (dim(PathwayDF)[1] == 1)
        break
      Pathway_accessions <- strsplit(PathwayDF$intersection, ",") 
      for (accession in Pathway_accessions)
      {
        inx <- which(grepl(accession, PathwayDF$intersection) == T)
        inx <- inx[-i] 
        if (length(inx) > 1)
        {
          Source <- rep(PathwayDF$term_name[i], length(inx))
          Target <- PathwayDF$term_name[c(inx)]
          PathwayNetwork <- rbind(PathwayNetwork , cbind(Source, Target, accession))
        }
      }
      PathwayDF <- PathwayDF[-i,]
      CoEnrichment <- setNames(aggregate(PathwayNetwork$accession, by = list(PathwayNetwork$Source, PathwayNetwork$Target),
                                         paste, collapse=","), c("Source", "Target", "Accesion"))
      
      CoEnrichment$ProteinCount <- str_count(CoEnrichment$Accesion, ",")
      CoEnrichment <- CoEnrichment[CoEnrichment$ProteinCount >= threshold,]
    }
    return(CoEnrichment)
  }
  
  output$path_enri_visu_prot <- renderCyjShiny({
    
    print("visualization")
    df_path_enri_id_prot()
    Enrich <- gost(df_path_enri_id_prot()[,1],evcodes = T, sources = c('KEGG', 'REAC'))
    Pathway <- Construct.COPathway(Enrich, input$overlap_min_prot)
    nodes_tot <- c(unique(Pathway[,1],unique(Pathway[,2])))
    
    
    path_enri.nodes <- data.frame(id=nodes_tot,
                                  type=nodes_tot,
                                  stringsAsFactors=FALSE)
    
    path_enri.edges <- data.frame(source=Pathway[,1],
                                  target=Pathway[,2],
                                  interaction=Pathway[,1],
                                  stringsAsFactors=FALSE)
    
    graph.json <- dataFramesToJSON(path_enri.edges, path_enri.nodes)
    cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
    
  })
  
  observeEvent(input$edge_wt_gene, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_gene)
    
    print(overlap_wt())
    setEdgeAttributes(session, attributeName="wt", sourceNodes=new_source_var(),
                      targetNodes=new_target_var(),
                      interactions=new_target_var(),
                      values=overlap_wt())
  })
  
  observeEvent(input$edge_wt_prot, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_prot)
    
    print(overlap_wt())
    setEdgeAttributes(session, attributeName="wt", sourceNodes=new_source_var(),
                      targetNodes=new_target_var(),
                      interactions=new_target_var(),
                      values=overlap_wt())
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  output$help_text_protein_set <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Please Upload or enter the Uniprot Accession Numbers.
          </b>
        </p>
      </center>
    ")
  })
  observeEvent(input$submit_protein_set,{
    if (is.null(input$file_protein_set) && is.null(input$text_protein_set)) {
      return(NULL)
    }else if(!is.null(input$file_protein_set)){
      parts <- strsplit(input$file_protein_set$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_protein_set$datapath)
      Accessions <- na.omit(Accessions)[,1]
      Accessions <- unique(Accessions)
      Accessions <- trimws(Accessions)
    }else{
      Accessions<-input$text_protein_set
      print(Accessions)
    }
    
    updateTextInput(session, "text_uniprot", value = paste(Accessions))
    updateTextInput(session, "text_prot_Int", value = paste(Accessions))
    updateTextInput(session, "text_prot_func", value = paste(Accessions))
    updateTextInput(session, "text_prot_local", value = paste(Accessions))
    updateTextInput(session, "text_prot_domain", value = paste(Accessions))
    updateTextInput(session, "text_prot_seq", value = paste(Accessions))
    updateTextInput(session, "text_prot_seq_evol", value = paste(Accessions))
    updateTextInput(session, "text_prot_seq_Patho", value = paste(Accessions))
    updateTextInput(session, "text_complex_prot", value = paste(Accessions))
    
    output$help_text_protein_set <- renderUI({
      HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Data Upload Complete. Please proceed with the analysis.
          </b>
        </p>
      </center>
    ")
    })
  })
  
  
  ###################################
  ###################################
  ##########   Uniprot    ###########
  ###################################
  ###################################
  
  df_uniprot <- reactive({
    print("running")
    if (is.null(input$file_uniprot) && is.null(input$text_uniprot)) {
      return(NULL)
    }
    else if(!is.null(input$file_uniprot)){
      parts <- strsplit(input$file_uniprot$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_uniprot$datapath)
      Accessions <- na.omit(Accessions)[,1]
      Accessions <- unique(Accessions)
      Accessions <- trimws(Accessions)
      print(Accessions)
    }else{
      Acessions<-strsplit(input$text_uniprot,",")
      print(length(Acessions[[1]]))
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_uniprot," ")
        
      }
      
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  plotUniprot <-  eventReactive(input$submit_uniprot, {
    
    Accessions <- df_uniprot()
    print(Accessions)
    hide("help_text_bio_pr")
    #print("Please Wait... Fetching Taxa Object. It may take a while")
    #TaxaObj <- GetNamesTaxa(Accessions)
    print("Please Wait... Fetching Gene Ontology Object. It may take a while")
    GeneOntologyObj <- GetProteinGOInfo(Accessions) 
    print("Done") 
    return(GeneOntologyObj)
  })
  
  output$help_text_bio_pr <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Gene Ontology (GO) terms 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  plotCE <- function() {
    # get data
    GO_df <- plotUniprot()
    return(Plot.GOSubCellular(GO_df,20))
    # ggplotly(bar_plot, tooltip = c("text"))
  }
  
  output$download_cell_plot <- downloadHandler(
    filename = function(){paste("Cellular-Component",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plotCE())
    }
  )
  
  plotBIO <- function() {
    # get data
    GO_df <- plotUniprot()
    return(PlotGOBiological(GO_df,20))
    # ggplotly(bar_plot, tooltip = c("text"))
  }
  
  output$download_bio_plot <- downloadHandler(
    filename = function(){paste("Biological-Process",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plotBIO())
    }
  )
  
  plotMol <- function() {
    # get data
    GO_df <- plotUniprot()
    return(Plot.GOMolecular(GO_df, 20))
    # ggplotly(bar_plot, tooltip = c("text"))
  }
  
  output$download_mole_plot <- downloadHandler(
    filename = function(){paste("Molecular-Function",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plotMol())
    }
  )
  
  output$uniprot_celplot <- renderPlot({
    GO_df <- plotUniprot()
    Plot.GOSubCellular(GO_df,20)
    ##ggplotly(plotCE(), tooltip = c("text"))
  })
  
  output$uniprotbioplot <- renderPlot({
    #ggplotly(plotBIO(), tooltip = c("text"))
    GO_df <- plotUniprot()
    PlotGOBiological(GO_df,20)
  })
  
  output$uniprot_molcplot <- renderPlot({
    GO_df <- plotUniprot()
    Plot.GOMolecular(GO_df,20)
    #ggplotly(plotMol(), tooltip = c("text"))
  })
  
  download_cel_table <- NULL
  
  output$uniprot_celtable <- shiny::renderDataTable({
    
    GO_df <- plotUniprot()
    CellularDF <- Goparse(GO_df, 5)
    CellularDF <- na.omit(CellularDF)
    download_cel_table <- CellularDF
    CellularDF
    
    
  })
  
  output$download_cell_comp <- downloadHandler(
    filename = function() {
      paste("Cellular-Component", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_cel_table, file, row.names = FALSE)
    }
  )
  
  download_bio_table <- NULL
  
  output$uniprot_biotable <- shiny::renderDataTable({
    
    GO_df <- plotUniprot()
    BiologicalDF <- Goparse(GO_df, 3)
    BiologicalDF <- na.omit(BiologicalDF)
    download_bio_table <- BiologicalDF
    BiologicalDF
  })
  
  output$download_bio_pro <- downloadHandler(
    filename = function() {
      paste("Biological-Process", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_bio_table, file, row.names = FALSE)
    }
  )
  
  download_mol_table <- NULL
  
  output$uniprot_molctable <- shiny::renderDataTable({
    
    GO_df <- plotUniprot()
    MolecularDF <- Goparse(GO_df, 4)
    MolecularDF <- na.omit(MolecularDF)
    download_mol_table <- MolecularDF
    MolecularDF
    
  })
  
  output$download_mole_func <- downloadHandler(
    filename = function() {
      paste("Molecular-Function", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_mol_table, file, row.names = FALSE)
    }
  )
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  #######  P-P Interactions   #######
  ###################################
  ###################################
  
  
  ############ Initializing Variables ###############
  
  df_interaction <- reactiveVal(0)
  df_names <- reactiveVal(0)
  
  # tbl.nodes <- data.frame(id=c("A", "B", "C"),
  #                       type=c("kinase", "TF", "glycoprotein"),
  #                       lfc=c(1, 1, 1),
  #                       count=c(0, 0, 0),
  #                       stringsAsFactors=FALSE)
  
  # tbl.edges <- data.frame(source=c("A", "B", "C"),
  #                       target=c("B", "C", "A"),
  #                       interaction=c("phosphorylates", "synthetic lethal", "unknown"),
  #                       stringsAsFactors=FALSE)
  
  ####################################################
  
  df_prot_Int <- reactive({
    print("running")
    if (is.null(input$file_prot_Int)&& is.null(input$text_prot_Int)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_Int)){
      parts <- strsplit(input$file_prot_Int$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_prot_Int$datapath)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }else{
      Acessions<-strsplit(input$text_prot_Int,",")
      print(length(Acessions[[1]]))
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_prot_Int," ")
        
      }
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    
    return(Accessions)
    
  })
  
  
  observeEvent(input$fit, ignoreInit=TRUE, {
    fit(session, 80)
  })
  
  
  observeEvent(input$loadStyleFile,  ignoreInit=TRUE, {
    if(input$loadStyleFile != ""){
      tryCatch({
        loadStyleFile(input$loadStyleFile)
      }, error=function(e) {
        msg <- sprintf("ERROR in stylesheet file '%s': %s", input$loadStyleFile, e$message)
        showNotification(msg, duration=NULL, type="error")
      })
      later(function() {updateSelectInput(session, "loadStyleFile", selected=character(0))}, 0.5)
    }
  })
  
  
  observeEvent(input$doLayout,  ignoreInit=TRUE,{
    if(input$doLayout != ""){
      strategy <- input$doLayout
      doLayout(session, strategy)
      later(function() {updateSelectInput(session, "doLayout", selected=character(0))}, 1)
    }
  })
  
  
  observeEvent(input$selectName,  ignoreInit=TRUE,{
    selectNodes(session, input$selectName)
  })
  
  
  observeEvent(input$sfn,  ignoreInit=TRUE,{
    selectFirstNeighbors(session)
  })
  
  
  observeEvent(input$fitSelected,  ignoreInit=TRUE,{
    fitSelected(session, 100)
  })
  
  
  observeEvent(input$getSelectedNodes, ignoreInit=TRUE, {
    output$selectedNodesDisplay <- renderText({" "})
    getSelectedNodes(session)
  })
  
  
  observeEvent(input$clearSelection,  ignoreInit=TRUE, {
    clearSelection(session)
  })  
  
  
  observeEvent(input$removeGraphButton, ignoreInit=TRUE, {
    removeGraph(session)
  })
  
  
  observeEvent(input$addRandomGraphFromDataFramesButton, ignoreInit=TRUE, {
    source.nodes <-  LETTERS[sample(1:5, 5)]
    target.nodes <-  LETTERS[sample(1:5, 5)]
    tbl.edges <- data.frame(source=source.nodes,
                            target=target.nodes,
                            interaction=rep("generic", length(source.nodes)),
                            stringsAsFactors=FALSE)
    all.nodes <- sort(unique(c(source.nodes, target.nodes, "orphan")))
    tbl.nodes <- data.frame(id=all.nodes,
                            type=rep("unspecified", length(all.nodes)),
                            stringsAsFactors=FALSE)
    addGraphFromDataFrame(session, tbl.edges, tbl.nodes)
  })
  
  
  observeEvent(input$selectedNodes, {
    newNodes <- input$selectedNodes;
    output$selectedNodesDisplay <- renderText({
      paste(newNodes)
    })
  })
  
  
  output$cyjShiny <- renderCyjShiny({
    print(" renderCyjShiny invoked")
    print("graph.json:")
    
    
    print("running...")
    
    
    # tryCatch({
    
    Accessions <- df_prot_int_id()
    print("Please Wait... Fetching interaction data. It may take a while")
    protein_interaction_df <- getInteraction(Accessions)
    df_interaction(protein_interaction_df)
    print("Fetched...")
    
    #migrating rowId to first colunm 
    # protein_interaction_df <- cbind(ID = rownames(protein_interaction_df),protein_interaction_df)
    # rownames(protein_interaction_df) <- 1:nrow(protein_interaction_df)
    
    #making nodes
    nodes <- as.character(protein_interaction_df[,1])
    for (i in 1:nrow(protein_interaction_df))
    {
      if(!(is.na(protein_interaction_df[i,2])))
      {
        data_df <- strsplit(as.character(protein_interaction_df[i,2]),"; ")
        for(j in data_df)
        {
          nodes <- c(nodes,j)
        }
      }
    }
    
    print(nodes)
    
    print("Please Wait... Fetching Gene Names. It may take a while")
    protein_gene_name <- getGeneNames(nodes)
    df_names(protein_gene_name)
    print("........................")
    print(as.character(protein_gene_name[,1]))
    print("Fetched...")
    edge_source <- character()
    edge_target <- character()
    
    for (i in 1:nrow(protein_interaction_df))
    {
      if(!(is.na(protein_interaction_df[i,2])))
      {
        data_df <- strsplit(as.character(protein_interaction_df[i,2]),"; ")
        for(j in data_df)
        {
          edge_source <- c(edge_source,rep(as.character(protein_gene_name[as.character(protein_interaction_df[i,1]),1]),length(j)))
          print(as.character(protein_gene_name[j,1]))
          edge_target <- c(edge_target,as.character(protein_gene_name[j,1]))
        }
      }
    }
    
    tbl.nodes <- data.frame(id=as.character(protein_gene_name[,1]),
                            type=as.character(protein_gene_name[,1]),
                            stringsAsFactors=FALSE)
    
    
    tbl.edges <- data.frame(source=edge_source,
                            target=edge_target,
                            interaction=edge_target,
                            stringsAsFactors=FALSE)
    
    # }, error = function(error_condition) {
    #   print("using defauslt value")
    # })
    
    graph.json <- dataFramesToJSON(tbl.edges, tbl.nodes)
    
    print(fromJSON(graph.json))
    cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
  })
  
  
  # observeEvent(input$submit_prot_Int, {
  
  #   print("running...")
  #   Accessions <- df_prot_Int()
  #   print("Please Wait... Fetching interaction data. It may take a while")
  #   protein_interaction_df <- getInteraction(Accessions)
  #   df_interaction(protein_interaction_df)
  #   print("Fetched...")
  
  #migrating rowId to first colunm 
  # protein_interaction_df <- cbind(ID = rownames(protein_interaction_df),protein_interaction_df)
  # rownames(protein_interaction_df) <- 1:nrow(protein_interaction_df)
  
  #making nodes
  # nodes <- as.character(protein_interaction_df[,1])
  # for (i in 1:nrow(protein_interaction_df))
  # {
  #   if(!(is.na(protein_interaction_df[i,2])))
  #   {
  #     data_df <- strsplit(protein_interaction_df[i,2],"; ")
  #     for(j in data_df)
  #     {
  #       nodes <- c(nodes,j)
  #     }
  #   }
  # }
  
  # print("Please Wait... Fetching Gene Names. It may take a while")
  # protein_gene_name <- getGeneNames(nodes)
  # df_names(protein_gene_name)
  # print("Fetched...")
  
  # # print("Rendering Visualization using Cytoscape")
  
  # # g <- graphNEL(as.character(protein_gene_name[,1]), edgemode="undirected")
  # for (i in 1:nrow(protein_interaction_df))
  # {
  #   if(!(is.na(protein_interaction_df[i,2])))
  #   {
  #     data_df <- strsplit(protein_interaction_df[i,2],"; ")
  #     for(j in data_df)
  #     {
  #       g <- graph::addEdge(as.character(protein_gene_name[as.character(protein_interaction_df[i,1]),1]), as.character(protein_gene_name[j,1]), g)
  #     }
  #   }
  # }
  
  # nodeDataDefaults(g, attr="label") <- "undefined"
  # nodeDataDefaults(g, attr="type") <- "undefined"
  # nodeDataDefaults(g, attr="flux") <- 0
  # edgeDataDefaults(g, attr="edgeType") <- "undefined"
  
  # rcy <- RCyjs(title="RCyjs vignette")
  # setGraph(rcy, g)
  # print("set")
  # strategies <- getLayoutStrategies(rcy)
  # print(strategies)
  
  # RCyjs::layout(rcy, "grid")
  # # print("lay")
  # fit(rcy, padding=200)
  # print("fit")
  # setDefaultStyle(rcy)
  
  # })
  
  df_prot_int_id <- eventReactive(input$submit_prot_Int, {
    hide("help_text_p_inte")
    Accessions <- df_prot_Int()
    return(Accessions)
  })
  
  output$help_text_p_inte <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Protein-Protein Interactions 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  getInteraction <- function(ProteinAccList) {
    
    if(!has_internet())
    {
      message("Please connect to the internet as the package requires internect connection.")
      return()
    }
    protein_interaction_df = data.frame()
    baseUrl <- "http://www.uniprot.org/uniprot/"
    Colnames = "interactor"
    for (ProteinAcc in ProteinAccList)
    {
      #to see if Request == 200 or not
      Request <- tryCatch(
        {
          GET(paste0(baseUrl , ProteinAcc,".xml") , timeout(60))
        },error = function(cond)
        {
          message("Internet connection problem occurs and the function will return the original error")
          message(cond)
        }
      )
      #this link return information in tab formate (format = tab)
      ProteinName_url <- paste0("?query=accession:",ProteinAcc,"&format=tab&columns=",Colnames)
      RequestUrl <- paste0(baseUrl , ProteinName_url)
      RequestUrl <- URLencode(RequestUrl)
      
      if (Request$status_code == 200){
        # parse the information in DataFrame
        ProteinDataTable <- tryCatch(read.csv(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
        if (!is.null(ProteinDataTable))
        {
          ProteinDataTable <- ProteinDataTable[1,]
          ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
          # add Dataframes together if more than one accession
          protein_interaction_df <- rbind(protein_interaction_df, ProteinInfoParsed)
          print(paste0(ProteinAcc," interactions Fetched.."))
        }
      }else {
        HandleBadRequests(Request$status_code)
      }
    }
    
    protein_interaction_df <- cbind(ID = rownames(protein_interaction_df),protein_interaction_df)
    rownames(protein_interaction_df) <- 1:nrow(protein_interaction_df)
    
    return(protein_interaction_df)
    
  }
  
  getGeneNames <- function(ProteinAccList) {
    
    # baseUrl <- "http://www.uniprot.org/uniprot/"
    # Colnames = "genes(PREFERRED)"
    
    # protein_gene_name = data.frame()
    # for (ProteinAcc in ProteinAccList)
    # {
    #   #to see if Request == 200 or not
    #   Request <- tryCatch(
    #     {
    #       GET(paste0(baseUrl , ProteinAcc,".xml") , timeout(10))
    #     },error = function(cond)
    #     {
    #       message("Internet connection problem occurs and the function will return the original error")
    #       message(cond)
    #     }
    #   ) 
    #   #this link return information in tab formate (format = tab)
    #   ProteinName_url <- paste0("?query=accession:",ProteinAcc,"&format=tab&columns=",Colnames)
    #   RequestUrl <- paste0(baseUrl , ProteinName_url)
    #   RequestUrl <- URLencode(RequestUrl)
    #   if (Request$status_code == 200){
    #     # parse the information in DataFrame
    #     ProteinDataTable <- tryCatch(read.csv(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
    #     if (!is.null(ProteinDataTable))
    #     {
    #       ProteinDataTable <- ProteinDataTable[1,]
    #       ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
    #       # add Dataframes together if more than one accession
    #       protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    #       print(paste0(ProteinAcc," name fetched"))
    #     }  else
    #   {
    #     ProteinDataTable <- as.character(ProteinAcc)
    #     ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
    #     print(ProteinInfoParsed)
    #     # add Dataframes together if more than one accession
    #     protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    #   }
    
    #   }else {
    #     HandleBadRequests(Request$status_code)
    
    #       ProteinDataTable <- as.character(ProteinAcc)
    #       ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
    #       # add Dataframes together if more than one accession
    #       protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    #   }
    # }
    
    # return(protein_gene_name)
    
    protein_gene_name = data.frame()
    # print(gene_names)
    # gene_names_df <- data.frame(
    #   key = gene_names[,1],
    #   pair = gene_names[,2]
    # )
    for (ProteinAcc in ProteinAccList)
    {
      ProteinDataTable <- as.character(lookup(ProteinAcc, as.data.frame(id_to_name), missing=ProteinAcc))
      ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
      # add Dataframes together if more than one accession
      protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    }
    
    return(protein_gene_name)
    
    
  }
  
  output$prot_int_table <- DT::renderDataTable({
    
    
    protein_interaction_df <- df_interaction()
    protein_gene_name <- df_names()
    print(protein_interaction_df)
    print("here")
    print(class(protein_interaction_df))
    if(df_names() == 0)
    {
      
      p_int_formatted <- data.frame()
      
    } else {
      
      protein_interaction_df[,1] <- as.character(protein_interaction_df[,1])
      
      p_int_formatted <- data.frame()
      count = 0
      n = 1
      for ( id in protein_interaction_df[,1])
      {
        count = count + 1
        if(!is.null(protein_interaction_df[,2]))
        {
          a = strsplit(as.character(protein_interaction_df[,2]),"; ")
          
          for(int_with in a[[count]])
          {
            p_int_row <- data.frame(id = as.character(paste0(as.character(lookup(id, as.data.frame(id_to_name), missing="Not found"))," ( ", id," )")),
                                    Interacts_With = as.character(paste0(as.character(lookup(int_with, as.data.frame(id_to_name), missing="Not found"))," ( ", int_with," )")),
                                    row.names = n)
            p_int_formatted <- rbind(p_int_formatted,p_int_row)
            n = n + 1
          }
        }
      }
      
      # for(i in 1:nrow(protein_interaction_df))
      # {
      #     protein_interaction_df[i,1] <- paste0(protein_interaction_df[i,1],
      #                               ' (',
      #                               protein_gene_name[protein_interaction_df[i,1],1],
      #                               ')')
      # }
      # print(protein_interaction_df)
      # colnames(protein_interaction_df)[2] <- "Interacts With"
      
    }
    
    p_int_formatted
    
  })
  
  output$prot_name_table <- DT::renderDataTable({
    protein_gene_name <- df_names()
    if(protein_gene_name == 0)
    {
      protein_gene_name <- data.frame()
    } else {
      
      protein_gene_name <- cbind(ID = rownames(protein_gene_name),protein_gene_name)
      rownames(protein_gene_name) <- 1:nrow(protein_gene_name)
      colnames(protein_gene_name)[2] <- "Names"
      
    } 
    protein_gene_name
    
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  #########  Gemne Mania  ###########
  ###################################
  ###################################
  
  df_genemania <- reactive({
    print("running")
    if (is.null(input$file_gene)&&is.null(input$text_gene)) {
      return(NULL)
    }else if(!is.null(input$file_gene)){
      parts <- strsplit(input$file_gene$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      gene_names <- read.csv(input$file_gene$datapath)
      gene_names <- na.omit(gene_names)
      gene_names <- gene_names[!duplicated(gene_names[, 1]), ]
    }else{
      gene_name<-strsplit(input$text_prot_Int," ")
      gene_names <- data.frame(gene_name[[1]][1])
      for (x in 2:length(gene_name[[1]])) {
        gene_names<-rbind(gene_names,gene_name[[1]][x])
      }
      
      print(gene_names)
      gene_names <- na.omit(gene_names)
      gene_names <- gene_names[!duplicated(gene_names[, 1]), ]
    }
    
    return(gene_names)
    
  })
  
  observeEvent(input$genemania_submit, {
    
    hide("help_text_gene_mania")
    print("running...")
    organism_id <- input$organismID
    gene_names <- df_genemania()
    base_url <- "http://genemania.org/search/"
    
    url <- paste0(base_url,organism_id)
    for ( names in as.character(gene_names) )
    {
      url <- paste0(url,"/",names)
    }
    # print(gene_mania_link())
    print(url)
    gene_mania_link(url)
    shinyjs::toggle("hide_link")
    
  })
  
  output$linkCo <- renderUI({
    tags$a(href = gene_mania_link(), "here", inline =TRUE)
  })
  
  output$help_text_gene_mania <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page submits a given gene list to 
          <a href ='http://genemania.org/'>GeneMania.org</a> to retrieve the Co-expression 
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ######  Protein Sequences  ########
  ###################################
  ###################################
  Proteins <- NULL
  df_prot_seq <- eventReactive(input$submit_prot_Seq, {
    print("running")
    if (is.null(input$file_prot_seq)&& is.null(input$text_prot_seq)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_seq)){
      parts <- strsplit(input$file_prot_seq$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      protein_Id <- unique(as.character(na.omit(read.csv(input$file_prot_seq$datapath)[,1])))
      Proteins <<- protein_Id
    }
    else{
      Acessions<-strsplit(input$text_prot_seq,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_prot_seq," ")
        
      }
      Proteins <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Proteins<-rbind(Proteins,Acessions[[1]][x])
      }
      
      print(Proteins)
      Proteins <- na.omit(Proteins)
      Proteins <- Proteins[!duplicated(Proteins[, 1]), ]
      
    }
    
    shinyjs::show("downloadData")
    return(Proteins)
    
  })
  
  Seqdata <- NULL
  
  output$help_text_prot_seq <- renderUI({
    HTML("<br>
    <br>
      <center>
        <p>
          <b>This page retrieves the full protein sequences from <a href ='https://www.uniprot.org/'>UniProt.org</a> of a given set of UniProt accessions, Please upload accessions to start analysis.
          </b>
        </p>
      </center>
    ")
  })
  
  output$help_text_prot_seq_evol <- renderUI({
    HTML("<br>
    <br>
      <center>
        <p>
          <b>
          This page performs Evolutionary analysis of protein sequences retrieved from <a href ='https://www.uniprot.org/'>UniProt.org</a>, Please upload accessions to start analysis.
          </b>
        </p>
      </center>
    ")
  })
  
  output$help_text_prot_seq_Patho <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves protein's pathological information from <a href ='https://www.uniprot.org/'>UniProt.org</a> of a given set of UniProt accessions, Please upload accessions to start analysis.
          </b>
        </p>
      </center>
    ")
  })
  output$SequencePlot <- renderPlot(
    {
      if (!is.null(df_prot_seq()))
      {
        hide("help_text_prot_seq")
        if(is.null(Seqdata))
        {
          Proteins <- df_prot_seq()
          Seqdata <<- GetSequences(Proteins)
        }
        PlotPhysicochemical(Seqdata)
      }
    }
    
  )
  output$GravyPlot <- renderPlot({
    if (!is.null(df_prot_seq()))
    {
      hide("help_text_prot_seq")
      if(is.null(Seqdata))
      {
        Proteins <- df_prot_seq()
        Seqdata <<- GetSequences(Proteins)
      }
      PlotGravy(Seqdata)
    }
    
  })
  output$ChargePlot <- renderPlot({
    if (!is.null(df_prot_seq()))
    {
      hide("help_text_prot_seq")
      if(is.null(Seqdata))
      {
        Proteins <- df_prot_seq()
        Seqdata <<- GetSequences(Proteins)
      }
      PlotCharge(Seqdata)
    }
    
  })
  output$AcidityPlot <- renderPlot({
    if (!is.null(df_prot_seq()))
    {
      hide("help_text_prot_seq")
      if(is.null(Seqdata))
      {
        Proteins <- df_prot_seq()
        Seqdata <<- GetSequences(Proteins)
      }
      PlotAcidity(Seqdata)
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("Sequences", ".FASTA")
    },    
    content = function(file) {
      
      Accessions <- df_prot_seq()
      for (Acc in Accessions)
      {
        Request <- tryCatch(
          {
            GET(paste0("https://www.uniprot.org/uniprot/" , Acc , ".Fasta") , timeout(10))
          },error = function(cond)
          {
            message("Internet connection problem occurs and the function will return the original error")
            message(cond)
          }
        )
        if (Request$status_code == 200)
        {
          OutNumber <<- OutNumber + 1
          Fastadata <- read.csv(paste0("https://www.uniprot.org/uniprot/" , Acc , ".Fasta") , header = F , sep = "\t")
          Sequences <- paste0(as.character(unlist(Fastadata)) , collapse = "\n")
          write.table(x = Sequences , file = paste0(FileName ,".fasta") , quote = F , row.names = F , col.names = F, append = T)
        }
        
      }
    }
  )
  
  ##################################
  ######## protein revolution ######
  
  df_prot_seq_evol <- eventReactive(input$submit_prot_seq_evol,{
    print("running")
    if (is.null(input$file_prot_seq_evol)&& is.null(input$text_prot_seq_evol)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_seq_evol)){
      parts <- strsplit(input$file_prot_seq_evol$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      protein_Id <- unique(as.character(na.omit(read.csv(input$file_prot_seq_evol$datapath)[,1])))
      Proteins <- protein_Id
    }
    else{
      Acessions<-strsplit(input$text_prot_seq_evol,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_prot_seq_evol," ")
        
      }
      Proteins <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Proteins<-rbind(Proteins,Acessions[[1]][x])
      }
      
      print(Proteins)
      Proteins <- na.omit(Proteins)
      Proteins <- Proteins[!duplicated(Proteins[, 1]), ]
      
    }
    return(Proteins)
    
  })
  
  GenesObj <- NULL
  
  output$GenePlot <- renderRadialNetwork(
    {
      if (!is.null(df_prot_seq_evol()))
      {
        if (is.null(GenesObj))
        {
          Proteins <- df_prot_seq_evol()
          GenesObj <- GetNamesTaxa(Proteins)
        }
        ConstructGenes(GenesObj)
      }
    }
  )
  
  output$Chromo <- renderPlot(
    if (!is.null(df_prot_seq_evol()))
    {
      if(is.null(GenesObj))
      {
        Proteins <- df_prot_seq_evol()
        GenesObj <- GetNamesTaxa(Proteins)
      }
      PlotChromosomeInfo(GenesObj)
    }
  )
  
  output$Phylogenetic <- renderPlot(
    {
      if(!is.null(df_prot_seq_evol()))
        if(is.null(Seqdata))
        {
          Proteins <- df_prot_seq_evol()
          Seqdata <<- GetSequences(Proteins)
        }
      ConstructPhylogeny(Seqdata)
    }
  )
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  #Pathogens
  df_prot_seq_Patho <- eventReactive(input$submit_prot_seq_Patho,{
    print("running")
    if (is.null(input$file_prot_seq_Patho)&& is.null(input$text_prot_seq_Patho)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_seq_Patho)){
      parts <- strsplit(input$file_prot_seq_Patho$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      type <- tolower(type)
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      protein_Id <- unique(as.character(na.omit(read.csv(input$file_prot_seq_Patho$datapath)[,1])))
      Proteins <- protein_Id
    }
    else{
      Acessions<-strsplit(input$text_prot_seq_Patho,",")
      if(length(Acessions[[1]])==1){
        Acessions<-strsplit(input$text_prot_seq_Patho," ")
        
      }
      Proteins <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Proteins<-rbind(Proteins,Acessions[[1]][x])
      }
      
      print(Proteins)
      Proteins <- na.omit(Proteins)
      Proteins <- Proteins[!duplicated(Proteins[, 1]), ]
      
    }
    
    return(Proteins)
    
  })
  
  Pathodata <- NULL
  DiseaseTable <- NULL
  
  output$DisaeseTable <- renderDataTable({
    if(!is.null(df_prot_seq_Patho()))
    {
      Proteins <- df_prot_seq_Patho()
      Pathodata <- GetPathology_Biotech(Proteins)
      DiseaseTable <- Get.diseases(Pathodata) 
    }
  }, escape = F)
  
  output$DiseasePlot <- renderBubbles({
    if(!is.null(df_prot_seq_Patho()))
    {
      if(!is.null(DiseaseTable))
      {
        Plot.NDiseases(DiseaseTable)
      }
      else {
        Proteins <- df_prot_seq_Patho()
        Pathodata <- GetPathology_Biotech(Proteins)
        DiseaseTable <- Get.diseases(Pathodata)
        Plot.NDiseases(DiseaseTable)
      }
    }
  })
  
  ########################################
  ##### Increases the Upload Limit #######
  ########################################
  
  options(shiny.maxRequestSize=10000*1024^2)
  
  ####___________________________####
  # RDS DATA ANALYSIS HERE
  
  rv <- reactiveValues(aggregate = NULL, genes = NULL, reductions = NULL, meta_nums = NULL, meta_cats = NULL, mysplitbydefault = NULL, pcs = NULL, use.pcs = NULL)
  
  observeEvent(input$rds_file, {
    req(input$rds_file)
    
    # Load the uploaded RDS file
    aggregate <- readRDS(input$rds_file$datapath)
    rv$aggregate <- aggregate
    
    # Extract and save the required values into reactiveValues
    rv$genes <- aggregate@assays$RNA
    rv$reductions <- attributes(aggregate@reductions)
    rv$meta_nums <- colnames(select_if(aggregate@meta.data, is.numeric))
    rv$meta_cats <- c(colnames(select_if(aggregate@meta.data, is.character)), colnames(select_if(aggregate@meta.data, is.factor)), colnames(select_if(aggregate@meta.data, is.logical)))
    rv$meta_cats <- rv$meta_cats[rv$meta_cats != "orig.ident"]
    rv$mysplitbydefault <- rv$meta_cats
    rv$pcs <- c('PC_1','PC_2','PC_3','PC_4','PC_5','PC_6','PC_7','PC_8','PC_9')
    rv$use.pcs <- 1:50
    
    # update values based on input from ui
    outVar_double = reactive({
      if (input$dataset == 'Genes'){mydata=rownames(rv$genes)}
      else if (input$dataset == 'Numeric Metadata') {mydata=rv$meta_nums}
      else if (input$dataset == 'PCs') {mydata=rv$pcs}
      mydata
    })
    
    # update values based on input from ui
    outVar_single = reactive({
      if (input$dataset_single == 'Genes'){mydata=rownames(rv$genes)}
      else if (input$dataset_single == 'Numeric Metadata') {mydata=rv$meta_nums}
      else if (input$dataset_single == 'PCs') {mydata=rv$pcs}
      mydata
    })
    
    # update values based on input from ui
    outVar_seperated = reactive({
      if (input$dataset_seperated == 'Genes'){mydata=rownames(rv$genes)}
      else if (input$dataset_seperated == 'Numeric Metadata') {mydata=rv$meta_nums}
      else if (input$dataset_seperated == 'PCs') {mydata=rv$pcs}
      mydata
    })
    
    getResChoices = reactive({
      req(rv$aggregate)
      req(input$identity_table)
      mydata = levels(eval(call("$", rv$aggregate, input$identity_table)))
      mydata
    })
    
    # Reduction Type for the Single Marker Plot
    observe({
      updateSelectInput(session, "reduction_single",
                        choices = rv$reductions
      )})
    
    # Reduction Type for the Double Marker Plot
    observe({
      updateSelectInput(session, "reduction_double",
                        choices = rv$reductions
      )})
    
    # Primary numeric value in the double marker plot
    observe({
      updateSelectInput(session, "numeric",
                        choices = outVar_double()
      )})
    
    # Secondary numeric value in the double marker plot
    observe({
      updateSelectInput(session, "numeric2",
                        choices = outVar_double(), 
                        selected = outVar_double()[2]
      )})
    
    # Only numeric input for the single marker plot
    observe({
      updateSelectInput(session, "numeric_single",
                        choices = outVar_single()
      )})
    
    # Double Marker identity
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "categorical",
                        choices = rv$meta_cats
      )})
    
    # Single Marker identity
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "categorical_single",
                        choices = rv$meta_cats
      )})
    
    # Multiple Feature Plot identity
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "multiple_feature_categorical_plot",
                        choices = rv$meta_cats
      )})
    
    # Multiple Feature Plot reductions
    observe({
      req(rv$reductions)
      updateSelectInput(session, "multiple_feature_reduction",
                        choices = rv$reductions
      )})
    
    # Cluster Tree identity
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "identity_tree",
                        choices = rv$meta_cats
      )})
    
    # Seperated Identity
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "identity_seperated",
                        choices = rv$meta_cats
      )})
    
    # Seperated Identity 2
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "identity_seperated2",
                        choices = rv$meta_cats
      )})
    
    # Seperated Numeric
    observe({
      updateSelectInput(session, "numeric_seperated",
                        choices = outVar_seperated()
      )})
    
    # Seperated Reduction
    observe({
      updateSelectInput(session, "reduction_seperated",
                        choices = rv$reductions
      )})
    
    # Seperated categorical Identity
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "identity_seperated_categorical",
                        choices = rv$meta_cats
      )})
    
    # Seperated categorical identity2
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "identity2_seperated_categorical",
                        choices = rv$meta_cats,
                        selected = rv$meta_cats[2]
      )})
    
    # Seperated categorical Reduction
    observe({
      req(rv$reductions)
      updateSelectInput(session, "reduction_seperated_categorical",
                        choices = rv$reductions
      )})
    
    # Seperated categorical split
    observe({
      req(rv$mysplitbydefault)
      updateSelectInput(session, "split_seperated_categorical",
                        choices = rv$mysplitbydefault
      )})
    
    # Multiple Feature Plot
    observe({
      updateSelectizeInput(session, "multiple_feature_list",
                           choices = rownames(rv$genes), server = TRUE
      )})
    
    # Table Identity
    observe({
      req(rv$meta_cats)
      updateSelectInput(session, "identity_table",
                        choices = rv$meta_cats
                        
      )})
    
    
    # Table Marker
    observe({
      req(getResChoices())
      updateSelectInput(session, "compare_table",
                        choices = getResChoices()
      )})
    
    # Table Compare
    observe({
      req(getResChoices())
      updateSelectInput(session, "markers_table",
                        choices = getResChoices()
      )})
    
    
    # Documentation
    output$markdown <- renderUI({
      includeMarkdown("Next step : trajectory analysis to be added")
    }) 
    
    ####___________________________________________####
    #### PLOTS FOR DOUBLE MARKER TAB ###
    # Marker Plot Double
    output$MarkerGenePlot <- renderPlot({
      req(rv$aggregate)
      req(input$numeric)
      req(input$numeric2)
      req(input$reduction_double)
      
      FeaturePlot(
        rv$aggregate,
        c(input$numeric, input$numeric2),
        blend=TRUE,
        reduction=input$reduction_double
      )
    })
    
    # Double Feature Violin Plot
    output$ViolinPlot <- renderPlot({
      req(input$categorical)
      req(rv$aggregate)
      req(input$numeric)
      req(input$numeric2)
      
      Idents(rv$aggregate) <- input$categorical
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      VlnPlot(object =  rv$aggregate, features = c(input$numeric, input$numeric2), pt.size = 0.05)
    })
    
    # Double Feature Categorical Feature Plot
    output$CategoricalPlot <- renderPlot({
      req(input$categorical)
      req(rv$aggregate)
      req(input$reduction_double)
      
      Idents(rv$aggregate) <- input$categorical
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      DimPlot(object = rv$aggregate, pt.size=0.5, reduction = input$reduction_double, label = T)
    })
    
    observeEvent(input$add_doublemarker, {
      hide("error_text_report")
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          fluidRow(
            column(width = 2),
            column(width = 8,
                   h3("Feature Plot (Double Marker)", align = "center"), 
                   plotOutput("MarkerGenePlotReport", height = 500),
                   br(),
                   br(),
                   h3("Violin Plot (Double Marker)", align = "center"), 
                   plotOutput("ViolinPlotReport", height = 500),
                   br(),
                   br(),
                   h3("Categorical Plot (Double Marker)", align = "center"), 
                   plotOutput("CategoricalPlotReport", height = 500),
                   br(),
                   br()
                   )
            )
        ),
      )
      
      output$MarkerGenePlotReport <- renderPlot({
        req(rv$aggregate)
        req(input$numeric)
        req(input$numeric2)
        req(input$reduction_double)
        
        FeaturePlot(
          rv$aggregate,
          c(input$numeric, input$numeric2),
          blend=TRUE,
          reduction=input$reduction_double
        )
      })
      
      output$ViolinPlotReport <- renderPlot({
        req(input$categorical)
        req(rv$aggregate)
        req(input$numeric)
        req(input$numeric2)
        
        Idents(rv$aggregate) <- input$categorical
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        VlnPlot(object =  rv$aggregate, features = c(input$numeric, input$numeric2), pt.size = 0.05)
      })
      
      output$CategoricalPlotReport <- renderPlot({
        req(input$categorical)
        req(rv$aggregate)
        req(input$reduction_double)
        
        Idents(rv$aggregate) <- input$categorical
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        DimPlot(object = rv$aggregate, pt.size=0.5, reduction = input$reduction_double, label = T)
      })
    })
    
    ####___________________________________________####
    #### PLOTS FOR SINGLE MARKER TAB ###
    # Marker Plot Single
    output$MarkerGenePlotSingle <- renderPlot({
      req(rv$aggregate)
      req(input$numeric_single)
      req(input$reduction_single)
      
      FeaturePlot(
        rv$aggregate,
        c(input$numeric_single),
        reduction=input$reduction_single
      )
    })
    
    # Single Feature Violin Plot
    output$ViolinPlotSingle <- renderPlot({
      req(input$categorical_single)
      req(rv$aggregate)
      req(input$numeric_single)
      
      Idents(rv$aggregate) <- input$categorical_single
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      VlnPlot(object =  rv$aggregate, features = c(input$numeric_single), pt.size = 0.05)
    })
    
    # Single Feature Categorical Feature Plot
    output$CategoricalPlotSingle <- renderPlot({
      req(input$categorical_single)
      req(rv$aggregate)
      req(input$categorical_single)
      req(input$reduction_single)
      
      Idents(rv$aggregate) <- input$categorical_single
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      DimPlot(object = rv$aggregate, group.by=input$categorical_single, pt.size=0.5, reduction = input$reduction_single, label = T)
    })
    
    observeEvent(input$add_singlemarker, {
      hide("error_text_report")
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          fluidRow(
            column(width = 2),
            column(width = 8,
                   h3("Feature Plot (Single Marker)", align = "center"), 
                   plotOutput("MarkerGeneSingleReport", height = 500),
                   br(),
                   br(),
                   h3("Violin Plot (Single Marker)", align = "center"), 
                   plotOutput("ViolinPlotSingleReport", height = 500),
                   br(),
                   br(),
                   h3("Categorical Plot (Single Marker)", align = "center"), 
                   plotOutput("CategoricalPlotSingleReport", height = 500),
                   br(),
                   br()
            )
          )
        ),
      )
      
      output$MarkerGeneSingleReport <- renderPlot({
        req(rv$aggregate)
        req(input$numeric_single)
        req(input$reduction_single)
        
        FeaturePlot(
          rv$aggregate,
          c(input$numeric_single),
          reduction=input$reduction_single
        )
      })
      
      output$ViolinPlotSingleReport <- renderPlot({
        req(input$categorical_single)
        req(rv$aggregate)
        req(input$numeric_single)
        
        Idents(rv$aggregate) <- input$categorical_single
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        VlnPlot(object =  rv$aggregate, features = c(input$numeric_single), pt.size = 0.05)
      })
      
      output$CategoricalPlotSingleReport <- renderPlot({
        req(input$categorical_single)
        req(rv$aggregate)
        req(input$categorical_single)
        req(input$reduction_single)
        
        Idents(rv$aggregate) <- input$categorical_single
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        DimPlot(object = rv$aggregate, group.by=input$categorical_single, pt.size=0.5, reduction = input$reduction_single, label = T)
      })
    })
    
    ####___________________________________________####
    #### PLOTS FOR MULTIPLE FEATURE PLOT TAB ###
    # Multiple Feature Categorical Plot
    output$MultipleFeatureCategoricalPlot <- renderPlot({
      req(input$multiple_feature_categorical_plot)
      req(rv$aggregate)
      req(input$multiple_feature_reduction)
      
      Idents(rv$aggregate) <- input$multiple_feature_categorical_plot
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      DimPlot(object = rv$aggregate, group.by=input$multiple_feature_categorical_plot, pt.size=0.5, reduction = input$multiple_feature_reduction, label = T)
    })
    
    # Multiple Feature Plot
    output$MultipleFeaturePlot <- renderPlot({
      req(rv$aggregate)
      req(input$multiple_feature_list)
      req(input$multiple_feature_reduction)
      
      FeaturePlot(
        rv$aggregate,
        input$multiple_feature_list,
        blend=FALSE,
        reduction=input$multiple_feature_reduction,
        ncol=4
      )
    })
    
    observeEvent(input$add_multiplefeature, {
      hide("error_text_report")
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          fluidRow(
            column(width = 2),
            column(width = 8,
                   h3("Multiple Feature Categorical Plot", align = "center"), 
                   plotOutput("MultipleFeatureCategoricalReport", height = 500),
                   br(),
                   br(),
                   h3("Multiple Feature Plot", align = "center"), 
                   plotOutput("MultipleFeaturePlotReport", height = 500),
                   br(),
                   br()
            )
          )
        ),
      )

      output$MultipleFeatureCategoricalReport <- renderPlot({
        req(input$multiple_feature_categorical_plot)
        req(rv$aggregate)
        req(input$multiple_feature_reduction)
        
        Idents(rv$aggregate) <- input$multiple_feature_categorical_plot
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        DimPlot(object = rv$aggregate, group.by=input$multiple_feature_categorical_plot, pt.size=0.5, reduction = input$multiple_feature_reduction, label = T)
      })
      
      output$MultipleFeaturePlotReport <- renderPlot({
        req(rv$aggregate)
        req(input$multiple_feature_list)
        req(input$multiple_feature_reduction)
        
        FeaturePlot(
          rv$aggregate,
          input$multiple_feature_list,
          blend=FALSE,
          reduction=input$multiple_feature_reduction,
          ncol=4
        )
      })
    })
    
    ####___________________________________________####
    #### PLOTS FOR CLUSTER TREE TAB ###
    # Cluster Tree Plot
    output$ClusterTree <- renderPlot({
      req(input$identity_tree)
      req(rv$aggregate)
      req(rv$use.pcs)
      
      Idents(rv$aggregate) <- input$identity_tree
      rv$aggregate <- BuildClusterTree(
        rv$aggregate, dims = rv$use.pcs)
      PlotClusterTree(rv$aggregate)
    })
    
    observeEvent(input$add_clustertree, {
      hide("error_text_report")
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          fluidRow(
            column(width = 2),
            column(width = 8,
                   h3("Cluster Tree", align = "center"), 
                   plotOutput("ClusterTreeReport", height = 500),
                   br(),
                   br()
            )
          )
        ),
      )
      
      output$ClusterTreeReport <- renderPlot({
        req(input$identity_tree)
        req(rv$aggregate)
        req(rv$use.pcs)
        
        Idents(rv$aggregate) <- input$identity_tree
        rv$aggregate <- BuildClusterTree(
          rv$aggregate, dims = rv$use.pcs)
        PlotClusterTree(rv$aggregate)
      })
    })
    
    ####___________________________________________####
    #### PLOTS FOR SEPARATE IDENTITY CATEGORICAL TAB ###
    
    # Seperated Identity Categorical Plot
    output$SeperatedIdentityCategorical <- renderPlot({
      req(input$identity_seperated_categorical)
      req(input$reduction_seperated_categorical)
      req(rv$aggregate)
      req(input$split_seperated_categorical)
      
      Idents(rv$aggregate) <- input$identity_seperated_categorical
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      DimPlot(rv$aggregate, reduction=input$reduction_seperated_categorical,
              split.by = input$split_seperated_categorical, ncol=4
      )
    })
    
    # Seperated Identity 2 Categorical Plot
    output$SeperatedIdentity2Categorical <- renderPlot({
      req(input$identity2_seperated_categorical)
      req(input$reduction_seperated_categorical)
      req(rv$aggregate)
      req(input$split_seperated_categorical)
      
      aggregate2 <- rv$aggregate
      Idents(aggregate2) <- input$identity2_seperated_categorical
      order <- sort(levels(aggregate2))
      levels(aggregate2) <- order
      DimPlot(aggregate2, reduction=input$reduction_seperated_categorical,
              split.by = input$split_seperated_categorical, ncol=4
      )
    })
    
    # Seperated Categorical table
    output$SeperatedCountsCategorical <- renderPlot({
      req(rv$aggregate)
      req(input$identity_seperated_categorical)
      req(input$identity2_seperated_categorical)
      
      length_data = as.data.frame(prop.table(table(eval(call('$', rv$aggregate[[]], input$identity_seperated_categorical)), 
                                                   eval(call('$', rv$aggregate[[]], input$identity2_seperated_categorical))),1))
      colnames(length_data) = c(input$identity_seperated_categorical, input$identity2_seperated_categorical, 'Freq')
      mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
      ggplot(length_data, aes_string(x=input$identity_seperated_categorical, y=input$identity2_seperated_categorical, fill='Freq')) + geom_tile() + scale_fill_gradientn(colours = mycol)
    })
    
    observeEvent(input$add_separatecategorical, {
      hide("error_text_report")
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          fluidRow(
            column(width = 2),
            column(width = 8,
                   h3("Separated Identity Categorical Plot", align = "center"), 
                   plotOutput("SeperatedIdentityCategoricalReport", height = 500),
                   br(),
                   br(), 
                   h3("Separated Identity Categorical Plot", align = "center"), 
                   plotOutput("SeperatedIdentity2CategoricalReport", height = 500),
                   br(),
                   br(), 
                   h3("Separated Counts Categorical Plot", align = "center"), 
                   plotOutput("SeperatedCountsCategoricalReport", height = 500)
            )
          )
        ),
      )
      
      output$SeperatedIdentityCategoricalReport <- renderPlot({
        req(input$identity_seperated_categorical)
        req(input$reduction_seperated_categorical)
        req(rv$aggregate)
        req(input$split_seperated_categorical)
        
        Idents(rv$aggregate) <- input$identity_seperated_categorical
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        DimPlot(rv$aggregate, reduction=input$reduction_seperated_categorical,
                split.by = input$split_seperated_categorical, ncol=4
        )
      })
      
      # Seperated Identity 2 Categorical Plot
      output$SeperatedIdentity2CategoricalReport <- renderPlot({
        req(input$identity2_seperated_categorical)
        req(input$reduction_seperated_categorical)
        req(rv$aggregate)
        req(input$split_seperated_categorical)
        
        aggregate2 <- rv$aggregate
        Idents(aggregate2) <- input$identity2_seperated_categorical
        order <- sort(levels(aggregate2))
        levels(aggregate2) <- order
        DimPlot(aggregate2, reduction=input$reduction_seperated_categorical,
                split.by = input$split_seperated_categorical, ncol=4
        )
      })
      
      # Seperated Categorical table
      output$SeperatedCountsCategoricalReport <- renderPlot({
        req(rv$aggregate)
        req(input$identity_seperated_categorical)
        req(input$identity2_seperated_categorical)
        
        length_data = as.data.frame(prop.table(table(eval(call('$', rv$aggregate[[]], input$identity_seperated_categorical)), 
                                                     eval(call('$', rv$aggregate[[]], input$identity2_seperated_categorical))),1))
        colnames(length_data) = c(input$identity_seperated_categorical, input$identity2_seperated_categorical, 'Freq')
        mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
        ggplot(length_data, aes_string(x=input$identity_seperated_categorical, y=input$identity2_seperated_categorical, fill='Freq')) + geom_tile() + scale_fill_gradientn(colours = mycol)
      })
      
    })
    
    ####___________________________________________####
    #### PLOTS FOR SEPARATE FEATURE TAB ###
    
    # Seperated Feature Plot
    output$SeperatedFeature <- renderPlot({
      req(input$identity_seperated)
      req(rv$aggregate)
      req(input$numeric_seperated)
      req(input$reduction_seperated)
      req(input$identity_seperated2)
      
      Idents(rv$aggregate) <- input$identity_seperated
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      FeaturePlot(rv$aggregate, c(input$numeric_seperated), reduction=input$reduction_seperated,
                  split.by = input$identity_seperated2, ncol=4
      )
    })
    
    # Seperated Dim Plot
    output$SeperatedDim <- renderPlot({
      req(input$identity_seperated)
      req(rv$aggregate)
      req(input$reduction_seperated)
      req(input$identity_seperated2)
      
      Idents(rv$aggregate) <- input$identity_seperated
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      DimPlot(rv$aggregate, reduction=input$reduction_seperated,
              split.by = input$identity_seperated2, ncol=4
      )
    })
    
    # Seperated Violin Plot
    output$SeperatedViolin <- renderPlot({
      req(input$identity_seperated)
      req(rv$aggregate)
      req(input$numeric_seperated)
      req(input$identity_seperated2)
      
      Idents(rv$aggregate) <- input$identity_seperated
      order <- sort(levels(rv$aggregate))
      levels(rv$aggregate) <- order
      VlnPlot(rv$aggregate, c(input$numeric_seperated), group.by = input$identity_seperated, split.by = input$identity_seperated2, ncol=4)
    })
    
    # Seperated Counts table
    output$SeperatedCounts <- renderTable({
      
      req(input$identity_seperated)
      req(rv$aggregate)
      req(input$identity_seperated2)
      
      marker = c(input$numeric_seperated)
      Idents(rv$aggregate) <- input$identity_seperated
      
      if(input$dataset_seperated == 'Numeric Metadata'){
        nm <- data.frame(matrix(unlist(eval(call('$', rv$aggregate, marker[1]))), nrow=length(eval(call('$', rv$aggregate, marker[1]))), byrow=T))
        colnames(nm) = marker
        rownames(nm) = labels(eval(call('$', rv$aggregate, marker[1])))
        widedat <- nm
      }
      else{widedat <- FetchData(rv$aggregate, marker)}
      
      widedat$Cluster <- Idents(rv$aggregate)
      widedat[["CellType"]] = eval(call("$", rv$aggregate, input$identity_seperated2))
      widedat$final = paste(widedat[["CellType"]], widedat$Cluster, sep="_")
      final_object = (aggregate(widedat[, 1:2], list(widedat$final), mean)[1:2])
      lab_list = widedat[["CellType"]]
      identities = widedat$Cluster
      
      num_list = widedat[[marker]]
      
      # df needs to be fixed
      tmp_df = data.frame(identities, num_list, lab_list)
      df = as.data.frame(pivot_wider(aggregate(tmp_df[2], list(tmp_df$identities, tmp_df$lab_list), mean), names_from = Group.2, values_from = num_list))
      df[is.na(df)] <- 0
      rownames(df) = df$Group.1
      drops <- c("Group.1")
      df = df[ , !(names(df) %in% drops)]
      
      df_p = as.data.frame.matrix(prop.table((table(eval(call("$", rv$aggregate, input$identity_seperated)), eval(call("$", rv$aggregate, input$identity_seperated2)))),2))
      df_p=df_p/colSums(df_p)
      
      merged_final = as.data.frame.matrix(merge(df, df_p, by.x = 'row.names', by.y = 'row.names', suffixes = c(".AvgExpression",".Proportion")))
      merged_final
    }, width = "100%", colnames=TRUE, rownames=TRUE, digits=4)
    
    observeEvent(input$add_separatedfeature, {
      hide("error_text_report")
      insertUI(
        selector = '#placeholder',
        ui = tagList(
          fluidRow(
            column(width = 2),
            column(width = 8,
                   h3("Separated Feature Plot", align = "center"), 
                   plotOutput("SeperatedFeatureReport", height = 500),
                   br(),
                   br(), 
                   h3("Separated Dimplot", align = "center"), 
                   plotOutput("SeperatedDimReport", height = 500),
                   br(),
                   br(), 
                   h3("Separated Counts Categorical Plot", align = "center"), 
                   plotOutput("SeperatedViolinReport", height = 500)
            )
          )
        ),
      )
      
      output$SeperatedFeatureReport <- renderPlot({
        req(input$identity_seperated)
        req(rv$aggregate)
        req(input$numeric_seperated)
        req(input$reduction_seperated)
        req(input$identity_seperated2)
        
        Idents(rv$aggregate) <- input$identity_seperated
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        FeaturePlot(rv$aggregate, c(input$numeric_seperated), reduction=input$reduction_seperated,
                    split.by = input$identity_seperated2, ncol=4
        )
      })
      
      # Seperated Dim Plot
      output$SeperatedDimReport <- renderPlot({
        req(input$identity_seperated)
        req(rv$aggregate)
        req(input$reduction_seperated)
        req(input$identity_seperated2)
        
        Idents(rv$aggregate) <- input$identity_seperated
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        DimPlot(rv$aggregate, reduction=input$reduction_seperated,
                split.by = input$identity_seperated2, ncol=4
        )
      })
      
      # Seperated Violin Plot
      output$SeperatedViolinReport <- renderPlot({
        req(input$identity_seperated)
        req(rv$aggregate)
        req(input$numeric_seperated)
        req(input$identity_seperated2)
        
        Idents(rv$aggregate) <- input$identity_seperated
        order <- sort(levels(rv$aggregate))
        levels(rv$aggregate) <- order
        VlnPlot(rv$aggregate, c(input$numeric_seperated), group.by = input$identity_seperated, split.by = input$identity_seperated2, ncol=4)
      })

    })
    
    # Marker Table
    output$markers <- renderTable({
      req(input$identity_table)
      req(rv$aggregate)
      req(input$markers_table)
      req(input$compare_table)
      
      Idents(rv$aggregate) <- input$identity_table
      if (as.logical(length(c(input$compare_table)))){FindMarkers(rv$aggregate, ident.1=input$markers_table, ident.2=input$compare_table)}
      else {FindMarkers(rv$aggregate, ident.1=input$markers_table)}
    }, rownames = TRUE, colnames = TRUE, width = "100%", digits=-5)
  })
  
  ####_________________________________________________________________________________####
  # INTEGRATED ANALYSIS HERE
  
  ####___________________________####
  # Non integration SelectizeInput here
  
  output$vlnplot16_selectize <- renderUI({
    selectizeInput("vlnplot16_selected", "Select Genes (vlnplot16)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$vlnplot17_selectize <- renderUI({
    selectizeInput("vlnplot17_selected", "Select Genes (vlnplot17)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$featureplot9_selectize <- renderUI({
    selectizeInput("featureplot9_selected", "Select Genes (featureplot9)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  # CCA Integration SelectizeInput here
  output$vlnplot3_selectize <- renderUI({
    selectizeInput("vlnplot3_selected", "Select Genes (vlnplot3)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$vlnplot4_selectize <- renderUI({
    selectizeInput("vlnplot4_selected", "Select Genes (vlnplot4)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$vlnplot5_selectize <- renderUI({
    selectizeInput("vlnplot5_selected", "Select Genes (vlnplot5)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$featureplot1_selectize <- renderUI({
    selectizeInput("featureplot1_selected", "Select Genes (featureplot1)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$featureplot2_selectize <- renderUI({
    selectizeInput("featureplot2_selected", "Select Genes (featureplot2)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  # Auto Annotation SelectizeInput here
  output$vlnplot7_selectize <- renderUI({
    selectizeInput("vlnplot7_selected", "Select Genes (vlnplot7)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$vlnplot8_selectize <- renderUI({
    selectizeInput("vlnplot8_selected", "Select Genes (vlnplot8)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$vlnplot9_selectize <- renderUI({
    selectizeInput("vlnplot9_selected", "Select Genes (vlnplot9)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$vlnplot10_selectize <- renderUI({
    selectizeInput("vlnplot10_selected", "Select Genes (vlnplot10)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$featureplot5_selectize <- renderUI({
    selectizeInput("featureplot5_selected", "Select Genes (featureplot5)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$ridgeplot2_selectize <- renderUI({
    selectizeInput("ridgeplot2_selected", "Select Genes (ridgeplot2)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  # Differential Expression SelectizeInput here
  output$de_selectize <- renderUI({
    selectizeInput("de_selected", "Select Genes (v1)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$idents1_selectize <- renderUI({
    selectizeInput("idents1_selected", "Select Genes (Idents)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  # Pseudobulk SelectizeInput here
  output$pseudobulk_selectize <- renderUI({
    selectizeInput("pseudobulk_selected", "Select Genes (vb1)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$idents2_selectize <- renderUI({
    selectizeInput("idents2_selected", "Select Genes (Idents)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  # Compare between SCDE and PSCDE Tab
  
  output$vlnplot14_selectize <- renderUI({
    selectizeInput("vlnplot14_selected", "Select Genes (vlnplot14)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  output$vlnplot15_selectize <- renderUI({
    selectizeInput("vlnplot15_selected", "Select Genes (vlnplot15)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  # Monocle3 Workflow Tab
  output$monocle_selectize <- renderUI({
    selectizeInput("monocle_selected", "Select Genes (monocle_ident_selection)", 
                   choices = NULL, 
                   multiple = TRUE)
  })
  
  ####___________________________####
  options(Seurat.object.assay.version = "v5")
  
  # Reactive values to save plot and df
  plots <- reactiveValues()
  df <- reactiveValues()
  
  df$complete <- FALSE
  
  # Track which tab has been opened
  count <- reactiveValues(
    nonintegration_tab = 0,
    CCA_tab = 0,
    autoannotation_tab = 0,
    diff_expression_tab = 0,
    pseudobulk_tab = 0, 
    compare_tab = 0, 
    monocle3_tab = 0
  )
  
  # Cleaner Tab
  cleaned_data <- reactive({
    req(input$file1)
    df <- read.csv(input$file1$datapath)
    df_no_duplicates <- df[!duplicated(df[, 1]), ]
    return(df_no_duplicates)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {paste0("clean_", input$file1$name)},
    content = function(file) {
      write.csv(cleaned_data(), file, row.names = FALSE)
    }
  )
  
  #### Read Raw count/expression Matrix Tab ####
  observeEvent(input$file2, {
    req(input$file2)
    
    df$mode <- "single"
    df$complete <- FALSE
    
    # Track which tab has been opened
    count <- reactiveValues(
      CCA_tab = 0,
      autoannotation_tab = 0,
      diff_expression_tab = 0,
      pseudobulk_tab = 0, 
      compare_tab = 0
    )
    
    merged_seurat_filtered.sub <- read.csv(input$file2$datapath, row.names = 1, check.names = FALSE)
    
    # Set idents
    df$idents <- sub("_.*", "", input$file2$name) %>% unique()
    
    #### ___________________________________________________________________________ ####
    ##THIS IS NON INTEGRATED (NEW TAB)
    # ___create Seurat Object with count data --------
    # include only genes that are are expressed in 3 or more cells and cells with complexity of 200 genes or more
    merged_seurat_filtered <- CreateSeuratObject(counts = merged_seurat_filtered.sub, project = "BII", min.cells = 3, min.features = 200)
    # str(merged_seurat_filtered) #
    
    # count matrix
    # merged_seurat_filtered@assays$RNA@counts[1:10,1:10] #
    
    # 2. QC --------
    merged_seurat_filtered[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_filtered, pattern = "^MT-")
    # str(merged_seurat_filtered) #
    
    # Show QC metrics for the first 5 cells
    # head(merged_seurat_filtered@meta.data, 5) #
    
    # We filter cells that have unique feature counts over 2,500 or less than 200
    # We filter cells that have >5% mitochondrial counts
    
    # ___Visualize QC metrics as a violin plot -------
    plots$vlnplot6 <- VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # ___feature-feature or gene-gene relationship --------
    plots$plot1 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plots$plot2 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    
    merged_seurat_filtered <- subset(merged_seurat_filtered, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
    # str(merged_seurat_filtered) #
    
    # unique(merged_seurat_filtered@meta.data$singleR.labels) #
    ##LINE BELOW NOT REQUIRED UNLESS YOU ARE DOING TRAJECTORY
    ### Idents(merged_seurat_filtered) <- merged_seurat_filtered$singleR.labels ###
    
    # 3. Normalization ----------
    merged_seurat_filtered <- NormalizeData(merged_seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
    # str(merged_seurat_filtered)#
    
    # ___identification of highly variable features ---------
    merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered, selection.method = "vst", nfeatures = 2000)
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(merged_seurat_filtered), 10)
    # top10 #
    
    # ___plot variable features with and without labels ---------
    plots$plot3_sub <- VariableFeaturePlot(merged_seurat_filtered)
    plots$plot4_sub <- LabelPoints(plot = plots$plot3_sub, points = top10, repel = TRUE)
    
    # 4. scaling the data (performed prior to linear dim reduction) ---------
    df$all.genes <- rownames(merged_seurat_filtered)
    merged_seurat_filtered <- ScaleData(merged_seurat_filtered, features = df$all.genes)
    # str(merged_seurat_filtered) #
    
    # 5. Linear Dimensionality Reduction ----------
    merged_seurat_filtered <- RunPCA(merged_seurat_filtered, features = VariableFeatures(object = merged_seurat_filtered))
    
    # ___Examine and visualize PCA results a few different ways -------
    # print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5) #
    
    # ___plot-1 --------
    plots$vizdimloadings1 <- VizDimLoadings(merged_seurat_filtered, dims = 1:2, reduction = "pca")
    
    # ___plot-2 --------
    plots$dimplot3 <- DimPlot(merged_seurat_filtered, reduction = "pca")
    
    # ___plot-3 heatmap -------
    # allows for easy exploration of the primary sources of heterogeneity in a dataset
    # and can be useful when trying to decide which PCs to include for further downstream analyses
    df$dimheatmap1 <- merged_seurat_filtered
    df$dimheatmap2 <- merged_seurat_filtered
    
    # ___to dertermine "dimensionality" of the dataset -------
    # essentially determine how many PCs to consider - we would ideally want to consider PCs that show maximum variations
    
    # JackStraw Procedure!
    # identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
    # NOTE: This process can take a long time for big datasets, comment out for expediency. More
    # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
    # computation time
    
    merged_seurat_filtered <- JackStraw(merged_seurat_filtered, num.replicate = 100)
    merged_seurat_filtered <- ScoreJackStraw(merged_seurat_filtered, dims = 1:20)
    
    plots$jackstrawplot2 <- JackStrawPlot(merged_seurat_filtered, dims = 1:20)
    # The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 
    # ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).
    
    # An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function).
    plots$elbowplot3 <- ElbowPlot(merged_seurat_filtered)
    # from the plot, it looks like majority of true signal is captured in the first 15 PCs.
    # PCs to consider = 15
    
    # 6. Cluster cells --------
    merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
    
    # The FindClusters() function contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
    # We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
    # Optimal resolution often increases for larger datasets. 
    merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 0.5)
    
    # Look at cluster IDs of the first 5 cells
    # head(Idents(merged_seurat_filtered), 5) #
    
    # 7. Run non-linear dimensional reduction (UMAP/tSNE) ---------
    merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)
    
    # note that you can set `label = TRUE` or use the LabelClusters function to help label
    # individual clusters
    plots$dimplot_umap <- DimPlot(merged_seurat_filtered, reduction = "umap")
    # plots$b1 <- DimPlot(merged_seurat_filtered, reduction = "umap.cca", label = TRUE) #
    # plots$b2 <- DimPlot(merged_seurat_filtered, reduction = "umap.cca", group.by = "Location", label = TRUE) #
    # plots$b3 <- DimPlot(merged_seurat_filtered, reduction = "umap.cca", group.by = "singleR.labels", label = TRUE) #
    
    df$merged_seurat_filtered <- merged_seurat_filtered
    
  })
  
  ####______________________________________________________________________####
  #Plots for Non integration 
  output$vlnplot6 <- renderPlot({ plots$vlnplot6 })
  
  output$combineplot8 <- renderPlot({ plots$plot1 + plots$plot2 })
  
  output$combineplot9 <- renderPlot({ plots$plot3_sub + plots$plot4_sub })
  
  output$vizdimloadings1 <- renderPlot({ plots$vizdimloadings1 })
  
  output$dimplot3 <- renderPlot({ plots$dimplot3 })
  
  output$dimheatmap1 <- renderPlot({ DimHeatmap(df$dimheatmap1, dims = 1, cells = 500, balanced = TRUE) })
  
  output$dimheatmap2 <- renderPlot({ DimHeatmap(df$dimheatmap2, dims = 1:5, cells = 500, balanced = TRUE) })
  
  output$jackstrawplot2 <- renderPlot({ plots$jackstrawplot2 })
  
  output$elbowplot3 <- renderPlot({ plots$elbowplot3 })
  
  output$dimplot_umap <- renderPlot({ plots$dimplot_umap })
  
  output$vlnplot16 <- renderPlot({ 
    req(input$vlnplot16_selected)
    req(df$merged_seurat_filtered)
    
    VlnPlot(df$merged_seurat_filtered, features = input$vlnplot16_selected, group.by = 'seurat_clusters') 
  })
  
  output$vlnplot17 <- renderPlot({ 
    req(input$vlnplot17_selected)
    req(df$merged_seurat_filtered)
    
    VlnPlot(df$merged_seurat_filtered, features = input$vlnplot17_selected, slot = "counts", log = TRUE)
  })
  
  output$featureplot9 <- renderPlot({
    req(input$featureplot9_selected)
    req(df$merged_seurat_filtered)
    
    if(length(input$featureplot9_selected) == 2) {
      FeaturePlot(df$merged_seurat_filtered, features = input$featureplot9_selected, blend = TRUE) 
    }
  })
  
  # Add to report button
  observeEvent(input$add_nonintegration, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Violin Plot (Non-integration)", align = "center"), 
                 plotOutput("vlnplot6Report", height = 500), br(), br(), 
                 h3("Combined Plot (Non-integration)", align = "center"), 
                 plotOutput("combineplot8Report", height = 500), br(), br(), 
                 h3("Combined Plot (Non-integration)", align = "center"), 
                 plotOutput("combineplot9Report", height = 500), br(), br(), 
                 h3("Visualise Dimensional Reduction Genes (Non-integration)", align = "center"), 
                 plotOutput("vizdimloadings1Report", height = 500), br(), br(), 
                 h3("Dimplot (Non-integration)", align = "center"), 
                 plotOutput("dimplot3Report", height = 500), br(), br(), 
                 h3("Dim Heatmap (Non-integration)", align = "center"), 
                 plotOutput("dimheatmap1Report", height = 500), br(), br(), 
                 h3("Dim Heatmap (Non-integration)", align = "center"), 
                 plotOutput("dimheatmap2Report", height = 500), br(), br(), 
                 h3("Jack Straw Plot (Non-integration)", align = "center"), 
                 plotOutput("jackstrawplot2Report", height = 500), br(), br(), 
                 h3("Elbow Plot (Non-integration)", align = "center"), 
                 plotOutput("elbowplot3Report", height = 500), br(), br(), 
                 h3("Dimplot UMAP (Non-integration)", align = "center"), 
                 plotOutput("dimplot_umapReport", height = 500), br(), br(), 
                 h3("Violin Plot (Non-integration)", align = "center"), 
                 plotOutput("vlnplot16Report", height = 500), br(), br(), 
                 h3("Violin Plot (Non-integration)", align = "center"), 
                 plotOutput("vlnplot17Report", height = 500), br(), br(),
                 h3("Feature Plot (Non-integration)", align = "center"), 
                 plotOutput("featureplot9Report", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$vlnplot6Report <- renderPlot({ plots$vlnplot6 })
    
    output$combineplot8Report <- renderPlot({ plots$plot1 + plots$plot2 })
    
    output$combineplot9Report <- renderPlot({ plots$plot3_sub + plots$plot4_sub })
    
    output$vizdimloadings1Report <- renderPlot({ plots$vizdimloadings1 })
    
    output$dimplot3Report <- renderPlot({ plots$dimplot3 })
    
    output$dimheatmap1Report <- renderPlot({ DimHeatmap(df$dimheatmap1, dims = 1, cells = 500, balanced = TRUE) })
    
    output$dimheatmap2Report <- renderPlot({ DimHeatmap(df$dimheatmap2, dims = 1:5, cells = 500, balanced = TRUE) })
    
    output$jackstrawplot2Report <- renderPlot({ plots$jackstrawplot2 })
    
    output$elbowplot3Report <- renderPlot({ plots$elbowplot3 })
    
    output$dimplot_umapReport <- renderPlot({ plots$dimplot_umap })
    
    output$vlnplot16Report <- renderPlot({ 
      req(input$vlnplot16_selected)
      req(df$merged_seurat_filtered)
      
      VlnPlot(df$merged_seurat_filtered, features = input$vlnplot16_selected, group.by = 'seurat_clusters') 
    })
    
    output$vlnplot17Report <- renderPlot({ 
      req(input$vlnplot17_selected)
      req(df$merged_seurat_filtered)
      
      VlnPlot(df$merged_seurat_filtered, features = input$vlnplot17_selected, slot = "counts", log = TRUE)
    })
    
    output$featureplot9Report <- renderPlot({
      req(input$featureplot9_selected)
      req(df$merged_seurat_filtered)
      
      if(length(input$featureplot9_selected) == 2) {
        FeaturePlot(df$merged_seurat_filtered, features = input$featureplot9_selected, blend = TRUE) 
      }
    })
  })
  
  
  # Update Gene SelectizeInput
  observeEvent(input$navbar, {
    req(input$navbar)
    req(df$all.genes)
    
    if(input$navbar == "nonintegration_tab" & count$nonintegration_tab == 0 & df$complete) {
      updateSelectizeInput(session, "vlnplot16_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "vlnplot17_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "featureplot9_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1:2], 
                           server = TRUE)
    }
  })
  
  ####___________________________________________________________________________####
  #DO YOU WANT TO INTEGRATE? THEN KEEP FOLLOWING,ELSE SCROLL DOWN TO THE merged_seurat_filtered 
  #CCA Integration (NEW TAB)
  observeEvent(input$file3, {
    req(input$file3)
    
    df$mode <- "integrated"
    df$complete <- FALSE
    
    # Track which tab has been opened
    count <- reactiveValues(
      nonintegration_tab = 0,
      CCA_tab = 0,
      autoannotation_tab = 0,
      diff_expression_tab = 0,
      pseudobulk_tab = 0, 
      compare_tab = 0
    )
    
    # Set idents
    df$idents <- sub("_.*", "", input$file3$name) %>% unique()
    
    files <- input$file3
    dirs <- files$name
    
    # Create a named list to store Seurat objects
    seurat_objects <- list()
    
    for (x in dirs) {
      # Read the CSV file
      cts <- read.csv(files$datapath[which(files$name == x)], row.names = 1, check.names = FALSE)
      
      # Create Seurat object
      seurat_objects[[x]] <- CreateSeuratObject(counts = cts)
    }
    
    # merge datasets
    merged_seurat <- merge(x = seurat_objects[[1]], y = seurat_objects[-1], 
                           add.cell.ids = dirs, project = 'OVA')
    
    #If there are breaks in the index then run this
    # Assuming ls() returns a character vector with the specified names
    
    # merged_seurat #
    
    # QC & filtering -----------------------
    
    # View(merged_seurat@meta.data) #
    # create a sample column
    merged_seurat$sample <- rownames(merged_seurat@meta.data)
    
    # Split sample column
    merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', 
                                        into = c('Location', 'Type', 'CD45State/CAN', 'Barcode'), sep = '_')
    
    ####______________________________________________________________________####
    #HarmonyIntegration (NEW TAB) (Skip to line 1087 for continuation to CCA integration)
    
    merged_seurat$mito.percent <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')
    # View(merged_seurat@meta.data) # 
    # explore QC
    
    # filter
    # merged_seurat #
    merged_seurat.filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                       nFeature_RNA > 200 & 
                                       mito.percent < 5)
    
    # standard workflow steps
    merged_seurat.filtered <- NormalizeData(merged_seurat.filtered)
    merged_seurat.filtered <- FindVariableFeatures(merged_seurat.filtered)
    merged_seurat.filtered <- ScaleData(merged_seurat.filtered)
    merged_seurat.filtered <- RunPCA(merged_seurat.filtered)
    plots$elbowplot1 <- ElbowPlot(merged_seurat.filtered)
    merged_seurat.filtered <- RunUMAP(merged_seurat.filtered, dims = 1:20, reduction = 'pca')
    
    plots$before <- DimPlot(merged_seurat.filtered, reduction = 'umap', group.by = 'CD45State/CAN')
    
    # run Harmony ----------- (Only when 'CD45State/CAN' has more than 1 level)
    if(length(levels(factor(merged_seurat.filtered@meta.data$`CD45State/CAN`))) > 1) {
      merged_seurat.harmony <- merged_seurat.filtered %>%
        RunHarmony(group.by.vars = 'CD45State/CAN', plot_convergence = FALSE)
      
      merged_seurat.harmony.embed <- Embeddings(merged_seurat.harmony, "harmony")
      # merged_seurat.harmony.embed[1:10,1:10] #
      
      # Do UMAP and clustering using ** Harmony embeddings instead of PCA **
      merged_seurat.harmony <- merged_seurat.harmony %>%
        RunUMAP(reduction = 'harmony', dims = 1:20) %>%
        FindNeighbors(reduction = "harmony", dims = 1:20) %>%
        FindClusters(resolution = 0.5)
      
      # visualize 
      plots$after <- DimPlot(merged_seurat.harmony, reduction = 'umap', group.by = 'CD45State/CAN')
    }
    
    ####______________________________________________________________________####
    
    # calculate mitochondrial percentage
    merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')
    
    # explore QC
    plots$vlnplot <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
    
    # filtering
    merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 500 &
                                       nFeature_RNA > 200 &
                                       mitoPercent < 5)
    
    df$all.genes <- rownames(merged_seurat_filtered)
    
    # Update selectizeInput immediately with new choices and default value once the gene names are determined
    ## CCA Integration
    if(count$CCA_tab == 0) {
      updateSelectizeInput(session, "vlnplot3_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "vlnplot4_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "vlnplot5_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "featureplot1_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "featureplot2_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      count$CCA_tab <- 1
    }
    
    # view(merged_seurat_filtered@meta.data) #
    
    # merged_seurat #
    
    # perform standard workflow steps to figure out if we see any batch effects --------
    merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
    merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
    
    top10 <- head(VariableFeatures(merged_seurat_filtered), 10)
    # top10 #
    
    # ___plot variable features with and without labels ---------
    plots$plot3 <- VariableFeaturePlot(merged_seurat_filtered)
    plots$plot4 <- LabelPoints(plot = plots$plot3, points = top10, repel = TRUE)
    
    merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
    merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
    plots$elbowplot2 <- ElbowPlot(merged_seurat_filtered)
    merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
    merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
    merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)
    
    # plot
    plots$p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Location')
    plots$p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
                        cols = c('red','green','blue'))
    plots$P6 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'CD45State/CAN')
    
    # perform integration to correct for batch effects ------
    #for V5
    merged_seurat_filtered <- IntegrateLayers(
      object = merged_seurat_filtered, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.cca",
      verbose = FALSE
    )
    
    #FOR V3
    #obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Location')
    #for(i in 1:length(obj.list)){
    #  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
    ##  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
    #}
    
    
    # select integration features
    #features <- SelectIntegrationFeatures(object.list = obj.list)
    
    # find integration anchors (CCA)
    #anchors <- FindIntegrationAnchors(object.list = obj.list,
    # anchor.features = features)
    
    # Check the number of dimensions available in the integrated.cca reduction
    
    # integrate data
    #seurat.integrated <- IntegrateData(anchorset = anchors)
    merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, reduction = "integrated.cca", dims = 1:20)
    merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 2, cluster.name = "cca_clusters")
    merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
    plots$p5 <- DimPlot(
      merged_seurat_filtered,
      reduction = "umap.cca",
      group.by = c("Location", "Type", "CD45State/CAN"),
      combine = FALSE, label.size = 2
    )
    
    # Scale data, run PCA and UMAP and visualize integrated data
    merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
    merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
    merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)
    
    plots$p3 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Location')
    plots$p4 <- DimPlot(merged_seurat_filtered, reduction = 'umap.cca', group.by = 'Type',
                        cols = c('red','green','blue'))
    
    plots$p8 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'seurat_clusters')
    # plots$p7 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'seurat_clusters_manual') #
    
    ####______________________________________________________####
    #Manual Annotation (IGNORE THIS)
    
    # Create a data frame with the cluster mapping
    cluster_mapping <- data.frame(
      Cluster = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29),
      Desc = c(
        "Epithelial_Cells", "T_Cells", "T_Cells", "T_Cells", "T_Cells", "NK Cells", 
        "T_Cells", "T_Cells", "T_Cells", "T_Cells", "NK Cells", 
        "NK Cells", "T_Cells", "Epithelial_Cells", "T_Cells", 
        "Monocyte", "T_Cells", "T_Cells", "B_Cells", "Monocyte", 
        "Monocyte", "Epithelial_Cells", "T_Cells", 
        "Epithelial_Cells", "Monocyte", "Fibroblast", 
        "B_Cells", "Epithelial_Cells", "Epithelial_Cells", 
        "Monocyte"
      )
    )
    
    # Merge the Seurat object with the cluster mapping
    merged_seurat_filtered$seurat_clusters_manual <- factor(
      merged_seurat_filtered$seurat_clusters,
      levels = cluster_mapping$Cluster,
      labels = cluster_mapping$Desc
    )
    
    # Plot the updated clusters
    plots$dimplot1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'seurat_clusters_manual')
    #________________________________________________
    
    df$vlnplot3 <- merged_seurat_filtered
    # view(merged_seurat_filtered@meta.data) #
    #______________________________________________________
    
    merged_seurat_filtered <- JackStraw(merged_seurat_filtered, num.replicate = 100)
    merged_seurat_filtered <- ScoreJackStraw(merged_seurat_filtered, dims = 1:20)
    
    plots$jackstrawplot <- JackStrawPlot(merged_seurat_filtered, dims = 1:20)
    
    #_______________________________________________________________________________________
    #before going for clusters 
    
    merged_seurat_filtered <- JoinLayers(merged_seurat_filtered)
    cluster1.markers <- FindMarkers(merged_seurat_filtered, ident.1 = 2, min.pct = 0.25)
    
    #Find all markers
    merged_seurat_filtered.markers <- FindAllMarkers(merged_seurat_filtered,
                                                     logfc.threshold = 0.25,
                                                     min.pct = 0.25,
                                                     only.pos = TRUE,
                                                     test.use = 'DESeq2',
                                                     slot = 'counts')
    
    features <- merged_seurat_filtered@commands$RunPCA.RNA$features
    df$vlnplot4 <-merged_seurat_filtered
    
    # ___VlnPlot() - you can plot raw counts as well ---------
    df$vlnplot5 <- merged_seurat_filtered
    
    # plot <- FeaturePlot(merged_seurat_filtered, features = c("CD45RA")) #
    # plots$plot <- HoverLocator(plot = plot, #
    #                          information = FetchData(merged_seurat_filtered, vars = c("ident", "PC_1", "nFeature_RNA"))) #
    
    df$featureplot1 <- merged_seurat_filtered
    # top10 <- merged_seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) #
    # plots$heatmap1 <- DoHeatmap(merged_seurat_filtered, features = top10$gene) #
    
    markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                         "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                         "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
    plots$dotplot1 <- DotPlot(merged_seurat_filtered, features = markers.to.plot, cols = c("blue", "red", "green", "yellow", "purple"), dot.scale = 8, split.by = "Location") +
      RotatedAxis() 
    plots$dimplot2 <- DimPlot(merged_seurat_filtered, reduction = "umap", split.by = "seurat_clusters_manual")
    
    df$featureplot2 <- merged_seurat_filtered
    
    # merged_seurat_filtered #
    
    # plots$featureplot3 <- FeaturePlot(merged_seurat_filtered, #
    #                                  features = c("GZMK", "CD19", "MUC1", "FAP", "CD14", "NCR1"), #
    #                                  split.by = "seurat_clusters_manual", #
    #                                  max.cutoff = 3, #
    #                                  cols = c("grey","red"), #
    #                                  reduction = "umap") #
    
    df$merged_seurat_filtered <- merged_seurat_filtered
    
  })
  
  ####______________________________________________________________________####
  #Plots for Harmony integration 
  output$elbowplot1 <- renderPlot({ plots$elbowplot1 })
  
  output$before <- renderPlot({ plots$before })
  
  output$after <- renderPlot({ plots$after })
  
  observeEvent(input$add_harmony, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Elbow Plot (Harmony Integration)", align = "center"), 
                 plotOutput("elbowplot1Report", height = 500), br(), br(), 
                 h3("Before (Harmony Integration)", align = "center"), 
                 plotOutput("beforeReport", height = 500), br(), br(), 
                 h3("After (Harmony Integration)", align = "center"), 
                 plotOutput("afterReport", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$elbowplot1Report <- renderPlot({ plots$elbowplot1 })
    
    output$beforeReport <- renderPlot({ plots$before })
    
    output$afterReport <- renderPlot({ plots$after })
  })
  
  ####______________________________________________________________________####
  #Plots for CCA integration
  
  output$vlnplot <- renderPlot({ plots$vlnplot })
  
  output$combineplot1 <- renderPlot({ plots$plot3 + plots$plot4 })
  
  output$elbowplot2 <- renderPlot({ plots$elbowplot2 })
  
  output$combineplot2 <- renderPlot({ grid.arrange(plots$p1, plots$p2, plots$P6, ncol = 2, nrow = 2) })
  
  output$combineplot3 <- renderPlot({ grid.arrange(plots$p1, plots$p2, plots$p3, plots$p4, ncol = 3, nrow = 3) })
  
  output$combineplot4 <- renderPlot({ plots$p1 + plots$p2 + plots$p3 + plots$p4 + plots$p5 })
  
  output$combineplot5 <- renderPlot({ plots$p3 + plots$p4 })
  
  output$combineplot6 <- renderPlot({ plots$p4 + plots$p5 })
  
  output$combineplot7 <- renderPlot({ plots$p3 }) #or p7 or p8
  
  output$dimplot1 <- renderPlot({ plots$dimplot1 })
  
  output$vlnplot3 <- renderPlot({ 
    req(input$vlnplot3_selected)
    req(df$vlnplot3)
    
    VlnPlot(df$vlnplot3, features = input$vlnplot3_selected, group.by = 'seurat_clusters_manual') 
  })
  
  output$jackstrawplot <- renderPlot({ plots$jackstrawplot })
  
  output$vlnplot4 <- renderPlot({ 
    req(input$vlnplot4_selected)
    req(df$vlnplot4)
    
    VlnPlot(df$vlnplot4, features = input$vlnplot4_selected, group.by = 'seurat_clusters_manual') 
  })
  
  output$vlnplot5 <- renderPlot({ 
    req(input$vlnplot5_selected)
    req(df$vlnplot5)
    
    VlnPlot(df$vlnplot5, features = input$vlnplot5_selected, slot = "counts", log = TRUE)
  })
  
  output$featureplot1 <- renderPlot({
    req(input$featureplot1_selected)
    req(df$featureplot1)
    
    if(length(input$featureplot1_selected) == 2) {
      FeaturePlot(df$featureplot1, features = input$featureplot1_selected, blend = TRUE) 
    }
  })
  
  output$dotplot1 <- renderPlot({ plots$dotplot1 })
  
  output$dimplot2 <- renderPlot({ plots$dimplot2 })
  
  output$featureplot2 <- renderPlot({ 
    req(input$featureplot2_selected)
    req(df$featureplot2)
    
    FeaturePlot(df$featureplot2, 
                features = input$featureplot2_selected, 
                split.by = "seurat_clusters_manual", 
                max.cutoff = 3, 
                cols = c("grey","red"), 
                reduction = "umap")  
  })
  
  #output$featureplot3 <- renderPlot({ plots$featureplot3 })
  
  observeEvent(input$add_cca, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Violin Plot (CCA Integration)", align = "center"), 
                 plotOutput("vlnplotReport", height = 500), br(), br(), 
                 h3("Combined Plot (CCA Integration)", align = "center"), 
                 plotOutput("combineplot1Report", height = 500), br(), br(), 
                 h3("Elbow Plot (CCA Integration)", align = "center"), 
                 plotOutput("elbowplot2Report", height = 500), br(), br(), 
                 h3("Combined Plot (CCA Integration)", align = "center"), 
                 plotOutput("combineplot2Report", height = 500), br(), br(), 
                 h3("Combined Plot (CCA Integration)", align = "center"), 
                 plotOutput("combineplot3Report", height = 500), br(), br(), 
                 h3("Combined Plot (CCA Integration)", align = "center"), 
                 plotOutput("combineplot4Report", height = 500), br(), br(), 
                 h3("Combined Plot (CCA Integration)", align = "center"), 
                 plotOutput("combineplot5Report", height = 500), br(), br(), 
                 h3("Combined Plot (CCA Integration)", align = "center"), 
                 plotOutput("combineplot6Report", height = 500), br(), br(), 
                 h3("Combined Plot (CCA Integration)", align = "center"), 
                 plotOutput("combineplot7Report", height = 500), br(), br(), 
                 h3("Dimplot (CCA Integration)", align = "center"), 
                 plotOutput("dimplot1Report", height = 500), br(), br(), 
                 h3("Violin Plot (CCA Integration)", align = "center"), 
                 plotOutput("vlnplot3Report", height = 500), br(), br(),
                 h3("Jack Straw Plot (CCA Integration)", align = "center"), 
                 plotOutput("jackstrawplotReport", height = 500), br(), br(),
                 h3("Violin Plot (CCA Integration)", align = "center"), 
                 plotOutput("vlnplot4Report", height = 500), br(), br(),
                 h3("Violin Plot (CCA Integration)", align = "center"), 
                 plotOutput("vlnplot5Report", height = 500), br(), br(),
                 h3("Feature Plot (CCA Integration)", align = "center"), 
                 plotOutput("featureplot1Report", height = 500), br(), br(),
                 h3("Dot Plot (CCA Integration)", align = "center"), 
                 plotOutput("dotplot1Report", height = 500), br(), br(),
                 h3("Dimplot (CCA Integration)", align = "center"), 
                 plotOutput("dimplot2Report", height = 500), br(), br(),
                 h3("Feature Plot (CCA Integration)", align = "center"), 
                 plotOutput("featureplot2Report", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$vlnplotReport <- renderPlot({ plots$vlnplot })
    
    output$combineplot1Report <- renderPlot({ plots$plot3 + plots$plot4 })
    
    output$elbowplot2Report <- renderPlot({ plots$elbowplot2 })
    
    output$combineplot2Report <- renderPlot({ grid.arrange(plots$p1, plots$p2, plots$P6, ncol = 2, nrow = 2) })
    
    output$combineplot3Report <- renderPlot({ grid.arrange(plots$p1, plots$p2, plots$p3, plots$p4, ncol = 3, nrow = 3) })
    
    output$combineplot4Report <- renderPlot({ plots$p1 + plots$p2 + plots$p3 + plots$p4 + plots$p5 })
    
    output$combineplot5Report <- renderPlot({ plots$p3 + plots$p4 })
    
    output$combineplot6Report <- renderPlot({ plots$p4 + plots$p5 })
    
    output$combineplot7Report <- renderPlot({ plots$p3 }) #or p7 or p8
    
    output$dimplot1Report <- renderPlot({ plots$dimplot1 })
    
    output$vlnplot3Report <- renderPlot({ 
      req(input$vlnplot3_selected)
      req(df$vlnplot3)
      
      VlnPlot(df$vlnplot3, features = input$vlnplot3_selected, group.by = 'seurat_clusters_manual') 
    })
    
    output$jackstrawplotReport <- renderPlot({ plots$jackstrawplot })
    
    output$vlnplot4Report <- renderPlot({ 
      req(input$vlnplot4_selected)
      req(df$vlnplot4)
      
      VlnPlot(df$vlnplot4, features = input$vlnplot4_selected, group.by = 'seurat_clusters_manual') 
    })
    
    output$vlnplot5Report <- renderPlot({ 
      req(input$vlnplot5_selected)
      req(df$vlnplot5)
      
      VlnPlot(df$vlnplot5, features = input$vlnplot5_selected, slot = "counts", log = TRUE)
    })
    
    output$featureplot1Report <- renderPlot({
      req(input$featureplot1_selected)
      req(df$featureplot1)
      
      if(length(input$featureplot1_selected) == 2) {
        FeaturePlot(df$featureplot1, features = input$featureplot1_selected, blend = TRUE) 
      }
    })
    
    output$dotplot1Report <- renderPlot({ plots$dotplot1 })
    
    output$dimplot2Report <- renderPlot({ plots$dimplot2 })
    
    output$featureplot2Report <- renderPlot({ 
      req(input$featureplot2_selected)
      req(df$featureplot2)
      
      FeaturePlot(df$featureplot2, 
                  features = input$featureplot2_selected, 
                  split.by = "seurat_clusters_manual", 
                  max.cutoff = 3, 
                  cols = c("grey","red"), 
                  reduction = "umap")  
    })
  })
  
  #_____________________________________________________________________
  #singleCellHaystack WIP (IGNORED, NOT INCLUDED)
  
  ####_____________________________________________________________________####
  # SECOND HALF OF PIPELINE
  
  observeEvent(df$merged_seurat_filtered, {
    merged_seurat_filtered <- df$merged_seurat_filtered
    print(-1)
    ####_____________________________________________________________________________________________####
    #AUTO ANNOTATION (NEW TAB)
    # get reference data -----------
    ref <- celldex::HumanPrimaryCellAtlasData()
    
    # View(as.data.frame(colData(ref))) #
    # View(merged_seurat_filtered@meta.data) #
    
    # expression values are log counts (log normalized counts)
    
    # run SingleR (default mode) ---------
    # default for SingleR is to perform annotation of each individual cell in the test dataset
    
    pbmc_counts <- GetAssayData(merged_seurat_filtered, layer = 'counts')
    
    pred <- SingleR(test = pbmc_counts,
                    ref = ref,
                    labels = ref$label.main)
    # pred #
    
    merged_seurat_filtered$singleR.labels <- pred$labels[match(rownames(merged_seurat_filtered@meta.data), rownames(pred))]
    if(df$mode == "single") {
      plots$a2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'singleR.labels')
    } else {
      plots$a2 <- DimPlot(merged_seurat_filtered, reduction = 'umap.cca', group.by = 'singleR.labels')
    }
    
    # Annotation diagnostics ----------
    # ...Based on the scores within cells -----------
    # pred #
    # pred$scores #
    
    plots$scoreheatmap1 <- plotScoreHeatmap(pred)
    
    # ...Based on deltas across cells ----------
    plots$deltadistribution <- plotDeltaDistribution(pred)
    
    # ...Comparing to unsupervised clustering ------------
    tab <- table(Assigned = pred$labels, Clusters = merged_seurat_filtered$seurat_clusters)
    plots$pheatmap1 <- pheatmap(log10(tab + 10), color = colorRampPalette(c('white','blue'))(20))
    #_____________________________________________________________________OR______
    #____________________________________________________________________________--
    #___________________________________________________________________________
    
    # run SingleR with multiple reference datasets (default mode) ---------
    
    # for pbmc data, we will use two datasets
    hpca <- celldex::HumanPrimaryCellAtlasData()
    dice <- celldex::DatabaseImmuneCellExpressionData()
    
    # ...1. Strategy 1: Using reference-specific labels ----------
    # hpca$label.main #
    # dice$label.main #
    
    # adding ref info to labels
    hpca$label.main <- paste0('HPCA.', hpca$label.main)
    dice$label.main <- paste0('DICE.', dice$label.main)
    
    # create a combined ref based on shared genes
    shared <- intersect(rownames(hpca), rownames(dice))
    combined <- cbind(hpca[shared,], dice[shared,])
    # combined #
    # combined$label.main #
    
    # run singleR using combined ref
    # savings counts into a separate object
    pbmc_counts <- GetAssayData(merged_seurat_filtered, layer = 'counts')
    
    com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
    # table(com.res1$labels) #
    
    merged_seurat_filtered$com.res1.labels <- com.res1[match(rownames(merged_seurat_filtered@meta.data), rownames(com.res1)), 'labels']
    # View(merged_seurat_filtered@meta.data) #
    
    plots$dimplot4 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)
    
    # ...2. Strategy 2: Comparing scores across references ----------
    # hpca$label.main #
    # dice$label.main #
    hpca$label.main <- gsub('HPCA\\.','', hpca$label.main)
    dice$label.main <- gsub('DICE\\.','', dice$label.main)
    
    com.res2 <- SingleR(test = pbmc_counts, 
                        ref = list(HPCA = hpca, DICE = dice),
                        labels = list(hpca$label.main, dice$label.main))
    
    # Check the final label from the combined assignment.
    # table(com.res2$labels) #
    
    # which reference scored best for which label?
    grouping <- paste0(com.res2$labels,'.', com.res2$reference)
    best_ref <- as.data.frame(split(com.res2, grouping))
    # view(best_ref) #
    
    # get de. genes from each individual references
    # metadata(com.res2$orig.results$HPCA)$de.genes #
    # metadata(com.res2$orig.results$DICE)$de.genes #
    
    merged_seurat_filtered$com.res2.labels <- com.res2[match(rownames(merged_seurat_filtered@meta.data), rownames(com.res2)), 'labels']
    # View(merged_seurat_filtered@meta.data) #
    
    # Combined diagnostics
    plots$scoreheatmap2 <- plotScoreHeatmap(com.res2)
    plots$dimplot5 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'com.res2.labels', label = TRUE)
    
    # ...3. Strategy 3: Using Harmonized Labels ----------
    
    hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')
    dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')
    
    # Using the same sets of genes:
    shared <- intersect(rownames(hpca.ont), rownames(dice.ont))
    hpca.ont <- hpca.ont[shared,]
    dice.ont <- dice.ont[shared,]
    
    # Showing the top 10 most frequent terms:
    # tail(sort(table(hpca.ont$label.ont)),10) #
    # tail(sort(table(dice.ont$label.ont)), 10) #
    
    # using label.ont instead on label.main while running SingleR
    
    com.res3 <- SingleR(test = pbmc_counts,
                        ref = list(HPCA = hpca.ont, DICE = dice.ont),
                        labels = list(hpca.ont$label.ont, dice.ont$label.ont))
    
    # table(com.res3$labels) #
    
    merged_seurat_filtered$com.res3.labels <- com.res3[match(rownames(merged_seurat_filtered@meta.data), rownames(com.res3)), 'labels']
    # View(merged_seurat_filtered@meta.data) #
    plots$dimplot6 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'com.res3.labels', label = TRUE)
    
    # How to map cell ontology terms? ----------------
    
    # colData(hpca.ont) #
    # colData(dice.ont) #
    print(0)
    hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")
    hpca.mapping <- read.delim(hpca.fle, header = F)
    
    #___________________________________________________________________________
    
    # 8. Finding differentially expressed features (cluster biomarkers) ---------
    # Seurat can help you find markers that define clusters via differential expression. 
    
    # ___find all markers of cluster 1 --------
    # cluster1.markers <- FindMarkers(merged_seurat_filtered, ident.1 = 2, min.pct = 0.25) #
    # head (cluster1.markers, n = 5) #
    
    # cluster1.markers_sorted <- cluster1.markers %>% #
    #  arrange(desc(avg_log2FC)) #
    # head (cluster1.markers_sorted, n = 5) #
    
    # ___find all markers distinguishing cluster 5 from clusters 0 and 3 --------
    # cluster5.markers <- FindMarkers(merged_seurat_filtered, ident.1 = 1, ident.2 = c(0, 3), min.pct = 0.25) #
    # head(cluster5.markers, n = 5) #
    
    # cluster5.markers_sorted <- cluster5.markers %>% #
    #  arrange(desc(avg_log2FC)) #
    # head(cluster5.markers_sorted, n = 20) #
    
    # ___find markers for every cluster compared to all remaining cells, report only the positive ones ---------
    merged_seurat_filtered.markers <- FindAllMarkers(merged_seurat_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    # merged_seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) #
    
    #Find all markers
    #FindAllMarkers(merged_seurat_filtered, #
    #               logfc.threshold = 0.25, #
    #               min.pct = 0.1, #
    #               only.pos = TRUE, #
    #               test.use = 'DESeq2', #
    #               slot = 'counts') #
    
    #Find Conserved Marker
    
    # markers_cluster3 <- FindConservedMarkers(merged_seurat_filtered, #
    #                                         ident.1 = 1, #
    #                                         grouping.var = 'Location') #
    
    # head(markers_cluster3) #
    # write.csv(markers_cluster3, file = 'conservedmarkers.csv', row.names = TRUE) #
    
    #_________________________________________________________
    
    #________________________________________________________
    # Assuming 'singleR.labels' is the column you're using for grouping
    # Assuming 'CD5L' is the feature you want to plot
    
    # Count the number of cells in each group
    group_counts <- table(merged_seurat_filtered$singleR.labels)
    
    # Identify groups with at least 10 cells
    valid_groups <- names(group_counts[group_counts >= 40])
    
    # Filter Seurat object to include only valid groups
    filtered_seurat <- subset(merged_seurat_filtered, subset = singleR.labels %in% valid_groups)
    
    # Create the violin plot
    df$vlnplot7 <- filtered_seurat
    #__________________________________________________________
    # 9. Visualization ---------
    # VlnPlot() (shows expression probability distributions across clusters)
    # FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. 
    # RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.
    
    # str(merged_seurat_filtered) #
    features <- merged_seurat_filtered@commands$RunPCA.RNA$features
    df$vlnplot8 <- merged_seurat_filtered
    df$vlnplot9 <- merged_seurat_filtered
    
    # ___VlnPlot() - you can plot raw counts as well ---------
    df$vlnplot10 <- merged_seurat_filtered
    
    
    # ___FeaturePlot()- visualize feature expression in low-dimensional space ---------
    # plots$featureplot4 <- FeaturePlot(merged_seurat_filtered, features = features[1:5]) #
    df$featureplot5 <- merged_seurat_filtered
    
    # Visualize co-expression of two features simultaneously
    #plots$featureplot6 <- FeaturePlot(merged_seurat_filtered, features = c("GZMK", "CCR7"), blend = TRUE)
    
    # ___interactive plots --------
    # Include additional data to display alongside cell names by passing in a data frame of
    # information Works well when using FetchData
    # works only with one feature
    # plot <- FeaturePlot(merged_seurat_filtered, features = c("CD5L")) #
    # plots$plot2 <- HoverLocator(plot = plot, #
    #  information = FetchData(merged_seurat_filtered, vars = c("ident", "PC_1", "nFeature_RNA", "seurat_clusters_manual"))) #
    
    # ___doHeatmap() --------
    top10 <- merged_seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    plots$doheatmap1 <- DoHeatmap(merged_seurat_filtered, features = top10$gene) + NoLegend()
    
    # ___RidgePlot() - Visualize single cell expression distribution in each cluster -------
    # plots$ridgeplot1 <- RidgePlot(merged_seurat_filtered, features = features[1:5], ncol=2) #
    df$ridgeplot2 <- merged_seurat_filtered
    # ___Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature --------
    # in each cluster. The color represents the average expression level
    # plots$dotplot2 <- DotPlot(merged_seurat_filtered, features = features[1:5]) + RotatedAxis() #
    
    # ___Single cell heatmap of feature expression -------
    plots$doheatmap2 <- DoHeatmap(subset(merged_seurat_filtered, downsample = 100), features = features[1:5], size = 3)
    
    # assigning cell type identity to clusters
    # new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
    #                      "NK", "DC", "Platelet")
    # names(new.cluster.ids) <- levels(pbmc)
    # pbmc <- RenameIdents(pbmc, new.cluster.ids)
    # DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    
    ## let's visualize top features
    #FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')
    
    
    # rename cluster 3 ident
    #Idents(ifnb_harmony)
    #ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')
    
    #DimPlot(ifnb_harmony, reduction = 'umap', label = T)
    
    # cells already have annotations provided in the metadata
    #View(ifnb_harmony@meta.data)
    
    ####______________________________________________________________________####
    #Differential Expression (NEW TAB)
    
    # Normalize the data
    if(df$mode == "integrated") {
      Idents(merged_seurat_filtered) <- merged_seurat_filtered$Location
    }
    
    # Normalize data
    merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)
    
    # Find DE features between
    if(length(df$idents) >= 2) {
      T1.de.markers <- FindMarkers(merged_seurat_filtered, ident.1 = df$idents[1], ident.2 =  df$idents[2])
      
      T1.de.markers <- T1.de.markers[order(-T1.de.markers$avg_log2FC), ]
      T1.de.markers.asc <- T1.de.markers[order(+T1.de.markers$avg_log2FC), ]
      
      df$T1.de.markers <- T1.de.markers
      df$T1.de.markers.asc <- T1.de.markers.asc
    }
    
    # View results
    # head(T1.de.markers) #
    # head(T1.de.markers.asc) #
    df$v1 <- merged_seurat_filtered
    
    ####______________________________________________________________________####
    ##Pseudobulk ((NEW TAB))
    print(1)
    if(df$mode == "integrated" & length(df$idents) >= 2) {
      pseudo_merged_seurat_filtered <- AggregateExpression(merged_seurat_filtered, assays = "RNA", return.seurat = T, group.by = c("Location", "seurat_clusters_manual"))
      
      # tail(Cells(pseudo_merged_seurat_filtered)) #
      
      Idents(pseudo_merged_seurat_filtered) <- pseudo_merged_seurat_filtered$Location
      df$pseudo_merged_seurat <- pseudo_merged_seurat_filtered
      
      bulk.T1.de <- FindMarkers(object = pseudo_merged_seurat_filtered, 
                                ident.1 = df$idents[1], 
                                ident.2 = df$idents[2],
                                test.use = "DESeq2")
      
      bulk.T1.de <- bulk.T1.de[order(-bulk.T1.de$avg_log2FC), ]
      bulk.T1.de.asc <- bulk.T1.de[order(+bulk.T1.de$avg_log2FC), ]
      
      df$bulk.T1.de <- bulk.T1.de
      df$bulk.T1.de.asc <- bulk.T1.de.asc
      
      # head(bulk.T1.de, n = 15) #
      # head(bulk.T1.de.asc, n = 15) #
      
      df$vb1 <- merged_seurat_filtered
      #change gene names as required
    }
    
    ####______________________________________________________________________####
    # Compare between SCDE and PSCDE ((NEW TAB))
    
    if(length(df$idents) >= 2) {
      names(bulk.T1.de) <- paste0(names(bulk.T1.de), ".bulk")
      bulk.T1.de$gene <- rownames(bulk.T1.de)
      
      names(T1.de.markers) <- paste0(names(T1.de.markers), ".sc")
      T1.de.markers$gene <- rownames(T1.de.markers)
      
      merge_dat <- merge(T1.de.markers, bulk.T1.de, by = "gene")
      merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]
      
      # Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
      common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                       merge_dat$p_val.sc < 0.05)]
      only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                        merge_dat$p_val.sc < 0.05)]
      only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                          merge_dat$p_val.sc > 0.05)]
      # print(paste0('# Common: ',length(common))) #
      # print(paste0('# DEG Only present single-cell: ',length(only_sc))) #
      # print(paste0('# DEG Only in bulk: ',length(only_bulk))) #
    }
    
    if(df$mode == "integrated" & length(df$idents >= 2)) {
      # create a new column to annotate sample-condition-celltype in the single-cell dataset
      merged_seurat_filtered$donor_id.stim <- paste0(merged_seurat_filtered$Location, "-", merged_seurat_filtered$seurat_clusters_manual)
      
      # generate violin plot 
      Idents(merged_seurat_filtered) <- merged_seurat_filtered$Location
      # print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc','p_val.bulk')]) #
      
      plots$vlnplot11 <- VlnPlot(merged_seurat_filtered, features = common[1:5], idents = df$idents, group.by = "Location") 
      
      plots$vlnplot12 <- VlnPlot(merged_seurat_filtered, features = common[1:2], idents = df$idents, group.by = "donor_id.stim", ncol = 1) 
      
      plots$vlnplot13 <- VlnPlot(merged_seurat_filtered, features = only_sc[1:5], idents = df$idents, group.by = "Location") 
    }
    
    # print(merge_dat[merge_dat$gene%in%c('KLF2','CCND3'),c('gene','p_val.sc','p_val.bulk')]) #
    
    df$vlnplot14 <- merged_seurat_filtered
    df$vlnplot15 <- merged_seurat_filtered
    
    ####_________________________________________________________________####
    ##Monocle3 Workflow (NEW TAB)
    
    Idents(merged_seurat_filtered) <- merged_seurat_filtered$singleR.labels
    df$e.seu_idents <- unique(levels(Idents(merged_seurat_filtered)))
    
    df$e.seu <- merged_seurat_filtered
    
    # Indicate the process is completed
    df$complete <- TRUE
  })
  
  observeEvent(input$monocle_selected, {
    ####_________________________________________________________________####
    ##Monocle3 Workflow (NEW TAB)
    
    req(df$e.seu)
    req(input$monocle_selected)
    
    merged_seurat_filtered <- df$e.seu
    # unique(merged_seurat_filtered@meta.data$singleR.labels) #
    
    e.seu <- subset(merged_seurat_filtered, idents = input$monocle_selected)
    
    # e.seu #
    # unique(e.seu@meta.data$singleR.labels) #
    
    # pre-processing using seurat
    e.seu <- NormalizeData(e.seu)
    e.seu <- FindVariableFeatures(e.seu)
    e.seu <- ScaleData(e.seu)
    e.seu <- RunPCA(e.seu)
    e.seu <- FindNeighbors(e.seu, dims = 1:30)
    e.seu <- FindClusters(e.seu, resolution = 0.9)
    e.seu <- RunUMAP(e.seu, dims = 1:30, n.neighbors = 50)
    
    plots$a1 <- DimPlot(e.seu, reduction = 'umap', group.by = 'singleR.labels', label = T) 
    plots$a2 <- DimPlot(e.seu, reduction = 'umap', group.by = 'seurat_clusters', label = T)
    
    # ...1 Convert to cell_data_set object ------------------------
    cds <- as.cell_data_set(e.seu)
    # cds #
    
    # to get cell metadata
    # colData(cds) #
    # to gene metdata
    # fData(cds) #
    # rownames(fData(cds))[1:10] #
    
    # since it misses the gene_short_name column, let's add it
    fData(cds)$gene_short_name <- rownames(fData(cds))
    
    # to get counts
    # counts(cds) #
    
    # ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
    # let's use the clustering information have
    
    # assign paritions
    reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
    names(reacreate.partition) <- cds@colData@rownames
    reacreate.partition <- as.factor(reacreate.partition)
    
    cds@clusters$UMAP$partitions <- reacreate.partition
    
    # Assign the cluster info 
    
    list_cluster <- e.seu@active.ident
    cds@clusters$UMAP$clusters <- list_cluster
    
    # Assign UMAP coordinate - cell embeddings
    
    cds@int_colData@listData$reducedDims$UMAP <- e.seu@reductions$umap@cell.embeddings
    
    # plot
    plots$cluster.before.trajectory <- plot_cells(cds,
                                                  color_cells_by = 'cluster',
                                                  label_groups_by_cluster = FALSE,
                                                  group_label_size = 5) +
      theme(legend.position = "right")
    
    plots$cluster.names <- plot_cells(cds,
                                      color_cells_by = "singleR.labels",
                                      label_groups_by_cluster = FALSE,
                                      group_label_size = 5) +
      scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
      theme(legend.position = "right")
    
    # ...3. Learn trajectory graph ------------------------
    cds <- learn_graph(cds, use_partition = FALSE)
    
    plots$plot_cell <- plot_cells(cds,
                                  color_cells_by = 'singleR.labels',
                                  label_groups_by_cluster = FALSE,
                                  label_branch_points = FALSE,
                                  label_roots = FALSE,
                                  label_leaves = FALSE,
                                  group_label_size = 5)
    
    
    # ...4. Order the cells in pseudotime -------------------
    
    # cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 0])) #
    
    # plots$plot_cell2 <- plot_cells(cds, #
    #                                color_cells_by = 'pseudotime', #
    #                                label_groups_by_cluster = FALSE, #
    #                                label_branch_points = FALSE, #
    #                                label_roots = FALSE, #
    #                                label_leaves = FALSE) #
    
    # cells ordered by monocle3 pseudotime
    
    # pseudotime(cds) #
    # cds$monocle3_pseudotime <- pseudotime(cds) #
    # data.pseudo <- as.data.frame(colData(cds)) #
    ##CHANGE THE seurat_clusters if you have actual clusters
    # plots$pseudotime <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) + #
    #  geom_boxplot() #
    
    # ...5. Finding genes that change as a function of pseudotime --------------------
    # deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4) #
    
    #deg_bcells %>% 
    #  arrange(q_value) %>% 
    #  filter(status == 'OK') %>% 
    #  head()
    
    # plots$featureplots7 <- FeaturePlot(e.seu, features = c('MALAT1', 'CCNB2', 'RPS4X', 'RPS14', 'RPLP1', 'RPS2')) #
    
    # visualizing pseudotime in seurat
    
    # e.seu$pseudotime <- pseudotime(cds) #
    # Idents(e.seu) <- e.seu$seurat_clusters #
    # plots$featureplots8 <- FeaturePlot(e.seu, features = "pseudotime", label = T) #
    
    ####_________________________________________________________________####
    #SAVE (NEW TAB)
    
    # 10. saving processed data ------
    df$merged_seurat_filtered_final <- merged_seurat_filtered
    
  })
  
  #### ______________________________________________________________________ ####
  # Plots for Auto Annotation
  output$combineplot11 <- renderPlot({ plots$P6 + plots$a2 }) 
  
  output$combineplot12 <- renderPlot({ plots$p1 + plots$a3 }) 
  
  output$combineplot13 <- renderPlot({ plots$a2 })
  
  output$scoreheatmap1 <- renderPlot(( plots$scoreheatmap1 ))
  
  output$deltadistribution <- renderPlot(( plots$deltadistribution ))
  
  output$pheatmap1 <- renderPlot(( plots$pheatmap1 ))
  
  output$dimplot4 <- renderPlot({ plots$dimplot4 })
  
  output$scoreheatmap2 <- renderPlot(( plots$scoreheatmap2 ))
  
  output$dimplot5 <- renderPlot({ plots$dimplot5 })
  
  output$dimplot6 <- renderPlot({ plots$dimplot6 })
  
  output$vlnplot7 <- renderPlot({ 
    req(input$vlnplot7_selected)
    req(df$vlnplot7)
    
    VlnPlot(df$vlnplot7, features = input$vlnplot7_selected, group.by = 'singleR.labels')
  })
  
  output$vlnplot8 <- renderPlot({
    req(input$vlnplot8_selected)
    req(df$vlnplot8)
    
    if(df$mode == "integrated") { VlnPlot(df$vlnplot8, features = input$vlnplot8_selected, group.by = 'seurat_clusters_manual') }
  })
  
  output$vlnplot9 <- renderPlot({ 
    req(input$vlnplot9_selected)
    req(df$vlnplot9)
    
    if(df$mode == "integrated") { VlnPlot(df$vlnplot9, features = input$vlnplot9_selected, group.by = 'Location') }
  })
  
  output$vlnplot10 <- renderPlot({ 
    req(input$vlnplot10_selected)
    req(df$vlnplot10)
    
    VlnPlot(df$vlnplot10, features = input$vlnplot10_selected, slot = "counts", log = TRUE)
  })
  
  output$featureplot5 <- renderPlot({ 
    req(input$featureplot5_selected)
    req(df$featureplot5)
    
    if(length(input$featureplot5_selected) == 2) { FeaturePlot(df$featureplot5, features = input$featureplot5_selected, blend = TRUE) }
  })
  
  #output$featureplot6 <- renderPlot({ plots$featureplot6 })
  
  output$doheatmap1 <- renderPlot({ plots$doheatmap1 })
  
  output$ridgeplot2 <- renderPlot({ 
    req(input$ridgeplot2_selected)
    req(df$ridgeplot2)
    
    RidgePlot(df$ridgeplot2, features = input$ridgeplot2_selected, ncol=2)
  })
  
  output$doheatmap2 <- renderPlot({ plots$doheatmap2 })
  
  observeEvent(input$add_auto, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Combined Plot (Auto Annotation)", align = "center"), 
                 plotOutput("combineplot11Report", height = 500), br(), br(), 
                 h3("Combined Plot (Auto Annotation)", align = "center"), 
                 plotOutput("combineplot12Report", height = 500), br(), br(), 
                 h3("Combined Plot (Auto Annotation)", align = "center"), 
                 plotOutput("combineplot13Report", height = 500), br(), br(), 
                 h3("Score Heat Map (Auto Annotation)", align = "center"), 
                 plotOutput("scoreheatmap1Report", height = 500), br(), br(), 
                 h3("Delta Distribution (Auto Annotation)", align = "center"), 
                 plotOutput("deltadistributionReport", height = 500), br(), br(), 
                 h3("Heat Map (Auto Annotation)", align = "center"), 
                 plotOutput("pheatmap1Report", height = 500), br(), br(), 
                 h3("Dimplot (Auto Annotation)", align = "center"), 
                 plotOutput("dimplot4Report", height = 500), br(), br(), 
                 h3("Score Heat Map (Auto Annotation)", align = "center"), 
                 plotOutput("scoreheatmap2Report", height = 500), br(), br(), 
                 h3("Dimplot (Auto Annotation)", align = "center"), 
                 plotOutput("dimplot5Report", height = 500), br(), br(), 
                 h3("Dimplot (Auto Annotation)", align = "center"), 
                 plotOutput("dimplot6Report", height = 500), br(), br(), 
                 h3("Violin Plot (Auto Annotation)", align = "center"), 
                 plotOutput("vlnplot7Report", height = 500), br(), br(), 
                 h3("Violin Plot (Auto Annotation)", align = "center"), 
                 plotOutput("vlnplot8Report", height = 500), br(), br(), 
                 h3("Violin Plot (Auto Annotation)", align = "center"), 
                 plotOutput("vlnplot9Report", height = 500), br(), br(), 
                 h3("Violin Plot (Auto Annotation)", align = "center"), 
                 plotOutput("vlnplot10Report", height = 500), br(), br(), 
                 h3("Feature Plot (Auto Annotation)", align = "center"), 
                 plotOutput("featureplot5Report", height = 500), br(), br(), 
                 h3("Heat Map (Auto Annotation)", align = "center"), 
                 plotOutput("doheatmap1Report", height = 500), br(), br(), 
                 h3("Ridge Plot (Auto Annotation)", align = "center"), 
                 plotOutput("ridgeplot2Report", height = 500), br(), br(), 
                 h3("Heat Map (Auto Annotation)", align = "center"), 
                 plotOutput("doheatmap2Report", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$combineplot11Report <- renderPlot({ plots$P6 + plots$a2 }) 
    
    output$combineplot12Report <- renderPlot({ plots$p1 + plots$a3 }) 
    
    output$combineplot13Report <- renderPlot({ plots$a2 })
    
    output$scoreheatmap1Report <- renderPlot(( plots$scoreheatmap1 ))
    
    output$deltadistributionReport <- renderPlot(( plots$deltadistribution ))
    
    output$pheatmap1Report <- renderPlot(( plots$pheatmap1 ))
    
    output$dimplot4Report <- renderPlot({ plots$dimplot4 })
    
    output$scoreheatmap2Report <- renderPlot(( plots$scoreheatmap2 ))
    
    output$dimplot5Report <- renderPlot({ plots$dimplot5 })
    
    output$dimplot6Report <- renderPlot({ plots$dimplot6 })
    
    output$vlnplot7Report <- renderPlot({ 
      req(input$vlnplot7_selected)
      req(df$vlnplot7)
      
      VlnPlot(df$vlnplot7, features = input$vlnplot7_selected, group.by = 'singleR.labels')
    })
    
    output$vlnplot8Report <- renderPlot({
      req(input$vlnplot8_selected)
      req(df$vlnplot8)
      
      if(df$mode == "integrated") { VlnPlot(df$vlnplot8, features = input$vlnplot8_selected, group.by = 'seurat_clusters_manual') }
    })
    
    output$vlnplot9Report <- renderPlot({ 
      req(input$vlnplot9_selected)
      req(df$vlnplot9)
      
      if(df$mode == "integrated") { VlnPlot(df$vlnplot9, features = input$vlnplot9_selected, group.by = 'Location') }
    })
    
    output$vlnplot10Report <- renderPlot({ 
      req(input$vlnplot10_selected)
      req(df$vlnplot10)
      
      VlnPlot(df$vlnplot10, features = input$vlnplot10_selected, slot = "counts", log = TRUE)
    })
    
    output$featureplot5Report <- renderPlot({ 
      req(input$featureplot5_selected)
      req(df$featureplot5)
      
      if(length(input$featureplot5_selected) == 2) { FeaturePlot(df$featureplot5, features = input$featureplot5_selected, blend = TRUE) }
    })
    
    #output$featureplot6 <- renderPlot({ plots$featureplot6 })
    
    output$doheatmap1Report <- renderPlot({ plots$doheatmap1 })
    
    output$ridgeplot2Report <- renderPlot({ 
      req(input$ridgeplot2_selected)
      req(df$ridgeplot2)
      
      RidgePlot(df$ridgeplot2, features = input$ridgeplot2_selected, ncol=2)
    })
    
    output$doheatmap2Report <- renderPlot({ plots$doheatmap2 })
  })
  
  ####______________________________________________________________________####
  # Plots for Differential Expression
  
  output$v1 <- renderPlot({ 
    req(input$de_selected)
    req(df$v1)
    
    if(df$mode == "integrated") { VlnPlot(df$v1, features = input$de_selected, idents = df$idents, group.by = "Location") }
  }) 
  
  observeEvent(input$idents1_selected, {
    req(input$idents1_selected)  # Ensure the input is available
    req(df$v1)
    
    if(length(input$idents1_selected) == 2) {
      # Replace df$markers with the value of input$idents1_selected
      T1.de.markers <- FindMarkers(df$v1, ident.1 = input$idents1_selected[1], ident.2 = input$idents1_selected[2])
      
      T1.de.markers <- T1.de.markers[order(-T1.de.markers$avg_log2FC), ]
      T1.de.markers.asc <- T1.de.markers[order(+T1.de.markers$avg_log2FC), ]
      
      df$T1.de.markers <- T1.de.markers
      df$T1.de.markers.asc <- T1.de.markers.asc
    }
  })
  
  observeEvent(input$add_diffexp, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Violin Plot (Diffenrential Expression)", align = "center"), 
                 plotOutput("v1Report", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$v1Report <- renderPlot({ 
      req(input$de_selected)
      req(df$v1)
      
      if(df$mode == "integrated") { VlnPlot(df$v1, features = input$de_selected, idents = df$idents, group.by = "Location") }
    }) 
  })
  
  ####______________________________________________________________________####
  # Plots for Pseudobulk
  output$vb1 <- renderPlot({ 
    req(input$pseudobulk_selected)
    req(df$vb1)
    
    if(df$mode == "integrated") { VlnPlot(df$vb1, features = input$pseudobulk_selected, idents = df$idents, group.by = "Location") }
  })
  
  observeEvent(input$idents2_selected, {
    req(input$idents2_selected)  # Ensure the input is available
    
    if(length(input$idents2_selected) == 2) {
      bulk.T1.de <- FindMarkers(object = df$pseudo_merged_seurat, 
                                ident.1 = input$idents2_selected[1], 
                                ident.2 = input$idents2_selected[2],
                                test.use = "DESeq2")
      
      bulk.T1.de <- bulk.T1.de[order(-bulk.T1.de$avg_log2FC), ]
      bulk.T1.de.asc <- bulk.T1.de[order(+bulk.T1.de$avg_log2FC), ]
      
      df$bulk.T1.de <- bulk.T1.de
      df$bulk.T1.de.asc <- bulk.T1.de.asc
    }
  })
  
  observeEvent(input$add_pseudobulk, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Violin Plot (Pseudobulk)", align = "center"), 
                 plotOutput("vb1Report", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$vb1Report <- renderPlot({ 
      req(input$pseudobulk_selected)
      req(df$vb1)
      
      if(df$mode == "integrated") { VlnPlot(df$vb1, features = input$pseudobulk_selected, idents = df$idents, group.by = "Location") }
    })
  })
  
  ####______________________________________________________________________####
  # Plots for Compare between SCDE and PSCDE
  output$vlnplot11 <- renderPlot({ plots$vlnplot11 })
  
  output$vlnplot12 <- renderPlot({ plots$vlnplot12 })
  
  output$vlnplot13 <- renderPlot({ plots$vlnplot13 })
  
  output$vlnplot14 <- renderPlot({ 
    req(input$vlnplot14_selected)
    req(df$vlnplot14)
    
    if(df$mode == "integrated") { VlnPlot(df$vlnplot14, features = input$vlnplot14_selected, idents = df$idents, group.by = "Location") }
  })
  
  output$vlnplot15 <- renderPlot({ 
    req(input$vlnplot15_selected)
    req(df$vlnplot15)
    
    if(df$mode == "integrated") { VlnPlot(df$vlnplot15, features = input$vlnplot15_selected, idents = df$idents, group.by = "donor_id.stim", ncol = 1) }
  })
  
  observeEvent(input$add_comparison, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Violin Plot (Compare between SCDE and PSCDE)", align = "center"), 
                 plotOutput("vlnplot11Report", height = 500), br(), br(), 
                 h3("Violin Plot (Compare between SCDE and PSCDE)", align = "center"), 
                 plotOutput("vlnplot12Report", height = 500), br(), br(), 
                 h3("Violin Plot (Compare between SCDE and PSCDE)", align = "center"), 
                 plotOutput("vlnplot13Report", height = 500), br(), br(), 
                 h3("Violin Plot (Compare between SCDE and PSCDE)", align = "center"), 
                 plotOutput("vlnplot14Report", height = 500), br(), br(), 
                 h3("Violin Plot (Compare between SCDE and PSCDE)", align = "center"), 
                 plotOutput("vlnplot15Report", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$vlnplot11Report <- renderPlot({ plots$vlnplot11 })
    
    output$vlnplot12Report <- renderPlot({ plots$vlnplot12 })
    
    output$vlnplot13Report <- renderPlot({ plots$vlnplot13 })
    
    output$vlnplot14Report <- renderPlot({ 
      req(input$vlnplot14_selected)
      req(df$vlnplot14)
      
      if(df$mode == "integrated") { VlnPlot(df$vlnplot14, features = input$vlnplot14_selected, idents = df$idents, group.by = "Location") }
    })
    
    output$vlnplot15Report <- renderPlot({ 
      req(input$vlnplot15_selected)
      req(df$vlnplot15)
      
      if(df$mode == "integrated") { VlnPlot(df$vlnplot15, features = input$vlnplot15_selected, idents = df$idents, group.by = "donor_id.stim", ncol = 1) }
    })
  })
  
  ####______________________________________________________________________####
  # Plots for Monocle3 workflow
  
  output$a1_monocle <- renderPlot({ plots$a1 }) 
  
  output$a2_monocle <- renderPlot({ plots$a2 }) 
  
  output$clusters <- renderPlot({ plots$cluster.before.trajectory})  # or cluster.names
  
  output$plot_cell <- renderPlot({ plots$plot_cell })
  
  output$plot_cell2 <- renderPlot({ plots$plot_cell2 })
  
  output$pseudotime <- renderPlot({ plots$pseudotime })
  
  output$featureplot7 <- renderPlot({ plots$featureplot7 })
  
  output$featureplot8 <- renderPlot({ plots$featureplot8 })
  
  observeEvent(input$add_monocle3, {
    hide("error_text_report")
    insertUI(
      selector = '#placeholder',
      ui = tagList(
        fluidRow(
          column(width = 2),
          column(width = 8,
                 h3("Dimplot (Monocle3 Workflow)", align = "center"), 
                 plotOutput("a1_monocleReport", height = 500), br(), br(),
                 h3("Dimplot (Monocle3 Workflow)", align = "center"), 
                 plotOutput("a2_monocleReport", height = 500), br(), br(),
                 h3("Cluster before Trajectory (Monocle3 Workflow)", align = "center"), 
                 plotOutput("clustersReport", height = 500), br(), br(),
                 h3("Plot Cells (Monocle3 Workflow)", align = "center"), 
                 plotOutput("plot_cellReport", height = 500), br(), br(),
                 h3("Plot Cells (Monocle3 Workflow)", align = "center"), 
                 plotOutput("plot_cell2Report", height = 500), br(), br(),
                 h3("Pseudotime (Monocle3 Workflow)", align = "center"), 
                 plotOutput("pseudotimeReport", height = 500), br(), br(),
                 h3("Feature Plot (Monocle3 Workflow)", align = "center"), 
                 plotOutput("featureplot7Report", height = 500), br(), br(),
                 h3("Feature Plot (Monocle3 Workflow)", align = "center"), 
                 plotOutput("featureplot8Report", height = 500), br(), br()
          )
        )
      ),
    )
    
    output$a1_monocleReport <- renderPlot({ plots$a1 }) 
    
    output$a2_monocleReport <- renderPlot({ plots$a2 }) 
    
    output$clustersReport <- renderPlot({ plots$cluster.before.trajectory})  # or cluster.names
    
    output$plot_cellReport <- renderPlot({ plots$plot_cell })
    
    output$plot_cell2Report <- renderPlot({ plots$plot_cell2 })
    
    output$pseudotimeReport <- renderPlot({ plots$pseudotime })
    
    output$featureplot7Report <- renderPlot({ plots$featureplot7 })
    
    output$featureplot8Report <- renderPlot({ plots$featureplot8 })
  })
  
  ####______________________________________________________________________####
  # Update selectizeInput as you enter the tab
  ## Auto Annotation
  observeEvent(input$navbar, {
    req(input$navbar)
    if(input$navbar == "autoannotation_tab" & count$autoannotation_tab == 0 & df$complete) {
      updateSelectizeInput(session, "vlnplot7_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "vlnplot8_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "vlnplot9_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "vlnplot10_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "featureplot5_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1:2], 
                           server = TRUE)
      
      updateSelectizeInput(session, "ridgeplot2_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      count$autoannotation_tab <- 1
    }
  })
  
  ## Differential expression 
  observeEvent(input$navbar, {
    req(input$navbar)
    if(input$navbar == "diff_expression_tab" & count$diff_expression_tab == 0 & df$complete) {
      updateSelectizeInput(session, "de_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "idents1_selected", 
                           choices = df$idents, 
                           selected = df$idents[1:2], 
                           server = TRUE)
      
      count$diff_expression_tab <- 1
    }
  })
  
  ## Pseudobulk
  observeEvent(input$navbar, {
    req(input$navbar)
    if(input$navbar == "pseudobulk_tab" & count$pseudobulk_tab == 0 & df$complete) {
      updateSelectizeInput(session, "pseudobulk_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "idents2_selected", 
                           choices = df$idents, 
                           selected = df$idents[1:2], 
                           server = TRUE)
      
      count$pseudobulk_tab <- 1
    }
  })
  
  
  ## Compare Between SCDE and PSCDE
  observeEvent(input$navbar, {
    req(input$navbar)
    if(input$navbar == "compare_tab" & count$compare_tab == 0 & df$complete) {
      updateSelectizeInput(session, "vlnplot14_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      updateSelectizeInput(session, "vlnplot15_selected", 
                           choices = df$all.genes, 
                           selected = df$all.genes[1], 
                           server = TRUE)
      
      count$compare_tab <- 1
    }
  })
  
  observeEvent(input$navbar, {
    req(input$navbar)
    if(input$navbar == "monocle3_tab" & count$monocle3_tab == 0 & df$complete) {
      updateSelectizeInput(session, "monocle_selected", 
                           choices = df$e.seu_idents, 
                           selected = df$e.seu_idents, 
                           server = TRUE)
      
      count$monocle3_tab <- 1
    }
  })
  
  ####______________________________________________________________________####
  # Download rds
  output$downloadReport <- downloadHandler(
    filename = function() { "merged_seurat_filtered_final.rds" },
    content = function(file) {
      saveRDS(df$merged_seurat_filtered_final, file)
    }
  )
  # session$onSessionEnded(stopApp)
}

app <- shinyApp(ui = ui, server = server)
app
