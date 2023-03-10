#' ---
#' title: Course5633 - Bioinformatics 
#' author: Qiao Yang
#' output:
#'   html_document:
#'     keep_tex: true
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     highlight: zenburn
#'     fig_width: 10
#'     fig_height: 9
#'     fig_dpi: 300
#' ---
#'
#'
#' # Course5633 Assignment 1
#+ message=FALSE

#### loading library ####
#install.packages("gtsummary")
library(dplyr)
library(magrittr)
library(tidyverse)
library(ggpubr)

#### Task 1 - Literature ####

#' ## Task 1 - Literature
#' ### Read the research article   
#' 
#' 
#'  The paper is "A single-cell and spatially resolved atlas of human breast cancers", published in 2021 on *Nature Genetics* journal. 
#'  This research investigated the tumor microenvironment of 26 primary breast cancer covering 11 ER+, 5 HER2+ and 10 TNBC subtypes, using single-cell RNA-seq and spatial transcriptomics analysis, and verifying some protein expression and interaction by CITE-seq. 
#' 
#' ### Answer the following questions   
#' 
#' **a. What is the medically relevant insight from the article?**  
#' 
#'  It discovers medical insights as below:  
#' 
#' 
#' - Developing single-cell RNA-seq subtype (scSubtype) method to call the intrinsic subtype on single-cell level and identify the heterogeneity of breast cancer in intrinsic subtyping on cellular level. For example, in some luminal and HER2E tumors, scSubtype predicted some basal-like cells, which is validated by IHC as well. This intrinsic subtype heterogeneity may predict the innate resistance to therapy and early relapse following therapy.   
#' 
#' - They constructed 7 robust gene modules (GMs) using driving genes of neoplastic cells, to represent the heterogeneity on cellular and spatial level. They further found that cancer cells manifest heterogeneous phenotypes in most tumors. With help of spatial transcriptomics data, GMs were spatially mapped to tumor regions, and they found that a mesenchymal-like state and proliferation are incompatible at cellular resolution in breast cancer.    
#' 
#' - Using scRNA-seq to depict immune milieu of breast cancer, and two novel macrophage population were found and further discovered their juxtaposition to PD-1+ lymphocytes and clinical relevance of their gene signatures in a larger cohort.  
#' 
#' - CAFs, PVL and endothelial cells were clustered and identified their pseudo-time trajectory state in the stromal compartment. They provide evidence that differentiation can drive transition between CAF subsets. They also suggest that signals from environment may control the differentiation or migration of mesenchymal cells. This would pave the way for aiming to control stromal-immune relevant bio-chemistry therapy.   
#' 
#' - Stromal-immune niches were spatially deconvoluted and they found distinct cancer phenotypes occur in mutually exclusive regions of breast cancers via spatially assigning GMs to breast tumors.   
#' 
#' - They also estimate the cellular proportions of bulk RNA-seq data, using gene signatures identified by scRNA-seq data instead. They found not only immune cells, major lineages (epithelial, immune and stromal) may also contribute to the primary breast cancer stratifying, which were partially associated with intrinsic subtype and genomic classifiers. 
#' 
#' **b. Which genomics technology/ technologies were used?**
#' 
#' - This study includes four bioinformatic technologies, covering single-cell RNA-seq (scRNA-seq), spatial transcriptomics, CITE-seq and bulk RNA-seq. It includes scRNA-seq data of 26 fresh surgical primary breast cancer, using Chromium Single-cell c2 3' and 5' Chemistry library on Nextseq 5000 platform and bulk RNA-seq data of matching 24 FFPE primary breast cancer, using High Pure Paraffin Kit on Hiseq 2500 platform. Four matching samples were stained with 10X Chromium 3' capture compatible TotalSeq-Antibodies (processing CITE-seq simultaneously) . Six matching samples were also performed spatial-resolved transcriptomics, using Visium Spatial Gene Expression Kit (10X Genomics).  
#' 
#' 
#' - This study also includes one immunostaining technology. Tumor samples were fixed and then processed for paraffin embedding to perform standard histological analysis. Immunohistochemistry (IHC) were performed on tumor blocks to estimate CK5 and ER. Tumor slides were also imaged using the Aperio CS2 Digital Pathology Slid Scanner and processed using Qupath (v0.2.0).  
#' 
#' ### Further related research questions 
#' 
#' **a. List and explain at least three questions/ hypotheses you can think of that extend the analysis presented in the paper**
#' 
#' 
#' - Here, it only studied primary breast cancers, but patients with breast cancer suffered from high probability of recurrence and metastasis, causing low survival rate (Riggio et al., 2020). It may extend the tumor tissues to recurrence tissue or metastasis tissue. In this way, it may investigate cellular proportions and cell interactions in recurrence and metastasis tumor, using gene signatures of scRNA-seq to deconvolute the cellular and intrinsic subtype heterogeneity during breast cancer progression.   
#' 
#' 
#' - Breast cancer have different subtypes, Luminal A/B, Her2-enriched, basal-like and normal-like intrinsic subtypes, mentioned in this paper. There are many methods to identify intrinsic subtypes of breast cancer, but they only mentioned one method in this paper, parker et al., (2009) and the test result only showed 66% agreement. Perhaps they can try to use other PAM50 calling algorithms or cohort with Prosigna Breast Cancer Prognostic Gene signature Assay, which is considered as the most accurate PAM50 test, to improve the accuracy of scSubtype.   
#' 
#' 
#' - Some researchers discovered that the presence of ERBB2 on the plasma membrane of TILs at time of diagnosis have associations with clinical efficacy (Suzuki et al., (2015)). In this study, it quantified 157 protein expression and transcriptomics data using CITE-seq and maintained single-cell information (B, T, myeloid and mesenchymal cells). Perhaps they can further confirm this phenomenon and find more robust controllable protein involved immune activation to give new insight to cancer therapy.   
#' 
#' 
#' **b. [Optional] Devise a computational analysis strategy for (some of) the listed questions under 3a**. 
#'  
#'  
#'  

# Improving PAM50 classification methods 

## function, Find the subtype identified by moste PAM50 classfication methos and verified by prosigna test.  
## Gene_matrix : gene expression matrix, gene symbol as rownames  
## methods: PAM50 classification method and their function calling method.  
## prosigna: it contains prosigna classification results, sample as rows.  

PAM50_improved = function( Gene_matrix, methods, prosigna ){
  
  res_pam50 = as.data.frame(matrix(ncol=0, nrow=ncol(Gene_matrix)))
  rownames(res_pam50) = colnames(Gene_matrix)
  
  ## call method to classify PAM0 results
  j =1
  for (i in methods) {
    res_pam50_temp = methods_function[i](Gene_matrix)
    res_pam50 = res_pam50_temp$results
    colnames(res_pam50)[j] = i
    j=j+1
  }

  ## integrate prosigna results
  res_pam50 = cbind(res_pam50, prosigna)
  
  ## Choose one PAM50 subtype identified by most methods and verified by prosigna
 Find_agreement_pm50subtype =  function(x) { 
    res= as.data.frame(table(as.factor(x[-length(x)]))) 
    subtype = which.max(res); 
    if (subtype == x[length(x)] )  
      subtype_final = subtype 
    else 
      subtype_final = ""  
    return(subtype)}
  
  res_agreement = apply(res_pam50, 2, Find_agreement_pm50subtype(x) )
  
  ##integrate the final results
  rownames(res_agreement) = rownames(res_pam50)
  res_pam50 = cbind(res_pam50, res_agreement)
  
  return(res_pam50)
}

#'  
#'  
#' 

#### Task 2 - Git repositories and R Markdown ####

#' ## Task 2 - Git repositories and R Markdown
#' 
#' ### Start a new Github repository
#' 
#' Start a new project in a Gitlab, Github or Figshare repository. Check with your doctoral supervisor if you can start a project in the repository of your lab or if you have to start your own repository.  
#' 
#' - I created a new repository "yqkiuo/Course5633Rrecord" to save documents of this course. 
#' 
#' ### Assignment 1 requirements  
#' 
#' All documentation of the Assignment 1 has to be provided in your Git/ Figshare repository as a (formatted) text document.  
#' 
#' - Finished. 
#' 
#' ### Assignment 2 requirements
#' 
#' All documentation of your hands-on work during Week 2 also has to be provided in your repository as R Markdown document(s).  
#' 
#' - I will do that. 




#### Task 3 - Introduction to R and online R course ####

#' 
#' ## Task 3 - Introduction to R and online R course
#' 
#' ### R and R studio installation
#' 
#' This part focus on R and R studio installation as well as R language learning and R coding tips. I basically finished them. 
#' 
#' - Install the most recent version of the R software on your computer by following the instructions provided at the R software website.  
#' 
#' - Install the most recent version of RStudio Desktop (Open Source version) on your computer.  
#' 
#' - Bioconductor is an add-on package for R providing tools for the analysis and comprehension of high-throughput genomic data.  
#' 
#' - Install the most recent version of the Bioconductor package on your computer.   
#' 
#' - Tidyverse is a toolbox for streamlining data import, modeling, transformation, curation, and visualization. Tidyverse tools enable you to make your scripts more reader friendly and overall more neat and efficient. https://www.tidyverse.org/packages/. 
#' 
#' - R online course   
#' 
#' - R-notebook for assignment 1  
#' 
#' - Reference for R language. 
#' 
#' - More information on R Markdown and the Markdown cheatsheet
#' 
#' 

#### Task 4 - R basic operations ####

#' ## Task 4 - R basic operations
#'         
#'        
#' **1. What is the square root of 10?**
message( "The answer : sqrt(10) = ",sqrt(10) )
message( "The answer (with two decimal places) : round(sqrt(10),2) = ", round(sqrt(10),2) )
#'    
#'      
#'       
#'         
  
#' **2. What is the logarithm of 32 to the base 2?**
message( "The answer : ", log2(32))
#'    
#'      
#'         
#'     
#' **3. What is the sum of the numbers from 1 to 1000? **
#'  
#' 
sum =0
for ( i in seq(1:1000) ){
  sum = sum+ i 
}
message( "The answer : ", sum)

#'     
#' **4. What is the sum of all even numbers from 2 to 1000?**

sum =0
for ( i in seq(1:1000) ){
  if ( (i%%2) == 0) {
  sum = sum + i 
}}
message( "The answer : ", sum)

#' **5. How many pairwise comparisons are there for 100 genes?**

#' This is special mathematical combination question.  

#' The combination formula is C(n,m).  

#' The function in r is choose().  


#This is special mathematical combination question
#The combination formula is C(n,m)
#The function in r is choose()
message( " The answer : ", choose(100,2))



#' **6. And how many ways to arrange 100 genes in triples?**


#' This is special mathematical arrangement question.  

#' The  formula is A(n,m) = n!/(n-m)!  

#' The useful function in r is prod().  

 
#This is special mathematical arrangement question
#The  formula is A(n,m) = n!/(n-m)!
#The useful function in r is prod()
message( " The answer : ", prod(1:100) / prod(1:(100-3) ) )

#### Task 5 - Using R example datasets ####

#' ## Task 5 - Using R example datasets
#' 
#' 
#' **1. Use the R internal CO2 dataset ("data(CO2)")** 

#loading inhouse data
data("CO2")

#' 
#' 
#' **2. Describe briefly the content of the CO2 dataset using the help function** 
#' 

#' Check basic data frame structure first
#' 

#Check basic data frame structure
help("CO2")
dim(CO2)
head(CO2)

#' If it contains missing data?
#If it contains missing data?
anyNA(CO2)

#' 

#' Hence this dataset contains intact CO2 concentration levels of the cold tolerance of the grass species *Echinochloa crus-galli*, covering six plants from Quebec and six plants from Mississippi. Half the plants of each type were chilled overnight and half were not as control.  
#' 
#' 
#' 
#' 
#' 
#' **3. What is the average and median CO2 uptake of the plants from Quebec and Mississippi?** 
#' 
#' 
#' 
#' In the second question, we checked the data frame.
#' The answers : 
CO2 %>%
  summarise(mean = mean(uptake), median = median(uptake))

#' 
#' **4. [Optional] In the "airway" example data from Bioconductor, how many genes are expressed in each sample? How many genes are not expressed in any sample?** 
#' 
#' We can check this dataset [here](http://bioconductor.org/packages/release/data/experiment/html/airway.html), 

#' First we install and load data
#+ message=FALSE
# BiocManager::install("airway")
library("airway")
data("airway")
#' We then check this data frame structure and content. 
help("airway")
airway 
#' 
#' This dataset includes read counts in genes for an RNA-Seq experiment on four human airway smooth muscle cell lines treated with dexamethasone.  
#' 
#' The answer of how many genes are expressed in each sample : 
colSums( assay(airway) != 0 )

#' How many genes are not expressed in any sample : 
nrow(assay(airway)[ rowSums( assay(airway) ) == 0,] )


#### Task 6 - R Functions ####

#' ## Task 6 - R Functions
#' 
#' **1. Write a function that calculates the ratio of the mean and the median of a given vector.**
#'  
#' x must be a vector
# x must be a vector
Ratio_MM = function(x) {
  if (!(is.null(x))) { #To check if the vector is null.
    ratio = mean(x, na.rm = TRUE)/median(x,na.rm = TRUE)
    return(ratio)
  }
}

#' This is an example of usage : 
#For example , x = c(1,2,3,10,5,10,7,8)
x = c(1,2,3,10,5,10,7,8)
ratio = Ratio_MM(x)
message( "The answer of the example : ", ratio)

#' **2. Write a function that ignores the lowest and the highest value from a given vector and calculate the mean.**
#' 
#' 
myMean = function(x) {
  if (!(is.null(x))) { #To check if the vector is null.
    new_x = x[c(-which.min(x), - which.max(x))]
    res = mean(new_x)
    return(res)
  }
}

#' This is an example of usage : 
x = c(1,2,3,10,5,10,7,8)
myMean = myMean(x)
message( "The answer of the example : ", myMean)

#' **3. Why, How and When to use pipe?**  
#'
#'
#' Pipe in R language is a consice and expressive method to integrate multiple operations together.   
#' The reason why we use pipe is it is useful 1) to prettify your code, 2) to make your code understood, 3) to avoid unnecessary local parameters and function definitions  
#' The pipe is offered by magrittr-package and triggered by %>%, which will take left-hand side values as input, and forward them into right-hand side expression.  
#' When the codes are sequential and the former one is the input of the latter one, but you want to minimize your code, then pipe would be helpful.   
#' 
#' 
#' **4. Write a short explanation of why apply-family could be useful in your work.**  
#'
#'
#' Apply-family includes apply(), lapply(), sapply(), mapply() and tapply(), efficiently avoiding loop functions, saving CPU consumption and reducing job running time. 
#' 


#### Task 7 - Basic visualization with R ####

#'
#' ## Task 7 - Basic visualization with R
#' ### 'magic_guys.csv' dataset  
#'
#' **1. Compare the distributions of the body heights of the two species graphically**  
#'  
#' Load dataset 'magic_guys.csv'
#' Check dataset
path= "/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/1 KI course/6 VT23 spring /Course5633/Course5633Rrecord/data/"
magic = read.csv(paste(path,"magic_guys.csv",sep = "") )
dim(magic)
head(magic)
table(magic$species)
#' **a. using the basic 'hist' function as well as 'ggplot' and 'geom_histogram' functions from the ggplot2 package. Optimize the plots for example by trying several different 'breaks'. Note that ggplot2-based functions give you many more options for changing the visualization parameters, try some of them.** 
#' 
#' 
#' With default breaks in hist()
p1 = magic %>%
  dplyr::filter(species == "jedi" ) %>%
  dplyr::select(length) %$%
  hist(length,plot = F)

p2 = magic %>%
  dplyr::filter(species == "sith" ) %>%
  dplyr::select(length) %$%
  hist(length, plot = F)

{
plot(p1, col = rgb(0,0,1,1/4), xlim=c(50,300), ylim = c(0,20), ) #"jedi"
plot(p2,col=rgb(1,0,0,1/4),xlim=c(50,300), ylim = c(0,20),add=T) #"sith"
}

#' With different breaks in hist()
p1 = 
  magic %>%
  dplyr::filter(species == "jedi" ) %>%
  dplyr::select(length) %$%
  hist(length,breaks = 15,plot = F)

p2 = 
  magic %>%
  dplyr::filter(species == "sith" ) %>%
  dplyr::select(length) %$%
  hist(length,breaks = 6,plot = F)

{
  plot(p1, col = rgb(0,0,1,1/4), xlim=c(50,300), ylim=c(0,20) ) #"jedi"
  plot(p2,col=rgb(1,0,0,1/4),ylim=c(0,20),add=T) #"sith"
  }

#' Using ggplot + geom_histogram()  
#' Computed variables: count
magic %>%
  dplyr::select(species,length) %>%
  ggplot(aes(length, fill= species)) +
  geom_histogram(alpha= 0.5, position = 'identity',bins = 10)+
  theme_classic()

#' Computed variables: density
magic %>%
  dplyr::select(species,length) %>%
  ggplot(aes(length, fill= species)) +
  geom_histogram(alpha= 0.5, aes(y = ..density..),position = 'identity',bins = 10)+
  theme_classic()

#'
#' **b. Do the same comparison as in a. but with boxplots. If you want to use the ggplot2-package, use the functions 'ggplot' and 'geom_boxplot'.** 
#' 
#' With boxplot()
boxplot( length ~ species, data = magic)

#' With ggplot + geom_boxplot

magic %>%
  ggplot(aes(x = species, y = length, fill = species)) +
  geom_boxplot() +
  theme_bw()

#' 
#' 
#' **c. Save the plots with the 'png', 'pdf', and 'svg' formats. In which situation would you use which file format?**  
#' 
#' PNG
# PNG example
path= "/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/1 KI course/6 VT23 spring /Course5633/Course5633Rrecord/data_plots/"
png(paste(path,"7_1_a_PNG.png",sep = "") )

# Code
boxplot( length ~ species, data = magic)

# Close device
dev.off()

#' PDF 
# PDF example
path= "/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/1 KI course/6 VT23 spring /Course5633/Course5633Rrecord/data_plots/"
pdf(paste(path,"7_1_a_PDF.pdf",sep = "") )

# Code
boxplot( length ~ species, data = magic)

# Close device
dev.off()


#' SVG 
# SVG example
path= "/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/1 KI course/6 VT23 spring /Course5633/Course5633Rrecord/data_plots/"
svg(paste(path,"7_1_a_SVG.svg",sep = "") )

# Code
boxplot( length ~ species, data = magic)

# Close device
dev.off()
#'
#' **In which situation would you use which file format?**  
#' 
#' 
#' PDF offers high quality with less file size, and can be viewed, printed or electronically transmitted with high figure quality. It is preferable for printing and sharing.   
#' PNG is an image format with pixels, easy to display but the quality depends on the resolution. It is better for figures with more color and details with high resolution.        
#' SVG supports animation and can be resized. But it suits for less complex elements, like less color.     
#'
#' 
#' 
#' 
#' ### 'microarray_data.csv' dataset  
#' 
#' 
#' 
#' **2. Load the gene expression data matrix from the 'microarray_data.tab' dataset**  
#' 
path= "/Users/qiaoyang/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Karolinska Ins/1 KI course/6 VT23 spring /Course5633/Course5633Rrecord/data/"
tab = read.table( paste(path, "microarray_data.tab",sep = ""), sep = "\t",
                  header = 1) 
#'  
#'  
#' **a. How big is the matrix in terms of rows and columns?**  
#' Chech this data
dim(tab)
#' it has 553 rows and 1000 columns.    
#' 
#' 
#' **b. Count the missing values per gene and visualize this result.**  
#' 

as.data.frame( colSums(is.na(tab)) )%>%
  dplyr::rename( `The number of Missing`= `colSums(is.na(tab))`)%>% 
  knitr::kable() %>%
  kableExtra::kable_styling("striped") %>%
  kableExtra::scroll_box(width = "50%",height = "200px")

#' **c. Find the genes for which there are more than X% (X=10%, 20%, 50%) missing values**  

#'   
#' Genes with missing data more than 10%
as.data.frame( colSums(is.na(tab)) )%>%
  dplyr::rename( `The number of Missing`= `colSums(is.na(tab))`)%>% 
  dplyr::mutate(`The percentage of Missing(%)` = `The number of Missing`/500*100 ) %>%
  dplyr::filter(`The percentage of Missing(%)` >10) %>%
  knitr::kable() %>%
  kableExtra::kable_styling("striped") %>%
  kableExtra::scroll_box(width = "50%",height = "200px")

#'   
#' Genes with missing data more than 20%
as.data.frame( colSums(is.na(tab)) )%>%
  dplyr::rename( `The number of Missing`= `colSums(is.na(tab))`)%>% 
  dplyr::mutate(`The percentage of Missing(%)` = `The number of Missing`/500*100 ) %>%
  dplyr::filter(`The percentage of Missing(%)` >20) %>%
  knitr::kable() %>%
  kableExtra::kable_styling("striped") %>%
  kableExtra::scroll_box(width = "50%",height = "200px")


#'   
#' Genes with missing data more than 50%
as.data.frame( colSums(is.na(tab)) )%>%
  dplyr::rename( `The number of Missing`= `colSums(is.na(tab))`)%>% 
  dplyr::mutate(`The percentage of Missing(%)` = `The number of Missing`/500*100 ) %>%
  dplyr::filter(`The percentage of Missing(%)` >50) %>%
  knitr::kable() %>%
  kableExtra::kable_styling("striped") %>%
  kableExtra::scroll_box(width = "50%",height = "200px")


#' **d. Replace the missing values by the average expression value for the particular gene.**
#' 
#' 
myReplace = function(x) { x[is.na(x)] = mean(x,na.rm = TRUE); return(x)}
myTab = apply(tab, 2, myReplace)

#' The number of missing values after replacing step  
#' Note: several genes are not expressed  
#'  
as.data.frame( colSums(is.na(myTab)) )%>%
  dplyr::rename( `The number of Missing`= `colSums(is.na(myTab))`)%>% 
  knitr::kable() %>%
  kableExtra::kable_styling("striped") %>%
  kableExtra::scroll_box(width = "50%",height = "200px")

#' 
#' **3. Visualize the data in the CO2 dataset in a way that gives you a deeper understanding of the data.**   
#'  

#' First, we could investigate the uptake difference between Quebec and Mississippi.  
#' 
#' 
#' 
#Boxplot by ggplot
CO2 %>%
  ggplot() +
  geom_boxplot(aes(x = Plant, y = uptake, fill = Treatment)) +
  theme_classic()
#' The level of CO2 uptake in Quebec is higher than that in Mississippi. And in Mississippi, chilled treatment highly reduced the uptakeof CO2, while slightly influenced in Quebec.  
#'  
#'  

CO2 %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(conc), y = uptake, fill = Treatment )) +
  facet_grid(Type~. ) +
  theme_classic()
#' The level of CO2 uptake is improving  with the increased level of ambient CO2 concentration in both Quebec and Mississippi. The tolerance ability could be improving with higher CO2 concentration, but would be adversely influenced by cold treatment, especially plants in Mississippi.    
#'  
#'  
#'  

#### Task 8 - Bioinformatics practice ####

#' ## Task 8 - Bioinformatics practice

#' ### Tidybiology
#' Install packages and chec datat  
#' 
#devtools::install_github("hirscheylab/tidybiology")
library(tidybiology)
data(package = "tidybiology")

#' Loading data  
#' 
data("chromosome")
dim(chromosome)

#' **a. Summary statistics for variables in chromosome dataset: variations, protein coding genes, and miRNAs**  
#' 
#' 
#' 
chromosome %>%
  dplyr::select(variations,protein_codinggenes,mi_rna) %>%
  summarise_all(  list(mean, median, max ) )

#' **b. How does the chromosome size distribute?**  
#' 

chromosome %>%
  ggplot() +
  geom_histogram( aes(length_mm), bins = 8) +
  theme_classic()
#' Most chromosomes have length in 50mm and there are less long chromosomes over 50mm.  
#'
#' **c. Does the number of protein coding genes or miRNAs correlate with the length of the chromosome?**  
#' 
#'  

chromosome %>%
  ggplot(aes(x = length_mm, y = `protein_codinggenes`))  +
  geom_point()  +
  geom_smooth(method = lm) +
  stat_cor(method = "pearson") +
  xlab("Chromosome length_mm") +
  ylab("The number of coding genes") + 
  theme_classic()

chromosome %>%
  ggplot(aes(x = length_mm, y = mi_rna))  +
  geom_point()  +
  geom_smooth(method = lm) +
  stat_cor(method = "pearson") +
  xlab("Chromosome length_mm") +
  ylab("The number of miRNAs") + 
  theme_classic()
#' The number of coding genes is significantly related to the length of chromosomes (R = 0.61 & p < 0.05), and the same as between the number of miRNAs and chromosome length (R =0.74, p < 0.05). Obviously, the number of miRNAs is more associated to chromosome length.  
#' 
#' **D. Perform similar summary and visualization for "proteins" data **
#' 
#' 
#' First, we load data and check it.
#' 
 
data("proteins")
dim(proteins)

#' Summary statistics for length and mass  
#' 
#' 
proteins %>%
  dplyr::select(length, mass) %>%
  summarise_all(list(mean, median, max))

#' Relationship between length and mass in protein  
#' 
#' 
proteins %>%
  dplyr::select(length, mass) %>%
  dplyr::filter(length < 10000) %>% ## removing two outliers with super high length
  ggplot(aes(x = length, y = mass)) +
  geom_point() +
  geom_smooth(method = lm) +
  stat_cor(method = "pearson") +
  xlab("Length of protien") +
  ylab("Mass of protein") + 
  labs(title = "The relationship of length and mass in proteins") +
  theme_bw()+
  theme(axis.text= element_text(size = 10, face = "bold"), 
        axis.title = element_text(size =15, face ="bold"),
        title = element_text(size =15, face ="bold"))

#' The mass of a protien is highly related to the length of protein from this figure (R = 1 & p < 0.01).   
