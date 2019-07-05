##############################
## PCAdapt outlier analysis
##
## Matt Brachmann (PhDMattyB)
##
## 2019-07-05
##
##############################

setwd('~/PhD/SNP Demographic modelling/Outliers_directory/')

library(tidyverse)
library(wesanderson)
library(patchwork)
library(janitor)
library(devtools)
library(skimr)
library(tvthemes)
library(rsed)
library(data.table)

theme_set(theme_bw())

# Other packages to load

library(vegan)
library(psych)
library(adegenet)
library(ggman)
library(qvalue)
library(pcadapt)
library(OutFLANK)
library(LEA)
library(diveRsity)

## MAP file ####
## load in the map file so we know where outliers are within the 
## genome
#MAP = read_tsv('icelandic_pops_plink.map')
MAP = read_tsv('Feb202019_Poly_Plink_input.map')
MAP$`#Chromosome`
chr_names = c('AC01', 'AC02', 'AC03', 'AC04p', 'AC04q.1:29', 'AC04q.2', 'AC05', 'AC06.1', 'AC06.2', 'AC07', 'AC08', 'AC09', 'AC10', 'AC11', 'AC12', 'AC13', 'AC14', 'AC15', 'AC16', 'AC17', 'AC18', 'AC19', 'AC20', 'AC21', 'AC22', 'AC23', 'AC24', 'AC25', 'AC26', 'AC27', 'AC28', 'AC30', 'AC31', 'AC32', 'AC33', 'AC34', 'AC35', 'AC36', 'AC37', 'contigs')

t(chr_names)
length(chr_names)
MAP = mutate(.data = MAP,
             CHROME = as.factor(case_when(
               `#Chromosome` == '1' ~ 'AC01',
               `#Chromosome` == '2' ~ 'AC02',
               `#Chromosome` == '3' ~ 'AC03',
               `#Chromosome` == '4' ~ 'AC04p',
               `#Chromosome` == '5' ~ 'AC04.1:29',
               `#Chromosome` == '6' ~ 'AC04q.2',
               `#Chromosome` == '7' ~ 'AC05',
               `#Chromosome` == '8' ~ 'AC06.1',
               `#Chromosome` == '9' ~ 'AC06.2',
               `#Chromosome` == '10' ~ 'AC07',
               `#Chromosome` == '11' ~ 'AC08',
               `#Chromosome` == '12' ~ 'AC09',
               `#Chromosome` == '13' ~ 'AC10',
               `#Chromosome` == '14' ~ 'AC11',
               `#Chromosome` == '15' ~ 'AC12',
               `#Chromosome` == '16' ~ 'AC13',
               `#Chromosome` == '17' ~ 'AC14',
               `#Chromosome` == '18' ~ 'AC15',
               `#Chromosome` == '19' ~ 'AC16',
               `#Chromosome` == '20' ~ 'AC17',
               `#Chromosome` == '21' ~ 'AC18',
               `#Chromosome` == '22' ~ 'AC19',
               `#Chromosome` == '23' ~ 'AC20',
               `#Chromosome` == '24' ~ 'AC21',
               `#Chromosome` == '25' ~ 'AC22',
               `#Chromosome` == '26' ~ 'AC23',
               `#Chromosome` == '27' ~ 'AC24',
               `#Chromosome` == '28' ~ 'AC25',
               `#Chromosome` == '29' ~ 'AC26',
               `#Chromosome` == '30' ~ 'AC27',
               `#Chromosome` == '31' ~ 'AC28',
               `#Chromosome` == '32' ~ 'AC30',
               `#Chromosome` == '33' ~ 'AC31',
               `#Chromosome` == '34' ~ 'AC32',
               `#Chromosome` == '35' ~ 'AC33',
               `#Chromosome` == '36' ~ 'AC34',
               `#Chromosome` == '37' ~ 'AC35',
               `#Chromosome` == '38' ~ 'AC36',
               `#Chromosome` == '39' ~ 'AC37',
               `#Chromosome` > '39' ~ 'Contigs')))

is.na(MAP$CHROME)
MAP$CHROME[is.na(MAP$CHROME)] = 'Contigs'
MAP$CHROME
#View(MAP3$CHROME3)

## Need to load the ped file  
## Once loaded we can split the genotype file up into each populations
PED = read_tsv('Feb202019_Poly_Plink_input.ped')

## POPN PED FILES #####
## This allows us to create ped files for each population!!
## Just need to put the Family id numbers for the popualtion in

# Pop_PED = PED %>% 
#   filter(`#FamilyID` %in% c('8', '4', '5'))

# write_tsv(Pop_PED, '', col_names = T)

## POPS ####
## We will use this ped file to delinate populations

POPS = PED %>% select(`#FamilyID`) %>% 
filter(`#FamilyID` %in% c())

## This allows us to name the populations 
POPS = mutate(.data = POPS,
              POP_name = as.factor(case_when(
                `#FamilyID` == "1" ~ "T.LGB",
                `#FamilyID` == '2' ~ 'V.BR',
                `#FamilyID` == '3' ~ 'V.SIL',
                `#FamilyID` == '8' ~ 'S.LGB',
                `#FamilyID` == '4'~ 'S.PL',
                `#FamilyID` == '5' ~ 'S.PI',
                `#FamilyID` == '6' ~ 'T.PL',
                `#FamilyID` == '7' ~ 'T.SB',
                `#FamilyID` == '9' ~ 'G.SB',
                `#FamilyID` == '10' ~ 'G.PI'))) 

## PCAdapt ####
## Load in the bed file for the popualtion you want to analyse
## Make the bed file in plink using the --make-bed flag
CHARR_POP = read.pcadapt('.bed', 
                         type = 'bed')

## Exploratory PCA to find the number of genetic groups
PCA = pcadapt(CHARR_POP, K = 15, method="mahalanobis", 
              min.maf = 0.01)
summary(PCA)
## Screeplot of the principal components to determine the 
## number of genetic groups
plot(PCA, option = "screeplot", K = 15)

## Re-load the data into R with the optimal k-value based on the screeplot
FISHY = pcadapt(CHARR_POP, K = 9, 
                method="mahalanobis", 
                min.maf = 0.01)
summary(FISHY)
FISHY$singular.values

## Exploratory graphs
# plot(FISHY , option = "manhattan")
# plot(FISHY, option = "qqplot")
# hist(FISHY$pvalues, xlab = "p-values", 
#      main = NULL, breaks = 50, col = "Grey")
# plot(FISHY, option = "stat.distribution")

## use this to get the PCA data that we will graph later on
# ggplot_FISHY = FISHY$scores
# ggplot_FISHY = as_tibble(ggplot_FISHY) %>%
#   rename(PC1 = 1,
#          PC2 = 2, 
#          PC3 = 3, 
#          PC4 = 4, 
#          PC5 = 5, 
#          PC6 = 6, 
#          PC7 = 7, 
#          PC8 = 8, 
#          PC9 = 9)
## Write a .csv file so we have the data outside of the script
# write_csv(ggplot_FISHY, '.csv')

## PCA PLOT ######
## Load the data
PCA_PLOT = read_csv('.csv')

## Making some better names for the popualtions we're working with
POPS = mutate(.data = POPS,
              Population = as.factor(case_when(
                `#FamilyID` == "1" ~ "Thingvallavatn - Large benthic",
                `#FamilyID` == '2' ~ 'Vatnshlidarvatn - Brown',
                `#FamilyID` == '3' ~ 'Vatnshlidarvatn - Silver',
                `#FamilyID` == '8' ~ 'Svinavatn - Large benthic',
                `#FamilyID` == '4'~ 'Svinavatn - Planktivorous',
                `#FamilyID` == '5' ~ 'Svinavatn - Piscivorous',
                `#FamilyID` == '6' ~ 'Thingvallavatn - Planktivorous',
                `#FamilyID` == '7' ~ 'Thingvallavatn - Small benthic',
                `#FamilyID` == '9' ~ 'Galtabol - Small benthic',
                `#FamilyID` == '10' ~ 'Galtabol - Piscivorous'))) 

Full_PCA_data = bind_cols(POPS, PCA_PLOT)

## Making a colour palete for the PCA
## use as many colors as there are genetic groups
PCAdapt_palete = c('#30C75D', '#94FBD0',
                   '#FA8656', '#94846D',
                   '#947666', '#C78930',
                   '#C76026', '#FA5548',
                   '#87FB96', '#44C726')

ggplot(data = Full_PCA_data, aes(x = PC1, y = PC2))+
  geom_point(aes(col = Population))+
  scale_color_manual(values = PCAdapt_palete)+
  labs(x = 'PC 1', 
       y = 'PC 2')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


## Outliers ####
## Outliers based on the PCA q-values 
## Make sure the the qvalue package is loaded 

## p-values for the snps based on the pca
Q = qvalue(FISHY$pvalues)
## calculating qvalues from the Q matrix
Qvalue = Q$qvalues
## Alpha value for statistical testing
## We're using a pretty stringent alpha value
alpha = 0.01
## Numbers of snps that are statistically significant <0.05
Outlier = which(Qvalue < alpha)
length(Outlier)
## Figure out with principle component the outliers are on
pc_outliers = get.pc(FISHY, Outlier)
pc1 = filter(pc_outliers, PC == 1)
pc2 = filter(pc_outliers, PC == 2)

## Find the pvalues
PVAL = FISHY$pvalues
# 
# ## bind the pvalues to the map file
PVALMAP = as.data.frame(cbind(MAP, PVAL))
PVALMAP = as_tibble(cbind(MAP, PVAL)) %>% 
  rename(CHROME = `#Chromosome`,
         MARKER_ID = `Marker ID`,
         GDIST = `Genetic distance`,
         PDIST = `Physical position`)

Outliers_pop = PVALMAP[Outlier,]

## write a .csv file so we have the outliers
# write_csv(Outliers_pop, '.csv')

## Manhattan plot #####
## read in the pcadapt outliers csv from earlier
Outliers_pop = read_csv('.csv')

non_contigs = Outliers_pop %>% 
  filter(CHROME3 != 'Contigs') %>% 
  arrange(CHROME3)
non_contigs$CHROME3

man_plot_palete = c('#0F203B', '#3E89FA')

man_plot = ggman(PVALMAP, 
                 chrom = 'CHROME3', 
                 pvalue = 'PVAL',
                 snp = 'MARKER_ID',
                 bp = 'PDIST',
                 pointSize = 2, 
                 title = 'CHARR MANHATTAN PLOT',
                 xlabel = 'Chromosome',
                 sigLine = -log10(max(Outliers_pop$PVAL)), 
                 ymax = 20, 
                 lineColour = 'black')

man_plot + 
  scale_color_manual(values = man_plot_palete)+
  theme(axis.text.x =
          element_text(size  = 8,
                       angle = 90,
                       hjust = 1,
                       vjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
