# ====================================================================
# PCA_analysis.R  |  Figure S4 (sub-cluster survival Cox / PCA)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : AML-sample-infor.csv, Combo_AML_enriched_AML_gene_by_lineage_Sample_fraction.csv, Combo_AML_gene_by_lineage_Sample_fraction.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript PCA_analysis.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

rm(list=ls())
suppressPackageStartupMessages({
  library(tidyverse)
})
library(ggplot2)
library(reshape2)
library(dplyr)
library("FactoMineR")
library("factoextra")
#https://cran.r-project.org/web/packages/factoextra/readme/README.html
setwd("~/Downloads/ALSF_AML_plot")

df_RSS <- read.table("AML_Cell_30_RSS_Leukemia_Enriched_Cell_Type_Sample.csv",sep = ",", header = TRUE)
df_RSS[2:ncol(df_RSS)]<- lapply(df_RSS[2:ncol(df_RSS)], FUN = function(y){as.numeric(y)})
colnames(df_RSS)
colnames(df_RSS)[1] <- 'Sample'
df_RSS<-df_RSS %>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))
colnames(df_RSS)<-gsub("_...","_(+)",colnames(df_RSS))
topreg<-read.table("topreg_enriched Leukemia_with_normal_BM.csv",sep = ",", header = TRUE)
df_RSS_f<-subset(df_RSS, select=c('Sample',topreg$Regulon))
df_RSS_f<-df_RSS[,-2]                
df_clic <- read.delim2("AML-sample-infor.csv", sep = ",", header = TRUE)
df_clic <-df_clic %>%mutate(Prognosis=0)%>%mutate(Prognosis=ifelse(Alive.deceased=='Alive',0,Prognosis))%>%
  mutate(Prognosis=ifelse(Alive.deceased=='Deceased',1,Prognosis))
df_clic$MRD<-df_clic$MRD...
df_clic$Blast<-df_clic$Blast..
df_clic <-merge(df_clic, df_Cell_Blast, by="Sample")
df <-merge(df_clic,df_RSS_f, by="Sample", all = T)

df<-df_RSS_f[,-1]
rownames(df)<-df_RSS_f$Sample

P <- PCA(df, scale.unit=TRUE, ncp=10, graph=F) 

var <- get_pca_var(P)
fviz_pca_var(res.pca, col.var = "black")
fviz_ca_biplot(P, repel = TRUE)
fviz_pca_ind(P, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(P, repel = TRUE, # Avoid text overlapping (doesn't scale to large datasets)
                col.var = "red", # Variables color
                col.ind = "black")
fviz_eig(P, addlabels = TRUE)

fviz_contrib(P, choice = "var", axes = 1, top = 10)
fviz_contrib(P, choice = "var", axes = 2, top = 10)
fviz_contrib(P, choice = "ind", axes = 1:2)
fviz_contrib(P, choice = "ind", axes = 3:4)
fviz_contrib(P, choice = "var", axes = 3, top = 6)
fviz_contrib(P, choice = "var", axes = 4, top = 6)
fviz_contrib(P, choice = "var", axes = 5, top = 6)
fviz_contrib(P, choice = "var", axes = 6, top = 6)

fviz_pca_ind(P, label="none", habillage=df$Alive.deceased,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2")


set.seed(123)
km.res <- kmeans(scale(df), 9, nstart = 25)
fviz_cluster(km.res, data = df,
             palette = c('#d62728',
                         '#aa40fc',
                         '#8c564b',
                         '#e377c2',
                         '#b5bd61',
                         '#17becf',
                         '#aec7e8',
                         '#ffbb78',
                         '#98df8a',
                         '#ff9896'),
             ggtheme = theme_minimal(),
             repel = TRUE,
             pointsize = 2, labelsize = 12,
             main = "Partitioning Clustering Plot"
)

fviz_nbclust(df, kmeans, method = "gap_stat")
res <- hcut(df, k = 4, stand = TRUE)
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))


df_gene<-read.delim2("Combo_AML_enriched_AML_gene_by_lineage_Sample_fraction.csv", sep = ",", header = TRUE)
#df_gene<-read.delim2("Combo_AML_gene_by_lineage_Sample_fraction.csv", sep = ",", header = TRUE)

#colnames(df_filtered)

genelist<-c('NPR3', 'ANGPT1', 'MECOM', 'ATP8B4', 'SELENOP', 'CD34', 'HOPX', 'BAALC', 'HPGDS',
            'SLC45A3', 'RHEX', 'SPNS2', 'MYCN', 'SPINK2', 'PROM1', 'HOXB6', 'CLEC3B', 'CDH4', 
            'CDH9', 'MALL', 'PARD3B', 'PDZRN4', 'FUT6', 'PRSS2', 'RHBDF1', 'DLK1', 'MFAP4', 
            'ITGA2', 'TM4SF1', 'TMEM98', 'JCAD', 'MMP7', 'MEIS1', 'MDK', 'NPDC1', 'ARHGEF17', 
            'RBPMS', 'SPARC', 'C1QTNF4', 'PLAU', 'SHANK3', 'SVOPL', 'SH3D21', 'TCTEX1D1', 'OBSL1',
            'ZNF521', 'KRT18', 'SV2A', 'COL24A1', 'HOXA3', 'NOG', 'INPP4B', 'SCHIP1', 'GCSAML', 
            'DYTN', 'KIAA1549', 'SLC2A5', 'MMRN1', 'NBL1', 'MN1', 'MPL', 'DOCK1', 'RAMP1', 
            'SLC18A2', 'MRC2', 'AIF1L', 'ITGA9', 'HOXB.AS3', 'HOXB3', 'CERCAM', 'C9orf43', 'AC016735.1',
            'RAB27B', 'TIE1', 'HTRA3', 'CAPN11', 'RXFP1', 'PADI3', 'FGF13', 'TIMP3', 'DEPTOR',
            'ZNF385D', 'HOXA9', 'BCL6B', 'WFDC1', 'LOX', 'TNNT3', 'EFNA1', 'VWA1', 'PHLDB1', 
            'HOXB7', 'HOXA6', 'LRP6', 'CCDC188', 'RHPN1', 'NPW', 'ANGPTL6', 'SHF', 'LTC4S',
            'TTC24', 'ATL1', 'ZBTB8A', 'ROBO3', 'TEC', 'SMIM24', 'IFIT5', 'MAP1A', 'EVA1B',
            'FGFBP3', 'ANGPT2', 'TOX', 'NECTIN2', 'ABHD17C', 'MEGF6', 'DOCK5', 'MLLT11', 'FHL1', 
            'C1orf21', 'PIK3R6', 'FUT7', 'ARID5B', 'PRKCH', 'NAB2', 'AGAP2', 'SIAE', 'GOT1', 
            'SAMD10', 'TRAF5', 'LKAAEAR1', 'AL034397.3', 'PGPEP1', 'KCNQ1OT1', 'GBP5', 'CCDC102A',
            'ESAM', 'AL590226.1', 'GIPC1', 'DUSP14', 'TOR4A', 'RGS3', 'EVPL', 'RTP4', 'B3GNTL1',
            'GNB5', 'ARSD', 'SVIL', 'FBXO41', 'PROK2', 'SCARF1', 'PDGFC', 'GPSM1', 'SLC9A3R2', 
            'BEND5', 'PLEKHA5', 'BANK1', 'LATS2', 'CXorf21', 'VSIG10', 'TBC1D16', 'LINC01268', 
            'RTN4R', 'ARAP2', 'ZNF66', 'C9orf139', 'ITGAV', 'ITGA5', 'GDF11', 'ACY3', 'IER5L', 
            'FAN1', 'GBP4', 'AC002454.1', 'GUCY1A1', 'TENT5A', 'PPM1H', 'GNG11', 'BMP1', 'AR',
            'MROH6', 'FAR2', 'SRGAP3', 'TRIM73', 'GPR174', 'THUMPD2', 'PPFIBP1', 'KCNQ5', 'CD69',
            'LTB4R2', 'FOSL1', 'GPX7', 'CTSF', 'TMIGD2', 'CRELD1', 'NR1H3', 'AGL', 'SDSL', 
            'SPEF2', 'HIST2H2BE', 'ENDOD1', 'IKZF2', 'ZNF254', 'AC090152.1', 'MZB1', 'NOTCH1', 
            'CHST14', 'MIR155HG', 'LINC01684', 'CLIP2', 'INKA1', 'HSPG2', 'MAPK12', 'ZBTB16', 
            'ATP1B1', 'TAPT1.AS1', 'AC092279.1', 'ZNF587B', 'ANKRD33B', 'SERPING1', 'TESPA1', 
            'DST', 'SORBS3', 'SMYD3', 'ATF3', 'GSEC', 'STXBP5', 'RTCA.AS1', 'DNASE2', 'ZNF737',
            'RUNX2', 'SYNJ2BP', 'ZNF660', 'HERC6', 'CUL7', 'MATK', 'CRIM1', 'GOLGA8A', 'KLHL13', 
            'AC246787.2', 'GCNT1', 'EMILIN1', 'DUSP10', 'MICU3', 'ALCAM', 'CLDN15', 'PALM', 
            'ABHD4', 'ATP8A1', 'EPB41L2', 'B4GALT6', 'EIF3J.DT', 'DDX60', 'SFXN3', 'CRYL1',
            'IFT81', 'MYO18A', 'HAGHL', 'TRAF3IP2', 'CASC15', 'KCNMB4', 'AEBP1', 'MSRB3', 
            'AC243960.1', 'CYSLTR1', 'TCEAL1', 'PLCB1', 'TMEM44', 'CDIP1', 'EPHB6', 'DENND1B',
            'EMP1', 'AHDC1', 'CD244', 'OAS2', 'AC002310.1', 'F2R', 'NAALADL1', 'NCK2', 'PTP4A3',
            'PHOSPHO2', 'TFPI', 'STMN3', 'ITGA6', 'S100Z', 'FAIM', 'CD200', 'KCNK17', 'IDH1', 
            'RTKN', 'COL9A2', 'SLC37A1', 'CLIP3', 'PTK7', 'NFIX', 'DENND3', 'UBASH3B', 'GMDS.DT',
            'SOCS5', 'CCL28', 'RBMS2', 'CPXM1', 'CBX1', 'APOL6', 'IFITM1', 'AL138724.1', 'ZEB1',
            'MED12L', 'STK32B', 'CCDC171', 'PLTP', 'MEX3B', 'GNAI1', 'CYB561', 'FAM30A', 'PDE6G', 
            'SLC4A7', 'CLU', 'CD109', 'MMP28', 'ACCS', 'LINC00963', 'ST3GAL4', 'CCDC189', 'CLSTN3',
            'MPI', 'JUP', 'KRBA2', 'ERG', 'TRIM6', 'LNCAROD', 'TSPAN7', 'GPR12', 'NR1I2', 'NPM2', 
            'HAAO', 'GOLGA8N', 'ACTA2', 'GTSF1', 'RTN2', 'IGSF10', 'DPP10', 'ECM1', 'NRXN2',
            'PTPRCAP', 'HSF4', 'KIF26B', 'RARRES2', 'PRX', 'MRC1', 'SOBP', 'DSC2', 'WASIR2', 
            'LINC02175', 'IL1RAP', 'JAG1', 'SLCO5A1', 'AC099489.1', 'RAB44', 'NT5DC3', 'TMTC2', 
            'TRIM46', 'ADGRE1', 'WDR49', 'TRPM4', 'SLC15A2', 'SNAI3', 'SUCNR1', 'LAMB2', 'NRG4',
            'DGKG', 'RGS18', 'DNAH1', 'LTK', 'CRIP2', 'AC011446.2', 'STAC3', 'FGFR1', 'AFF2', 
            'AC127502.2', 'CASS4', 'AK5', 'FBXL19', 'LGALS12', 'SNHG4', 'LINC01341', 'RBKS', 
            'TTC28', 'NTNG2', 'ZG16B', 'NDST1', 'AMN', 'SAP25', 'RAB7B', 'RFX8', 'TPSAB1', 
            'THSD7A', 'LINC00539', 'EGF', 'DLGAP2', 'KAZN', 'GJA4', 'AGAP1', 'GP9', 'GOLGA8O',
            'GOLGA8R', 'GALNT5', 'AC009315.1', 'GPC1', 'MYO5C', 'MS4A4E', 'CACNA2D4', 'CFD', 
            'CD14', 'VSTM1', 'CSF1R', 'S100A9', 'MNDA', 'PTAFR', 'STAB1', 'OLIG1', 'VCAN',
            'HMOX1', 'LRRC25', 'CXCL16', 'FPR1', 'CDA', 'CD33', 'CLEC4E', 'RETN', 'S100A12', 
            'CLEC4G', 'HCAR2', 'CLEC5A', 'CXCL1', 'CETP', 'CLEC12B', 'HAL', 'SERPINB2', 'MSLN', 
            'CEBPE', 'HLX', 'ANPEP', 'PTGDS', 'CD96', 'CD7', 'LAT', 'CCL3', 'GNLY', 'PRF1', 
            'CD3E', 'BCL11B', 'TRBC1', 'PBX4', 'TRGC1', 'NMUR1', 'NLRP2', 'FEZ1', 'TRBV28', 
            'AP005482.1', 'TSPAN2', 'TRGC2', 'CD3D', 'SPON2', 'AC108134.3', 'GALNT3', 'PHLDA1', 
            'LINC00649', 'TNFRSF4', 'LAX1', 'PGGHG', 'TRGV9', 'NT5E', 'CST7', 'ADRB2', 'ZAP70',
            'PDE3B', 'SH2D1A', 'CCL5', 'TRBC2', 'RUNX3', 'ECHDC2', 'CD9', 'CD19', 'FCER2', 'DNTT',
            'PAX5', 'CSRP2', 'CLEC14A', 'UMODL1', 'ZDHHC19', 'VAT1L', 'CD180', 'LAMA5', 'MILR1', 
            'PPP1R14A', 'CHD7', 'GRAP', 'IGF2BP3', 'STAP1', 'LINC00926', 'PIK3C2B', 'IGHM', 
            'CDKN2C', 'CDKN2A', 'KYNU', 'RASGRP3', 'BIK', 'C16orf74', 'GALNT14', 'PXDN', 
            'SLC24A3', 'CLC', 'PRKAR2B', 'KLF1', 'GATA1', 'XK', 'TFR2', 'VWF', 'NTRK1', 'TPSD1',
            'FAM171A1', 'DEPDC7', 'LINC01515', 'PRCD', 'LTBP1', 'DLC1', 'P2RY1', 'AQP3', 'GDPD1',
            'ADD2', 'PCAT18', 'PTH2R', 'PLA2G4A', 'TMEM156', 'SELENBP1', 'PLD1', 'MYL4', 
            'MICALL2', 'MFSD3', 'TSGA10', 'TRIM66', 'CYP2E1', 'ITGA2B', 'HILPDA', 'CNRIP1', 
            'TMEM102', 'TNFRSF25', 'PRKAR1B', 'RAB6B', 'CSTF1', 'BPI', 'NETO2', 'CPA3', 'HGF', 
            'MPO', 'AZU1', 'TSPOAP1', 'TRH', 'PRTN3', 'ELANE', 'PLPPR3', 'CTSG', 'IGFBP2', 'STAR',
            'CKM', 'ADGRG3', 'CPNE7', 'TMEM221', 'PRRT4', 'ARTN', 'PIWIL4', 'ITGA7', 'BCL2L10', 
            'CCL23', 'ILDR2', 'KCNE5', 'GFI1', 'AC009961.4', 'CPT1B', 'MMP9', 'STX1A', 'ABTB2',
            'IL10', 'CBX2', 'SLITRK4', 'RNF175', 'BGLAP', 'TRIM71', 'S100A1', 'NETO1', 'MMP19', 
            'COL2A1', 'PSD3', 'MT1G', 'FCGR1B', 'PRLR', 'RET', 'SIGLEC1', 'SLPI', 'PLB1', 
            'S100B', 'TPSB2', 'IRX1', 'GAS2L1', 'RNASE3', 'CLECL1', 'IL5RA', 'TNFSF9', 'COL23A1',
            'MS4A3', 'FCGR1A', 'FBLN2', 'HOMER3', 'WNT5B', 'GPR34', 'TBL1Y', 'GPR153', 'SHOX2',
            'VASN', 'LRRC7', 'HES1', 'FST', 'GCNT4', 'RELN', 'MYO7B', 'EGR2', 'CLEC4D', 'FRMD3',
            'ATP8A2', 'AC139493.2', 'PPARG', 'DEFB1', 'AC097637.3', 'CILP2', 'AGBL4', 'FAM81B',
            'IRGM', 'APBA1', 'COL14A1', 'PHACTR3', 'PLCH1', 'LBX2', 'CCNA1', 'AC125603.2', 
            'AC244502.1', 'PITX1', 'RUNX1T1', 'ST18', 'GGT5', 'MYO18B', 'KRT17', 'AGXT', 'IL9R',
            'GJA1', 'TTLL10', 'IL1R1', 'LINC01114', 'LCT', 'PLA2R1', 'SCN2A', 'AC104088.3', 
            'COL5A2', 'RFTN2', 'AC020550.2', 'EGOT', 'AC134508.2', 'KCNMB2', 'LINC01091',
            'ADCY2', 'FAT1', 'ADAMTS12', 'FGF10', 'AC008945.1', 'LINC02147', 'SNCAIP', 'FOXC1', 
            'FOXCUT', 'MDFI', 'HOXA2', 'HOXA11', 'TRGV6', 'INHBA', 'OPN1SW', 'SCARA3', 'WT1',
            'PMCH', 'IRX3', 'IRX5', 'CCL1', 'ODF4', 'KRT23', 'AC103702.1', 'HOXB8', 'PIEZO2', 
            'SLC52A3', 'TBX1', 'TMPRSS3', 'BMX', 'TSIX', 'COL4A5', 'LINC01645', 'PTPN14', 
            'SLC27A6', 'SLC7A2', 'SYT4', 'AC027601.6', 'DUXAP8', 'GPC4', 'AC013267.2', 
            'LINC01833', 'SIX3', 'AL596442.3', 'AGR2', 'AC104232.1', 'AL392086.3', 'PDZD7',
            'ADRB1', 'AP005273.2', 'PPP1R27', 'PTGER1', 'APOC2', 'AL139220.2', 'HORMAD1', 'SHE', 
            'TDRD10', 'OLFML2B', 'NR5A2', 'SYT2', 'ITGB5', 'SPSB4', 'UNC5C', 'NKD2', 'AC116351.1',
            'PDGFRB', 'UNCX', 'HOXA13', 'CAV2', 'TRBVB', 'ADAMDEC1', 'OLFML2A', 'EBF3', 'KITLG',
            'DCN', 'POU4F1', 'CD276', 'CSPG4', 'SEZ6L2', 'ASPHD1', 'NTN1', 'HOXB9', 'ABCA9', 
            'ABCA6', 'DOCK6', 'SIGLEC12', 'BMP2', 'AL118508.1', 'SPP1', 'TSLP', 'MAMDC2', 'MKRN3',
            'FENDRR', 'FOXF1', 'SIPA1L2', 'AC083855.2', 'AC044893.1', 'VLDLR', 'POSTN', 'TUSC8', 
            'GSX1', 'AL121672.3', 'PRAME', 'NKD1', 'LINC01731', 'HIST2H2AA4', 'DUSP27', 
            'AL591848.4', 'TPO', 'IL1RL1', 'NRP2', 'CCR3', 'AC111000.4', 'UGT2B11', 'FLT4',
            'SYCP2L', 'Z98745.2', 'NFE4', 'FBXL13', 'CLEC2L', 'KLRG2', 'GLIS3', 'ADGRD2', 'LHX6', 
            'IL2RA', 'PCDH15', 'LINC00958', 'ITPKA', 'PAK6', 'CACNA1H', 'PTX4', 'AL031710.1', 
            'IRX6', 'FOXL1', 'SOX15', 'RBFOX3', 'USP41', 'CFAP47', 'MYH11', 'GREM1', 'MEIS2', 
            'S100A16', 'MYBPH', 'CDK15', 'MAB21L4', 'CHDH', 'IL17RB', 'ENPP3', 'LY6K', 'SPINK4', 
            'SLC28A3', 'TCN1', 'MYRF', 'AC061961.1', 'KCNJ3', 'TIMP4', 'EPHA3', 'LINC02506', 
            'PCDHGB2', 'GDA', 'GPR152', 'SLITRK6', 'ALPK3', 'HS3ST4', 'SULT1A2', 'AC092127.2', 
            'SSC5D', 'FAM43B', 'PTPRU', 'TEKT2', 'DAB1', 'LINC02609', 'LINC01357', 'PDZK1', 
            'SYCP1', 'AQP10', 'ESRRG', 'LYPD8', 'GREB1', 'EGR4', 'SLC9A2', 'GLI2', 'ABCB11', 
            'AC012499.1', 'ASIC4', 'SUSD5', 'FBXL2', 'CCK', 'MAGI1', 'ROBO2', 'IGSF11', 
            'LINC02010', 'KDR', 'CDH12', 'SLC13A3', 'CTCFL', 'REG4', 'DMRTA1', 'SNORC', 'MAPK10', 
            'CNTN5', 'LINC01538', 'CEACAM6', 'WASIR1', 'PGBD5', 'NNMT', 'NPIPB2', 'MEG8', 'GABRE',
            'SERPINI2', 'MYT1L', 'AC093390.2', 'TBX18', 'KIAA1217', 'RFLNA', 'AC008403.3', 
            'SCUBE1', 'PTGFR', 'CSN1S1', 'LIFR', 'LGSN', 'LZTS1', 'AC090015.1', 'PGF', 'APLN', 
            'AC004556.3', 'EXOC3L2', 'PAQR5', 'MGP', 'TMEM52B', 'EDIL3', 'MYLK4', 'ANKRD18B',
            'HPSE2', 'H19', 'PYGM', 'SP7', 'LHX9', 'CNTN4', 'FOXL2', 'HMX3', 'ARC', 'LINC01436',
            'RASL10A', 'GNG4', 'BEGAIN', 'DOC2A', 'AC091138.1', 'PDGFA', 'SPAG6', 'JPH4', 'CCDC3',
            'LINC02284', 'ANO7', 'AL355597.1', 'MATN4', 'KCNN2', 'PENK', 'SDC2', 'TGM5', 'RHOXF2B',
            'MEGF10', 'SOSTDC1', 'TACC2', 'LRP4', 'ADGRL2', 'SATB2', 'PRDM5', 'NDNF', 'AC010378.2', 
            'MLXIPL', 'FOXD4', 'AL449043.1', 'AL450311.1', 'HMX2', 'LINC02381', 'NAV3', 'C12orf56', 
            'PRAC2', 'AC079466.1', 'NXF3', 'CT45A1', 'CT45A10', 'MAGEA4',
            'SLC44A3', 'MNX1', 'BDNF', 'C2CD4B', 'BHLHE23', 'LPAR4', 'CAV1', 'NDST3', 'SELP',
            'THY1','IL3RA','CD33','FCGR2A','CD36','TFRC','GYPA','GYPB',
            'IL7R','CD19','MME','PDCD1','CD274','CTLA4','HAVCR2','PTPRC','CD3D','CD3E','CD3G')
df_filtered<-df_gene%>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))
df_filtered<-subset(df_filtered, select=c('Sample',genelist))
df_filtered<-df_gene[,-2]
df_filtered[2:ncol(df_filtered)] = lapply(df_filtered[2:ncol(df_filtered)], FUN = function(y){as.numeric(y)})

df<-df_filtered[,-1]
rownames(df)<-df_filtered$Sample

P <- PCA(df[,-1], scale.unit=TRUE, ncp=10, graph=F) 

var <- get_pca_var(P)
fviz_pca_var(res.pca, col.var = "black")
fviz_ca_biplot(P, repel = TRUE)
fviz_pca_ind(P, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(P, repel = TRUE, # Avoid text overlapping (doesn't scale to large datasets)
                col.var = "red", # Variables color
                col.ind = "black")
fviz_eig(P, addlabels = TRUE)

fviz_contrib(P, choice = "var", axes = 1, top = 6)
fviz_contrib(P, choice = "var", axes = 2, top = 6)
fviz_contrib(P, choice = "var", axes = 3, top = 6)
fviz_contrib(P, choice = "var", axes = 4, top = 6)
fviz_contrib(P, choice = "var", axes = 5, top = 6)
fviz_contrib(P, choice = "var", axes = 6, top = 6)

fviz_pca_ind(P, label="none", habillage=df$Alive.deceased,
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2")

df<-df[,-1]

rownames(df)<-df_RSS_f$Sample

df[is.na(df)] <- 0
set.seed(123)
km.res <- kmeans(scale(df), 6, nstart = 25)
fviz_cluster(km.res, data = df,
             palette = c('#d62728',
                         '#aa40fc',
                         '#8c564b',
                         '#e377c2',
                         '#b5bd61',
                         '#17becf',
                         '#aec7e8',
                         '#ffbb78',
                         '#98df8a',
                         '#ff9896'),
             ggtheme = theme_minimal(),
             repel = TRUE,
             main = "Partitioning Clustering Plot"
)

fviz_nbclust(df, kmeans, method = "gap_stat")
res <- hcut(df, k = 4, stand = TRUE)
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))