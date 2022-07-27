#fast version of variant aggregation for sparseq v1
#lauren c, jeff wrana's lab
#June 7 2021

#install.packages("dplyr")
#install.packages("openxlsx")
library(dplyr)
library(openxlsx)

#move to workig directory of the sparseq run
setwd("/Volumes/PromisePegasus/SPAR_SEQ/clinical_runs/runcl29")

#import counts
#whatever you named the count table
counts<-read.xlsx("runcl29_countTable.xlsx", sheet = 1)
srbd<-read.csv("SparSeq_Srbd_top1_Var_summary_table.txt", sep = ",", stringsAsFactors = F,)
head(counts)
head(srbd)

srbd<-srbd[,c("sample", "Srbd_total_read.", "Srbd_top1_read.",   "aa_var")]
#aggregate srbd 
srbd<-as.data.frame(srbd %>% 
     group_by(sample) %>% 
     mutate(aa_var = paste(aa_var, collapse = "; ")) )

srbd<-distinct(srbd)

#combine
combined<-merge(counts, srbd, by = "sample", all=T)
combined$type<-"Control"
#W ids are until the end of May, June samples will have X
#keep both because old may samples may show up still
combined[grepl("^W", combined$sample, ignore.case=TRUE), ]$type<-"Sample"
combined[grepl("^X", counts$sample, ignore.case=TRUE), ]$type<-"Sample"
#they are sometimes naming the negative human samples with negative-W#....
combined[grepl("^Neg", combined$sample, ignore.case=TRUE), ]$type<-"Sample"


#split controls and samples
controlsOnly<-subset(combined, combined$type == "Control")
samplesOnly<-subset(combined, combined$type == "Sample")

#cutoff v1 - nmads calculation for NEW PASS FAIL CUTOFF
#Take SRBD of All H2O, HEK, noRTs
#use bowtie counts 
head(controlsOnly)
controlsOnly<-controlsOnly[,c("sample", "Ngene", "ACTB", "ACTG", "PPIB", "Rdrp",  "Spoly", "Srbd")]
#-take median of those values
srbd_median<-median(controlsOnly$Srbd)
#-absolute value of difference between each value and the median value
controlsOnly$medianDiff<-controlsOnly$Srbd - srbd_median
controlsOnly$medianDiff<-abs(controlsOnly$medianDiff)
#-then take the median of that to get MAD
srbd_diff_median<-median(controlsOnly$medianDiff)
#-use original median of the controls + (2xMAD)
#-that value is now the threshold for pass or fail for the QC check
combinedmedians<-srbd_median + (2*srbd_diff_median) #
if (combinedmedians < 10) {combinedmedians<-10} #this is  now a minimum threshold

#import variant details/selected variants
#the actual file is important_vars_detailed_V3.txt
#it is much easier to import this after converting to xlsx 
vardet<-read.xlsx("run29_v1_selectedVariants.xlsx",  sheet=1)
vardet<-vardet[,c( "sample","Srbd_tophit","Spbs_tophit", "Srbd_tophit_percentage", "Srbd_WT_percentage", "Srbd_total_read_count", 
                   "Spbs_tophit_percentage", "Spbs_total_read_count")]
#locate samples where there is weak n501y
# if it has weak n501y with a lot of WT plus pro681arg - may actually be just WT instead of n501y
#due to some issue with the pro681arg blocking primers and causing false positives
#if less than 50% we have to mark it 
# vardet<-subset(vardet, grepl("N501Y", vardet$Srbd_tophit, ignore.case=TRUE))
# vardet<-subset(vardet, vardet$Srbd_tophit_percentage < 50.00)
# srbdTopHitWarnings<-vardet$sample
#use this list to give warnings about low n501y percent
#june 10 - jeff asked to change it to marking rows with a difference > 15 between tophit% and wt% 
#srbd tophit percent and srbd wt percent
vardet$srbd_diff<-abs(vardet$Srbd_WT_percentage - vardet$Srbd_tophit_percentage)
srbdTopHitWarnings<-subset(vardet, ((vardet$srbd_diff < 15) & (vardet$Srbd_tophit != "WT") ) )


#compile results of samples for E484K/Glu484Lys and N501Y/Asn501Tyr
samplesOnly$E484K<-"X"
samplesOnly$N501Y<-"X"
samplesOnly$QC_Check_MinReads<-"X"
samplesOnly$QC_Check_N501Y_E484K_Coverage<-"X"
samplesOnly$QC_Check<-"X" #this is an overall check of the other two checks
#fill in the date of the tests
samplesOnly$Date<-"Processed June 30 2021"

samplesOnly[grepl("Asn501Tyr", samplesOnly$aa_var), ]$N501Y<-"+"
samplesOnly[!(grepl("Asn501Tyr", samplesOnly$aa_var)), ]$N501Y<-"-"
samplesOnly[grepl("Glu484Lys", samplesOnly$aa_var), ]$E484K<-"+"
samplesOnly[!(grepl("Glu484Lys", samplesOnly$aa_var)), ]$E484K<-"-"
samplesOnly[samplesOnly$Srbd > combinedmedians,]$QC_Check_MinReads<-"Pass"
samplesOnly[samplesOnly$Srbd < combinedmedians,]$QC_Check_MinReads<-"Fail"
samplesOnly[samplesOnly$Srbd == combinedmedians,]$QC_Check_MinReads<-"Pass"

#check list of fails 
samplesOnly[samplesOnly$sample %in% srbdTopHitWarnings$sample,]$QC_Check_N501Y_E484K_Coverage<-"Fail"
samplesOnly[!(samplesOnly$sample %in% srbdTopHitWarnings$sample),]$QC_Check_N501Y_E484K_Coverage<-"Pass"
#
samplesOnly$QC_Check <- with(samplesOnly, ifelse(QC_Check_MinReads == "Pass" & QC_Check_N501Y_E484K_Coverage == "Pass", "Pass", "Indeterminate"))

#cut extra columns 
head(samplesOnly)
finalReport<-samplesOnly[,c("sample", "E484K", "N501Y", "QC_Check_MinReads", "QC_Check_N501Y_E484K_Coverage", "QC_Check", "Date")]
head(finalReport)

write.xlsx(finalReport, file="sparseq_report_CLrun29v1.xlsx", sheetName = "CLRun 29v1", row.names = F)
#write.table(finalReport, file="sparseq_report_runXX.txt", quote=F, row.names=F, sep = "\t")


#end


