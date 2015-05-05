#' Supplementary Matrical for MetaPathways v2.5: Quantitative functional, taxonomic, and usability improvements
#' ---------------------------------------------------------------------------------
#' Kishori M. Konwar, Niels W. Hanson, Maya P. Bhatia, Dongjae Kim, Shang-Ju Wu, Aria S. Hahn, Connor Morgan-Lang, Hiu Kan Cheung, and Steven J. Hallam
#'
#' The MetaPathways v2.5 release is available from GitHub with binaries compatible for 
#' various Linux and Mac OSX Operating Systems:
#' 
#' * <http://github.com/hallamlab/metapathways2/releases>

#' * Load required libraries

library(ggplot2)
library(reshape2)
library(plyr)

#' ## Weighted Taxonomic Distance (WTD)
#' 
#' * Demonstrate the classification of pathways by the WTD. 
#' * Predicted MetaCyc pathways with WTD vallues can be found in the `<sample>/results/pgdbs/` folder

wd <- "~/Dropbox/manuscripts/metapathways_2_5/analysis/data/"
setwd(wd)
df1 <- read.table("01_pathways_4093112.pwy.txt", header=T, sep="\t")
last <- read.table("02_last_blast_test_last_vs_cog_mscore_40.txt", comment.char="#", sep="\t", header=FALSE)
blast <- read.table("02_last_blast_test_blast_vs_cog_mscore_40.txt", comment.char="#", sep="\t", header=FALSE)
ORFs <- read.table("03_rpkm_master_table_ORF_count.txt", sep="\t", header=T)
RPKMs <- read.delim("03_rpkm_master_table_rpkm.txt", sep="\t", header=T)

df1$hazard_class = "None"

# calculate summary statistics of negative WTD scores and set hazard classes
wtd_stats <- summary(df1$WTD[df1$WTD < 0])
df1$hazard_class[df1$WTD > wtd_stats["Median"]] = "Low"
df1$hazard_class[df1$WTD > wtd_stats["1st Qu."] & df1$WTD <= wtd_stats["Median"]] = "Medium"
df1$hazard_class[df1$WTD <= wtd_stats["1st Qu."]] = "High"
df1$hazard_class[df1$WTD >= 0] = "None"

# factor new levels
df1$hazard_class <- factor(df1$hazard_class, levels=c("None", "Low", "Medium", "High"))

# create Figure 1a
p1 <- ggplot(df1, aes(x=WTD)) + geom_histogram(aes(fill=hazard_class), binwidth=0.2)
p1 <- p1 + scale_fill_manual(values=c("#CCCCCC", "#AFE591", "#FAA982", "#DF737D")) 
p1 <- p1 + theme_bw(base_family="Gill Sans")
p1 <- p1 + ylab("Pathway Frequency")
p1

#' ## LAST Bit-score and E-value
#'
#' * related to the bitscore and evalue

last$algo = "LAST"
blast$algo = "BLAST"

df <- rbind(last[-13], blast)
df <- cbind(paste(df$V1, df$V2, sep = "_"),df)

header <- c("hit", "query", "subject", "identity", "al_length", "missmatch", "gap_open", 
            "query_start", "query_end", "subject_start", "subject_end",
            "eval", "bit_score", "algo")

colnames(df) <- header
df.m <- melt(df)

## Fit Linear Model
df.m.eval <- dcast(subset(df.m, variable=="eval"), hit~algo, max)
df.m.eval[df.m.eval == "-Inf"] = NA
df.m.eval <- df.m.eval[complete.cases(df.m.eval),]
eval_fit <- lm(log(LAST) ~ log(BLAST), data=df.m.eval)

summary(eval_fit)

p2 <-ggplot(df.m.eval, aes(x=BLAST, y=LAST))
p2 <- p2 + geom_point(size=3, alpha=0.5)
p2 <- p2 + theme_bw(base_family="Gill Sans")
p2 <- p2 + stat_smooth(formula="y ~ x", method="lm", se = FALSE, size=2)
p2 <- p2 + scale_x_log10() + scale_y_log10() 
p2 <- p2 + xlab("BLAST E-value") + ylab("LAST E-value")
p2

## ORF Counts vs. Reads per kilobase per-million mapped (RPKM)

ORFs$measure <- "ORFs"
RPKMs$measure <- "RPKMs"

ORFs.m <- melt(ORFs, id.vars=c("pwy_short", "pwy_long", "measure"))
RPKMs.m <- melt(RPKMs, id.vars=c("pwy_short", "pwy_long", "measure"))
df <- rbind(ORFs.m, RPKMs.m)
df.measure <- dcast(df, pwy_short*variable~measure)
df.measure$ORFs <- as.numeric(df.measure$ORFs)

# remove those with both zero
df.measure.sub <- subset(df.measure, ORFs != 0 & RPKMs != 0)
p3 <- ggplot(df.measure.sub, aes(ORFs, RPKMs)) 
p3 <- p3 + geom_point(alpha=0.3) 
p3 <- p3 + scale_x_log10(breaks=c(1,10,100,1000)) 
p3 <- p3 + scale_y_log10(breaks=c(1,10,100,1000,10000)) 
p3 <- p3 + stat_smooth(method="lm", size=2) 
p3 <- p3+ theme_bw(base_family="Gill Sans")

rpkm_fit <- lm(log(RPKMs)~log(ORFs),data=df.measure.sub)
summary(rpkm_fit)

