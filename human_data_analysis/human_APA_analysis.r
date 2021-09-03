library(DESeq2)
library(dplyr)

##################### DESeq2 analysis ##########################

input_file = "~/Downloads/KRAS_Derry_QAPA_TPM.csv"
derry = read.csv(input_file, header = T)[, 1:11]

input_file = "~/Downloads/quant.combined.sf"
counts.data = read.table(file = input_file, 
                         header = T) %>%
  select(Name, contains("NumReads")) %>%
  mutate_if(is.numeric, as.integer)
row.names(counts.data) = counts.data$Name
counts.data = select(counts.data, -Name)

# [Optional] Keep genes with detected reads in at least one sample.
total = rowSums(counts.data)
subset = counts.data[total > 0, ]

# Prepare samples for DESeq
samples = colnames(subset)
column.data = data.frame(samples) %>%
  mutate(samples = gsub("NumReads_Derry_[0-9]_|NumReads_Derry_1[0-9]_", 
                        "", 
                        samples),
         samples = gsub("_S[0-9]_R1_001|_S1[0-9]_R1_001", 
                        "", 
                        samples),
         lysate = gsub("17R|17N|152R|152N", 
                       "", 
                       samples),
         phenotype = gsub("R[1-3]", "R", samples),
         phenotype = gsub("N[1-3]", "N", phenotype),
         kras = gsub("R", "", phenotype),
         kras = gsub("N", "", kras),
         mutation = gsub("17", "", phenotype),
         mutation = gsub("152", "", mutation)) %>%
  mutate(mutation = as.factor(mutation),
         kras = as.factor(kras),
         lysate = as.factor(lysate))
row.names(column.data) = samples

# Create DESeqDataSet object from matrix with counts.
# One could do the analysis by first selecting a week and
# then comparing two conditions. Besides, one could do
# the analysis by introducing additional factor for a week.
# For the help consider documentation at ?DESeq2::results
# Levels: (N, R) and (152, 17)
dds = DESeq2::DESeqDataSetFromMatrix(countData = subset,
                                     colData = column.data,
                                     design =~ kras + mutation + kras:mutation)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds = DESeq2::DESeq(dds)
resultsNames(dds)

res_152 = DESeq2::results(dds, list( c("mutation_R_vs_N") ), alpha = 0.05)
res_152 = data.frame(res_152) %>%
  mutate(Name = row.names(res_152))
res_17 = DESeq2::results(dds, list( c("mutation_R_vs_N", "kras17.mutationR") ), alpha = 0.05)
res_17 = data.frame(res_17) %>%
  mutate(Name = row.names(res_17))

res = full_join(res_152, res_17, by = "Name", suffix = c("_152", "_17"))
res = select(res, Name, contains("152"), contains("17"))

write.csv(x = data.frame(res),
          file = "~/Downloads/differential_APA.csv",
          row.names = F)


##################    APA patterns   #######################


input_file = "~/Downloads/KRAS_Derry_QAPA_TPM.csv"
derry = read.csv(input_file, header = T)[, 1:11]

derry_deseq = read.csv("~/Downloads/differential_APA.csv",
                       header = T) %>%
  mutate(Name = sapply(stringr::str_split(Name, ",", ), last)) %>%
  tidyr::separate(Name, c("Transcript", "Gene_Name", "Genome", 
                          "Chr", "LastExon.Start", "LastExon.End",
                          "Strand", "Region", "UTR3.Start", "UTR3.End"), 
                  sep = "_") %>%
  select(-Region, -Genome, -Transcript) %>%
  mutate(LastExon.Start = as.integer(LastExon.Start),
         LastExon.End = as.integer(LastExon.End),
         UTR3.Start = as.integer(UTR3.Start),
         UTR3.End = as.integer(UTR3.End)) %>%
  left_join(derry) %>%
  filter(!is.na(APA_ID))

lfc_threshold = 0.5849625 # FC = 1.5
switches_152_deseq = group_by(derry_deseq, Gene) %>%
  filter(any(log2FoldChange_152 > lfc_threshold & padj_152 < 0.05) & 
           any(log2FoldChange_152 < -lfc_threshold & padj_152 < 0.05)) %>%
  filter(n() == 2) %>%
  arrange(Gene) %>%
  mutate(APA_type = "switch") %>%
  select(-contains("quants"), -contains("17"))

switches_17_deseq = group_by(derry_deseq, Gene) %>%
  filter(any(log2FoldChange_17 > lfc_threshold & padj_17 < 0.05) & 
           any(log2FoldChange_17 < -lfc_threshold & padj_17 < 0.05)) %>%
  filter(n() == 2) %>%
  arrange(Gene) %>%
  mutate(APA_type = "switch") %>%
  select(-contains("quants"), -contains("152"))

proximal_152_deseq = group_by(derry_deseq, Gene) %>%
  mutate(proximal = as.integer(Length == min(Length))) %>%
  filter(abs(log2FoldChange_152) > lfc_threshold & padj_152 < 0.05) %>%
  filter(n() == 1 & proximal == 1) %>%
  ungroup() %>%
  mutate(APA_type = "proximal") %>%
  arrange(Gene) %>%
  select(-contains("quants"), -contains("17"), -proximal)

proximal_17_deseq = group_by(derry_deseq, Gene) %>%
  mutate(proximal = as.integer(Length == min(Length))) %>%
  filter(abs(log2FoldChange_17) > lfc_threshold & padj_17 < 0.05) %>%
  filter(n() == 1 & proximal == 1) %>%
  ungroup() %>%
  mutate(APA_type = "proximal") %>%
  arrange(Gene) %>%
  select(-contains("quants"), -contains("152"), -proximal)

distal_152_deseq = group_by(derry_deseq, Gene) %>%
  mutate(proximal = as.integer(Length == min(Length))) %>%
  filter(abs(log2FoldChange_152) > lfc_threshold & padj_152 < 0.05) %>%
  filter(n() == 1 & proximal != 1) %>%
  ungroup() %>%
  mutate(APA_type = "distal") %>%
  arrange(Gene) %>%
  select(-contains("quants"), -contains("17"), -proximal)

distal_17_deseq = group_by(derry_deseq, Gene) %>%
  mutate(proximal = as.integer(Length == min(Length))) %>%
  filter(abs(log2FoldChange_17) > lfc_threshold & padj_17 < 0.05) %>%
  filter(n() == 1 & proximal != 1) %>%
  ungroup() %>%
  mutate(APA_type = "distal") %>%
  arrange(Gene) %>%
  select(-contains("quants"), -contains("152"), -proximal)

APA_17_deseq = data.frame(bind_rows(switches_17_deseq, distal_17_deseq, proximal_17_deseq))
APA_152_deseq = data.frame(bind_rows(switches_152_deseq, distal_152_deseq, proximal_152_deseq))

APA_17_deseq %>% 
  group_by(APA_type) %>% 
  summarise(n = length(unique(Gene_Name)))
APA_152_deseq %>% 
  group_by(APA_type) %>% 
  summarise(n = length(unique(Gene_Name)))

write.csv(APA_17_deseq,
          "~/Downloads/tandem_APA_17_deseq.csv",
          row.names = F)
write.csv(APA_152_deseq,
          "~/Downloads/tandem_APA_152_deseq.csv",
          row.names = F)