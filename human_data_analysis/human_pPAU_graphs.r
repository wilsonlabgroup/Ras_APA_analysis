library(dplyr)
library(ggplot2)

# Download tables with TPM and PAU values for isoforms. Enter paths to downloaded files
tpm_full = readxl::read_xlsx("~/Downloads/abh0562_Suppl. Excel_seq4_v1.xlsx") %>%
           filter(Gene_Name != "NA",
                  !grepl("4[3-4][0-9]", Gene_Name))
tpm_tbl = reshape2::melt(tpm_full %>%
                           group_by(Gene_Name) %>% 
                           filter(n()>1) %>%
                           select(Gene_Name, contains("mean"))) %>%
          mutate(variable = gsub("_mean", "", variable)) %>%
          tidyr::separate(variable, c("X1", "X2", "genotype", "KD_condition")) %>%
          select(-X1, -X2) %>%
          mutate(sample = paste(genotype, KD_condition, sep = "_")) %>%
          group_by(Gene_Name, sample) %>%
          summarise(TPM = sum(value, na.rm = T))
pau_full = readxl::read_xlsx("~/Downloads/abh0562_Suppl. Excel_seq6_v1.xlsx") %>%
           filter(Gene_Name != "NA",
                  !grepl("44[0-9]", Gene_Name)) %>%
           mutate_at(vars(matches("pPAU")), as.numeric)
pau_tbl = reshape2::melt(pau_full %>%
                          select(Gene_Name, contains("pPAU"))) %>%
          mutate(variable = gsub("[0-9]_pPAU", "", variable)) %>%
          tidyr::separate(variable, c("X1", "X2", "genotype", "KD_condition")) %>%
          select(-X1, -X2) %>%
          mutate(sample = paste(genotype, KD_condition, sep = "_")) %>%
          group_by(Gene_Name, sample) %>%
          summarise(value_mean = mean(value, na.rm = T),
                    value_sd = sd(value, na.rm = T)/sqrt(3))

PAU_diff = pau_tbl %>% 
           left_join(tpm_tbl, by = c("Gene_Name", "sample")) %>%
           tidyr::separate(sample, c("genotype", "knockdown"), "_") %>%
           group_by(Gene_Name, genotype) %>% 
           filter(n() == 2) %>%
           summarise(PAU_ctrl = value_mean[knockdown == "siRLUC"],
                     TPM = (TPM[knockdown == "siCFIm25"] + TPM[knockdown == "siRLUC"])/2,
                     diff = value_mean[knockdown == "siCFIm25"] - value_mean[knockdown == "siRLUC"],
                     Z = diff/sqrt(value_sd[knockdown == "siCFIm25"]**2 + value_sd[knockdown == "siRLUC"]**2),
                     Z = abs(Z)) %>%
           arrange(Z) %>%
           filter(!is.na(Z) & Z < 1e6) %>%
           mutate(Z = ifelse(Z > 50, 50, Z))

ggplot(data = pau_tbl,
       aes(x = sample,
           y = value_mean)) +
  geom_boxplot(notch = T, fill = "green", alpha = 0.5) +
  ggthemes::theme_base() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = PAU_diff,
       aes(x = log(TPM, 10),
           y = diff,
           color = Z)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  ggthemes::theme_base() +
  scale_color_distiller(palette = "Greens") +
  xlab(latex2exp::TeX("log_{10}(TPM)")) +
  ylab(latex2exp::TeX("pPAU_{KD} - pPAU_{CTRL}")) +
  facet_wrap(~genotype) +
  coord_cartesian(ylim = c(-100, 100))

ggplot(data = PAU_diff %>%
         filter(Z > 2),
       aes(x = diff)) +
  geom_density(notch = T) +
  geom_vline(xintercept = 0) +
  ggthemes::theme_base() +
  facet_wrap(~genotype, ncol = 1) +
  xlab(latex2exp::TeX("pPAU_{KD} - pPAU_{CTRL}")) +
  coord_cartesian(xlim = c(-50, 50))



