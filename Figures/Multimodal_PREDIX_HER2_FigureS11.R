# FigS10b
# selected features for All
library(data.table);library(ggplotify);library(ggpubr)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/consensus_signatures.csv") 
df=df[df$Scenario=="Global",]
df$Modality=factor(df$Modality,levels = c("Clin","DNA","RNA","Prot","WSI"))
group_colors <- c(
  Clin = "#434279", 
  DNA = "#f2a104", 
  RNA = "#72a2c0",
  Prot= "#9467bd",
  WSI = "#00743f"
)
df=df[order(df$Modality,-df$mean_imp),]
df$Feature <- factor(df$Feature, levels = rev(unique(df$Feature)))

p <- ggplot(df, aes(x = mean_imp, y = Feature, fill = Modality)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, alpha = 0.9) +
  scale_fill_manual(values = group_colors) + 
  labs(x = "Mean Importance", y = "") +     
  theme_classic() + 
  theme(
    axis.text.x = element_text(color = "black"), 
    axis.text.y = element_text(color = "black"), 
    legend.position = "top", 
    legend.title = element_text(face = "bold")
  )


print(p)
Feature=df$Feature
ggsave("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS11/Imp_Plot.pdf", plot = p, width = 7, height = 10, dpi = 300)


library(IOBR);library(data.table);library(tableone)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
df$Clin_prolifvalu[df$Clin_prolifvalu=="Unknown"]=NA
df$Clin_prolifvalu=as.numeric(df$Clin_prolifvalu)
df$DNA_HRD=as.numeric(df$DNA_HRD)
Feature=setdiff(Feature,c("Clin_ER","Clin_ANYNODES","Clin_TUMSIZE","DNA_coding_mutation_TP53_oncokb","Clin_Arm"))
str(df[,Feature]) 
res=batch_wilcoxon(
  df,
  target = "pCR",
  feature =Feature
)
res=res[res$p.value<0.05,]

df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/clin_multiomics_curated_metrics_PREDIX_HER2.txt")%>%as.data.frame()
cat=c("Clin_ER","Clin_ANYNODES","Clin_TUMSIZE","DNA_coding_mutation_TP53_oncokb","Clin_Arm")
tableOne <- CreateTableOne(vars = cat, strata = c("pCR"), data = df,
                           factorVars = cat)
print(tableOne)


# FigS10c

# FigS10d AUC modality vs arm
library(data.table);library(ggplotify);library(ggpubr)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure7/consensus_performance.csv") 
df$Model <- factor(df$Model, levels = c("Clin", "DNA", "RNA", "Prot", "WSI","Fused_ElasticNet"))

# 同样，可以固定 Scenario 的颜色顺序
df$Scenario <- factor(df$Scenario, levels = c("DHP", "T-DM1","Global"))

# --- 2. 使用 ggplot2 绘图 ---
# --- 2. Plotting ---
p <- ggplot(df, aes(x = Model, y = Pooled_OOF_AUROC, color = Scenario)) +
  
  # Add a horizontal reference line at y = 0.5 (representing a random guess baseline)
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  
  # Add error bars for the 95% Confidence Intervals
  geom_errorbar(aes(ymin = CI_lower_95, ymax = CI_upper_95),
                position = position_dodge(width = 0.6), # Dodge to prevent overlapping of Scenarios
                width = 0.2,                            # Width of the error bar caps
                linewidth = 0.8) + 
  
  # Add points for the AUROC mean values
  geom_point(position = position_dodge(width = 0.6), 
             size = 3) +
  
  # Manually map the specific colors to each Scenario
  scale_color_manual(values = c("Global" = "#e5c06e", 
                                "DHP" = "#8491B4FF", 
                                "T-DM1" = "#91D1C2FF")) +
  
  # Set plot title and axis labels
  labs(title = "Model Performance across Different Modalities",
       x = "",
       y = "Pooled OOF AUROC (95% CI)",
       color = "Scenario") +
  
  # Apply theme and customize appearance details
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    
    # Bold x-axis text and rotate 45 degrees to prevent long names from overlapping
    axis.text.x = element_text(size = 11, face = "bold", angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 10),
    
    # Remove minor grid lines for a cleaner look
    panel.grid.minor = element_blank(),                    
    
    # Add faint vertical grid lines to separate modalities clearly
    panel.grid.major.x = element_line(linetype = "dotted", color = "grey80"),
    
    # Position the legend at the top
    legend.position = "top", 
    legend.title = element_text(face = "bold")
  )

print(p)
ggsave("E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS11/AUC_Modality_Plot.pdf", plot = p, width = 7, height = 5, dpi = 300)





