config <- list(
  traits = c("INC", "SEV", "DON", "FDK"),
  seed = 123L,
  folds = 5L,
  repetitions = 100L,
  maf = 0.05,
  heterozygosity = 0.10,
  output = "results/pipeline",
  unreplicated_2025 = c("Q Qual AYN", "Topcross", "DH"),
  historical = "data/FHB_Project_Training_Data.csv",
  testing_2025 = "data/FHB_Project_Testing_Data_2025.csv",
  locations = c(
    UNL = "data/FHB_Project_Testing_Data_UNL_2026.csv",
    UIUC = "data/FHB_Project_Testing_Data_UIUC_2026.csv",
    SDSU = "data/FHB_Project_Testing_Data_SDSU_2026.csv"
  ),
  markers = "data/fhb_analysis_2026_production_final.vcf.gz"
)
