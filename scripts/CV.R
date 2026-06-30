cv_gblup <- function(pheno_df, GRM, k = 5, seed = 123) {
  set.seed(seed)
  
  pheno_df <- pheno_df %>%
    dplyr::filter(FullSampleName %in% rownames(GRM))
  
  folds <- sample(rep(1:k, length.out = nrow(pheno_df)))
  fold_acc <- numeric(k)
  
  for (fold in 1:k) {
    train <- pheno_df[folds != fold, ]
    test  <- pheno_df[folds == fold, ]
    
    ids_train <- train$FullSampleName
    ids_test  <- test$FullSampleName
    
    GRM_train <- GRM[ids_train, ids_train]
    GRM_test  <- GRM[ids_test, ids_train]
    
    fit <- rrBLUP::mixed.solve(
      y = train$y,
      K = GRM_train,
      X = matrix(1, nrow = nrow(train), ncol = 1)
    )
    
    gebv_test <- GRM_test %*% fit$u
    
    fold_acc[fold] <- cor(gebv_test, test$y, use = "complete.obs")
  }
  
  list(
    fold_accuracies = fold_acc,
    mean_accuracy = mean(fold_acc)
  )
}

traits <- c("INC", "SEV", "DON")

cv_all_traits <- lapply(traits, function(t) {
  message("Running CV for trait: ", t)
  
  ph <- merged_data %>%
    dplyr::filter(TRAIT == t) %>%
    dplyr::select(FullSampleName, y) %>%
    dplyr::filter(!is.na(y))
  
  res <- cv_gblup(ph, GRM_aligned, k = 5)
  
  data.frame(
    TRAIT = t,
    mean_accuracy = res$mean_accuracy,
    fold1 = res$fold_accuracies[1],
    fold2 = res$fold_accuracies[2],
    fold3 = res$fold_accuracies[3],
    fold4 = res$fold_accuracies[4],
    fold5 = res$fold_accuracies[5]
  )
})

cv_results_df <- dplyr::bind_rows(cv_all_traits)
cv_results_df

library(ggplot2)

ggplot(cv_results_df, aes(x = TRAIT, y = mean_accuracy)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Cross-Validation Accuracy by Trait",
       y = "Mean Accuracy",
       x = "Trait")

