#### PACKAGES ####
library("ggplot2")
library("psych")
library("corrplot")
library("vegan")
library("mvabund")
library("car")
library("eulerr")
library("reshape2")
library("DHARMa")  
library("iNEXT")
library("lme4")
library("lmerTest")
library("MuMIn")
library("effects")
library("emmeans")

#################


#### FUNCTIONS ####
stacked.bars <- function(taxon_level, var_list, number, input_community, input_tax, input_variable, tax_disp, X_levels = NULL, trans, bars = T) {
  library("reshape2")  
  mycol_red2  <- c("white","#fffae5","#ffe158","#ffd412","#ffc616","#ffbd18","#ffb41a","#ffa61e","#ff9415","#ff860f","#ff7909","#ff6600","#f36206","#ea600a","#e25d0e","#d65a14","#d65112","#d5410e","#d52d0a","#d51b06","#d40902","#d40000","#c50000","#920000","#6b0000","#3c0000","#120000","#000000")
  # taxon_level = taxonomic level at which to make the plot (e.g., order, family, genus, otu)
  # var_list = factors by which to cluster (e.g., treatment1, treatment 2)
  # input_community = the community matrix (otus in cols, samples in rows)
  # input_variable = a data.frame containing the variables in columns and samples in rows
  # input_tax = the taxonomy matrix, with [, c(kingdom:otu)]
  # number = number of taxa to display
  # trans = transformation to scale the gradient
  for(i in 1:ncol(input_variable)) droplevels(as.factor(input_variable[,i])) -> input_variable[,i]
  
  # express counts in proportions (relative to sample)
  comp_df <- input_community/rowSums(input_community)
  row.names(comp_df) <- row.names(input_community)
  
  # bring in variables
  var_comp_df <- data.frame(merge(input_variable,comp_df,by="row.names"),row.names = 1)
  
  # average (aggregate) by variables of interest
  if(length(var_list) == 1){
    agg_vars_df <- aggregate(var_comp_df[,-c(1:ncol(input_variable))],
                             by=list(var_comp_df[,as.factor(var_list)]),
                             data=var_comp_df,
                             FUN="mean")
    
    # simplify the variable info to a single column and transpose for next step
    simpl_df <- agg_vars_df[, -c(1:length(var_list))]
    row.names(simpl_df) <- agg_vars_df$var_list
  } else {
    agg_vars_df <- aggregate(var_comp_df[, -c(1:ncol(input_variable))],
                             by=var_comp_df[, c(var_list)],
                             data=var_comp_df,
                             FUN = "mean")
    
    # simplify the variable info to a single column and transpose for next step
    simpl_df <- agg_vars_df[,-c(1:length(var_list))]
    row.names(simpl_df) <- apply(agg_vars_df[, var_list], MARGIN = 1, paste, collapse = "_" )
  }
  
  t_simpl_df <- t(simpl_df)
  
  # bring in the taxonomy info, organize cols and name the otu col
  if(isTRUE(which(colnames(input_tax) == "metabolism") == 2)) {
    lowest_taxon_level <- 2
  } else 
  {
    if(isTRUE(which(colnames(input_tax) == "species") > 0)) {
      lowest_taxon_level <- which(colnames(input_tax) == "species") + 1
    } else
    { 
      lowest_taxon_level <- which(colnames(input_tax) == "otu") + 1
    }}
  
  tax_comp <- merge(input_tax, t_simpl_df, by = "row.names")[, c(2:lowest_taxon_level, (lowest_taxon_level + 1):(ncol(input_tax) + ncol(t_simpl_df) + 1))]
  
  #names(tax_comp)[names(tax_comp) == "Row.names"] <- "otu"
  
  # get col number of taxlevel
  tax <- which(colnames(tax_comp) == taxon_level)
  
  # collapse (aggregate) to desired taxonomy level 
  if(tax == 1) 
  {
    agg_vars_df <- aggregate(tax_comp[, -c(1:(ncol(input_tax) + 1))],
                             by = list(tax_comp[, 1:tax]),
                             FUN ="sum")
  } else {
    agg_vars_df <- aggregate(tax_comp[, -c(1:(ncol(input_tax) + 1))],
                             by = tax_comp[, 1:tax],
                             FUN = "sum")
  }
  

  
  # add taxon overall prop abundance (after averaging over groups) and abundance rank
  abr_vars_df <- data.frame(agg_vars_df[, 1:tax],
                            "abundance" = rowMeans(agg_vars_df[,-(1:tax)]),
                            "rank" = rank(1/rowMeans(agg_vars_df[,-(1:tax)]),ties.method = "random"),
                            agg_vars_df[,-(1:tax)])
  
  # reduce the data to the desired number of taxa to display
  reduced_df <- abr_vars_df[abr_vars_df$rank <= number, ]
  
  # stacking the data.frame to get taxon abundances all in one column (abundances are group)
  stacked_df <- melt(reduced_df,id.vars = 1:(tax+2))
  
  # organizing plot axes
  if(tax == 1) {
    stacked_df$Y <- stacked_df[, 1]
  } else 
  {
    stacked_df$Y <- apply(stacked_df[, tax_disp], MARGIN = 1, paste, collapse = "_" )
  }
  
  Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(kingdom, phylum, class, order, family, genus, otu, rank)), ]$Y))}, error = function(e) {
    Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(kingdom, phylum, class, order, family, genus, rank)), ]$Y))}, error = function(e) {
      Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(kingdom, phylum, class, order, family, rank)), ]$Y))}, error = function(e) {
        Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(kingdom, phylum, class, order, rank)), ]$Y))}, error = function(e) {
          Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(kingdom, phylum, class, rank)), ]$Y))}, error = function(e) {
            Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(kingdom, phylum, rank)), ]$Y))}, error = function(e) {
              Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(kingdom, rank)), ]$Y))}, error = function(e) {
                Y_levels <- tryCatch({rev(unique(stacked_df[with(stacked_df, order(rank)), ]$Y))})})})})})})})})
  
  
  if(is.null(X_levels)==T) {
    unique(stacked_df$variable) -> X_levels
  } else
  {
    X_levels -> X_levels}
  
  if(bars == T){
    ggplot(data = stacked_df, aes(fill = factor(Y,       levels = Y_levels),
                                  y = value,
                                  x = factor(variable,levels = X_levels))) +
      geom_bar(position="stack", stat="identity",width = 1) +
      #scale_fill_brewer(palette = "Paired") +
      scale_x_discrete(position = "top") +
      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.title = element_blank(),
            title = element_blank(),
            legend.title = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(size = 7),
            panel.spacing.y = unit(c(0,0),"mm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "darkgrey"),
            panel.border = element_rect(fill = "transparent",colour = "black", linewidth = 1),    
            axis.ticks.length.y = unit(1, "mm"),
            axis.ticks.y = element_line(colour="#333333",linewidth = .1),
            axis.ticks.length.x = unit(0, "mm"),
            plot.margin = unit(c(1, 1, 1, 1), "mm"))
  } else {
    ggplot(data = stacked_df, aes(x = factor(variable,levels = X_levels),
                                  y = factor(Y,       levels = Y_levels))) +
      geom_tile(colour="white", size=.1,aes(fill = value)) +
      scale_fill_gradientn(colors = mycol_red2, limits = c(min(stacked_df$value), max(stacked_df$value)),trans = trans, expand = 0) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.title = element_blank(),
            title = element_blank(),
            legend.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = "transparent",colour = "black",linewidth = 1),    
            axis.ticks.length = unit(1, "mm"),
            axis.ticks = element_line(colour="#333333", linewidth = .1),
            axis.ticks.length.x = unit(0, "mm"),
            plot.margin = unit(c(1,1,1,1), "mm"))
  }
}
occupancy <- function(community, taxonomy) {
  com <- community[, rownames(taxonomy)]
  com[com != 0] <- 1
  occupancy <- colSums(com)/nrow(com)
}
abundance <- function(community, taxonomy, proportional = T) {
  com <- community[, rownames(taxonomy)]
  if(proportional == F){
    abundance <- colSums(com)}
  else
    abundance <- colSums(com)/sum(colSums(com))
}
rarefaction <- function(x, sample, replicate, round = F) {
  library("progress")
  library("vegan")
  pb <- progress_bar$new(total = replicate)
  rar_list <- list()
  for (i in 1:replicate) {
    rar_list[[i]] <- rrarefy(x, sample)
    pb$tick()
    Sys.sleep(1 / 100)
  }
  rar_table <- Reduce(`+`, rar_list)/replicate
  
  if(round == T){
    rar_table <- round(rar_table)
    rar_table <- data.frame(rar_table[, colSums(rar_table) > 0])
    rm(rar_list)
    invisible(rar_table)
    } else {
      rar_table <- data.frame(rar_table[, colSums(rar_table) > 0])
      rm(rar_list)
      invisible(rar_table)
    }
  }
gimper      <- function(model, vars = NULL, ord = "occupancy", summary = NULL, p.adjust = "none", ci = .95){
  
  if(is.null(vars) == T) vars <- c(colnames(data.frame(t(model$coefficients)[, -1])))
  ci_level <- ifelse(ci == .99, 2.58, 1.96)
  vars <- gsub("\\(|)" , ".", vars)
  coefs <- data.frame(t(model$coefficients))[, vars]
  abundance <- colSums(model$y)
  occ <- model$y; occ[occ != 0] <- 1
  occupancy <- colSums(occ)/nrow(occ)
  stdrs <- data.frame(t(model$stderr.coefficients))[, vars]
  lowers<- coefs - ci_level * stdrs;dim(lowers)
  uppers<- coefs + ci_level * stdrs;dim(uppers)
  
  sumr <- list()
  for(i in c(vars)) {
    sumr[[i]] <- data.frame(row.names = rownames(coefs),
                            "coefs"  = coefs[, i],
                            "abundance" = abundance,
                            "occupancy" = occupancy,
                            "lowers" = lowers[, i],
                            "uppers" = uppers[, i])
  }
  
  if(is.null(summary) == F){
    for(i in c(vars)) {
      sumr[[i]][, "pval"] <- p.adjust(data.frame(summary$uni.p)[, i], method = p.adjust)
    }}
  
  sumr_negative <- list()
  sumr_positive <- list()
  sumr_total    <- list()
  
  if(is.null(summary)) {
    for(i in c(vars)) {
      sumr_negative[[i]] <- sumr[[i]][sumr[[i]]$coefs < 0 & sumr[[i]]$lowers < 0 & sumr[[i]]$uppers < 0, ]
      sumr_positive[[i]] <- sumr[[i]][sumr[[i]]$coefs > 0 & sumr[[i]]$lowers > 0 & sumr[[i]]$uppers > 0, ]
      if(ord == "coefs"){
        sumr_total[[i]] <- rbind(sumr_positive[[i]][order(sumr_positive[[i]]$coefs, decreasing = T), ],
                                 sumr_negative[[i]][order(sumr_negative[[i]]$coefs, decreasing = T), ])
      } else {
        sumr_total[[i]] <- rbind(sumr_positive[[i]][order(sumr_positive[[i]]$occupancy, decreasing = T), ],
                                 sumr_negative[[i]][order(sumr_negative[[i]]$occupancy, decreasing = T), ])}}} 
  else {
    for(i in c(vars)) {
      sumr_negative[[i]] <- sumr[[i]][sumr[[i]]$coefs < 0 & sumr[[i]]$lowers < 0 & sumr[[i]]$uppers < 0 & sumr[[i]]$pval < 0.05, ]
      sumr_positive[[i]] <- sumr[[i]][sumr[[i]]$coefs > 0 & sumr[[i]]$lowers > 0 & sumr[[i]]$uppers > 0 & sumr[[i]]$pval < 0.05, ]
      if(ord == "coefs"){
        sumr_total[[i]] <- rbind(sumr_positive[[i]][order(sumr_positive[[i]]$coefs, decreasing = T), ],
                                 sumr_negative[[i]][order(sumr_negative[[i]]$coefs, decreasing = T), ])
      } else {
        sumr_total[[i]] <- rbind(sumr_positive[[i]][order(sumr_positive[[i]]$occupancy, decreasing = T), ],
                                 sumr_negative[[i]][order(sumr_negative[[i]]$occupancy, decreasing = T), ])}}}
  
  overview <- data.frame(row.names = vars)
  for(i in c(vars)) {
    overview[i, "neg_responses"] <- nrow(sumr_negative[[i]])
    overview[i, "pos_responses"] <- nrow(sumr_positive[[i]])
    overview[i, "total"] <- nrow(sumr_total[[i]])
  }
  final <- list("tables" = sumr_total,
                "summary" = overview[order(overview$total, decreasing = T), ])
  invisible(final)
} 
###################


#### DATA PREPARATION ####
fish_var <- read.csv("C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/fish_var.csv", row.names = 1, header = T, stringsAsFactors = T);dim(fish_var)
fish_var$sample_id <- rownames(fish_var)
fish_otu <- read.csv("C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/fish_otu.csv", row.names = 1, header = T, stringsAsFactors = T);dim(fish_otu)
fish_otu <- fish_otu[rownames(fish_var), order(colSums(fish_otu), decreasing = T)];dim(fish_otu)
fish_otu <- fish_otu[, colSums(fish_otu) > 0];dim(fish_otu)
all.equal(rownames(fish_otu), rownames(fish_var))

fish_tax <- read.csv("C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/fish_tax.csv", header = T, row.names = 1, stringsAsFactors = T)[, c(1, 3, 5, 7, 9, 11)];dim(fish_tax)
fish_tax <- fish_tax[colnames(fish_otu), ];dim(fish_tax)
fish_tax$otu <- rownames(fish_tax)

#rarefaction to equalize seq_depth
#fish_rar <- rarefaction(fish_otu, sample = 3500, replicate = 100, round = T);dim(fish_rar)
#fish_rar <- fish_rar[, order(colSums(fish_rar), decreasing = T)];dim(fish_rar)
#save(fish_rar, file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_rar.Rdata")
load(file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_rar.Rdata")

fish_rar_tax <- fish_tax[colnames(fish_rar), ]
fish_rar_tax$abundance <- abundance(fish_rar, fish_rar_tax)
fish_rar_tax$occupancy <- occupancy(fish_rar, fish_rar_tax)

fish_var$ESN <- estimateD(t(fish_rar), q = 2, nboot = 0)$qD
fish_var$rich   <- specnumber(fish_rar)

#########################


#### DATA EXAMINATION ####

#seq_depth
fish_var$log_age <- log(fish_var$age_years)
fish_var$log_weight <- log(fish_var$weight_g)
fish_var$log2_grain <- log2(fish_var$grain_um)
fish_cor <- corr.test(fish_var[, c("condition_factor", "log_age", "log_weight", "length_cm", "log2_grain", "vms_SAR")],
                      method = "spearman",
                      adjust = "none", ci = F)

corrplot(fish_cor[[1]], p.mat = fish_cor[[4]],
         method = "color", order = "original", addCoef.col = "black",
         tl.col = "black", tl.srt = 45, tl.cex =1 ,
         number.cex = 1, addgrid.col = "white",
         outline = F, diag = F, type = "lower",
         cl.ratio = 0.1, insig = "blank", sig.level = 1)

#########################


#### STACKED BARS ####
th <- theme(axis.text.x = element_text(size = 8, angle = 90), 
            axis.text.y = element_text(size = 5),
            axis.title = element_blank(), 
            title = element_blank(),
            legend.title = element_blank(),
            panel.spacing.y = unit(c(0,0),"mm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "darkgrey"),
            panel.border = element_rect(fill = "transparent", colour = "black", linewidth = 1),
            axis.ticks.length.y = unit(1, "mm"), axis.ticks.y = element_line(colour = "#333333", linewidth = .1),
            axis.ticks.length.x = unit(0, "mm"), plot.margin = unit(c(1, 1, 1, 1), "mm"))

x_species <- paste(fish_var[order(fish_var$species, fish_var$sex, fish_var$sample_id), "species"],
                   fish_var[order(fish_var$species, fish_var$sex, fish_var$sample_id), "sex"],
                   fish_var[order(fish_var$species, fish_var$sex, fish_var$sample_id), "sample_id"],
                   sep = "_")


fish_var_bars_fam <- stacked.bars(input_variable = fish_var,
                                          input_community = fish_rar,
                                          input_tax = fish_rar_tax,
                                          taxon_level = "family",
                                          var_list = c("species", "sex", "sample_id"),
                                          number = 30,
                                          tax_disp = c("phylum", "class", "family"),
                                          X_levels = x_species,
                                          trans = "identity")

fish_var_bars_fam +
  scale_fill_manual(values = c("#b15928",
                               "#d2ebf2","#c0dce6","#aecddb","#98bccd","#85acc1","#729cb4","#5e8ca7","#4b7d9c","#386d8f","#145078",
                               "#cbe5d0","#94c1a6","#247852",
                               "#c6dfaf","#7ec06e","#33a02c",
                               "#ffe88d","#ffde5c","#ffcc00",
                               "#dccce4","#a384bf","#6a3d9a","#800080",
                               "#FDBF6F","#ff7f00",
                               "#fbcecd","#fb9a99","#e31a1c",
                               "#6f3819")) + th

#####################


#### PERMANOVA ####
#fish_permanova_bray <- adonis2(fish_rar ~ species * (sex + condition_factor + log(age_years) + log(weight_g) + log2(grain_um) + vms_SAR),
#                               data = fish_var, permutations = 9999, method = "bray", by = "terms")
#save(fish_permanova_bray, file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_permanova_bray.Rdata")
load(file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_permanova_bray.Rdata")
fish_permanova_bray

##################


#### DIFFERENTIAL ABUNDANCE ANALYSIS ####

#reduce dataset to most OTUs in at least 1% of the samples, and of 0.1% abundance
fish_red <- fish_rar[, fish_rar_tax$occupancy > 0.01 & fish_rar_tax$abundance > 0.001];dim(fish_red)
fish_red_tax <- fish_rar_tax[colnames(fish_red), ]

#relevel species factor for manyglm for each species
fish_var$species_lim <- relevel(fish_var$species, ref = "Limanda_limanda")
fish_var$species_ple <- relevel(fish_var$species, ref = "Pleuronectes_platessa")

#fish_mglm_bug <- manyglm(mvabund(fish_red) ~ species + log(age_years) + condition_factor + log2(grain_um) + vms_SAR, data = fish_var, theta.method = "ML")
#save(fish_mglm_bug, file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_mglm_bug.Rdata")
load(file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_mglm_bug.Rdata")

#fish_mglm_lim <- manyglm(mvabund(fish_red) ~ species_lim + log(age_years) + condition_factor + log2(grain_um) + vms_SAR, data = fish_var, theta.method = "ML")
#save(fish_mglm_lim, file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_mglm_lim.Rdata")
load(file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_mglm_lim.Rdata")

#fish_mglm_ple <- manyglm(mvabund(fish_red) ~ species_ple + log(age_years) + condition_factor + log2(grain_um) + vms_SAR, data = fish_var, theta.method = "ML")
#save(fish_mglm_ple, file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_mglm_ple.Rdata")
load(file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/fish_mglm_ple.Rdata")

fish_gimp_bug <- gimper(fish_mglm_bug, ci = .95)
fish_gimp_lim <- gimper(fish_mglm_lim, ci = .95)
fish_gimp_ple <- gimper(fish_mglm_ple, ci = .95)
  
fish_gimp <- data.frame("otu" =  c(colnames(fish_mglm_bug$y), colnames(fish_mglm_lim$y), colnames(fish_mglm_ple$y)), 
                        "host" = c(rep("bug", ncol(fish_mglm_bug$y)), rep("lim", ncol(fish_mglm_lim$y)), rep("ple", ncol(fish_mglm_ple$y))))
fish_gimp[, c(colnames(fish_red_tax))] <- fish_red_tax[fish_gimp$otu, ]
fish_gimp$pooled_coef <- -c(
  (fish_mglm_bug$coefficients[2, ] + fish_mglm_bug$coefficients[3, ])/2,
  (fish_mglm_lim$coefficients[2, ] + fish_mglm_lim$coefficients[3, ])/2,
  (fish_mglm_ple$coefficients[2, ] + fish_mglm_ple$coefficients[3, ])/2)
fish_gimp$pooled_stdr <- c(
  sqrt((fish_mglm_bug$stderr.coefficients[2, ]^2 + fish_mglm_bug$stderr.coefficients[3, ]^2)/2),
  sqrt((fish_mglm_lim$stderr.coefficients[2, ]^2 + fish_mglm_lim$stderr.coefficients[3, ]^2)/2),
  sqrt((fish_mglm_ple$stderr.coefficients[2, ]^2 + fish_mglm_ple$stderr.coefficients[3, ]^2)/2))
fish_gimp$pooled_lwr <- fish_gimp$pooled_coef - 1.96 * fish_gimp$pooled_stdr
fish_gimp$pooled_upr <- fish_gimp$pooled_coef + 1.96 * fish_gimp$pooled_stdr
#fish_gimp$pooled_p <- p.adjust(pt(fish_gimp$pooled_coef/fish_gimp$pooled_stdr, df = nrow(fish_var), lower.tail = T), method = "none")

fish_gimp$species <- as.factor(with(fish_gimp,
                                    ifelse(pooled_lwr > 0, host, 
                                    ifelse(pooled_upr < 0, paste("no", host, sep = "_"),  NA))))

fish_gimp_summary <- fish_gimp_bug$summary[c("vms_SAR", "log2.grain_um.", "condition_factor", "log.age_years."), ] 
fish_gimp_summary["species", ] <- c("", "", sum(summary(fish_gimp$species)[1:6]))
fish_gimp_summary["---lim", ]  <- c(summary(fish_gimp$species)["no_lim"], summary(fish_gimp$species)["lim"], sum(summary(fish_gimp$species)["no_lim"], summary(fish_gimp$species)["lim"]))
fish_gimp_summary["---ple", ]  <- c(summary(fish_gimp$species)["no_ple"], summary(fish_gimp$species)["ple"], sum(summary(fish_gimp$species)["no_ple"], summary(fish_gimp$species)["ple"]))
fish_gimp_summary["---bug", ]  <- c(summary(fish_gimp$species)["no_bug"], summary(fish_gimp$species)["bug"], sum(summary(fish_gimp$species)["no_bug"], summary(fish_gimp$species)["bug"]))

#final table showing differentially abundant and species specific species
fish_gimp_summary
########################################


#### FOREST PLOTS ####
fish_gimp_grain <- fish_gimp_bug$tables$log2.grain[order(fish_gimp_bug$tables$log2.grain$occupancy, decreasing = T), ]
fish_gimp_grain$taxon <- paste(rownames(fish_gimp_grain), fish_red_tax[rownames(fish_gimp_grain), "genus"])
gg_fish_grain <-  ggplot(fish_gimp_grain[order(fish_gimp_grain$abundance, decreasing = T), ][1:10, ],
                         aes(x = reorder(taxon, -coefs), y = coefs)) + 
  geom_errorbar(aes(ymin = lowers, ymax = uppers), col = "black", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefs), size = 2) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour = "black", linewidth = 1)) +
  ggtitle("fold change with with grainsize  + 95%CIs");gg_fish_grain

fish_gimp_vms_SAR <- fish_gimp_bug$tables$vms_SAR[order(fish_gimp_bug$tables$vms_SAR$occupancy, decreasing = T), ]
fish_gimp_vms_SAR$taxon <- paste(rownames(fish_gimp_vms_SAR), fish_red_tax[rownames(fish_gimp_vms_SAR), "genus"])
gg_fish_vms_SAR <-  ggplot(fish_gimp_vms_SAR[order(fish_gimp_vms_SAR$abundance, decreasing = T), ][1:10, ],
                           aes(x = reorder(taxon, -coefs), y = coefs)) + 
  geom_errorbar(aes(ymin = lowers, ymax = uppers), col = "black", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefs), size = 2) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour = "black", linewidth = 1)) +
  ggtitle("fold change with vms_SAR + 95%CIs");gg_fish_vms_SAR

fish_gimp_age <- fish_gimp_bug$tables$log.age_years[order(fish_gimp_bug$tables$log.age_years$occupancy, decreasing = T), ]
fish_gimp_age$taxon <- paste(rownames(fish_gimp_age), fish_red_tax[rownames(fish_gimp_age), "genus"])
gg_fish_age <-  ggplot(fish_gimp_age[order(fish_gimp_age$abundance, decreasing = T), ][1:10, ],
                       aes(x = reorder(taxon, -coefs), y = coefs)) + 
  geom_errorbar(aes(ymin = lowers, ymax = uppers), col = "#54abf4", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefs), size = 2) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour = "black", linewidth = 1)) +
  ggtitle("fold change with with age in years + 95%CIs");gg_fish_age

fish_gimp_condition <- fish_gimp_bug$tables$condition_factor[order(fish_gimp_bug$tables$condition_factor$occupancy, decreasing = T), ]
fish_gimp_condition$taxon <- paste(rownames(fish_gimp_condition), fish_red_tax[rownames(fish_gimp_condition), "genus"])
gg_fish_condition <-  ggplot(fish_gimp_condition[order(fish_gimp_condition$abundance, decreasing = T), ][1:10, ],
                       aes(x = reorder(taxon, -coefs), y = coefs)) + 
  geom_errorbar(aes(ymin = lowers, ymax = uppers), col = "#0500f4", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefs), size = 2) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour = "black", linewidth = 1)) +
  ggtitle("fold change with condition factor + 95%CIs");gg_fish_condition

gg_fish_species <-  ggplot(fish_gimp[fish_gimp$occupancy > 0 & fish_gimp$species %in% c("bug", "lim", "ple"), ],
                           aes(x = reorder(paste(otu, genus), pooled_coef), y = -pooled_coef)) + 
  facet_wrap(host ~ ., nrow = 3, scales = "free_y") +
  scale_color_manual(values = c("#187500", "#FFB100", "#E97777"), guide = "none") + 
  geom_errorbar(aes(ymin = pooled_lwr, ymax = pooled_upr, col = host), alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = pooled_coef), size = 2) + 
  ylim(0, 20) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour = "black", linewidth = 1)) +
  ggtitle("species specific OTUs + 95%CIs");gg_fish_species

#####################


#### VENN DIAGRAM ####

# venn differential abundances
venn_dif_abun <- c(
  "lim"         = as.numeric(fish_gimp_summary["---lim", 2]),      # Only A
  "ple"         = as.numeric(fish_gimp_summary["---ple", 2]),      # Only B
  "bug"         = as.numeric(fish_gimp_summary["---bug", 2]),      # Only C
  "lim&ple"     = as.numeric(fish_gimp_summary["---lim", 1]),    # Intersection of A and B
  "lim&bug"     = as.numeric(fish_gimp_summary["---ple", 1]),    # Intersection of A and C
  "ple&bug"     = as.numeric(fish_gimp_summary["---bug", 1]),    # Intersection of B and C
  "lim&ple&bug" = 0   # Intersection of A, B, and C
)

plot(venn(venn_dif_abun),
     fills = c("#187500", "#FFB100", "#E97777"),
     edges = TRUE, quantities = TRUE, main = "differentially abundant OTUs")



#####################


#### DIVERSITY ANALYSIS ####

# effective species number (ESN)
lmer_fish_ESN <- lmer(ESN ~ species * (sex + log2(grain_um) + vms_SAR + log(weight_g) +
                       log(age_years) + condition_factor) + (1|station),
                     na.action = na.fail, REML = F,
                     data = fish_var);AICc(lmer_fish_ESN)
plot(simulateResiduals(lmer_fish_ESN))
#dr_lmer_fish_ESN <- dredge(lmer_fish_ESN, trace = 2)
#save(dr_lmer_fish_ESN, file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/dr_lmer_fish_ESN.Rdata")
load(file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/dr_lmer_fish_ESN.Rdata")
dr_lmer_fish_ESN
lmer_fish_ESN_best <- get.models(dr_lmer_fish_ESN, subset = 1)[[1]]
summary(lmer_fish_ESN_best);
Anova(lmer_fish_ESN_best)
r.squaredGLMM(lmer_fish_ESN_best)

# OTU richness
lmer_fish_rich <- lmer(rich ~ species * 
                       (sex + log2(grain_um) + vms_SAR + log(weight_g) +
                          log(age_years) + condition_factor) + (1|station),
                      na.action = na.fail, REML = F,
                      data = fish_var);AICc(lmer_fish_rich)
plot(simulateResiduals(lmer_fish_rich))
#dr_lmer_fish_rich <- dredge(lmer_fish_rich, trace = 2)
#save(dr_lmer_fish_rich, file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/dr_lmer_fish_rich.Rdata")
load(file = "C:/Users/Bonthond/Documents/GitHub/flatfish_microbiota/Rdata/dr_lmer_fish_rich.Rdata")
dr_lmer_fish_rich
lmer_fish_rich_best <- get.models(dr_lmer_fish_rich, subset = 1)[[1]]
summary(lmer_fish_rich_best);
Anova(lmer_fish_rich_best)
r.squaredGLMM(lmer_fish_rich_best)
plot(effect(lmer_fish_rich_best, term = "species", resid = T), smooth.residuals = F)
plot(effect(lmer_fish_rich_best, term = "log2(grain_um)", resid = T), smooth.residuals = F)
plot(effect(lmer_fish_rich_best, term = "log2(grain_um):species", resid = T), smooth.residuals = F)
plot(effect(lmer_fish_rich_best, term = "log(weight_g) ", resid = T), smooth.residuals = F)
plot(effect(lmer_fish_rich_best, term = "log(age_years) ", resid = T), smooth.residuals = F)

#posthoc test to compare levels in species
pairs(emmeans(lmer_fish_rich_best, ~ species), adjust = "holm")
emtrends(lmer_fish_rich_best, ~ species, var = "log2(grain_um)", infer = T, adjust = "holm")
############################

