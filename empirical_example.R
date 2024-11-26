setwd("~/Documents/Ben's Stuff/0 KU/Dissertation/Simulations/Triangle_Plots_best_practices/empirical_example")

# Connect RStudio to GitHub
library(usethis)
use_git()
use_github()
usethis::use_readme_rmd()



library(knitr)
knit(input="README.rmd", output = "README.md")

#usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)
#usethis::use_git_remote("origin", url = "https://github.com/omys-omics/empirical_example", overwrite = TRUE)
#########


# Load packages#
library(triangulaR)
library(vcfR)
library(SNPfiltR)
#library(bgchm)



# Define new function
sim.hybrid.CIs <- function(vcfR = NULL, pm = NULL, p1 = NULL, p2 = NULL, reps = 1000, plot.samples = T,
                           colors = NULL, ind.labels = F, cex = 2, alpha = 1, jitter = 0, max.overlaps = 10) 
{
  if (any(is.na(pm$pop))) {
    stop("All individuals must be assigned to a population (no NAs in popmap)")
  }
  if (!(identical(sort(unique(colnames(vcfR@gt)[-1])), sort(unique(pm$id))))) {
    stop("There is at least one individual in the vcfR object that is not in the popmap, or vice versa")
  }
  m <- extract.gt(vcfR)
  m[m == "0|0"] <- 0
  m[m == "0|1"] <- 1
  m[m == "1|0"] <- 1
  m[m == "1|1"] <- 2
  m[m == "0/0"] <- 0
  m[m == "0/1"] <- 1
  m[m == "1/0"] <- 1
  m[m == "1/1"] <- 2
  p1.gts <- m[, pm[pm$pop == p1, ]$id]
  p2.gts <- m[, pm[pm$pop == p2, ]$id]
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  af_p1 <- (rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts == 
                                                                2, na.rm = TRUE)))/(2 * rowSums(!is.na(p1.gts)))
  af_p2 <- (rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts == 
                                                                2, na.rm = TRUE)))/(2 * rowSums(!is.na(p2.gts)))
  af_diff <- abs(af_p1 - af_p2)
  af <- data.frame(p1 = af_p1, p2 = af_p2)
  s <- nrow(af)
  
  
  p1.allele <- ifelse(af_p1 > af_p2, 2, 0)
  p2.allele <- ifelse(af_p2 > af_p1, 2, 0)
  g <- matrix(nrow = nrow(m), ncol = ncol(m))
  g[m == p1.allele] <- 0
  g[m == 1] <- 1
  g[m == p2.allele] <- 2
  g[is.na(m)] <- NA
  g[m == -9] <- NA
  colnames(g) <- colnames(m)
  rownames(g) <- rownames(m)
  
  af[af$p1 > 0.5, "p1"] <- 1 - af[af$p1 > 0.5, "p1"]
  af[af$p2 < 0.5, "p2"] <- 1 - af[af$p2 < 0.5, "p2"]
  
  # counts of observed 1 allele and total sampled alleles
  p1.gts <- g[, pm[pm$pop == p1, ]$id]
  p2.gts <- g[, pm[pm$pop == p2, ]$id]
  p1.gts[] <- sapply(p1.gts, as.numeric)
  p2.gts[] <- sapply(p2.gts, as.numeric)
  
  # observed for each parental pop
  obs1.p1 <- rowSums(p1.gts == 1, na.rm = TRUE) + (2 * rowSums(p1.gts ==2, na.rm = TRUE))
  obs1.p2 <- rowSums(p2.gts == 1, na.rm = TRUE) + (2 * rowSums(p2.gts ==2, na.rm = TRUE))
  #total number of allele sampled for each parental pop
  samp.p1 <- 2 * rowSums(!is.na(p1.gts))
  samp.p2 <- 2 * rowSums(!is.na(p2.gts))
  
  af <- data.frame(p1 = af$p1,
                   obs1.p1 = obs1.p1,
                   samp.p1 = samp.p1,
                   p2 = af$p2,
                   obs1.p2 = obs1.p2,
                   samp.p2 = samp.p2)
  
  # define parameters for beta posterior
  alpha.p1 <- obs1.p1 + 1
  beta.p1 <- samp.p1 - obs1.p1 + 1
  alpha.p2 <- obs1.p2 + 1
  beta.p2 <- samp.p2 - obs1.p2 + 1
  
  sim.genotypes <- matrix(ncol = 3, nrow = 0)
  for (i in 1:reps) {
    
    # simulate true allele frequencies in parental populations
    # sample from the beta distribution
    sim.freqs.p1 <- rbeta(length(alpha.p1), alpha.p1, beta.p1)
    sim.freqs.p2 <- rbeta(length(alpha.p2), alpha.p2, beta.p2)
    
    # simulate allele frequency in simulated F1 population
    sim.freqs.f1 <- (sim.freqs.p1 + sim.freqs.p2) / 2
    
    # simulate P1 alleles
    # sample from the binomial distribution
    sim.p1.allele.1 <- rbinom(length(sim.freqs.p1), size = 1, prob = sim.freqs.p1)
    sim.p1.allele.2 <- rbinom(length(sim.freqs.p1), size = 1, prob = sim.freqs.p1)
    
    # simulate P2 alleles
    # sample from the binomial distribution
    sim.p2.allele.1 <- rbinom(length(sim.freqs.p2), size = 1, prob = sim.freqs.p2)
    sim.p2.allele.2 <- rbinom(length(sim.freqs.p2), size = 1, prob = sim.freqs.p2)
    
    # simulate F1 alleles
    # sample from the binomial distribution
    sim.f1.allele.1 <- rbinom(length(sim.freqs.p1), size = 1, prob = sim.freqs.p1)
    sim.f1.allele.2 <- rbinom(length(sim.freqs.p2), size = 1, prob = sim.freqs.p2)
    
    # simulate F2 alleles
    # sample from the binomial distribution
    sim.f2.allele.1 <- rbinom(length(sim.freqs.f1), size = 1, prob = sim.freqs.f1)
    sim.f2.allele.2 <- rbinom(length(sim.freqs.f1), size = 1, prob = sim.freqs.f1)
    
    # simulate BC F1xP1 alleles
    # sample from the binomial distribution
    sim.bcp1.allele.1 <- rbinom(length(sim.freqs.p1), size = 1, prob = sim.freqs.p1)
    sim.bcp1.allele.2 <- rbinom(length(sim.freqs.f1), size = 1, prob = sim.freqs.f1)
    
    # simulate BC F1xP2 alleles
    # sample from the binomial distribution
    sim.bcp2.allele.1 <- rbinom(length(sim.freqs.p2), size = 1, prob = sim.freqs.p2)
    sim.bcp2.allele.2 <- rbinom(length(sim.freqs.f1), size = 1, prob = sim.freqs.f1)
    
    
    sim.data <- data.frame(p1 = af$p1,
                           sim.freqs.p1 = sim.freqs.p1,
                           p2 = af$p2,
                           sim.freqs.p2 = sim.freqs.p2,
                           sim.p1.allele.1 = sim.p1.allele.1,
                           sim.p1.allele.2 = sim.p1.allele.2,
                           sim.p1.genotype = sim.p1.allele.1 + sim.p1.allele.2,
                           sim.p2.allele.1 = sim.p2.allele.1,
                           sim.p2.allele.2 = sim.p2.allele.2,
                           sim.p2.genotype = sim.p2.allele.1 + sim.p2.allele.2,                           
                           sim.f1.allele.1 = sim.f1.allele.1,
                           sim.f1.allele.2 = sim.f1.allele.2,
                           sim.f1.genotype = sim.f1.allele.1 + sim.f1.allele.2,
                           sim.freqs.f1 = sim.freqs.f1,
                           sim.f2.allele.1 = sim.f2.allele.1,
                           sim.f2.allele.2 = sim.f2.allele.2,
                           sim.f2.genotype = sim.f2.allele.1 + sim.f2.allele.2,
                           sim.bcp1.allele.1 = sim.bcp1.allele.1,
                           sim.bcp1.allele.2 = sim.bcp1.allele.2,
                           sim.bcp1.genotype = sim.bcp1.allele.1 + sim.bcp1.allele.2,
                           sim.bcp2.allele.1 = sim.bcp2.allele.1,
                           sim.bcp2.allele.2 = sim.bcp2.allele.2,
                           sim.bcp2.genotype = sim.bcp2.allele.1 + sim.bcp2.allele.2)
    
    p1.hi <- sum(sim.data$sim.p1.genotype) / (2 * nrow(sim.data))
    p1.het <- sum(sim.data$sim.p1.genotype==1) / nrow(sim.data)    
    p2.hi <- sum(sim.data$sim.p2.genotype) / (2 * nrow(sim.data))
    p2.het <- sum(sim.data$sim.p2.genotype==1) / nrow(sim.data)    
    f1.hi <- sum(sim.data$sim.f1.genotype) / (2 * nrow(sim.data))
    f1.het <- sum(sim.data$sim.f1.genotype==1) / nrow(sim.data)
    f2.hi <- sum(sim.data$sim.f2.genotype) / (2 * nrow(sim.data))
    f2.het <- sum(sim.data$sim.f2.genotype==1) / nrow(sim.data)
    bcp1.hi <- sum(sim.data$sim.bcp1.genotype) / (2 * nrow(sim.data))
    bcp1.het <- sum(sim.data$sim.bcp1.genotype==1) / nrow(sim.data)
    bcp2.hi <- sum(sim.data$sim.bcp2.genotype) / (2 * nrow(sim.data))
    bcp2.het <- sum(sim.data$sim.bcp2.genotype==1) / nrow(sim.data)
    
    sim.genotypes <- rbind(sim.genotypes, c(p1.hi, p1.het, "P1"))
    sim.genotypes <- rbind(sim.genotypes, c(p2.hi, p2.het, "P2"))
    sim.genotypes <- rbind(sim.genotypes, c(f1.hi, f1.het, "F1"))
    sim.genotypes <- rbind(sim.genotypes, c(f2.hi, f2.het, "F2"))
    sim.genotypes <- rbind(sim.genotypes, c(bcp1.hi, bcp1.het, "P1xF1"))
    sim.genotypes <- rbind(sim.genotypes, c(bcp2.hi, bcp2.het, "P2xF1"))
    
  }
  
  sim.genotypes <- data.frame(sim.genotypes)
  colnames(sim.genotypes) <- c("hi", "het", "class")
  
  #  return(sim.genotypes)
  
  
  #  f1.hi <- quantile(sim.genotypes[,1], probs = c(0.025, 0.975))
  #  f1.het <- quantile(sim.genotypes[,2], probs = c(0.025, 0.975))
  #  
  #  CIs <- data.frame(hi.lower = f1.hi[1],
  #                    hi.upper = f1.hi[2],
  #                    het.lower = f1.het[1],
  #                    het.upper = f1.het[2])
  #  rownames(CIs) <- "F1"
  
  
  # Define white polygons outside the triangle
  masking_polygons <- list(
    data.frame(  # Left triangle
      x = c(0, 0.525, 0, 0),
      y = c(0, 1.05, 1.05, 0)
    ),
    data.frame(  # Right triangle
      x = c(1, 1, 0.475, 1),
      y = c(0, 1.05, 1.05, 0)
    ),
    data.frame(  # Left rectangle
      x = c(-0.05, 0, 0, -0.05, -0.05),
      y = c(0, 0, 0.5, 0.5, 0)
    ),
    data.frame(  # Right rectangle
      x = c(1, 1.05, 1.05, 1, 1),
      y = c(0, 0, 0.5, 0.5, 0)
    ),
    data.frame(  # Bottom rectangle
      x = c(-0.05, 1.05, 1.05, -0.05, -0.05),
      y = c(0, 0, -0.05,-0.05, 0)
    )
  )
  
  # Define ellipse centers
  ellipse_centers <- sim.genotypes %>%
    group_by(class) %>%
    summarize(
      center_x = mean(as.numeric(hi)),
      center_y = mean(as.numeric(het)),
      .groups = "drop"
    )
  
  if (is.null(colors)) {
    color_ramp <- colorRampPalette(c("orange", "blue", "green", "red1", "yellow", "purple"))
    colors <- color_ramp(length(unique(pm$pop)))
  }
  
  if (plot.samples) {
    hi <- hybridIndex(vcfR = vcfR, pm = pm, p1 = p1, p2 = p2)
    p <- ggplot(hi, aes(x = hybrid.index, y = heterozygosity, color = as.factor(pop))) +
      stat_ellipse(data = sim.genotypes, aes(x = as.numeric(hi), y = as.numeric(het), color = NULL, group = as.factor(class)), geom = "polygon", alpha = 0.1) +
      lapply(masking_polygons, function(poly) {
        geom_polygon(data = poly, aes(x = x, y = y), fill = "white", alpha = 1, color = NA)
      }) +
      geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
      geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
      stat_function(fun = function(hi) 2 * hi * (1 - hi), xlim = c(0, 1), color = "black", linetype = "dashed") +
      geom_text(
        data = ellipse_centers,
        aes(x = center_x, y = center_y, label = class), 
        color = "black",
        size = 5,
        fontface = "bold"
      ) +
      geom_jitter(cex = cex, alpha = alpha, width = jitter, height = jitter) +
      guides(shape = guide_legend(override.aes = list(size = 5), order = 2, label.theme = element_text(face = "italic"))) +
      xlab("Hybrid Index") +
      ylab("Interclass Heterozygosity") +
      scale_color_manual("pop", values = colors) +
      ylim(c(-0.05, 1.05)) +
      xlim(c(-0.05, 1.05)) +
      theme_classic()
    
    if (ind.labels) {
      p <- p + geom_label_repel(aes(label = id), size = 2, max.overlaps = max.overlaps)
    }
  }
  
  if (!plot.samples) {
    p <- ggplot() +
      stat_ellipse(data = sim.genotypes, aes(x = as.numeric(hi), y = as.numeric(het), color = NULL, group = as.factor(class)), geom = "polygon", alpha = 0.1) +
      lapply(masking_polygons, function(poly) {
        geom_polygon(data = poly, aes(x = x, y = y), fill = "white", alpha = 1, color = NA)
      }) +
      geom_segment(aes(x = 0.5, xend = 1, y = 1, yend = 0), color = "black") +
      geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 1), color = "black") +
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), color = "black") +
      stat_function(fun = function(hi) 2 * hi * (1 - hi), xlim = c(0, 1), color = "black", linetype = "dashed") +
      geom_text(
        data = ellipse_centers,
        aes(x = center_x, y = center_y, label = class), 
        color = "black",
        size = 5,
        fontface = "bold"
      ) +
      guides(shape = guide_legend(override.aes = list(size = 5), order = 2, label.theme = element_text(face = "italic"))) +
      xlab("Hybrid Index") +
      ylab("Interclass Heterozygosity") +
      ylim(c(-0.05, 1.05)) +
      xlim(c(-0.05, 1.05)) +
      theme_classic()
  }
  return(p)
}


# Read in empirical data
fox_sparrow_90 <- read.vcfR("snps_3_90.vcf.gz")
fox_sparrow_pm <- read.table("passerella_samples_3.txt", header = T)







#################################
##   Even parental sampling    ##
#################################

fox_sparrow_aims <- alleleFreqDiff(vcfR = fox_sparrow_90, pm = fox_sparrow_pm, 
                                   p1 = "ili", p2 = "una", difference = 1)
fox_sparrow_hi_het <-hybridIndex(vcfR = fox_sparrow_aims, pm = fox_sparrow_pm, p1 = "ili", p2 = "una")
triangle.plot(fox_sparrow_hi_het)

library(dplyr)
library(ggrepel)
sim.hybrid.CIs(vcfR = fox_sparrow_aims, pm = fox_sparrow_pm, p1 = "ili", p2 = "una", reps = 1000, 
               plot.samples = T, colors = NULL, ind.labels = F, max.overlaps = 10) 










