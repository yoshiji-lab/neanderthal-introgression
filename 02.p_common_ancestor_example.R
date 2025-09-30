.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))
library(dplyr)

# Define the function
p_common_ancestor <- function(length_bp,
                              recomb_cM_Mb,
                              t_split_yrs = 550000,
                              archaic_sample_age_yrs = 50000,
                              gen_time_yrs = 29) {
  # Recombination rate: Morgans per bp per generation
  r <- (recomb_cM_Mb / 100) / 1e6
  
  # Total branch length in generations
  t_tot_gen <- (t_split_yrs + (t_split_yrs - archaic_sample_age_yrs)) / gen_time_yrs
  
  # Expected tract length
  L_exp <- 1 / (r * t_tot_gen)
  
  # Scaled length
  x <- length_bp / L_exp
  
  # Gamma survival function (k = 2)
  return(exp(-x) * (1 + x))
}


# ---------- toy replication of Zeberg & Pääbo 2020 ----------
toy <- data.frame(
  region_id = c("core 49.4 kb", "full 333.8 kb"),
  length_bp = c(49400, 333800),
  recomb_cM_Mb = c(0.53, 0.53)
)

toy <- toy %>%
  rowwise() %>%
  mutate(P_common_ancestor = p_common_ancestor(length_bp, recomb_cM_Mb)) %>%
  ungroup()

print(toy)


# Test yourself ----
toy2 <- data.frame(
  region_id = c("core 49.4 kb", "full 333.8 kb"),
  length_bp = c(49400, 333800),
  recomb_cM_Mb = c(0.526, 0.526)
)

toy2 <- toy2 %>%
  rowwise() %>%
  mutate(P_common_ancestor = p_common_ancestor(length_bp, recomb_cM_Mb)) %>%
  ungroup()

print(toy2)

# Optional: test values like linspace in Python
# values <- seq(0.525, 0.527, length.out = 11)
# for (i in values) {
#   cat(i, "\n")
#   cat(p_common_ancestor(49400, i), "\n")
#   cat(p_common_ancestor(333800, i), "\n")
# }


p_common_ancestor(49400, 0.591216)
p_common_ancestor(333800, 0.591216)


p_common_ancestor(333800, 3.804968)
