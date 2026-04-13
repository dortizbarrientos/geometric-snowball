###############################################################################
# 0. Libraries
###############################################################################
library(ggplot2)
library(viridis)
library(cowplot)
library(latex2exp)  # for LaTeX-style labels (TeX(...))

###############################################################################
# 1. Single-locus Environmental Model
###############################################################################
environment_model <- function(s_a, s_b, phi) {
  if (s_a <= 0 || s_b <= 0) return(1)
  numerator   <- s_a * s_b * (1 + phi^2)
  denominator <- (s_a + phi*s_b) * (phi*s_a + s_b)
  if (denominator == 0) return(1)
  
  val <- numerator / denominator
  pmin(1, pmax(0, val))
}

make_pdme_data <- function(ratios = c(1,2,4,8), n_points = 200) {
  phi_vals <- seq(0, 1, length.out = n_points)
  df <- expand.grid(ratio = ratios, phi = phi_vals)
  df$pdme <- mapply(
    function(r, ph) environment_model(r, 1, ph),
    df$ratio, df$phi
  )
  df
}

###############################################################################
# 2. Multi-locus Polygenic Model
###############################################################################
pdmp_polygenic <- function(n, ratio, phi) {
  if (n %% 2 != 0) stop("n must be even for this derivation.")
  half <- n / 2
  
  # First half-loci
  s1_1 <- rep(ratio,       half)
  s2_1 <- rep(ratio * phi, half)
  
  # Second half-loci
  s1_2 <- rep(ratio * phi, half)
  s2_2 <- rep(ratio,       half)
  
  s1 <- c(s1_1, s1_2)
  s2 <- c(s2_1, s2_2)
  
  sigma_1 <- sum(s1)
  sigma_2 <- sum(s2)
  
  same_allele_sum <- sum((s1 / sigma_1) * (s2 / sigma_2))
  1 - same_allele_sum
}

make_pdmp_data_phi <- function(n_vals = c(2,4,6,8,10),
                               phi_vals = c(0,0.2,0.5,0.8,1)) {
  df <- expand.grid(n = n_vals, phi = phi_vals)
  # ratio=1 => symmetrical selection
  df$pdmp <- mapply(
    function(nn, ph) pdmp_polygenic(nn, 1, ph),
    df$n, df$phi
  )
  df
}

###############################################################################
# 3. Panel (A): P_{DM,E} vs. phi, multiple selection ratios
###############################################################################
dfA <- make_pdme_data(ratios = c(1,2,4,8), n_points = 200)
pA <- ggplot(dfA, aes(x=phi, y=pdme, colour=factor(ratio))) +
  geom_line(size=1) +
  scale_colour_viridis_d(
    name = TeX("$s_{A}/s_{B}$"),
    end  = 0.9,
    option = "D"
  ) +
  labs(
    x = TeX("$varphi$ (Environmental similarity)"),
    y = TeX("$P_{DM,E}$")
  ) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  ggtitle("A)") +
  theme_minimal(base_size=13, base_family="Helvetica") +
  theme(
    text            = element_text(face="plain"), 
    panel.border    = element_rect(colour="black", fill=NA, size=0.5),
    legend.position = "bottom"
  )

###############################################################################
# 4. Panel (B): P_{DM,P} vs. n, multiple phi lines (ratio=1)
###############################################################################
dfB <- make_pdmp_data_phi(
  n_vals   = c(2,4,6,8,10),
  phi_vals = c(0,0.2,0.5,0.8,1)
)

pB <- ggplot(dfB, aes(x=n, y=pdmp, colour=factor(phi))) +
  geom_line(size=1) +
  geom_point(size=2) +
  scale_colour_viridis_d(
    name   = TeX("$varphi$"),
    end    = 0.9,
    option = "C"
  ) +
  labs(
    x = TeX("$n$ (Number of loci)"),
    y = TeX("$P_{DM,P}$")
  ) +
  coord_cartesian(ylim=c(0,1)) +
  ggtitle("B)") +
  theme_minimal(base_size=13, base_family="Helvetica") +
  theme(
    text            = element_text(face="plain"),
    panel.border    = element_rect(colour="black", fill=NA, size=0.5),
    legend.position = "bottom"
  )

###############################################################################
# 5. Combine Panels
###############################################################################
final_plot <- cowplot::plot_grid(
  pA, pB,
  ncol       = 2,
  rel_widths = c(1,1),
  align      = "hv"
)

print(final_plot)

# Optionally save:
ggsave("PDME_PanelA__PDMP_PanelB.pdf", final_plot, width=10, height=5)
