# Fixtures carry noise deliberately. On a perfectly deterministic panel the
# solver finds an EXACT single-donor match, which makes the pre-period RMSPE
# zero (so the RMSPE ratio is undefined) and leaves one positive-weight
# donor (so the bias-correction regression is skipped). Neither is a bug in
# the estimator; both make the fixture test nothing.
toy <- function(effect = 5) {
  periods <- 2010:2017
  base <- outer(c(10, 12, 14, 16, 18, 13), rep(1, length(periods)))
  trend <- outer(rep(1, 6), seq_along(periods))
  Y <- base + trend +
    withr::with_seed(42, matrix(stats::rnorm(6 * length(periods), sd = 0.4),
                                nrow = 6))
  Y[6, 7:8] <- Y[6, 7:8] + effect
  dimnames(Y) <- list(c(paste0("d", 1:5), "treated"), periods)
  Y
}
