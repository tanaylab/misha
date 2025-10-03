library(misha)

# Create a simple test genomic database
gsetroot(GROOT)

# Create a simple PSSM matrix (AT-rich motif)
pssm <- matrix(c(
  0.8, 0.1, 0.05, 0.05,  # Position 1: mostly A
  0.05, 0.05, 0.05, 0.85,  # Position 2: mostly T
  0.8, 0.1, 0.05, 0.05,  # Position 3: mostly A
  0.05, 0.05, 0.05, 0.85   # Position 4: mostly T
), nrow=4, byrow=TRUE)
colnames(pssm) <- c("A", "C", "G", "T")

# Test 1: Create PWM vtrack without spatial weighting (backward compatibility)
gvtrack.create("test_pwm_nospatial", NULL, 
               func="pwm", 
               params=list(pssm=pssm, bidirect=TRUE))

# Test 2: Create PWM vtrack with spatial weighting
# Spatial factors weight different positions - higher in the middle
spat_factors <- c(0.5, 1.0, 2.0, 1.0, 0.5)
spat_bin <- 50L

gvtrack.create("test_pwm_spatial", NULL,
               func="pwm",
               params=list(
                 pssm=pssm, 
                 bidirect=TRUE,
                 spat_factor=spat_factors,
                 spat_bin=spat_bin
               ))

# Test 3: pwm.max with spatial
gvtrack.create("test_pwm_max_spatial", NULL,
               func="pwm.max",
               params=list(
                 pssm=pssm,
                 bidirect=TRUE,
                 spat_factor=spat_factors,
                 spat_bin=spat_bin
               ))

# Test 4: pwm.max.pos with spatial
gvtrack.create("test_pwm_maxpos_spatial", NULL,
               func="pwm.max.pos",
               params=list(
                 pssm=pssm,
                 bidirect=TRUE,
                 spat_factor=spat_factors,
                 spat_bin=spat_bin
               ))

# Run a simple extraction on a small interval
test_interval <- gintervals(1, 1000000, 1000500, 1)

cat("Testing basic PWM (no spatial):\n")
result1 <- gextract("test_pwm_nospatial", intervals=test_interval, iterator=test_interval)
print(head(result1))

cat("\nTesting spatial PWM:\n")
result2 <- gextract("test_pwm_spatial", intervals=test_interval, iterator=test_interval)
print(head(result2))

cat("\nTesting pwm.max with spatial:\n")
result3 <- gextract("test_pwm_max_spatial", intervals=test_interval, iterator=test_interval)
print(head(result3))

cat("\nTesting pwm.max.pos with spatial:\n")
result4 <- gextract("test_pwm_maxpos_spatial", intervals=test_interval, iterator=test_interval)
print(head(result4))

# Clean up
gvtrack.rm(c("test_pwm_nospatial", "test_pwm_spatial", "test_pwm_max_spatial", "test_pwm_maxpos_spatial"), force=TRUE)

cat("\nAll spatial PWM tests completed successfully!\n")
