scenarios <- data.frame(
  scenario_id = sprintf("R%02d", 1:12),
  scenario_name = c(
    "central_case2",
    "low_clock_weak_jump",
    "high_clock_strong_jump",
    "sigma2_hard_B07",
    "alpha_hard_H09",
    "low_alpha_high_clock_B04",
    "high_alpha_high_sigma2_H02",
    "high_beta_high_rho_B08",
    "high_beta_low_rho_H07",
    "interior_high_beta_H03",
    "low_beta_high_rho",
    "minimum_alpha_strong_jump"
  ),
  source = c(
    "locked_case2",
    "locked_low_clock",
    "locked_high_clock",
    "locked_B07",
    "locked_H09",
    "locked_B04",
    "locked_H02",
    "locked_B08",
    "locked_H07",
    "locked_H03",
    "existing_bounds_interaction",
    "existing_bounds_interaction"
  ),
  alpha = c(
    1.25, 1.45, 1.55, 1.15, 1.40, 1.10,
    1.58, 1.55, 1.55, 1.32, 1.25, 1.08
  ),
  sigma2 = c(
    1.00, 0.60, 1.40, 1.30, 0.70, 0.80,
    1.25, 0.90, 1.00, 0.95, 1.10, 1.30
  ),
  sigma1 = c(
    1.20, 1.60, 0.70, 1.50, 1.30, 0.80,
    0.85, 0.90, 1.60, 1.20, 0.70, 1.00
  ),
  rho = c(
    1.45, 0.90, 2.20, 1.10, 1.00, 1.80,
    1.05, 1.90, 1.20, 1.50, 1.90, 1.70
  ),
  anchor_class = c(
    "central",
    "low_clock;weak_jump;extreme_ecf_stress",
    "high_clock;strong_jump;propagation_stress",
    "sigma2_hard;low_alpha;strong_jump",
    "alpha_hard",
    "low_alpha;high_clock",
    "high_alpha;high_sigma2;low_beta;low_rho",
    "high_beta;high_rho",
    "high_beta;low_rho",
    "interior;high_beta",
    "low_beta;high_rho",
    "minimum_alpha;strong_jump;high_clock"
  ),
  long_span_subset = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
  ),
  stress_subset = c(
    FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
  ),
  stringsAsFactors = FALSE
)
scenarios$beta <- scenarios$sigma1 * scenarios$rho
scenarios$jump_strength <- cut(
  scenarios$sigma2,
  breaks = c(-Inf, 0.80, 1.15, Inf),
  labels = c("weak", "moderate", "strong")
)
scenarios$clock_speed <- cut(
  scenarios$rho,
  breaks = c(-Inf, 1.15, 1.75, Inf),
  labels = c("low", "moderate", "high")
)
scenarios$alpha_range <- cut(
  scenarios$alpha,
  breaks = c(-Inf, 1.20, 1.45, Inf),
  labels = c("low", "central", "high")
)

profile_path <- function(n_steps, terminal) {
  if (n_steps == 10000000L && terminal == 50) {
    return("finite_window_correction_numerical_profiles.csv")
  }
  file.path(
    "overarching_profiles",
    sprintf(
      "finite_window_profile_n%d_T%s.csv",
      as.integer(n_steps),
      format(terminal, scientific = FALSE, trim = TRUE)
    )
  )
}

rows <- list()
row_index <- 0L
add_run <- function(
    scenario_index, panel, experiment, stage, replication,
    n_steps, terminal, seed, priority) {
  row_index <<- row_index + 1L
  scenario <- scenarios[scenario_index, ]
  run_key <- paste(
    experiment,
    scenario$scenario_id,
    sprintf("rep%02d", replication),
    paste0("n", as.integer(n_steps)),
    paste0("T", format(terminal, scientific = FALSE, trim = TRUE)),
    sep = "__"
  )
  rows[[row_index]] <<- data.frame(
    run_key = run_key,
    case = 10000L + row_index,
    panel = panel,
    experiment = experiment,
    stage = stage,
    scenario_id = scenario$scenario_id,
    replication = replication,
    seed = as.integer(seed),
    n_steps = as.integer(n_steps),
    terminal = as.numeric(terminal),
    delta = terminal / n_steps,
    correction_profile = profile_path(n_steps, terminal),
    priority = as.integer(priority),
    frozen_status = "planned",
    stringsAsFactors = FALSE
  )
}

# Experiment A: eight development paths and one untouched path per regime.
for (scenario_index in seq_len(nrow(scenarios))) {
  for (replication in 1:8) {
    add_run(
      scenario_index = scenario_index,
      panel = "development",
      experiment = "production",
      stage = if (replication <= 3L) {
        "feasibility"
      } else {
        "primary_expansion"
      },
      replication = replication,
      n_steps = 10000000L,
      terminal = 50,
      seed = 110000000L + scenario_index * 1000L + replication,
      priority = if (replication <= 3L) 1L else 2L
    )
  }
  add_run(
    scenario_index = scenario_index,
    panel = "untouched",
    experiment = "production",
    stage = "untouched_final",
    replication = 9L,
    n_steps = 10000000L,
    terminal = 50,
    seed = 120000000L + scenario_index * 1000L + 1L,
    priority = 6L
  )
}

# Experiment B: lower-n characterization. The n=10M arm reuses production
# development replications 1:3 in analysis rather than duplicating paths.
infill_cells <- list(
  list(n = 1000000L, code = 210000000L),
  list(n = 4000000L, code = 220000000L)
)
for (cell in infill_cells) {
  for (scenario_index in seq_len(nrow(scenarios))) {
    for (replication in 1:3) {
      add_run(
        scenario_index = scenario_index,
        panel = "development",
        experiment = "infill",
        stage = "feasibility_characterization",
        replication = replication,
        n_steps = cell$n,
        terminal = 50,
        seed = cell$code + scenario_index * 1000L + replication,
        priority = 3L
      )
    }
  }
}

# Experiment C: fixed-Delta long-span ladder on the six frozen stress anchors.
long_cells <- list(
  list(n = 2500000L, terminal = 25, code = 310000000L),
  list(n = 5000000L, terminal = 50, code = 320000000L),
  list(n = 10000000L, terminal = 100, code = 330000000L)
)
long_indices <- which(scenarios$long_span_subset)
for (cell in long_cells) {
  for (scenario_index in long_indices) {
    for (replication in 1:3) {
      add_run(
        scenario_index = scenario_index,
        panel = "development",
        experiment = "long_span",
        stage = "feasibility_characterization",
        replication = replication,
        n_steps = cell$n,
        terminal = cell$terminal,
        seed = cell$code + scenario_index * 1000L + replication,
        priority = 4L
      )
    }
  }
}

# Experiment E: four additional production-design paths in five frozen
# high-risk regimes. The first three form the stress feasibility stage.
stress_indices <- which(scenarios$stress_subset)
for (scenario_index in stress_indices) {
  for (replication in 1:4) {
    add_run(
      scenario_index = scenario_index,
      panel = "development",
      experiment = "heavy_tail_stress",
      stage = if (replication <= 3L) {
        "stress_feasibility"
      } else {
        "stress_expansion"
      },
      replication = replication,
      n_steps = 10000000L,
      terminal = 50,
      seed = 410000000L + scenario_index * 1000L + replication,
      priority = 5L
    )
  }
}

seed_manifest <- do.call(rbind, rows)
if (anyDuplicated(seed_manifest$seed)) {
  stop("Seed blocks are not disjoint.")
}
if (anyDuplicated(seed_manifest$run_key)) {
  stop("Run keys are not unique.")
}
if (anyDuplicated(seed_manifest$case)) {
  stop("Run case identifiers are not unique.")
}

run_design <- merge(
  seed_manifest,
  scenarios,
  by = "scenario_id",
  all.x = TRUE,
  sort = FALSE
)
run_design <- run_design[
  match(seed_manifest$run_key, run_design$run_key),
  ,
  drop = FALSE
]
runner_design <- data.frame(
  panel = run_design$panel,
  run_id = run_design$run_key,
  case = run_design$case,
  seed = run_design$seed,
  alpha = run_design$alpha,
  sigma2 = run_design$sigma2,
  sigma1 = run_design$sigma1,
  rho = run_design$rho,
  n_steps = run_design$n_steps,
  terminal = run_design$terminal,
  stringsAsFactors = FALSE
)

write.csv(
  scenarios,
  "frozen_scenario_design.csv",
  row.names = FALSE
)
write.csv(
  seed_manifest,
  "frozen_seed_manifest.csv",
  row.names = FALSE
)
write.csv(
  runner_design,
  "frozen_run_design.csv",
  row.names = FALSE
)

cat("Frozen scenarios:", nrow(scenarios), "\n")
cat("Frozen runs:", nrow(seed_manifest), "\n")
print(table(seed_manifest$experiment, seed_manifest$panel))
print(table(seed_manifest$priority))
