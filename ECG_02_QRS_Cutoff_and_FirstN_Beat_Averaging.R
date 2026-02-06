###############################################################################
#
# Standardized / publishable pipeline (NO metadata required):
#  1) Select folder with "*_2min.csv" files (GUI)
#  2) Pool QRS, fit 2-component Gaussian mixture (mclust)
#  3) Compute candidate cutoffs + save comparison PDFs
#  4) Choose cutoff strategy (GUI) — simplified 2-step UI
#  5) Enter N beats in ONE textbox (default = 50)
#  6) For each file: find first row where QRS >= cutoff, extract next N beats
#     - If not enough rows to reach N beats, skip that file (explicit warning)
#  7) Compute mean of numeric columns for that segment
#  8) Save summary CSV + cleaned XLSX with renamed columns
#
# OUTPUT STRUCTURE:
#   <OUTPUT>/ECG_Run_First<N>_<CutoffLabel>_YYYYmmdd_HHMMSS/
#     01_QRS_Cutoffs/
#     02_FirstNBeats_Summary/
#
# Windows/RStudio Tk robustness:
#  - Avoids tkwait.visibility() (causes ".X was deleted before its visibility changed")
#  - Uses tkwait.variable(done) and handles window close (X)
#  - Forces listbox colors to avoid "black box" / invisible text
###############################################################################

suppressPackageStartupMessages({ library(tcltk) })

ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}
pkgs <- c("dplyr","tidyr","stringr","openxlsx","mclust","tibble")
invisible(lapply(pkgs, ensure_pkg))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(mclust)
  library(tibble)
})

# ---------------------------- GUI helpers ------------------------------------

pick_dir <- function(caption = "Select a folder") {
  p <- tcltk::tk_choose.dir(caption = caption)
  if (is.na(p) || !nzchar(p)) stop("No folder selected.")
  normalizePath(p, winslash = "/", mustWork = TRUE)
}

# Robust modal single-selection dialog (no tkwait.visibility)
gui_select_one <- function(title, choices, width_chars = 60, height_rows = 14) {
  if (!length(choices)) stop("No choices provided.")
  
  tt <- tktoplevel()
  tkwm.title(tt, title)
  tkwm.resizable(tt, 0, 0)
  
  res  <- tclVar("")
  done <- tclVar(0)
  
  frm <- tkframe(tt)
  
  lb <- tklistbox(
    frm,
    selectmode = "single",
    height = height_rows,
    width  = width_chars,
    exportselection = FALSE
  )
  
  # Force visible colors (avoids black/invisible listbox)
  tkconfigure(lb,
              background       = "white",
              foreground       = "black",
              selectbackground = "#cce5ff",
              selectforeground = "black")
  
  sb <- tkscrollbar(frm, orient = "vertical", command = function(...) tkset(lb, ...))
  tkconfigure(lb, yscrollcommand = function(...) tkset(sb, ...))
  
  for (ch in choices) tkinsert(lb, "end", ch)
  
  tkgrid(lb, sb, sticky = "nsew")
  tkgrid.configure(sb, sticky = "ns")
  
  btn_frm <- tkframe(tt)
  
  on_ok <- function() {
    sel <- as.integer(tkcurselection(lb))
    if (length(sel) && !is.na(sel)) tclvalue(res) <- choices[sel + 1] else tclvalue(res) <- ""
    tclvalue(done) <- 1
  }
  on_cancel <- function() {
    tclvalue(res) <- ""
    tclvalue(done) <- 1
  }
  
  ok_btn <- tkbutton(btn_frm, text = "OK", width = 10, command = on_ok)
  cancel_btn <- tkbutton(btn_frm, text = "Cancel", width = 10, command = on_cancel)
  
  tkbind(lb, "<Double-Button-1>", function() on_ok())
  tkbind(tt, "<Escape>", function() on_cancel())
  
  # Handle window close (X)
  tkbind(tt, "<Destroy>", function() {
    if (as.integer(tclvalue(done)) == 0) {
      tclvalue(res) <- ""
      tclvalue(done) <- 1
    }
  })
  
  tkgrid(frm, padx = 10, pady = c(10, 5))
  tkgrid(btn_frm, padx = 10, pady = c(0, 10), sticky = "e")
  tkgrid(ok_btn, cancel_btn, padx = 5)
  
  # Modal behavior
  tkfocus(lb)
  suppressWarnings(try(tkgrab.set(tt), silent = TRUE))
  tkwait.variable(done)
  suppressWarnings(try(tkgrab.release(tt), silent = TRUE))
  suppressWarnings(try(tkdestroy(tt), silent = TRUE))
  
  out <- tclvalue(res)
  if (!nzchar(out)) stop("No selection made.")
  out
}

# One textbox integer input (robust)
gui_textbox_integer <- function(title, label, default_value = 50) {
  tt <- tktoplevel()
  tkwm.title(tt, title)
  tkwm.resizable(tt, 0, 0)
  
  done <- tclVar(0)
  val  <- tclVar(as.character(default_value))
  
  lbl <- tklabel(tt, text = label)
  entry <- tkentry(tt, textvariable = val, width = 24)
  
  btn_frm <- tkframe(tt)
  
  on_ok <- function() tclvalue(done) <- 1
  on_cancel <- function() tclvalue(done) <- 2
  
  ok_btn <- tkbutton(btn_frm, text = "OK", width = 10, command = on_ok)
  cancel_btn <- tkbutton(btn_frm, text = "Cancel", width = 10, command = on_cancel)
  
  tkgrid(lbl, padx = 10, pady = c(10, 2), sticky = "w")
  tkgrid(entry, padx = 10, pady = c(0, 10), sticky = "we")
  tkgrid(btn_frm, padx = 10, pady = c(0, 10), sticky = "e")
  tkgrid(ok_btn, cancel_btn, padx = 5)
  
  tkfocus(entry)
  tkselection.range(entry, 0, "end")
  tkbind(entry, "<FocusIn>", function() tkselection.range(entry, 0, "end"))
  tkbind(entry, "<Return>", function() on_ok())
  tkbind(tt, "<Escape>", function() on_cancel())
  
  tkbind(tt, "<Destroy>", function() {
    if (as.integer(tclvalue(done)) == 0) tclvalue(done) <- 2
  })
  
  suppressWarnings(try(tkgrab.set(tt), silent = TRUE))
  tkwait.variable(done)
  suppressWarnings(try(tkgrab.release(tt), silent = TRUE))
  suppressWarnings(try(tkdestroy(tt), silent = TRUE))
  
  res <- as.integer(tclvalue(done))
  if (res != 1) stop("Cancelled.")
  
  out <- tclvalue(val)
  n <- suppressWarnings(as.integer(out))
  if (is.na(n) || !is.finite(n)) stop("Invalid integer for N beats.")
  if (n < 1) stop("N beats must be >= 1.")
  n
}

# ---------------------------- Robust file + ID handling -----------------------

FILE_PATTERN <- "_2min\\.csv$"

list_2min_files <- function(dir_injected) {
  list.files(dir_injected, pattern = FILE_PATTERN, full.names = TRUE)
}

mouse_id_from_filename <- function(file_path) {
  stem <- tools::file_path_sans_ext(basename(file_path))
  if (grepl("_", stem)) sub("_.*$", "", stem) else stem
}

# ---------------------------- Utilities --------------------------------------

to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

rename_cols_base <- function(df, old_to_new) {
  nm <- names(df)
  for (old in names(old_to_new)) {
    if (old %in% nm) nm[nm == old] <- old_to_new[[old]]
  }
  names(df) <- nm
  df
}

# ---------------------------- Defaults ----------------------------------------

qrs_col_default <- "QRS.Interval..s."
posterior_p2_default <- 0.10

# old -> new
rename_map_old_to_new <- c(
  "RR.Interval..s."           = "RR_Interval_sec",
  "Heart.Rate..BPM."          = "Heart_Rate_BPM",
  "PR.Interval..s."           = "PR_Interval_sec",
  "P.Duration..s."            = "P_Duration_sec",
  "QRS.Interval..s."          = "QRS_Interval_sec",
  "QT.Interval..s."           = "QT_Interval_sec",
  "QTc..s."                   = "QTc_sec",
  "JT.Interval..s."           = "JT_Interval_sec",
  "Tpeak.Tend.Interval..s."   = "Tpeak_Tend_Interval_sec",
  "P.Amplitude..mV."          = "P_Amplitude_mV",
  "Q.Amplitude..mV."          = "Q_Amplitude_mV",
  "R.Amplitude..mV."          = "R_Amplitude_mV",
  "S.Amplitude..mV."          = "S_Amplitude_mV",
  "ST.Height..mV."            = "ST_Height_mV",
  "T.Amplitude..mV."          = "T_Amplitude_mV"
)

# ---------------------------- 1) Select inputs --------------------------------

message("Step 1/3 — Select folder containing '*_2min.csv' files.")
dir_injected <- pick_dir("Select folder with *_2min.csv files")

message("Step 2/3 — Select OUTPUT folder (a run subfolder will be created inside).")
out_root <- pick_dir("Select OUTPUT folder")

files <- list_2min_files(dir_injected)
if (!length(files)) stop(sprintf("No '*_2min.csv' files found in:\n%s", dir_injected))

# ---------------------------- QRS column (easier) -----------------------------

# Auto-detect QRS column if default exists; otherwise ask once
df0 <- read.csv(files[1], stringsAsFactors = FALSE, dec = ",", check.names = FALSE)

if (qrs_col_default %in% names(df0)) {
  qrs_col <- qrs_col_default
} else {
  qrs_col <- gui_select_one(
    title   = "Pick the QRS column (not auto-detected)",
    choices = names(df0)
  )
}
message("Using QRS column: ", qrs_col)

# ---------------------------- 2) Pool QRS -------------------------------------

read_qrs_from_file <- function(f, qrs_col) {
  df <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE, dec = ",", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !(qrs_col %in% names(df))) return(numeric())
  nums <- to_num(df[[qrs_col]])
  nums[is.finite(nums)]
}

all_qrs <- unlist(lapply(files, read_qrs_from_file, qrs_col = qrs_col), use.names = FALSE)
if (length(all_qrs) < 20) stop("Too few QRS values to fit a mixture (need >= 20).")

# ---------------------------- 3) Fit mixture ----------------------------------

mix <- Mclust(all_qrs, G = 2, modelNames = "V")

mu   <- mix$parameters$mean
sig  <- sqrt(mix$parameters$variance$sigmasq)
pi_k <- mix$parameters$pro

ord <- order(mu)
mu <- mu[ord]; sig <- sig[ord]; pi_k <- pi_k[ord]

f1 <- function(x) pi_k[1] * dnorm(x, mu[1], sig[1])
f2 <- function(x) pi_k[2] * dnorm(x, mu[2], sig[2])
fmix <- function(x) f1(x) + f2(x)
post2 <- function(x) f2(x) / fmix(x)

xs_between <- seq(mu[1], mu[2], length.out = 5000)
cut_equal  <- uniroot(function(x) f1(x) - f2(x), lower = mu[1], upper = mu[2])$root
cut_valley <- xs_between[which.min(fmix(xs_between))]

posterior_p2 <- posterior_p2_default
cut_p2 <- uniroot(function(x) post2(x) - posterior_p2, lower = mu[1], upper = mu[2])$root

# Comp1 percentiles available
comp1_pcts <- c(60, 65, 70, 75, 80, 85, 90, 95)
comp1_cuts <- qnorm(comp1_pcts / 100, mean = mu[1], sd = sig[1])

# ---------------------------- 4) Choose cutoff strategy (simpler GUI) --------

message("Step 3/3 — Choose cutoff strategy.")
method <- gui_select_one(
  "Cutoff method",
  c("Comp1 percentile (recommended)",
    "Equal-density intersection",
    "Mixture valley",
    sprintf("Posterior P2 = %.0f%%", posterior_p2 * 100)),
  width_chars = 40,
  height_rows = 8
)

cutoff <- NA_real_
cut_label <- NA_character_

if (method == "Comp1 percentile (recommended)") {
  pct_choice <- gui_select_one(
    "Pick Comp1 percentile",
    paste0(comp1_pcts, "th percentile"),
    width_chars = 25,
    height_rows = length(comp1_pcts) + 2
  )
  pct_num <- as.numeric(sub("th percentile", "", pct_choice))
  cutoff <- qnorm(pct_num / 100, mean = mu[1], sd = sig[1])
  cut_label <- paste0("Comp1P", pct_num)
  
} else if (method == "Equal-density intersection") {
  cutoff <- cut_equal
  cut_label <- "EqualDens"
  
} else if (method == "Mixture valley") {
  cutoff <- cut_valley
  cut_label <- "Valley"
  
} else if (method == sprintf("Posterior P2 = %.0f%%", posterior_p2 * 100)) {
  cutoff <- cut_p2
  cut_label <- paste0("PostP2_", round(posterior_p2 * 100))
  
} else {
  stop("Unexpected selection.")
}

message(sprintf("Using cutoff = %.6f s (%s)", cutoff, cut_label))

# ---------------------------- 5) ONE textbox for N beats ----------------------

segment_n_beats <- gui_textbox_integer(
  title = "How many beats (rows) to average per mouse?",
  label = "Enter N beats (default 50):",
  default_value = 50
)

# ---------------------------- Create run subfolder ----------------------------

run_id  <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(out_root, paste0("ECG_Run_First", segment_n_beats, "_", cut_label, "_", run_id))
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

out_qrs <- file.path(run_dir, "01_QRS_Cutoffs")
out_sum <- file.path(run_dir, "02_FirstNBeats_Summary")
dir.create(out_qrs, showWarnings = FALSE, recursive = TRUE)
dir.create(out_sum, showWarnings = FALSE, recursive = TRUE)

# ---------------------------- Save parameters report ---------------------------

param_file <- file.path(out_qrs, "QRS_cutoff_parameters.txt")
cat(
  "--- Mixture parameters ---\n",
  sprintf("Comp1: mu=%.6f, sd=%.6f, pi=%.4f\n", mu[1], sig[1], pi_k[1]),
  sprintf("Comp2: mu=%.6f, sd=%.6f, pi=%.4f\n\n", mu[2], sig[2], pi_k[2]),
  "--- Candidate cutoffs (seconds) ---\n",
  sprintf("Equal-density intersection: %.6f\n", cut_equal),
  sprintf("Mixture valley:            %.6f\n", cut_valley),
  sprintf("Posterior P2=%.0f%%:         %.6f\n", posterior_p2 * 100, cut_p2),
  paste(sprintf("Comp1 %sth percentile: %.6f", comp1_pcts, comp1_cuts), collapse = "\n"),
  "\n--- Chosen cutoff ---\n",
  sprintf("Label: %s\n", cut_label),
  sprintf("Cutoff (s): %.6f\n", cutoff),
  sprintf("N beats: %d\n", segment_n_beats),
  "\n",
  file = param_file
)

# ---------------------------- Plots -------------------------------------------

pdf(file.path(out_qrs, "QRS_cutoff_comparison.pdf"), width = 9, height = 6)
hist(all_qrs, breaks = 80, freq = FALSE,
     main = "QRS intervals & candidate cutoffs",
     xlab = "QRS interval (s)", border = "grey")

xs <- seq(min(all_qrs), max(all_qrs), length.out = 2000)
lines(xs, pi_k[1] * dnorm(xs, mu[1], sig[1]), lwd = 2)
lines(xs, pi_k[2] * dnorm(xs, mu[2], sig[2]), lwd = 2, lty = 2)

abline(v = cut_equal,  col = "red",    lty = 2, lwd = 2)
abline(v = cut_valley, col = "blue",   lty = 2, lwd = 2)
abline(v = cut_p2,     col = "purple", lty = 2, lwd = 2)

cols_q <- rev(rainbow(length(comp1_cuts)))
for (i in seq_along(comp1_cuts)) abline(v = comp1_cuts[i], col = cols_q[i], lty = 3, lwd = 1.5)

legend("topright", bty = "n", cex = 0.85,
       legend = c(
         sprintf("equal-dens: %.6f", cut_equal),
         sprintf("valley:     %.6f", cut_valley),
         sprintf("P2=%.0f%%:    %.6f", posterior_p2 * 100, cut_p2),
         paste0("Comp1 ", comp1_pcts, "%=", sprintf("%.6f", comp1_cuts))
       ),
       col = c("red", "blue", "purple", cols_q),
       lty = c(2, 2, 2, rep(3, length(comp1_cuts))),
       lwd = c(2, 2, 2, rep(1.5, length(comp1_cuts))))
dev.off()

pdf(file.path(out_qrs, "QRS_chosen_cutoff.pdf"), width = 9, height = 6)
hist(all_qrs, breaks = 80, freq = FALSE,
     main = sprintf("QRS intervals (chosen cutoff = %.6f s)", cutoff),
     xlab = "QRS interval (s)", border = "grey")
lines(xs, pi_k[1] * dnorm(xs, mu[1], sig[1]), lwd = 2)
lines(xs, pi_k[2] * dnorm(xs, mu[2], sig[2]), lwd = 2, lty = 2)
abline(v = cutoff, col = "blue", lty = 2, lwd = 2)
dev.off()

# ---------------------------- First-N beats summary ----------------------------

summary_list <- list()
missing_ids <- character()
short_after_cutoff <- character()

for (file_path in files) {
  mouse_id <- mouse_id_from_filename(file_path)
  
  df <- tryCatch(
    read.csv(file_path, stringsAsFactors = FALSE, dec = ",", check.names = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(df) || !(qrs_col %in% names(df))) {
    warning(sprintf("[%s] Could not read or missing QRS column. Skipping.", mouse_id))
    next
  }
  
  df[[qrs_col]] <- to_num(df[[qrs_col]])
  idx <- which(df[[qrs_col]] >= cutoff)
  
  if (!length(idx)) {
    warning(sprintf("[%s] No QRS >= cutoff (%.6f).", mouse_id, cutoff))
    missing_ids <- c(missing_ids, mouse_id)
    next
  }
  
  start_idx <- idx[1]
  
  # Require FULL N beats; otherwise skip (explicit)
  if (start_idx + (segment_n_beats - 1) > nrow(df)) {
    warning(sprintf("[%s] Not enough rows after cutoff to extract %d beats. Skipping.", mouse_id, segment_n_beats))
    short_after_cutoff <- c(short_after_cutoff, mouse_id)
    next
  }
  
  end_idx <- start_idx + (segment_n_beats - 1)
  
  seg <- df[start_idx:end_idx, , drop = FALSE] %>%
    mutate(across(where(is.character), ~ to_num(.x)))
  
  stats <- seg %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(mouse_id = mouse_id) %>%
    select(mouse_id, everything())
  
  summary_list[[mouse_id]] <- stats
}

summary_df <- bind_rows(summary_list)

out_csv <- file.path(out_sum, sprintf("Summary_First%d_Averages_2min.csv", segment_n_beats))
write.csv(summary_df, out_csv, row.names = FALSE)

summary_clean <- rename_cols_base(summary_df, rename_map_old_to_new)

out_xlsx <- file.path(out_sum, sprintf("Summary_First%d_Final_2min.xlsx", segment_n_beats))
openxlsx::write.xlsx(summary_clean, out_xlsx, overwrite = TRUE)

if (length(missing_ids)) {
  message("Mouse IDs with no QRS >= cutoff: ", paste(missing_ids, collapse = ", "))
}
if (length(short_after_cutoff)) {
  message("Mouse IDs skipped (not enough rows to reach N beats after cutoff): ",
          paste(short_after_cutoff, collapse = ", "))
}

message("\nDone ✅ Standardized outputs are ready.")
message("Next: create metadata OUTSIDE this script, then run Script 02.")
message("Outputs written under: ", run_dir)
