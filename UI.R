# app.R
suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(lubridate)
  library(vars)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(readr)
  library(strucchange)
  library(sandwich)
  library(car)
})

# ---------- Helpers ----------
prep_df <- function(df) {
  df <- as.data.frame(df)
  if ("DATE" %in% names(df)) df <- df %>% mutate(Date = as.Date(DATE))
  else if ("Date" %in% names(df)) df <- df %>% mutate(Date = as.Date(Date))
  else stop("No DATE/Date column found.")
  df %>%
    arrange(Date) %>%
    mutate(
      SP500_Returns = suppressWarnings(as.numeric(SP500_Returns)),
      Real_GDP      = if ("Real_GDP" %in% names(.)) suppressWarnings(as.numeric(Real_GDP)) else NA_real_,
      M2            = if ("M2" %in% names(.))       suppressWarnings(as.numeric(M2))       else NA_real_
    )
}

roll_sum4  <- function(x){ n <- length(x); out <- rep(NA_real_, n); if (n>=4) for(i in 4:n){xs<-x[(i-3):i]; out[i]<-if(all(is.finite(xs))) sum(xs) else NA_real_}; out }
roll_mean4 <- function(x){ n <- length(x); out <- rep(NA_real_, n); if (n>=4) for(i in 4:n){xs<-x[(i-3):i]; out[i]<-if(all(is.finite(xs))) mean(xs) else NA_real_}; out }
roll_comp4 <- function(r_pct){
  n <- length(r_pct); out <- rep(NA_real_, n)
  if (n>=4){
    one_plus <- 1 + (r_pct/100); one_plus[!is.finite(one_plus)] <- NA_real_
    for (i in 4:n){ xs <- one_plus[(i-3):i]; out[i] <- if (all(is.finite(xs))) (prod(xs)-1)*100 else NA_real_ }
  }
  out
}

make_yearly <- function(dfq, returns_method = "sum4", overlap = TRUE, anchor_q = 4) {
  dfq <- prep_df(dfq) %>% mutate(q = quarter(Date), y = year(Date))
  ret_y <- switch(returns_method,
                  sum4 = roll_sum4(dfq$SP500_Returns),
                  compound4 = roll_comp4(dfq$SP500_Returns),
                  mean4 = roll_mean4(dfq$SP500_Returns),
                  roll_sum4(dfq$SP500_Returns))
  df_aug <- dfq %>%
    mutate(SP500_y = ret_y, GDP_y = roll_sum4(Real_GDP), M2_y = roll_sum4(M2))
  if (overlap) {
    df_aug %>%
      transmute(Date, SP500_Returns = SP500_y, Real_GDP = GDP_y, M2 = M2_y) %>%
      filter(is.finite(SP500_Returns), is.finite(Real_GDP), is.finite(M2))
  } else {
    df_aug %>%
      filter(q == anchor_q) %>%
      transmute(Date, SP500_Returns = SP500_y, Real_GDP = GDP_y, M2 = M2_y) %>%
      filter(is.finite(SP500_Returns), is.finite(Real_GDP), is.finite(M2))
  }
}

yq_date <- function(year, quarter) {
  m <- c(1,4,7,10)[pmax(1, pmin(4, quarter))]
  as.Date(sprintf("%04d-%02d-01", as.integer(year), m))
}

crit_label_to_key <- function(lbl) {
  switch(lbl, "AIC"="AIC(n)", "HQ"="HQ(n)", "BIC (Schwarz)"="SC(n)", "FPE"="FPE", "Override"="OVERRIDE", "SC(n)")
}

fit_var <- function(dat, lag_max = 8, crit_key = "SC(n)", const_type = "const", p_override = NULL) {
  type_val <- match.arg(const_type, c("const","none","trend"))
  if (!is.null(p_override)) {
    p <- as.integer(p_override)
    v <- vars::VAR(dat, p = p, type = type_val); v$call$type <- type_val
    return(list(model=v, p=p, selection=NULL))
  }
  sel <- vars::VARselect(dat, lag.max = lag_max, type = type_val)
  p <- as.integer(sel$selection[[crit_key]]); if (is.na(p) || p < 1) p <- 2
  v <- vars::VAR(dat, p = p, type = type_val); v$call$type <- type_val
  list(model = v, p = p, selection = sel)
}

build_var_ols <- function(dat_var, dep, p) {
  stopifnot(dep %in% colnames(dat_var))
  k  <- ncol(dat_var)
  nm <- colnames(dat_var)
  out <- dat_var
  for (j in seq_len(k)) {
    for (L in 1:p) {
      out[[paste0(nm[j], ".l", L)]] <- dplyr::lag(out[[nm[j]]], L)
    }
  }
  Xcols <- grep("\\.l[0-9]+$", colnames(out), value = TRUE)
  keep  <- stats::complete.cases(out[, Xcols, drop = FALSE]) &
    stats::complete.cases(out[, dep,  drop = FALSE])
  out2 <- out[keep, , drop = FALSE]
  frm  <- as.formula(paste(dep, "~", paste(Xcols, collapse = " + ")))
  list(data = out2, formula = frm, keep_idx = which(keep))
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Macro–Return VAR Explorer (Enhanced)"),
  tabsetPanel(
    id = "tabs",
    tabPanel(
      "VAR & IRFs",
      sidebarLayout(
        sidebarPanel(
          selectInput("dataset", "Dataset",
                      c("data_all", "data_all_sec", "data_all_yearly"), selected = "data_all"),
          uiOutput("year_quarter_ui"),
          checkboxGroupInput("predictors", "Include predictors",
                             choices = c("Real_GDP","M2"), selected = c("Real_GDP","M2")),
          conditionalPanel(
            condition = "input.dataset == 'data_all_yearly'",
            selectInput("yr_ret_method", "Yearly aggregation for SP500_Returns",
                        choices = c("Sum of 4 quarters (log returns)"="sum4",
                                    "Compounded 4 quarters (simple %)"="compound4",
                                    "Average of 4 quarters"="mean4"),
                        selected = "sum4"),
            radioButtons("yr_overlap", "Yearly overlap",
                         choices = c("Overlapping (rolling 4q)"="yes",
                                     "Non-overlapping (anchor quarter only)"="no"),
                         selected = "yes"),
            conditionalPanel(
              condition = "input.yr_overlap == 'no'",
              sliderInput("yr_anchor_q", "Anchor quarter (for non-overlap)", min=1, max=4, value=4, step=1)
            )
          ),
          sliderInput("lagmax", "Lag max (for selection)", min=1, max=12, value=8, step=1),
          selectInput("crit", "Lag selection criterion",
                      c("BIC (Schwarz)","AIC","HQ","FPE","Override"), selected="BIC (Schwarz)"),
          conditionalPanel(
            condition = "input.crit == 'Override'",
            sliderInput("plag", "Manual lag (p)", min=1, max=12, value=2, step=1)
          ),
          selectInput("const_type", "Deterministic term",
                      c("const","none","trend"), selected="const"),
          sliderInput("irf_h", "IRF horizon (quarters)", min=4, max=20, value=12, step=1),
          radioButtons("irf_kind", "IRF type",
                       choices = c("Orthogonal (Cholesky)"="ortho",
                                   "Generalized (order-invariant)"="girf"),
                       selected = "ortho"),
          actionButton("run", "Run VAR", class = "btn-primary"),
          br(), br(),
          helpText("Tip: Non-overlapping yearly = pick one anchor quarter (e.g., Q4) per year.")
        ),
        mainPanel(
          tags$h4("Status"), verbatimTextOutput("status"),
          tags$h4("Chosen lag & summary"), verbatimTextOutput("summary"),
          tags$h4("Diagnostics"), verbatimTextOutput("diagnostics"),
          tags$h4("Granger causality to SP500_Returns"), verbatimTextOutput("granger"),
          fluidRow(column(6, plotOutput("irf_gdp")), column(6, plotOutput("irf_m2"))),
          tags$h4("FEVD of SP500_Returns"), tableOutput("fevd_tbl")
        )
      )
    ),
    tabPanel(
      "Structural Breaks (Bai–Perron)",
      sidebarLayout(
        sidebarPanel(
          helpText("Detect multiple breaks in a single VAR equation (OLS, Bai–Perron). Uses current dataset/period and lag p."),
          selectInput("bp_dep", "Equation (dependent variable)",
                      choices = c("SP500_Returns","Real_GDP","M2"), selected = "SP500_Returns"),
          numericInput("bp_max_breaks", "Max breaks (K max)", value = 5, min = 0, max = 10, step = 1),
          sliderInput("bp_trim", "Trimming proportion (min segment size)", min = 0.10, max = 0.35, value = 0.15, step = 0.01),
          radioButtons("bp_vcov", "Covariance for inference",
                       choices = c("IID (default)"="iid", "HAC (Newey–West, lag=3)"="hac"),
                       selected = "hac"),
          checkboxInput("bp_force_k", "Force number of breaks (override BIC)?", value = FALSE),
          conditionalPanel(
            condition = "input.bp_force_k == true",
            numericInput("bp_k", "Number of breaks (k)", value = 1, min = 0, max = 10, step = 1)
          ),
          actionButton("bp_run", "Detect breaks", class = "btn-warning")
        ),
        mainPanel(
          tags$h4("Break selection & dates"),
          verbatimTextOutput("bp_summary"),
          tags$h4("BIC by number of breaks"),
          plotOutput("bp_bic_plot", height = 260),
          tags$h4("Segment-wise coefficients (selected k)"),
          tableOutput("bp_coef_tbl")
        )
      )
    ),
    # ---------- New Tab: Break date converter ----------
    tabPanel(
      "Break date converter",
      sidebarLayout(
        sidebarPanel(
          radioButtons("conv_mode", "Date source",
                       choices = c("Use current period (recommended)" = "current",
                                   "Custom quarterly sequence" = "custom"),
                       selected = "current"),
          conditionalPanel(
            condition = "input.conv_mode == 'custom'",
            numericInput("conv_start_year", "Start year", value = 1960, min = 1800, step = 1),
            sliderInput("conv_start_q", "Start quarter", min = 1, max = 4, value = 1, step = 1),
            numericInput("conv_n", "Number of observations (n)", value = 64, min = 1, step = 1)
          ),
          textInput("conv_fracs", "Fractional break(s)",
                    placeholder = "e.g. 0.125, 0.5 0.796875; separate by comma/space/semicolon"),
          actionButton("conv_run", "Convert", class = "btn-info"),
          helpText("Rule: index = round(fraction * n), clamped to [1, n]. n = length of the chosen date sequence.")
        ),
        mainPanel(
          tags$h4("Converted breaks"),
          tableOutput("conv_table")
        )
      )
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  # Prepare datasets (reactive)
  all_data <- reactive({
    returns_method <- if (is.null(input$yr_ret_method)) "sum4" else input$yr_ret_method
    overlap_flag   <- if (is.null(input$yr_overlap)) TRUE else input$yr_overlap == "yes"
    anchor_q       <- if (is.null(input$yr_anchor_q)) 4 else input$yr_anchor_q
    lst <- list(
      data_all        = if (exists("data_all"))     prep_df(get("data_all"))     else NULL,
      data_all_sec    = if (exists("data_all_sec")) prep_df(get("data_all_sec")) else NULL,
      data_all_yearly = if (exists("data_all"))     make_yearly(get("data_all"),
                                                                returns_method = returns_method,
                                                                overlap = overlap_flag,
                                                                anchor_q = anchor_q) else NULL
    )
    validate(need(!is.null(lst$data_all) || !is.null(lst$data_all_sec),
                  "Load data_all / data_all_sec into the R session first."))
    lst
  })
  
  # Year/quarter UI
  output$year_quarter_ui <- renderUI({
    dat <- all_data()[[input$dataset]]; req(dat)
    yrs <- year(range(dat$Date, na.rm = TRUE))
    tagList(
      sliderInput("years", "Year range", min = min(yrs), max = max(yrs),
                  value = c(min(yrs), max(yrs)), step = 1, sep = ""),
      fluidRow(
        column(6, sliderInput("start_q", "Start quarter", min = 1, max = 4, value = 1, step = 1)),
        column(6, sliderInput("end_q",   "End quarter",   min = 1, max = 4, value = 4, step = 1))
      )
    )
  })
  
  # Period slice
  period_slice <- reactive({
    dat <- all_data()[[input$dataset]]
    req(dat, input$years, input$start_q, input$end_q)
    y0 <- input$years[1]; y1 <- input$years[2]
    d_start <- yq_date(y0, input$start_q); d_end <- yq_date(y1, input$end_q)
    if (d_end < d_start) { tmp <- d_start; d_start <- d_end; d_end <- tmp }
    dat %>% filter(Date >= d_start, Date <= d_end)
  })
  
  # ---- VAR run ----
  run_res <- eventReactive(input$run, {
    dat <- period_slice()
    validate(need(!is.null(dat) && nrow(dat) > 10, "Selected period has too few observations."))
    vars_to_use <- c("SP500_Returns", intersect(input$predictors, c("Real_GDP","M2")))
    validate(need(length(vars_to_use) >= 2, "Select at least one predictor."))
    dat_var <- dat %>% dplyr::select(dplyr::all_of(vars_to_use)) %>% tidyr::drop_na()
    crit_key <- crit_label_to_key(input$crit)
    if (crit_key == "OVERRIDE") {
      p_man <- input$plag
      validate(need(nrow(dat_var) >= (p_man + 20),
                    "Not enough observations for the chosen manual lag p. Reduce p or widen the period."))
      res <- fit_var(dat_var, const_type = input$const_type, p_override = p_man)
    } else {
      validate(need(nrow(dat_var) >= (input$lagmax + 20),
                    "Not enough observations for the chosen lag.max. Reduce lag.max or widen the period."))
      res <- fit_var(dat_var, lag_max = input$lagmax, crit_key = crit_key, const_type = input$const_type)
    }
    serial <- vars::serial.test(res$model, lags.pt = 16, type = "PT.asymptotic")
    arch   <- vars::arch.test(res$model,  lags.multi = 12)
    normt  <- vars::normality.test(res$model)
    rootsM <- Mod(vars::roots(res$model))
    ortho_flag <- (input$irf_kind == "ortho")
    irf_gdp <- if ("Real_GDP" %in% colnames(dat_var))
      vars::irf(res$model, impulse="Real_GDP", response="SP500_Returns",
                n.ahead=input$irf_h, boot=TRUE, runs=500, ci=0.95, ortho=ortho_flag) else NULL
    irf_m2  <- if ("M2" %in% colnames(dat_var))
      vars::irf(res$model, impulse="M2", response="SP500_Returns",
                n.ahead=input$irf_h, boot=TRUE, runs=500, ci=0.95, ortho=ortho_flag) else NULL
    gc_txt <- capture.output({
      others <- setdiff(colnames(dat_var), "SP500_Returns")
      
      ## =========================================================
      ## 1) Equation-level Granger causality (full VAR, SP500 eq.)
      ## =========================================================
      cat("=== Equation-level Granger causality (full VAR, SP500 equation) ===\n")
      
      # SP500 equation from the full VAR
      eq_sp <- res$model$varresult$SP500_Returns
      
      # Helper: build and run linearHypothesis for one variable's lags
      test_in_sp500_eq <- function(varname) {
        # Only if that variable is in the VAR
        if (!varname %in% colnames(dat_var)) return(NULL)
        
        # Coefficient names for this equation
        cf_names <- names(coef(eq_sp))
        
        # All possible lag terms for this variable
        lag_terms <- paste0(varname, ".l", 1:res$p)
        lag_terms <- lag_terms[lag_terms %in% cf_names]
        
        if (length(lag_terms) == 0) return(NULL)
        
        hyp <- paste(lag_terms, "= 0")
        cat("\n--- H0: past", varname, "do NOT Granger-cause SP500_Returns (in full VAR) ---\n")
        print(car::linearHypothesis(eq_sp, hyp))
      }
      
      # Does M2 Granger-cause SP500 in the full 3-var VAR?
      test_in_sp500_eq("M2")
      
      # Does Real_GDP Granger-cause SP500 in the full 3-var VAR?
      test_in_sp500_eq("Real_GDP")
      
      ## =========================================================
      ## 2) Pairwise (bivariate) Granger tests (SP500 <-> others)
      ## =========================================================
      cat("\n\n=== Pairwise (bivariate) Granger causality ===\n")
      
      if (length(others) == 0) {
        cat("No other variables besides SP500_Returns in VAR.\n")
      } else {
        for (v in others) {
          cat("\n---------------------------------------------\n")
          cat("Bivariate VAR for:", v, "<-> SP500_Returns\n")
          cat("---------------------------------------------\n")
          
          # Build a bivariate VAR with SP500_Returns and the current variable v
          subdat <- dat_var[, c("SP500_Returns", v), drop = FALSE]
          
          subvar <- vars::VAR(
            subdat,
            p    = res$p,                  # same lag length as main VAR
            type = res$model$call$type     # same deterministic term
          )
          
          # 1) Does v (M2 or Real_GDP) Granger-cause SP500?
          cat("\n--- Does", v, "Granger-cause SP500_Returns? (bivariate) ---\n")
          print(vars::causality(subvar, cause = v))
          
          # 2) Does SP500 Granger-cause v?
          cat("\n--- Does SP500_Returns Granger-cause", v, "? (bivariate) ---\n")
          print(vars::causality(subvar, cause = "SP500_Returns"))
        }
      }
    })
    
    fe <- vars::fevd(res$model, n.ahead=input$irf_h)
    fe_tbl <- NULL
    if ("SP500_Returns" %in% names(fe)) {
      comp_names <- colnames(fe$SP500_Returns)
      fe_tbl <- tibble(horizon = 1:input$irf_h)
      for (nm in comp_names) fe_tbl[[paste0("from_", nm)]] <- fe$SP500_Returns[, nm]
    }
    list(
      dat_range = range(dat$Date, na.rm = TRUE),
      vars_used = colnames(dat_var),
      p = res$p, model = res$model,
      serial = serial, arch = arch, normt = normt, roots = rootsM,
      irf_gdp = irf_gdp, irf_m2 = irf_m2, fevd_tbl = fe_tbl,
      gc_text = paste(gc_txt, collapse = "\n"), crit_used = input$crit,
      dat_var = dat_var
    )
  }, ignoreInit = TRUE)
  
  # VAR outputs
  output$status <- renderText({
    if (!input$run) return("Pick options and press 'Run VAR'.")
    rr <- run_res(); req(rr)
    add_yearly <- if (input$dataset == "data_all_yearly")
      paste0("\nYearly aggregation: ",
             input$yr_ret_method %||% "sum4", "; Overlap: ",
             if ((input$yr_overlap %||% "yes") == "yes") "Yes" else paste0("No (Q", input$yr_anchor_q %||% 4, ")"))
    else ""
    paste0("Dates: ", paste(format(rr$dat_range, "%Y-%m-%d"), collapse = " → "),
           "\nVariables: ", paste(rr$vars_used, collapse = ", "),
           "\nChosen lag p (", rr$crit_used, "): ", rr$p,
           "\nIRF type: ", if (input$irf_kind == "ortho") "Orthogonal (Cholesky)" else "Generalized (order-invariant)",
           add_yearly)
  })
  output$summary <- renderPrint({ rr <- run_res(); req(rr); summary(rr$model) })
  output$diagnostics <- renderPrint({ rr <- run_res(); req(rr)
  cat("--- Serial correlation (Portmanteau) ---\n"); print(rr$serial)
  cat("\n--- ARCH effects ---\n"); print(rr$arch)
  cat("\n--- Normality ---\n"); print(rr$normt)
  cat("\n--- Stability roots (modulus) ---\n"); print(rr$roots)
  })
  output$granger <- renderPrint({ rr <- run_res(); req(rr); cat(rr$gc_text) })
  output$irf_gdp <- renderPlot({ rr <- run_res(); req(rr); validate(need(!is.null(rr$irf_gdp), "Real_GDP not included.")); plot(rr$irf_gdp) })
  output$irf_m2  <- renderPlot({ rr <- run_res(); req(rr); validate(need(!is.null(rr$irf_m2),  "M2 not included.")); plot(rr$irf_m2) })
  output$fevd_tbl <- renderTable({ rr <- run_res(); req(rr); rr$fevd_tbl })
  
  # ---------- Bai–Perron detection ----------
  bp_res <- eventReactive(input$bp_run, {
    rr <- isolate(run_res())
    validate(need(!is.null(rr), "Run the VAR first to set data and lag p."))
    dat_var <- rr$dat_var
    p       <- rr$p
    dep     <- input$bp_dep
    validate(need(dep %in% colnames(dat_var), "Chosen equation not in current model."))
    
    ols <- build_var_ols(dat_var, dep = dep, p = p)
    df  <- ols$data
    frm <- ols$formula
    
    n    <- nrow(df)
    h    <- max(1, floor(input$bp_trim * n))
    Kmax <- input$bp_max_breaks
    
    bp <- strucchange::breakpoints(frm, data = df, h = h, breaks = Kmax)
    if (isTRUE(input$bp_force_k)) {
      k_star <- as.integer(input$bp_k)
    } else {
      bic_vals <- stats::BIC(bp)
      k_star   <- which.min(bic_vals) - 1L
    }
    if (is.na(k_star)) k_star <- 0L
    
    bd <- try(strucchange::breakdates(bp, breaks = k_star), silent = TRUE)
    if (inherits(bd, "try-error")) bd <- NULL
    
    ci_df <- NULL
    ci_try <- try(strucchange::confint(bp, breaks = k_star), silent = TRUE)
    if (!inherits(ci_try, "try-error") && !is.null(ci_try$confint)) {
      ci_mat <- as.data.frame(ci_try$confint)
      if (ncol(ci_mat) == 3) {
        names(ci_mat) <- c("index","start","end")
        ci_df <- data.frame(
          break_id = seq_len(nrow(ci_mat)),
          index    = ci_mat$index,
          start    = ci_mat$start,
          end      = ci_mat$end,
          row.names = NULL
        )
      }
    }
    
    # Segment-wise coefficients (optional HAC SEs)
    vcov_fun <- if (input$bp_vcov == "hac")
      function(x) sandwich::NeweyWest(x, lag = 3, prewhite = FALSE) else NULL
    
    bp_idxs <- bp$breakpoints
    if (length(bp_idxs) == 0 || is.na(bp_idxs[1])) bp_idxs <- integer(0)
    if (k_star < length(bp_idxs)) bp_idxs <- bp_idxs[seq_len(k_star)]
    idx <- c(0L, bp_idxs, n)
    
    seg_coefs <- list()
    for (s in seq_len(length(idx) - 1L)) {
      i1 <- idx[s] + 1L; i2 <- idx[s + 1L]
      seg_fit <- stats::lm(frm, data = df[i1:i2, , drop = FALSE])
      seg_vc  <- if (is.null(vcov_fun)) stats::vcov(seg_fit) else vcov_fun(seg_fit)
      cf <- stats::coef(seg_fit)
      se <- sqrt(diag(seg_vc))
      seg_coefs[[s]] <- data.frame(
        segment   = s,
        term      = names(cf),
        estimate  = as.numeric(cf),
        std_error = as.numeric(se),
        row.names = NULL
      )
    }
    seg_tbl <- dplyr::bind_rows(seg_coefs)
    
    bic_vals <- stats::BIC(bp)
    bic_df   <- data.frame(k = seq_along(bic_vals) - 1L, BIC = as.numeric(bic_vals))
    
    list(
      dep = dep, p = p, n = n, h = h, Kmax = Kmax,
      k_star = k_star, breaks_idx = bp$breakpoints,
      breakdates = bd, ci_df = ci_df,
      seg_tbl = seg_tbl, bic_df = bic_df
    )
  }, ignoreInit = TRUE)
  
  output$bp_summary <- renderPrint({
    res <- bp_res(); req(res)
    cat("Equation:", res$dep, "| lag p:", res$p, "| n:", res$n, "\n")
    cat("Kmax:", res$Kmax, "| trimming h:", res$h, "\n")
    cat("Selected #breaks (BIC unless overridden):", res$k_star, "\n")
    if (!is.null(res$breakdates) && length(res$breakdates)) {
      cat("\nBreak dates (approx., fractional scale):\n"); print(res$breakdates)
    } else {
      cat("\nNo breaks selected or dates unavailable.\n")
    }
    if (!is.null(res$ci_df)) {
      cat("\nBreakpoint confidence intervals (indices):\n"); print(res$ci_df)
    }
  })
  
  output$bp_bic_plot <- renderPlot({
    res <- bp_res(); req(res)
    ggplot(res$bic_df, aes(k, BIC)) +
      geom_line() + geom_point() +
      geom_vline(xintercept = res$k_star, linetype = 2) +
      labs(x = "Number of breaks (k)", y = "BIC",
           title = "BIC across break counts (lower is better)") +
      theme_minimal(base_size = 12)
  })
  
  output$bp_coef_tbl <- renderTable({
    res <- bp_res(); req(res)
    res$seg_tbl
  })
  
  # ---------- Break date converter ----------
  parse_fracs <- function(txt) {
    if (is.null(txt) || !nzchar(txt)) return(numeric(0))
    parts <- unlist(strsplit(txt, "[,;\\s]+"))
    vals  <- suppressWarnings(as.numeric(parts))
    vals[is.finite(vals)]
  }
  qy_lab <- function(d) paste0("Q", lubridate::quarter(d), "/", lubridate::year(d))
  
  conv_tbl <- eventReactive(input$conv_run, {
    fracs <- parse_fracs(input$conv_fracs)
    validate(need(length(fracs) > 0, "Provide at least one fractional value."))
    
    if (input$conv_mode == "current") {
      dates <- period_slice()$Date
      validate(need(length(dates) > 0, "Current period has no dates."))
      n <- length(dates)
    } else {
      validate(need(!is.null(input$conv_n) && input$conv_n >= 1, "Set a positive n."))
      start <- yq_date(input$conv_start_year, input$conv_start_q)
      n <- as.integer(input$conv_n)
      dates <- start + months(3 * (0:(n-1)))
    }
    
    idx <- round(fracs * n)
    idx[idx < 1] <- 1L
    idx[idx > n] <- n
    
    data.frame(
      fraction = fracs,
      n_obs    = n,
      index    = idx,
      date     = dates[idx],
      qy       = qy_lab(dates[idx]),
      row.names = NULL
    )
  }, ignoreInit = TRUE)
  
  output$conv_table <- renderTable({
    conv_tbl()
  })
}

shinyApp(ui, server)

