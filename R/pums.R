#' ACS PUMS variables RAILS uses
#'
#' The Census API variable names fetched by [pums_fetch()], with what each is
#' used for. Supplied as a function rather than a constant so it can be printed,
#' subsetted, or extended before being passed back to [pums_fetch()].
#'
#' @param geography `"region"` (the default) or `"state"`. The state form adds
#'   `ST`, the state FIPS code.
#'
#' @return A named character vector: names are the ACS variable codes, values
#'   describe what each contributes.
#'
#' @seealso [pums_fetch()], [pums_recode()].
#'
#' @examples
#' pums_variables()
#'
#' @export
pums_variables <- function(geography = c("region", "state")) {
  geography <- match.arg(geography)
  v <- c(
    PWGTP   = "person weight",
    AGEP    = "age, used for the adult filter and age groups",
    SEX     = "sex",
    HINCP   = "household income, banded",
    SCHL    = "educational attainment, banded",
    TEN     = "tenure, banded into own/rent/other",
    HISP    = "Hispanic origin",
    RACWHT  = "race: White",
    RACBLK  = "race: Black",
    RACASN  = "race: Asian",
    RACAIAN = "race: American Indian or Alaska Native",
    RACNH   = "race: Native Hawaiian",
    RACPI   = "race: Pacific Islander",
    RACSOR  = "race: some other race",
    REGION  = "Census region"
  )
  if (geography == "state") v <- c(v, ST = "state FIPS code")
  v
}

#' State FIPS to USPS abbreviation
#'
#' @keywords internal
#' @noRd
pums_fips <- function() {
  c("01" = "AL", "02" = "AK", "04" = "AZ", "05" = "AR", "06" = "CA",
    "08" = "CO", "09" = "CT", "10" = "DE", "11" = "DC", "12" = "FL",
    "13" = "GA", "15" = "HI", "16" = "ID", "17" = "IL", "18" = "IN",
    "19" = "IA", "20" = "KS", "21" = "KY", "22" = "LA", "23" = "ME",
    "24" = "MD", "25" = "MA", "26" = "MI", "27" = "MN", "28" = "MS",
    "29" = "MO", "30" = "MT", "31" = "NE", "32" = "NV", "33" = "NH",
    "34" = "NJ", "35" = "NM", "36" = "NY", "37" = "NC", "38" = "ND",
    "39" = "OH", "40" = "OK", "41" = "OR", "42" = "PA", "44" = "RI",
    "45" = "SC", "46" = "SD", "47" = "TN", "48" = "TX", "49" = "UT",
    "50" = "VT", "51" = "VA", "53" = "WA", "54" = "WV", "55" = "WI",
    "56" = "WY")
}

#' Download ACS PUMS microdata from the Census API
#'
#' Fetches the person-level PUMS records RAILS uses as its reference sample.
#' The result is raw ACS codes; pass it to [pums_recode()] to get analysis
#' variables.
#'
#' @param year ACS 1-year vintage, e.g. `2022`.
#' @param geography `"region"` fetches the four Census regions in one request.
#'   `"state"` fetches the 50 states plus DC, one request each -- the API caps
#'   how many `ucgid` values a single call accepts, so asking for all 51 at once
#'   returns an error page rather than data.
#' @param variables Named character vector of ACS variable codes; see
#'   [pums_variables()].
#' @param key A Census API key. Free, from
#'   \url{https://api.census.gov/data/key_signup.html}. Defaults to the
#'   `CENSUS_API_KEY` environment variable -- set it in `~/.Renviron` rather
#'   than writing the key into a script, since scripts get committed.
#' @param cache Optional path to a CSV. If the file exists it is read and no
#'   request is made; otherwise the download is written there. Repeat runs of an
#'   analysis should always set this: the state form is 51 API calls.
#' @param verbose Report progress via [message()].
#'
#' @return A data frame of raw PUMS records, one row per person, with the
#'   requested variables as numeric columns.
#'
#' @section Network access:
#' This function requires an internet connection and the \pkg{jsonlite} package.
#' Nothing else in RAILS touches the network.
#'
#' @seealso [pums_recode()] to turn the result into analysis variables, then
#'   [pums_reference()] to get cells and margins.
#'
#' @examples
#' \dontrun{
#' # Sys.setenv(CENSUS_API_KEY = "...")   # better: put it in ~/.Renviron
#' raw <- pums_fetch(2022, cache = "PUMS_2022.csv")
#' ref <- pums_reference(pums_recode(raw))
#' }
#'
#' @export
pums_fetch <- function(year = 2022,
                       geography = c("region", "state"),
                       variables = NULL,
                       key = Sys.getenv("CENSUS_API_KEY"),
                       cache = NULL,
                       verbose = TRUE) {

  geography <- match.arg(geography)
  if (is.null(variables)) variables <- pums_variables(geography)

  if (!is.null(cache) && file.exists(cache)) {
    say(verbose, "reading cached PUMS from ", cache)
    return(utils::read.csv(cache, stringsAsFactors = FALSE))
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("pums_fetch() needs the jsonlite package. ",
         'install.packages("jsonlite")', call. = FALSE)
  }
  if (!nzchar(key)) {
    stop("No Census API key. Request a free one at ",
         "https://api.census.gov/data/key_signup.html and set it with ",
         'Sys.setenv(CENSUS_API_KEY = "..."), preferably from ~/.Renviron so ',
         "it never lands in a script.", call. = FALSE)
  }

  ucgids <- if (geography == "region") {
    paste0("0200000US", 1:4)
  } else {
    paste0("0400000US", names(pums_fips()))
  }

  ## The API accepts several ucgid values per call for regions, but not 51
  ## states, so states are fetched one at a time.
  batches <- if (geography == "region") list(ucgids) else as.list(ucgids)

  get_one <- function(ids) {
    url <- paste0(
      "https://api.census.gov/data/", year, "/acs/acs1/pums?get=",
      paste(names(variables), collapse = ","),
      "&ucgid=", paste(ids, collapse = ","),
      "&key=", key
    )
    parsed <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      stop("Census API request failed for ", paste(ids, collapse = ", "),
           ": ", conditionMessage(e),
           "\nA non-JSON response usually means a bad key, an unavailable ",
           "year, or too many geographies in one call.", call. = FALSE)
    })
    out <- as.data.frame(parsed[-1, , drop = FALSE], stringsAsFactors = FALSE)
    colnames(out) <- parsed[1, ]
    out
  }

  say(verbose, "fetching ", length(batches), " request(s) from the ", year,
      " ACS 1-year PUMS")
  parts <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    parts[[i]] <- get_one(batches[[i]])
    if (length(batches) > 1 && i %% 10 == 0) {
      say(verbose, "  ", i, "/", length(batches))
    }
  }

  dt <- do.call(rbind, parts)

  ## The API appends geography columns whose names differ only in case from the
  ## requested ones; keep the first of each so renaming has a single target.
  dt <- dt[, !duplicated(toupper(colnames(dt))), drop = FALSE]

  ## Everything arrives as zero-padded character. Coerce the requested variables
  ## now, so that a comparison like HISP != 1 cannot silently misclassify.
  for (v in intersect(names(variables), colnames(dt))) {
    dt[[v]] <- suppressWarnings(as.numeric(dt[[v]]))
  }

  if (!is.null(cache)) {
    utils::write.csv(dt, cache, row.names = FALSE)
    say(verbose, "cached to ", cache)
  }
  dt
}

#' Recode raw ACS PUMS into RAILS analysis variables
#'
#' Applies the banding and labelling used in the RAILS application: age groups,
#' sex, race and ethnicity, income bands, educational attainment, tenure, and
#' Census region. Every result is a factor with fixed levels, so cell tables
#' built from different extracts line up.
#'
#' @param data Raw PUMS records from [pums_fetch()], or a data frame with the
#'   same ACS column names.
#' @param year The ACS vintage `data` came from. Only used to check the codings
#'   below are the right ones; see the note.
#' @param adults_only Keep only records with `AGEP > 17`. `TRUE` for a
#'   reference sample of adults. Set `FALSE` when you need complete population
#'   counts, for instance to compute state population shares.
#'
#' @return A data frame with one row per retained person and columns `weight`,
#'   `age`, `agegroup`, `sex`, `race_eth`, `income`, `edu`, `homeown`, `region`,
#'   and, when `ST` was fetched, `state` as a two-letter USPS code. Demographic
#'   columns may be `NA` where a code fell outside the bands.
#'
#' @section Codings are vintage-specific:
#' The `SCHL`, `TEN`, `SEX`, `HISP`, `RAC*` and `HINCP` bandings follow the
#' **2022** ACS 1-year PUMS data dictionary. Category codes do change between
#' releases, so a different `year` raises a warning: check each mapping against
#' that year's dictionary at
#' \url{https://www.census.gov/programs-surveys/acs/microdata/documentation.html}
#' before trusting the output.
#'
#' Income bands are `<35k`, `35k-50k`, `50k-75k`, `75k-100k`, `>100k`; note that
#' `HINCP` can be negative, and values below -60000 are treated as missing rather
#' than as the lowest band.
#'
#' @seealso [pums_fetch()], [pums_reference()].
#'
#' @examples
#' # A handful of synthetic records in raw ACS coding:
#' raw <- data.frame(
#'   PWGTP = c(100, 250, 90), AGEP = c(30, 70, 16), SEX = c(1, 2, 2),
#'   HINCP = c(42000, 120000, 20000), SCHL = c(21, 16, 12), TEN = c(1, 3, 2),
#'   HISP = c(1, 5, 1), RACWHT = c(1, 0, 0), RACBLK = c(0, 0, 1),
#'   RACASN = c(0, 0, 0), REGION = c(1, 4, 3)
#' )
#' pums_recode(raw)
#'
#' @export
pums_recode <- function(data, year = 2022, adults_only = TRUE) {

  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!identical(as.integer(year), 2022L)) {
    warning("pums_recode(): the codings here follow the 2022 ACS PUMS data ",
            "dictionary, but year = ", year, " was given. Category codes ",
            "change between vintages -- verify SCHL, TEN, HISP and HINCP ",
            "against that year's dictionary before using the result.",
            call. = FALSE)
  }

  data <- data[, !duplicated(toupper(colnames(data))), drop = FALSE]

  need <- c("PWGTP", "AGEP", "SEX", "HINCP", "SCHL", "TEN", "HISP", "REGION")
  missing_cols <- setdiff(need, colnames(data))
  if (length(missing_cols)) {
    stop("`data` is missing the ACS variable(s): ",
         paste(missing_cols, collapse = ", "),
         ". Fetch them with pums_fetch() or see pums_variables().",
         call. = FALSE)
  }

  num <- function(v) if (v %in% colnames(data)) {
    suppressWarnings(as.numeric(data[[v]]))
  } else {
    rep(NA_real_, nrow(data))
  }

  age <- num("AGEP")
  keep <- if (isTRUE(adults_only)) !is.na(age) & age > 17 else rep(TRUE, nrow(data))

  pick <- function(v) num(v)[keep]
  age  <- age[keep]

  hisp   <- pick("HISP")
  white  <- pick("RACWHT")
  black  <- pick("RACBLK")
  asian  <- pick("RACASN")
  inc    <- pick("HINCP")
  schl   <- pick("SCHL")
  ten    <- pick("TEN")
  reg    <- pick("REGION")
  sexcode <- pick("SEX")

  ## Ethnicity, then race, then the combined variable. HISP == 1 is
  ## "Not Spanish/Hispanic/Latino"; every other code is Hispanic origin.
  eth  <- ifelse(is.na(hisp), NA_character_,
                 ifelse(hisp != 1, "Hispanic", "Non-Hispanic"))
  race <- ifelse(!is.na(white) & white == 1, "White",
          ifelse(!is.na(black) & black == 1, "Black",
          ifelse(!is.na(asian) & asian == 1, "Asian", "Others")))

  race_eth <- ifelse(is.na(eth), NA_character_,
              ifelse(eth == "Hispanic", "Hispanic",
              ifelse(race == "White", "NH White",
              ifelse(race == "Black", "NH Black",
              ifelse(race == "Asian", "NH Asian", "Others")))))

  sex <- ifelse(is.na(sexcode), NA_character_,
                ifelse(sexcode == 1, "Male",
                       ifelse(sexcode == 2, "Female", NA_character_)))

  income <- cut(inc,
                breaks = c(-60000, 35000, 50000, 75000, 100000, Inf),
                labels = c("<35k", "35k-50k", "50k-75k", "75k-100k", ">100k"),
                right = FALSE)

  edu <- ifelse(is.na(schl), NA_character_,
         ifelse(schl < 12, "Less than highschool",
         ifelse(schl < 16, "Some highschool",
         ifelse(schl <= 17, "Highschool graduate",
         ifelse(schl < 21, "Some college", "College graduate or advanced")))))

  homeown <- ifelse(is.na(ten), NA_character_,
             ifelse(ten %in% c(1, 2), "Own",
             ifelse(ten == 3, "Rent",
             ifelse(ten == 4, "Others", NA_character_))))

  region <- ifelse(is.na(reg), NA_character_,
                   c("Northeast", "Midwest", "South", "West")[reg])

  agegroup <- ifelse(age <= 24, "18-24",
              ifelse(age <= 44, "25-44",
              ifelse(age <= 64, "45-64",
              ifelse(age <= 74, "65-74", "75+"))))

  out <- data.frame(
    weight   = pick("PWGTP"),
    age      = age,
    agegroup = factor(agegroup, levels = c("18-24", "25-44", "45-64",
                                           "65-74", "75+")),
    sex      = factor(sex, levels = c("Female", "Male")),
    race_eth = factor(race_eth, levels = c("Hispanic", "NH Asian", "NH Black",
                                           "NH White", "Others")),
    income   = income,
    edu      = factor(edu, levels = c("Less than highschool", "Some highschool",
                                      "Highschool graduate", "Some college",
                                      "College graduate or advanced")),
    homeown  = factor(homeown, levels = c("Own", "Rent", "Others")),
    region   = factor(region, levels = c("Northeast", "Midwest", "South",
                                         "West")),
    stringsAsFactors = FALSE
  )

  if ("ST" %in% colnames(data)) {
    fips <- pums_fips()
    out$state <- unname(fips[sprintf("%02d", pick("ST"))])
    if (anyNA(out$state)) {
      warning("pums_recode(): ", sum(is.na(out$state)),
              " record(s) have a state code outside the 50 states and DC.",
              call. = FALSE)
    }
  }

  rownames(out) <- NULL
  out
}

#' Build a RAILS reference sample from recoded PUMS
#'
#' Takes recoded PUMS microdata and produces the three things RAILS needs from a
#' reference sample: cleaned person-level records, the aggregated cell table,
#' and the population margins.
#'
#' Records with a missing value in any of `vars` cannot enter a cell, so they are
#' dropped -- but dropping them would also lose their share of the population.
#' The weights are therefore rescaled so the retained records still sum to the
#' population total computed *before* the drop. This is the step most easily got
#' wrong by hand, and it is why the totals are worth taking from here rather than
#' from an ad-hoc aggregation.
#'
#' @param data Recoded PUMS from [pums_recode()], or any data frame with a
#'   `weight` column and the columns in `vars`.
#' @param vars Covariates to aggregate over. Defaults to the seven used in the
#'   application.
#' @param order Highest interaction order to compute margins for. Must be at
#'   least the `scope` you will pass to [rails()].
#' @param rescale Rescale the weights of the retained records to the population
#'   total before incomplete records were dropped. Leave `TRUE` unless you have a
#'   reason to want the reference sample to under-count.
#'
#' @return An object of class `pums_reference`: a list with `micro` (the
#'   retained person-level records), `cells` (the aggregated table, ready for
#'   `rails(..., aggregated = TRUE)`), `totals` (the named margin vector for
#'   `pop_totals`), `n_pop` (the population total the weights sum to), `n_dropped`
#'   and `vars`.
#'
#' @seealso [pums_fetch()], [pums_recode()], [rails()].
#'
#' @examples
#' # Standing in for real PUMS, with the same column names and levels:
#' pop <- rails_simulate(5000, seed = 1)
#' fake <- data.frame(weight = 1 / pop$ps, pop[, c("agegroup", "sex", "income")])
#'
#' ref <- pums_reference(fake, vars = c("agegroup", "sex", "income"), order = 2)
#' ref
#' head(ref$cells)
#'
#' @export
pums_reference <- function(data,
                           vars = c("agegroup", "sex", "edu", "homeown",
                                    "income", "race_eth", "region"),
                           order = 3,
                           rescale = TRUE) {

  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  check_cells(data, vars, "`data`")

  w <- as.numeric(data$weight)
  n_pop <- sum(w, na.rm = TRUE)

  complete <- stats::complete.cases(data[, vars, drop = FALSE]) & !is.na(w)
  micro <- data[complete, , drop = FALSE]
  n_dropped <- sum(!complete)

  if (!nrow(micro)) {
    stop("Every record has a missing value in `vars`; nothing to aggregate.",
         call. = FALSE)
  }

  if (isTRUE(rescale)) {
    micro$weight <- micro$weight / sum(micro$weight) * n_pop
  } else {
    n_pop <- sum(micro$weight)
  }

  cells  <- as.data.frame(rails_cells(micro, vars, weights = "weight"))
  totals <- rails_totals(cells, vars = vars, order = order)

  out <- list(micro = micro, cells = cells, totals = totals,
              n_pop = n_pop, n_dropped = n_dropped, vars = vars, order = order)
  class(out) <- "pums_reference"
  out
}

#' @export
print.pums_reference <- function(x, ...) {
  cat("PUMS reference sample\n")
  cat("  covariates : ", paste(x$vars, collapse = ", "), "\n", sep = "")
  cat("  records    : ", format(nrow(x$micro), big.mark = ","), sep = "")
  if (x$n_dropped) {
    cat(" (", format(x$n_dropped, big.mark = ","),
        " dropped for missing covariates; weights rescaled)", sep = "")
  }
  cat("\n")
  cat("  cells      : ", nrow(x$cells), "\n", sep = "")
  cat("  margins    : ", length(x$totals), " up to order ", x$order, "\n",
      sep = "")
  cat("  population : ", format(round(x$n_pop), big.mark = ","), "\n", sep = "")
  invisible(x)
}
