# Download ACS PUMS microdata from the Census API

Fetches the person-level PUMS records RAILS uses as its reference
sample. The result is raw ACS codes; pass it to
[`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md)
to get analysis variables.

## Usage

``` r
pums_fetch(
  year = 2022,
  geography = c("region", "state"),
  variables = NULL,
  key = Sys.getenv("CENSUS_API_KEY"),
  cache = NULL,
  verbose = TRUE
)
```

## Arguments

- year:

  ACS 1-year vintage, e.g. `2022`.

- geography:

  `"region"` fetches the four Census regions in one request. `"state"`
  fetches the 50 states plus DC, one request each – the API caps how
  many `ucgid` values a single call accepts, so asking for all 51 at
  once returns an error page rather than data.

- variables:

  Named character vector of ACS variable codes; see
  [`pums_variables()`](https://extrasane.github.io/RAILSpkg/reference/pums_variables.md).

- key:

  A Census API key. Free, from
  <https://api.census.gov/data/key_signup.html>. Defaults to the
  `CENSUS_API_KEY` environment variable – set it in `~/.Renviron` rather
  than writing the key into a script, since scripts get committed.

- cache:

  Optional path to a CSV. If the file exists it is read and no request
  is made; otherwise the download is written there. Repeat runs of an
  analysis should always set this: the state form is 51 API calls.

- verbose:

  Report progress via
  [`message()`](https://rdrr.io/r/base/message.html).

## Value

A data frame of raw PUMS records, one row per person, with the requested
variables as numeric columns.

## Network access

This function requires an internet connection and the jsonlite package.
Nothing else in RAILS touches the network.

## See also

[`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md)
to turn the result into analysis variables, then
[`pums_reference()`](https://extrasane.github.io/RAILSpkg/reference/pums_reference.md)
to get cells and margins.

## Examples

``` r
if (FALSE) { # \dontrun{
# Sys.setenv(CENSUS_API_KEY = "...")   # better: put it in ~/.Renviron
raw <- pums_fetch(2022, cache = "PUMS_2022.csv")
ref <- pums_reference(pums_recode(raw))
} # }
```
