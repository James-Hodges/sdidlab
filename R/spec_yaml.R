#' Read and write specifications as YAML
#'
#' `write_spec_yaml()` writes one specification; `write_specs_yaml()` writes
#' a *sweep* -- several specifications in one file, with the settings they
#' share hoisted into a `defaults:` block so the file shows exactly what
#' differs between them. `read_specs()` reads either form (or every YAML
#' file in a directory) and `read_spec_yaml()` reads a single-spec file.
#'
#' Writing then reading is an identity: the specifications that come back
#' have the same [spec_hash()] as those written.
#'
#' A sweep file looks like this:
#'
#' ```yaml
#' defaults:
#'   design: {treat_time: 1989}
#'   inference: {replications: 100}
#' specs:
#' - label: SDID
#' - label: Synthetic control
#'   estimator: sc
#' ```
#'
#' @param spec An [sdid_spec()].
#' @param specs A list of [sdid_spec()]s.
#' @param path File to write, or `NULL` to return the text instead. For
#'   `read_specs()`, a file or a directory.
#' @param title Comment line at the top of a sweep file.
#'
#' @return The writers return the path invisibly (or the YAML text when
#'   `path = NULL`); `read_specs()` returns a list of specifications and
#'   `read_spec_yaml()` a single one.
#'
#' @examples
#' specs <- list(
#'   sdid_spec(label = "SDID", design = list(treat_time = 1989)),
#'   sdid_spec(label = "SC", estimator = "sc", design = list(treat_time = 1989))
#' )
#' cat(write_specs_yaml(specs, NULL, title = "Prop 99 estimators"))
#'
#' f <- tempfile(fileext = ".yml")
#' write_specs_yaml(specs, f)
#' back <- read_specs(f)
#' identical(vapply(back, spec_hash, ""), vapply(specs, spec_hash, ""))
#' @name spec_yaml
NULL

#' @rdname spec_yaml
#' @export
write_spec_yaml <- function(spec, path) {
  body <- yaml::as.yaml(spec_to_plain(spec), indent = 2L)
  text <- paste0("# sdidlab specification\n# spec_hash: ", spec_hash(spec), "\n", body)
  if (is.null(path)) return(text)
  writeLines(text, path)
  invisible(path)
}

#' @rdname spec_yaml
#' @export
read_spec_yaml <- function(path) {
  if (!file.exists(path)) cli::cli_abort("No spec file at {.file {path}}.")
  raw <- yaml::read_yaml(path)
  if (!is.null(raw$specs)) {
    cli::cli_abort("{.file {basename(path)}} holds several specifications; use {.fn read_specs}.")
  }
  do.call(sdid_spec, raw[setdiff(names(raw), "schema_version")])
}

#' @rdname spec_yaml
#' @export
read_specs <- function(path) {
  if (dir.exists(path)) {
    files <- sort(list.files(path, pattern = "[.]ya?ml$", full.names = TRUE))
    return(unlist(lapply(files, read_specs), recursive = FALSE))
  }
  if (!file.exists(path)) cli::cli_abort("No spec file at {.file {path}}.")
  raw <- yaml::read_yaml(path)
  if (is.null(raw$specs)) return(list(read_spec_yaml(path)))

  defaults <- raw$defaults %||% list()
  specs <- lapply(seq_along(raw$specs), function(i) {
    merged <- merge_spec(merge_spec(spec_defaults(), defaults), raw$specs[[i]])
    tryCatch(
      do.call(sdid_spec, merged[setdiff(names(merged), "schema_version")]),
      error = function(e) {
        cli::cli_abort(c("Spec {i} in {.file {basename(path)}} is invalid.",
                         x = conditionMessage(e)), call = NULL)
      }
    )
  })
  labels <- vapply(specs, function(s) s$label, character(1))
  dupes <- unique(labels[duplicated(labels) & !is.na(labels)])
  if (length(dupes)) {
    cli::cli_abort(c("Duplicate labels in {.file {basename(path)}}.", x = "{.val {dupes}}"))
  }
  specs
}

#' @rdname spec_yaml
#' @export
write_specs_yaml <- function(specs, path, title = "sdidlab specifications") {
  if (inherits(specs, "sdid_spec")) specs <- list(specs)
  plain <- lapply(specs, spec_to_plain)
  base <- spec_to_plain(spec_defaults())
  shared <- shared_fields(plain, base)
  body <- list()
  if (length(shared)) body$defaults <- shared
  body$specs <- lapply(plain, function(p) strip_shared(p, shared, base))
  text <- paste0(
    "# ", title, "\n",
    "# Placebo standard errors depend on the RNG stream: fit these in file\n",
    "# order under one seed, or with run_specs(seed_each = TRUE).\n",
    yaml::as.yaml(body, indent = 2L)
  )
  if (is.null(path)) return(text)
  writeLines(text, path)
  invisible(path)
}

# YAML cannot tell a zero-length vector from an empty map, so empty vectors
# become list() consistently -- which is also what reading them back gives.
spec_to_plain <- function(spec) {
  plain <- unclass(spec)
  attr(plain, "spec_hash") <- NULL
  rapply(plain, how = "replace", f = function(x) {
    if (is.list(x)) return(x)
    if (!length(x)) return(list())
    x
  })
}

# A field is shared if every spec has the same value and that value differs
# from the default (a field left at its default is simply absent).
shared_fields <- function(plain, base) {
  out <- list()
  for (nm in names(base)) {
    if (nm %in% c("id", "label", "note", "tags", "schema_version")) next
    vals <- lapply(plain, function(p) p[[nm]])
    if (is.list(base[[nm]]) && !is.null(names(base[[nm]])) && !identical(nm, "options")) {
      sub <- shared_fields(lapply(vals, function(v) v %||% list()), base[[nm]])
      if (length(sub)) out[[nm]] <- sub
    } else {
      first <- vals[[1]]
      if (all(vapply(vals, function(v) identical(v, first), logical(1))) &&
          !identical(first, base[[nm]])) {
        out[[nm]] <- first
      }
    }
  }
  out
}

strip_shared <- function(p, shared, base) {
  keep <- list()
  for (nm in names(p)) {
    if (nm == "schema_version") next
    v <- p[[nm]]
    if (nm %in% c("id", "label", "note", "tags")) {
      if (length(v) && !identical(v, base[[nm]]) && !all(is.na(unlist(v)))) keep[[nm]] <- v
      next
    }
    s <- shared[[nm]]
    if (is.list(v) && !is.null(names(v)) && !identical(nm, "options")) {
      d <- base[[nm]]
      sub <- list()
      for (k in names(v)) {
        if (!is.null(s[[k]]) && identical(s[[k]], v[[k]])) next
        if (is.null(s[[k]]) && identical(d[[k]], v[[k]])) next
        sub[[k]] <- v[[k]]
      }
      if (length(sub)) keep[[nm]] <- sub
    } else {
      if (!is.null(s) && identical(s, v)) next
      if (is.null(s) && identical(v, base[[nm]])) next
      keep[[nm]] <- v
    }
  }
  keep
}
