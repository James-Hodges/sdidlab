#' A specification's stable identifier
#'
#' Twelve hex characters of a SHA-256 over the specification's canonical
#' JSON. Metadata (`id`, `label`, `note`, `tags`) and exclusion notes are
#' excluded, so relabelling a specification does not change its identity,
#' and any field still at its default is dropped before hashing, so adding a
#' new field with a default leaves every existing hash unchanged. Field
#' order never matters.
#'
#' @param spec An [sdid_spec()].
#'
#' @return A 12-character string.
#'
#' @examples
#' a <- sdid_spec(design = list(treat_time = 1989), label = "first name")
#' b <- sdid_spec(label = "second name", design = list(treat_time = 1989L))
#' identical(spec_hash(a), spec_hash(b))
#'
#' @export
spec_hash <- function(spec) {
  cached <- attr(spec, "spec_hash")
  if (!is.null(cached)) return(cached)
  spec_hash_impl(spec)
}

SPEC_HASH_EXCLUDE <- c("schema_version", "id", "label", "note", "tags")

spec_hash_impl <- function(spec) {
  payload <- unclass(spec)[setdiff(names(spec), SPEC_HASH_EXCLUDE)]
  attr(payload, "spec_hash") <- NULL
  payload$pool$exclusion_notes <- NULL

  defaults <- spec_defaults()[setdiff(names(spec_defaults()), SPEC_HASH_EXCLUDE)]
  defaults$pool$exclusion_notes <- NULL

  canonical <- canonicalise(drop_defaults(payload, defaults))
  json <- jsonlite::toJSON(canonical, auto_unbox = TRUE, digits = NA,
                           null = "null", na = "null")
  substr(digest::digest(as.character(json), algo = "sha256", serialize = FALSE), 1L, 12L)
}

# Recursively remove anything identical to its default. A section left
# entirely at its defaults disappears.
drop_defaults <- function(x, default) {
  if (!is.list(x) || !is.list(default)) {
    return(if (identical(x, default)) NULL else x)
  }
  out <- list()
  for (nm in names(x)) {
    kept <- if (nm %in% names(default)) drop_defaults(x[[nm]], default[[nm]]) else x[[nm]]
    if (!is.null(kept) && length(kept)) out[[nm]] <- kept
  }
  if (!length(out)) NULL else out
}

# Sort names at every level so construction order cannot leak into the hash.
canonicalise <- function(x) {
  if (!is.list(x)) return(x)
  x <- lapply(x, canonicalise)
  if (is.null(names(x))) return(x)
  x[order(names(x))]
}
