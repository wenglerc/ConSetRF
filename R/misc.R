#' Degenerative Partitionsfolge
#'
#' Erzeugt eine degenerative Folge von Partitionen, die immer feiner werden,
#' in dem die Mengen stets in eine fixe Anzahl an gleichgroßen Teilmengen
#' geteilt wird.
#'
#' @param S Grundmenge, die partitioniert werden soll.
#' @param pieces Anzahl der Teile, in die gesplittet werden soll.
#' @param max_iterations maximale Anzahl an Folgengliedern.
#'
#' @return Eine Liste je Folgenglied. Jede dieser Listen enthält die
#' immer feiner werdende Aufspaltungen der Menge S.
#'
#' @export
#'
#' @examples
#' partition_seq( seq(0, 1, length.out = 100), pieces = 3 )
#'
partition_seq <- function(S, pieces = 2, max_iterations = 1000) {
    iterations <- min(length(S), max_iterations)

    result <- lapply(1:(log(iterations, base = pieces) + 1), function(iter) {
        split(S, cut(seq_along(S),
                     # Intervalle in pieces Teile teilen
                     breaks = pieces ^ iter,
                     labels = FALSE))
    })

    return(result)
}

