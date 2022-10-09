#' Konfidenz-Obermenge eines Excursion Sets
#'
#' Berechnung einer Konfidenz-Obermenge der Menge
#' \eqn{{s: \mu(s) = level}} bzw. \eqn{{s: \mu(s) > level}} bzgl. eines Zufallsfelds
#' \eqn{X(s) = \mu(s) + \sigma(s)\epsilon(s)} zu einer kompakten Grundmenge S mit
#' Eulercharakteristik 1.
#' Die Berechnung basiert auf der Metode aus
#' "False Discovery Control for Random Fields" (2004) - M. Perone Pacifico,
#' C. Genovese, I. Verdinelli and L. Wasserman
#'
#' @param data Realisierungen des Zufallsfeldes.
#' @param S Grundmenge S (1D und EC(S) = 1).
#' @param partitions Partitionierung der Grundmenge (degenerierend und verfeinernd).
#' @param alpha Signifikanzniveau für die Konfidenzmenge.
#' @param level Grenzwert des Excursion Sets.
#' @param equal Wahrheitswert, ob einseitig (FALSE) oder zweiseitig (TRUE) getestet
#' werden soll. Beim einseitigen Test nimmt die Nullhypothese stets an, dass das
#' angegebene Level überschritten wird.
#' @param pmethod Approximationsmethode der p-Werte:
#' "tgkf" - Gaußsche Kinematische Formel für t-Verteilungen (default),
#' "mboot"- Multiplier Bootstrap.
#' @param mb.iter Anzahl der Iterationen für den Multiplier Bootstrap (default: 1000).
#'
#' @export
#'
#' @return Konfidenzmenge U mit \eqn{P(S_0 \subset U) \geq 1 - \alpha}
#'
ConfSet <- function(data,
                    S,
                    partitions,
                    alpha,
                    level,
                    equal = T,
                    pmethod = "tgkf",
                    mb.iter = 1000) {
    # Initialisierungen ----
    samplesize <- ncol(data)

    if (equal){
        t.statistic <-
            abs(sqrt(samplesize) * (rowSums(data) / samplesize - level) /
            apply(data, 1, stats::sd))
    } else {
        t.statistic <-
            -sqrt(samplesize) * (rowSums(data) / samplesize - level) /
            apply(data, 1, stats::sd)
    }

    if (pmethod == "mboot") mboot <- mBoot(data, iter = mb.iter) %>% t()

    # Berechnungen ----
    sets.per.niveau <-
        lapply(alpha, function(niveau) {
            U <- c()

            # Aussere Schleife: Loop über die Folge an Partitionen ----
            for (partS in partitions) {
                # Reduktion der aktuellen Partition auf jene Teilmengen,
                # die noch nicht in U enthalten sind: partS_n \ U_{n-1}
                partS <-
                    lapply(partS, function(vec) if (!all(vec %in% U)) vec) %>%
                    Filter(Negate(is.null), .)

                if (all(sapply(partS, is.null))) break

                N <- length(partS)

                # Innere Schleife: Berechnung von U_n nach Partition partS^*_n ----

                # Berechnung und Sortierung der realisierten Werte der Teststatistik
                teststatistics <- vapply(partS,
                                         function(vec)
                                             t.statistic[S %in% vec] %>% max(),
                                         numeric(1))
                sortorder <- order(teststatistics, decreasing = T)
                teststatistics <-
                    sort(teststatistics, decreasing = T)

                # Vereinigungen V_k ab dem k-ten Partitionselement (geordnet wie Teststatisik)
                sortedUnions <-
                    lapply(1:N, function(i) {
                        partS[sortorder[i:N]] %>% unlist(use.names = F) %>% sort()
                    })

                # Berechnung der p-Werte
                accept.index = 1
                if (pmethod == "tgkf") {
                    for (i in N:1) {
                        pvalue <- tGKF(sortedUnions[[i]],
                                       data[S %in% sortedUnions[[i]],],
                                       threshold = teststatistics[i])
                        # Akzeptiere V_k, wenn p(x_k, V_k) >= alpha
                        if (pvalue < alpha) {
                            accept.index <- ifelse(i == N, 0, i + 1)
                            break
                        }
                    }

                } else if (pmethod == "mboot") {
                    for (i in N:1) {
                        mboot.max <- mboot[, S %in% sortedUnions[[i]]] %>%
                            data.frame() %>% do.call(pmax, .)
                        pvalue <-
                            length(mboot.max[mboot.max >= teststatistics[i]]) /
                            length(mboot.max)
                        # Akzeptiere V_k, wenn p(x_k, V_k) >= alpha
                        if (pvalue < alpha) {
                            accept.index <- ifelse(i == N, 0, i + 1)
                            break
                        }
                    }
                }

                # Erweitere U_n = U_{n-1} + V_k
                if (accept.index > 0)
                    U <- base::union(U, sortedUnions[[accept.index]])
            }

            sort(U)
        })

    if(length(sets.per.niveau) == 1){
        return(sets.per.niveau[[1]])
    } else { return(sets.per.niveau) }

}
