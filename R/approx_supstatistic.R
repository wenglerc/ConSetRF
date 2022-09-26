# library(dplyr)
# library(stats)
# library(pracma)


#' Gaussian Kinematic Formular für eindimensionale t-verteilte Zufallsfelder
#'
#' Berechnung der erwarteteten Eulercharakteristik per GKF für das Excursion
#' Set \eqn{S_0 = {\mu(s) > level}} eines eindimensionalen Zufallsfeldes
#' zu einem fixen Level.
#' Die Grundmenge des Zufallsfeldes muss zusammenhängend sein,
#' also eine Eulercharakteristik von 1 besitzen
#'
#' @param Set Grundmenge des Zufallsfeldes.
#' @param X Realisierungen des Zufallsfeldes, wobei in den Zeilen die Werte
#'     pro Punkt der Grundmenge und in den Spalten die Werte pro Realisierung
#'     enthalten sind.
#' @param threshold Level des Exkursion Sets.
#'
#' @return Die erwartetete Eulercharakteristik des Excursion Sets \eqn{S_0}.
#'
#' @export
#'
tGKF <- function(Set, X, threshold) {

    # EC Dichten p_0 und p_1
    degree = ncol(X) - 1
    p0 <- 1 - stats::pt(threshold, degree)
    p1 <- (2 *pi)^(-1) * (1 + threshold ^ 2 / degree)^(-(degree - 1) / 2)

    # LKC_0
    L0 <- 1

    # LKC_1 gemaess Schaetzung
    residuals <-
        (X - rowSums(X) / ncol(X)) / apply(X, 1, sd)

    if (length(Set) > 1) {
        df.residuals <- apply(residuals, 2, function(values) {
            func <- stats::splinefun(Set, values, method = "natural")
            pracma::fderiv(func, Set, method = "central")
        }) %>% apply(1, stats::sd)
        L1 <- # Trapezregel
            sum(diff(Set) / 2 * (df.residuals[-length(df.residuals)] + df.residuals[-1]))
    } else {
        L1 <- 0
    }

    # Bestimme die EEC
    min(p0 * L0 + p1 * L1, 1)
}




#' Multiplier Bootstrap fuer das Supremum eines t-verteilten Zufallsfeldes
#'
#' Der Multiplier Bootstrap stellt die Datenbasis zur Approximation
#' der Verteilung des Supremums der T-Statistik eines eindimensionalen
#' Zufallsfeldes her.
#' Hierzu werden Residuen der Realisierungen des Zufallsfeldes genommen und mit
#' standardnormalverteilten Zufallsvariablen verrechnet, um die Verteilung
#' zu simulieren.
#'
#' @param X Realisierungen des Zufallsfeldes, wobei in den Zeilen die Werte
#'     pro Punkt der Grundmenge und in den Spalten die Werte pro Realisierung
#'     enthalten sind.
#' @param iter Anzahl der Iterationen zur Approximation der Verteilung.
#'
#' @return Matrix mit den gesampelten Daten, wobei in den Zeilen die Werte
#'     pro Punkt der Grundmenge und in den Spalten die Werte pro Iteration
#'     enthalten sind.
#'
#' @export
#'
mBoot <- function(X, iter = 1000) {

    S.size = nrow(X)
    samplesize = ncol(X)

    # Residuen
    residuals <-
        sqrt(samplesize / (samplesize - 1)) * (X - rowSums(X) / samplesize) %>%
        as.matrix()
    # Standardnormalverteilte Multiplier
    g <- stats::rnorm(samplesize * iter) %>%
        array(c(samplesize, iter, S.size)) %>%
        aperm(c(1, 3, 2))
    # Anwendung der Multiplier
    g.times.res <- array(residuals, dim = c(dim(residuals), iter)) %>%
        aperm(perm = c(2,1,3)) * g
    # Standardabweichung
    means <- g.times.res %>% colSums() %>% array(dim = c(dim(.), samplesize)) %>%
        aperm(perm = c(3,1,2)) /samplesize
    sd <- (g.times.res - means)^2 %>% colSums() %>% sqrt() %>%
        array(c(dim(.), samplesize)) %>% aperm(c(3, 1, 2)) * sqrt(samplesize - 1)^(-1)
    # Summen gewichteter Residuen über die Samples
    mboot <- ( g.times.res / sd / sqrt(samplesize)) %>% colSums()

    mboot
}

