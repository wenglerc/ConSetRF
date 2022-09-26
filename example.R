library(doParallel)     # Parallelisierung aufwaendiger Schleifen
library(plyr)
library(dplyr)          # Pipe-Operator
library(reshape2)       # Umformen von Dataframes

library(mvtnorm)        # Multivariate Normalverteilungen

library(ggplot2)        # Plotten
theme_set(theme_light())

library(tictoc)


# Generieren eines Zufallsprozesses X(s) = mu(s) + sigma(s)*epsilon(s) --------

mu1 <- function(S) {
    sin(S * 5 * pi) * 5 * S
} # Erwartungswertfunktion

sigma1 <- function(S){
    rep(1, length(S))
} # Varianzfunktion

eps1 <- function(S) {
    outer(S, S, FUN = function(x, y) exp( -( x - y )^2 / 0.0005 ))
} # Fehlerprozess

linearModell <- function(n, S, mu, sigma, eps) {
    mu(S) + sigma(S) * t(mvtnorm::rmvnorm(n, sigma = eps(S)))
} # Lineares Modell

# Initialisierung globaler Variablen -----------------------------------

randomseed <- 1
samplesize.max <- 500
iterations <- 500
gridpoints <- 200               # Anzahl equidistanter Gitterpunkte

alpha <- 0.05
level <- 0                      # Grenzwert aus der Null-Hypothese
S <- seq(0, 1, length.out = gridpoints)  # Grundmenge S = [0,1]
partitions <- partition_seq(S, pieces = 4)  # Folge von Partitionen (degenerierend, feiner werdend)


# Testdaten --------------------------------------------------

set.seed(randomseed)

# Realisierungen von X
cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)
data.all <- foreach(i = 1:iterations) %dopar% {
    data.frame(linearModell(samplesize.max, S, mu1, sigma1, eps1))
}
stopCluster(cluster)
rm(cluster)

data.mu <- data.frame(S = S, Erwartungswert = mu(S))
S0 = S[mu(S) > level]

# Darstellung von Realisierungen
data.plot <- reshape2::melt(
    data.frame(S = S, data.all[[1]])[, 1:10],
    id.vars = "S",
    variable.name = 'Sample',
    value.name = "Wert"
)
p.ex <- ggplot() +
    labs(x = "Grundmenge S", y = "", colour = "Realisierung") +
    theme(legend.position = "none")
p.ex + geom_line(data = data.plot, aes(S, Wert, colour = Sample), alpha = 0.5) +
    geom_line(data = data.mu, aes(S, Erwartungswert), colour = "red")

# Berechnungen --------------------------------------------------

# U <- ConfSet(data.all[[4]], S, partitions, alpha, level, pmethod = "mboot")
# cat("S0 nicht in U: ", length( S0[ !(S0 %in% U)] ) / length(S0),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0)] ) / length(U))
#
# U <- ConfSet(data.all[[4]], S, partitions, alpha, level, pmethod = "mboot", mb.iter = 500)
# cat("S0 nicht in U: ", length( S0[ !(S0 %in% U)] ) / length(S0),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0)] ) / length(U))
#
# mbm <- microbenchmark(
#     V1 = ConfSet(data.all[[1]], S, partitions, alpha, level, pmethod = "mboot"),
#     V2 = ConfSet(data.all[[1]], S, partitions, alpha, level, pmethod = "mboot", mb.iter = 500),
#     times = 10
# )
# autoplot(mbm)


## Berechnungen nach Samplegroesse und Signifikanzniveau --------------------

samplesize.list <- c(25, 50, 100, 250, 500)
names(samplesize.list) <- paste("N", samplesize.list, sep = "")

signif.list <- c(0.9, 0.95, 0.99)
names(signif.list) <- paste("A", signif.list, sep = "")


#### Ergebnisse GKF --------------------------------------------------

cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

# Obermengen je Samplegroesse
tic("GKF parallel")
results.gkf <-
    foreach(data = data.all, .packages = c("dplyr", "ConfSetRF")) %dopar% {

        l <- lapply(samplesize.list, function(N) {
            ConfSet(data[, 1:N], S, partitions, signif.list, level, pmethod = "tgkf")
        })
        unlist(l, recursive = F)

    }
stopCluster(cluster)
toc()


# Ueberdeckungswahrscheinlichkeit
covering.gkf <-
    sapply(results.gkf, function(samplelist) {
        sapply(samplelist, function(U)
            all(S0 %in% U))
    }) %>% t() %>% data.frame() %>%
    sapply(function(fullcov)
        length(fullcov[fullcov == T]) / iterations)
covering.gkf <- data.frame(
    value = covering.gkf,
    samplesize = rep(samplesize.list, each = length(signif.list)),
    niveau = rep(as.character(signif.list), times = length(samplesize.list))
)


#### Ergebnisse MB --------------------------------------------------

cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

# Obermengen je Samplegroesse
tic("MB parallel")
results.mb <-
    foreach(data = data.all, .packages = c("dplyr", "ConfSetRF")) %dopar% {

        l <- lapply(samplesize.list, function(N) {
            ConfSet(data[, 1:N], S, partitions, signif.list, level,
                    pmethod = "mboot", mb.iter = 500)
        })
        unlist(l, recursive = F)

    }
stopCluster(cluster)
toc()


# Ueberdeckungswahrscheinlichkeit
covering.mb <-
    sapply(results.mb, function(samplelist) {
        sapply(samplelist, function(U)
            all(S0 %in% U))
    }) %>% t() %>% data.frame() %>%
    sapply(function(fullcov)
        length(fullcov[fullcov == T]) / iterations)
covering.mb <- data.frame(
    value = covering.mb,
    samplesize = rep(samplesize.list, each = length(signif.list)),
    niveau = rep(as.character(signif.list), times = length(samplesize.list))
)


# Auswertung --------------------------------------------------

data.gkf <- data.frame(covering.gkf, method = "t-GKF")
data.mb <- data.frame(covering.mb, method = "Multiplier Bootstrap")
cov.plot <- rbind(data.gkf, data.mb)

p.samp <- ggplot(data = cov.plot,
                 aes(x = samplesize, y = value,
                     color = method, linetype = niveau)) +
    geom_line() + geom_point() +
    geom_hline(yintercept = 1 - alpha, col = "red4", size = 0.75,
               linetype = "dashed") +
    labs(x = "Anzahl an Samples", y = "Überdeckungsrate",
         color = "Methode", linetype = "Signifikanzniveau") +
    theme(legend.position = c(0.85, 0.3))
p.samp



