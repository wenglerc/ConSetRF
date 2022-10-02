library(doParallel)     # Parallelisierung aufwaendiger Schleifen
library(plyr)
library(dplyr)          # Pipe-Operator
library(reshape2)       # Umformen von Dataframes

library(ConfSetRF)

library(mvtnorm)        # Multivariate Normalverteilungen

library(ggplot2)        # Plotten
library(patchwork)
theme_set(theme_light(base_size = 14))
library(stringr)


library(tictoc)
library(microbenchmark)


# Funktionen Plot-Daten ---------------------------------------------------

cov.per.model <- function(results, covset, methodname){
    cov <-
        sapply(results, function(samplelist) {
            sapply(samplelist, function(U) if(all(covset %in% U)) 1 else 0)
        }) %>% t() %>% data.frame(model = gsub("\\.\\d+", "",rownames(.)))
    cov <- reshape2::melt(
        cov,
        id.vars = "model",
        variable.name = "samplesize"
    ) %>% mutate(method = methodname,
                 samplesize = as.numeric(gsub("N", "", samplesize)))

    cov
}

overcov.per.model <- function(results, covset, methodname){
    cov <-
        sapply(results, function(samplelist) {
            sapply(samplelist, function(U) length( U[!(U %in% covset)] ) / length(U))
        }) %>% t() %>% data.frame(model = gsub("\\.\\d+", "",rownames(.)))
    cov <- reshape2::melt(
        cov,
        id.vars = "model",
        variable.name = "samplesize"
    ) %>% mutate(method = methodname,
                 samplesize = as.numeric(gsub("N", "", samplesize)))

    cov
}

# Modell-Komponenten ----------------------------------------------------

mu1 <- function(S) {
    sin(S * 2 * pi)
} # Erwartungswertfunktion 1
mu2 <- function(S) {
    sin(S * 4 * pi) * 2 * S
} # Erwartungswertfunktion 2
mu3 <- function(S) {
    sin(S * 12 * pi) * (S - 0.5) * 2
} # Erwartungswertfunktion 3

sigma <- function(S) {
    rep(0.3, length(S))
} # Varianzfunktion

c1 <- function(S, h = 0.075) {
    outer(S, S, FUN = function(x, y) exp( -( x - y )^2 / (2*h^2) ))
} # Fehlerprozess 1

c2 <- function(S, a = 30, l = 0.5) {
    outer(S, S, FUN = function(x, y) (1 + ( x - y )^2 / 2*a*l^2)^(-a) )
} # Fehlerprozess 2


# Variablen zur Steuerung der Auswertung -----------------------------------

randomseed <- 1
samplesize.max <- 500
iterations <- 1000
gridpoints <- 200       # Anzahl equidistanter Gitterpunkte
alpha <- 0.05           # Signifikanzniveau

# Testdaten -----------------------------------

set.seed(randomseed)

S <- seq(0, 1, length.out = gridpoints)     # Grundmenge
partitions <- partition_seq(S, pieces = 4)  # Folge von Partitionen
level <- 0                                  # Grenzwert der Null-Hypothese

mu1.S <- mu1(S)
mu2.S <- mu2(S)
mu3.S <- mu3(S)
sigma.S <- sigma(S)
c1.S <- c1(S)
c2.S <- c2(S)

S0.1 <- S[mu1.S > level]
S0.2 <- S[mu2.S > level]
S0.3 <- S[mu3.S > level]

# Realisierungen von X
cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

tic("Daten generieren")
eps.all <- foreach(i = 1:(3 * iterations)) %dopar% {
    if (i <= iterations) {
        t(mvtnorm::rmvnorm(samplesize.max, sigma = c1.S))
    } else if (i <= (2 * iterations)) {
        t(mvtnorm::rmvt(samplesize.max, sigma = c1.S, df = 3))
    } else {
        t(mvtnorm::rmvt(samplesize.max, sigma = c2.S, df = 3))
    }
}

data1 <- foreach(i = 1:(2 * iterations)) %dopar% {
    if (i <= iterations)
        mu2.S + sigma.S * eps.all[[i]]
    else
        mu3.S + sigma.S * eps.all[[i - iterations]]
}
data2 <- foreach(eps = eps.all) %dopar% {
    mu1.S + sigma.S * eps
}
toc()

stopCluster(cluster)


names(eps.all) <- paste(rep(c("A.", "B.", "C."), each = iterations),
                        rep(1:iterations, times = 3),
                        sep = "")
names(data1) <- paste(rep(c("A2.", "A3."), each = iterations),
                      rep(1:iterations, times = 2),
                      sep = "")
names(data2) <- paste(rep(c("A1.", "B.", "C."), each = iterations),
                      rep(1:iterations, times = 3),
                      sep = "")

data.all <- c(data1, data2)

rm(data1, data2, cluster, c1.S, c2.S, sigma.S)
gc()


# Darstellung von Realisierungen der Modelle A,B,C

data.examp.list <- list(
    data.A1.examp = data.all[[which(grepl("A1", names(data.all)) == T)[1]]],
    data.B.examp = data.all[[which(grepl("B", names(data.all)) == T)[1]]],
    data.C.examp = data.all[[which(grepl("C", names(data.all)) == T)[1]]]
)
data.noise.list <- list(
    data.A1.examp = eps.all[[which(grepl("A", names(eps.all)) == T)[1]]],
    data.B.noise = eps.all[[which(grepl("B", names(eps.all)) == T)[1]]],
    data.C.noise = eps.all[[which(grepl("C", names(eps.all)) == T)[1]]]
)

max.examp <- max(sapply(data.examp.list, function(mat) max(mat[,1:10])))
min.examp <- min(sapply(data.examp.list, function(mat) min(mat[,1:10])))
max.noise <- max(sapply(data.noise.list, function(mat) max(mat[,1:10])))
min.noise <- min(sapply(data.noise.list, function(mat) min(mat[,1:10])))

p <- ggplot() + labs(x = "", y = "") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))

plots.examp <- vector( mode = "list", length(data.noise.list) )
plots.noise <- vector( mode = "list", length(data.noise.list) )
for ( i in 1:length(data.noise.list) ){
    data.examp <- reshape2::melt(
        data.frame(S = S, data.examp.list[[i]])[, 1:11],
        id.vars = "S", variable.name = 'Sample', value.name = "Wert"
    )
    data.noise <- reshape2::melt(
        data.frame(S = S, data.noise.list[[i]])[, 1:11],
        id.vars = "S", variable.name = 'Sample', value.name = "Wert"
    )

    plots.examp[[i]] <- p + ylim(min.examp, max.examp) +
        ggtitle(paste("Modell",
                      gsub(".*\\.(.+)\\..*", "\\1", names(data.examp.list)[i]))) +
        geom_line(data = data.examp, aes(S, Wert, colour = Sample, linetype = Sample),
                  alpha = 0.5) +
        geom_line(data = data.frame(S = S, Erwartungswert = mu1.S),
                  aes(S, Erwartungswert), colour = "red")

    plots.noise[[i]] <- p +  ylim(min.noise, max.noise) +
        ggtitle(paste("Fehlerprozess Modell",
                      gsub(".*\\.(.+)\\..*", "\\1", names(data.noise.list)[i]))) +
        geom_line(data = data.noise, aes(S, Wert, colour = Sample, linetype = Sample))
}

# (plots.examp[[1]] + plots.examp[[2]] + plots.examp[[3]]) /
#     (plots.noise[[1]] + plots.noise[[2]] + plots.noise[[3]])

# Darstellung von Realisierungen der Modelle A1,A2,A3

data.exampA.list <- list(
    data.A1.examp = data.all[[which(grepl("A1", names(data.all)) == T)[1]]],
    data.A2.examp = data.all[[which(grepl("A2", names(data.all)) == T)[1]]],
    data.A3.examp = data.all[[which(grepl("A3", names(data.all)) == T)[1]]]
)

max.exampA <- max(sapply(data.exampA.list, function(mat) max(mat[,1:10])))
min.exampA <- min(sapply(data.exampA.list, function(mat) min(mat[,1:10])))

plots.exampA <- vector( mode = "list", length(data.exampA.list) )
for ( i in 1:length(data.exampA.list) ){
    data.examp <- reshape2::melt(
        data.frame(S = S, data.exampA.list[[i]])[, 1:11],
        id.vars = "S", variable.name = 'Sample', value.name = "Wert"
    )

    plots.exampA[[i]] <- p + ylim(min.exampA, max.exampA) +
        ggtitle(paste("Modell",
                      gsub(".*\\.(.+)\\..*", "\\1", names(data.exampA.list)[i]))) +
        geom_line(data = data.examp, aes(S, Wert, colour = Sample, linetype = Sample),
                  alpha = 0.5)
}

# (plots.exampA[[1]] +
#         geom_line(data = data.frame(S = S, Erwartungswert = mu1.S),
#                   aes(S, Erwartungswert), colour = "red")) +
# (plots.exampA[[2]] +
#          geom_line(data = data.frame(S = S, Erwartungswert = mu2.S),
#                    aes(S, Erwartungswert), colour = "red")) +
# (plots.exampA[[3]] +
#          geom_line(data = data.frame(S = S, Erwartungswert = mu3.S),
#                    aes(S, Erwartungswert), colour = "red"))

# Speicherplatz frei machen
rm(
    eps.all, mu1.S, mu2.S, mu3.S,
    data.examp, data.noise, data.examp.list, data.noise.list, data.exampA.list,
    max.examp, min.examp, max.noise, min.noise, max.exampA, min.exampA,
    p, plots.examp, plots.noise, plots.exampA
)
gc()

# Berechnungen ----------------------------------------------------------------


# U <- ConfSet(data.all[[4]], S, partitions, alpha, level, pmethod = "mboot")
# cat("S0 nicht in U: ", length( S0[ !(S0 %in% U)] ) / length(S0),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0)] ) / length(U))
#
# U <- ConfSet(data.all[[4]], S, partitions, alpha, level, pmethod = "mboot", mb.iter = 500)
# cat("S0 nicht in U: ", length( S0[ !(S0 %in% U)] ) / length(S0),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0)] ) / length(U))
#
#
### Laufzeit -----------------------------------------------------------------
# mbm <- microbenchmark(
#     GKF = ConfSet(data.all[[1]][, 1:50], S, partitions, alpha, level,
#                   pmethod = "gkf"),
#     MB500 = ConfSet(data.all[[1]][, 1:50], S, partitions, alpha, level,
#                      pmethod = "mboot", mb.iter = 500),
#     MB1000 = ConfSet(data.all[[1]][, 1:50], S, partitions, alpha, level,
#                  pmethod = "mboot", mb.iter = 1000),
#     times = 100
# )
# autoplot(mbm)


## Berechnungen per Samplesize -----------------------------------------------

samplesize.list <- c(25, 50, 100, 250, 500)
names(samplesize.list) <- paste("N", samplesize.list, sep = "")

#### Ergebnisse GKF ---------------------------------------------------------

cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

# Obermengen je Samplegroesse
tic("GKF parallel")
results.gkf <-
    foreach(data = data.all, .packages = c("dplyr", "ConfSetRF")) %dopar% {

        lapply(samplesize.list, function(N) {
            ConfSet(data[, 1:N], S, partitions, alpha, level, pmethod = "tgkf")
        })
    }
stopCluster(cluster)
toc()

names(results.gkf) <- names(data.all)


# Coverage je Modell ABC
cov.gkf.ABC <-
    cov.per.model(results.gkf[-(1:(2*iterations))], S0.1, "t-GKF")
summarise.cov.gkf.ABC <- cov.gkf.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")
# Overcoverage je Modell ABC
overcov.gkf.ABC <-
    overcov.per.model(results.gkf[-(1:(2*iterations))], S0.1, "t-GKF")
summarise.overcov.gkf.ABC <- overcov.gkf.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")


# Coverage je Erwartungswertfunktion
cov.gkf.means <- rbind(
    cov.per.model(results.gkf[(2 * iterations + 1):(3 * iterations)], S0.1, "t-GKF"),
    cov.per.model(results.gkf[1:iterations], S0.2, "t-GKF"),
    cov.per.model(results.gkf[(iterations + 1):(2 * iterations)], S0.3, "t-GKF")
)
summarise.cov.gkf.means <- cov.gkf.means %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")
# Overcoverage je Erwartungswertfunktion
overcov.gkf.means <-rbind(
    overcov.per.model(results.gkf[(2 * iterations + 1):(3 * iterations)], S0.1, "t-GKF"),
    overcov.per.model(results.gkf[1:iterations], S0.2, "t-GKF"),
    overcov.per.model(results.gkf[(iterations + 1):(2 * iterations)], S0.3, "t-GKF")
)
summarise.overcov.gkf.means <- overcov.gkf.means %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")

rm(results.gkf) # Speicherplatz frei machen
gc()

#### Ergebnisse MB ---------------------------------------------------------

cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

# Obermengen je Samplegroesse
tic("MB parallel")
results.mb <-
    foreach(data = data.all[-(1:(2*iterations))], .packages = c("dplyr", "ConfSetRF")) %dopar% {

        lapply(samplesize.list, function(N) {
            ConfSet(data[, 1:N], S, partitions, alpha, level,
                    pmethod = "mboot", mb.iter = 500)
        })
        # unlist(l, recursive = F)

    }
stopCluster(cluster)
toc()

names(results.mb) <- names(data.all[-(1:(2*iterations))])

# Coverage je Modell ABC
cov.mb.ABC <-
    cov.per.model(results.mb, S0.1, "Multiplier Bootstrap")
summarise.cov.mb.ABC <- cov.mb.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")
# Overcoverage je Modell ABC
overcov.mb.ABC <-
    overcov.per.model(results.mb, S0.1, "Multiplier Bootstrap")
summarise.overcov.mb.ABC <- overcov.mb.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")

rm(results.mb) # Speicherplatz frei machen
gc()


# Auswertung ----------------------------------------------------------------
pdf(file = "Auswertungen.pdf", width = 6, height = 4)

## Coverage -----------------------------------------------------------------
cov.plot.ABC <- rbind(summarise.cov.gkf.ABC, summarise.cov.mb.ABC)

pcov <- ggplot(cov.plot.ABC, aes(x = samplesize)) +
    labs(x = "Anzahl an Samples", y = "Überdeckungsrate", color = "Modell",
         linetype = "Methode", shape = "Modell") +
    scale_linetype(labels = c("t-GKF", str_wrap("Multiplier Bootstrap", 5))) +
    theme(legend.justification = c(1,0.5), legend.position = c(0.999,0.5),
          legend.box.background = element_rect(fill = "white", color = "darkgray")) +
    geom_hline(yintercept = 1 - alpha, col = "red3", size = 0.75,
               linetype = "dashed")

# Coverage je Modell ABC
pcovABC <- pcov +
    geom_line(aes(y = value_mean, color = model,
               linetype = method)) +
    geom_point(aes(y = value_mean, color = model, shape = model))
pcovABC
# Coverage je Erwartungswertfunktion
pcovmean <- pcov +
    geom_line(data = summarise.cov.gkf.means,
              aes(y = value_mean, color = model)) +
    geom_point(data = summarise.cov.gkf.means,
               aes(y = value_mean, color = model, shape = model))
pcovmean

## Overcoverage -----------------------------------------------
overcov.plot.ABC <- rbind(summarise.overcov.gkf.ABC, summarise.overcov.mb.ABC)

povercov <- ggplot(overcov.plot.ABC, aes(x = samplesize)) +
    labs(x = "Anzahl an Samples",
         y = expression(paste("Anteil von ", U, "\\", S[0], " in ", U )),
         color = "Modell", linetype = "Methode", shape = "Modell") +
    scale_linetype(labels = c("t-GKF", str_wrap("Multiplier Bootstrap", 5))) +
    theme(legend.justification = c(1,1), legend.position = c(0.999,0.999),
          legend.box.background = element_rect(fill = "white", color = "darkgray"))+
    geom_hline(yintercept = alpha, col = "red3", size = 0.75,
               linetype = "dashed")
# Overcoverage je Modell ABC Durchschnittswerte
povercovABC <- povercov +
    geom_line(aes(y = value_mean, color = model, linetype = method)) +
    geom_point(aes(y = value_mean, color = model, shape = model))+
    ggtitle("Durchschnittliche Überschätzung") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovABC

# Overcoverage je Modell ABC GKF
povercovABC.gkf <- povercov +
    geom_boxplot(data = overcov.gkf.ABC,
                 aes(x = reorder(as.character(samplesize), samplesize),
                     y = value, color = model)) +
    ggtitle("Überschätzung mittels t-GKF") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovABC.gkf

# Overcoverage je Modell ABC MB
povercovABC.mb <- povercov +
    geom_boxplot(data = overcov.mb.ABC,
                 aes(x = reorder(as.character(samplesize), samplesize),
                     y = value, color = model)) +
    ggtitle("Überschätzung mittels Multiplier Bootstrap") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovABC.mb

# Overcoverage je Erwartungswertfunktion Durchschnittswerte
povercovmean <- povercov +
    geom_line(data = summarise.overcov.gkf.means,
              aes(y = value_mean, color = model)) +
    geom_point(data = summarise.overcov.gkf.means,
               aes(y = value_mean, color = model, shape = model))+
    ggtitle("Durchschnittliche Überschätzung") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovmean

povercovmean.gkf <- povercov +
    geom_boxplot(data = overcov.gkf.means,
                 aes(x = reorder(as.character(samplesize), samplesize),
                     y = value, color = model)) +
    ggtitle("Überschätzung mittels t-GKF") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovmean.gkf

dev.off()



