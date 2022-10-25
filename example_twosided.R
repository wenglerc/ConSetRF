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

overcov.per.model <- function(results, S, covset, methodname){
    cov <-
        sapply(results, function(samplelist) {
            sapply(samplelist, function(U) {
                length( U[!(U %in% covset)] ) / (length(S) - length(covset))
            })
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

mu1 <- function(S, scale = 10, length = 1/3, edgefunc = F) {
    left <- (1 - length) / 2
    right <- (1 + length) / 2
    if (edgefunc == F)
        sapply(S, function(s) {
            if (s <= left) (scale * (s - left) ^ 2)
            else if (s < right) 0
            else (-scale * (s - right) ^ 2)
        })
    else
        sapply(S, function(s) {
            if (s <= left) sin( s * 8*pi / (1 - length) + pi/2) - 1/left*s
            else if (s < right) 0
            else sin( s * 8*pi / (1 - length) + pi/2) + 1/right*s - 2
        })
} # Erwartungswertfunktion 1
mu2 <- function(S) {
    sapply(S, function(s){
        if (s <= 1/5) 15*(s - 1/5)^2
        else if (s < 2/5) 0
        else if (s <= 3/5) (sin(s*10*pi + pi/2) - 1)/5
        else if (s < 4/5) 0
        else -15*(s - 4/5)^2
    })
} # Erwartungswertfunktion 2
mu3 <- function(S) {
    sapply(S, function(s){
        if (s <= 1/9) 20*(s - 1/9)^2
        else if (s < 2/9) 0
        else if (s <= 3/9) (sin(s*18*pi + pi/2) - 1)/10
        else if (s < 4/9) 0
        else if (s <= 5/9) -(sin(s*18*pi + pi/2) - 1)/10
        else if (s < 6/9) 0
        else if (s <= 7/9) (sin(s*18*pi + pi/2) - 1)/10
        else if (s < 8/9) 0
        else -20*(s - 8/9)^2
    })
} # Erwartungswertfunktion 3

sigma <- function(S) {
    rep(0.2, length(S))
} # Varianzfunktion

c1 <- function(S, h = 0.1) {
    outer(S, S, FUN = function(x, y) exp( -( x - y )^2 / (2*h^2) ))
} # Fehlerprozess 1

c2 <- function(S, a = 30, l = 0.8) {
    outer(S, S, FUN = function(x, y) (1 + ( x - y )^2 / 2*a*l^2)^(-a) )
} # Fehlerprozess 2


# Variablen zur Steuerung der Auswertung -----------------------------------

randomseed <- 1
samplesize.max <- 400
iterations <- 1500
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
mu4.S <- mu1(S, scale = 1)    # near 0 at the border
mu5.S <- mu1(S, length = 1/50)  # shorter intervall
mu6.S <- mu1(S, length = 3/4)   # longer intervall
mu.all <- list(mu1.S, mu2.S, mu3.S, mu4.S, mu5.S, mu6.S)

rm(mu1.S, mu2.S, mu3.S, mu4.S, mu5.S, mu6.S)
gc()

sigma.S <- sigma(S)
c1.S <- c1(S)
c2.S <- c2(S)

S0 <- lapply(mu.all, function(mu) {S[mu == 0]})

# Realisierungen von X simulieren
cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

tic("Daten generieren")
eps.all <- foreach(i = 1:(3 * iterations)) %dopar% {
    if (i <= iterations) {
        t(mvtnorm::rmvt(samplesize.max, sigma = c1.S, df = 3))
    } else if (i <= (2 * iterations)) {
        t(mvtnorm::rmvt(samplesize.max, sigma = c2.S, df = 3))
    } else {
        t(mvtnorm::rmvnorm(samplesize.max, sigma = c1.S))
    }
}

data.A <- foreach(mu = mu.all) %:% foreach(eps = eps.all[-(1:(2*iterations))]
            ) %dopar% { mu + sigma.S * eps }
data.BC <- foreach(eps = eps.all[1:(2*iterations)]
            ) %dopar% { mu.all[[1]] + sigma.S * eps }

toc()

stopCluster(cluster)

# Benennung und Ordnung der Daten
names(eps.all) <- paste(rep(c("B.", "C.", "A."), each = iterations),
                        1:iterations, sep = "")

names(data.A) <- paste("A", 1:length(data.A), ".", sep = "")

names(data.BC) <- paste(rep(c("B.", "C."), each = iterations),
                        1:iterations, sep = "")

data.all <- c(data.BC, unlist(data.A, recursive = F))

rm(data.A, data.BC, cluster, c1.S, c2.S)
gc()


# Darstellung von Realisierungen der Modelle A,B,C

data.examp.list <- list(
    data.A1.examp = data.all[[which(grepl("A1", names(data.all)) == T)[1]]],
    data.B.examp = data.all[[which(grepl("B", names(data.all)) == T)[1]]],
    data.C.examp = data.all[[which(grepl("C", names(data.all)) == T)[1]]]
)
data.noise.list <- list(
    data.A1.examp = sigma.S * eps.all[[which(grepl("A", names(eps.all)) == T)[1]]],
    data.B.noise = sigma.S * eps.all[[which(grepl("B", names(eps.all)) == T)[1]]],
    data.C.noise = sigma.S * eps.all[[which(grepl("C", names(eps.all)) == T)[1]]]
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
                  alpha = 0.8) +
        geom_line(data = data.frame(S = S, Erwartungswert = mu.all[[1]]),
                  aes(S, Erwartungswert), colour = "red", size = 1)

    plots.noise[[i]] <- p +  ylim(min.noise, max.noise) +
        ggtitle(paste("Fehlerprozess",
                      gsub(".*\\.(.+)\\..*", "\\1", names(data.noise.list)[i]))) +
        geom_line(data = data.noise, aes(S, Wert, colour = Sample, linetype = Sample))
}

# (plots.examp[[1]] + plots.examp[[2]] + plots.examp[[3]]) /
#     (plots.noise[[1]] + plots.noise[[2]] + plots.noise[[3]])

# Darstellung von Realisierungen der Modelle A1,..., A6

data.exampA.list <- list(
    data.A1.examp = data.all[[which(grepl("A1", names(data.all)) == T)[1]]],
    data.A2.examp = data.all[[which(grepl("A2", names(data.all)) == T)[1]]],
    data.A3.examp = data.all[[which(grepl("A3", names(data.all)) == T)[1]]],
    data.A4.examp = data.all[[which(grepl("A4", names(data.all)) == T)[1]]],
    data.A5.examp = data.all[[which(grepl("A5", names(data.all)) == T)[1]]],
    data.A6.examp = data.all[[which(grepl("A6", names(data.all)) == T)[1]]]
)

# max.exampA <- max(sapply(data.exampA.list, function(mat) max(mat[,1:10])))
# min.exampA <- min(sapply(data.exampA.list, function(mat) min(mat[,1:10])))
max.exampA <- 1.5
min.exampA <- -1.5

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
                  alpha = 0.8) +
        geom_line(data = data.frame(S = S, Erwartungswert = mu.all[[i]]),
                  aes(S, Erwartungswert), colour = "red", size = 1)
}

# (plots.exampA[[1]] + plots.exampA[[2]] + plots.exampA[[3]]) /
#     (plots.exampA[[4]] + plots.exampA[[5]] + plots.exampA[[6]])

# Speicherplatz frei machen
rm(
    eps.all, mu.all, sigma.S,
    data.examp, data.noise, data.examp.list, data.noise.list, data.exampA.list,
    max.examp, min.examp, max.noise, min.noise, max.exampA, min.exampA,
    p, plots.examp, plots.noise, plots.exampA
)
gc()

# Berechnungen ----------------------------------------------------------------

# S0.2 = S0[[2]]
#
# U <- ConfSet(data.all[[1]], S, partitions, alpha, level, testcase = 1,
#              pmethod = "tGKF")
# cat("S0 nicht in U: ", length( S0.2[ !(S0.2 %in% U)] ) / length(S0.2),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0.2)] ) / length(U))
#
# U <- ConfSet(data.all[[1]], S, partitions, alpha, level, testcase = 1,
#              pmethod = "mboot", mb.iter = 500)
# cat("S0 nicht in U: ", length( S0.2[ !(S0.2 %in% U)] ) / length(S0.2),
#     "\nU nicht in S0: ", length( U[ !(U %in% S0.2)] ) / length(U))

#
### Laufzeit -----------------------------------------------------------------
# mbm <- microbenchmark(
#     GKF = ConfSet(data.all[[3]][, 1:50], S, partitions, alpha, level,
#                   testcase = 1, pmethod = "tgkf"),
#     MB500 = ConfSet(data.all[[3]][, 1:50], S, partitions, alpha, level,
#                     testcase = 1, pmethod = "mboot", mb.iter = 500),
#     MB1000 = ConfSet(data.all[[3]][, 1:50], S, partitions, alpha, level,
#                      testcase = 1, pmethod = "mboot", mb.iter = 1000),
#     MB2000 = ConfSet(data.all[[3]][, 1:50], S, partitions, alpha, level,
#                      testcase = 1, pmethod = "mboot", mb.iter = 2000),
#     times = 500
# )
# autoplot(mbm)


## Berechnungen per Samplesize -----------------------------------------------

samplesize.list <- c(20, 50, 100, 200, 400)
names(samplesize.list) <- paste("N", samplesize.list, sep = "")

#### Ergebnisse GKF ---------------------------------------------------------

cluster <- makeCluster(detectCores() - 1)
registerDoParallel(cluster)

# Obermengen je Samplegroesse
tic("GKF parallel")
results.gkf <-
    foreach(data = data.all, .packages = c("dplyr", "ConfSetRF")) %dopar% {

        lapply(samplesize.list, function(N) {
            ConfSet(data[, 1:N], S, partitions, alpha, level, testcase = 1,
                    pmethod = "tgkf")
        })
    }
stopCluster(cluster)
toc()

names(results.gkf) <- names(data.all)

# Coverage je Modell ABC
cov.gkf.ABC <-
    cov.per.model(results.gkf[1:(3*iterations)], S0[[1]], "t-GKF")
summarise.cov.gkf.ABC <- cov.gkf.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")
# Overcoverage je Modell ABC
overcov.gkf.ABC <-
    overcov.per.model(results.gkf[1:(3*iterations)], S, S0[[1]], "t-GKF")
summarise.overcov.gkf.ABC <- overcov.gkf.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")


# Coverage je Erwartungswertfunktion
cov.gkf.means <-
    lapply(1:length(S0), function(i){
        cov.per.model(results.gkf[((i+1)*iterations + 1):((i+2)*iterations)],
                      S0[[i]], "t-GKF")
        }) %>% do.call(rbind, .)
summarise.cov.gkf.means <- cov.gkf.means %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")
# Overcoverage je Erwartungswertfunktion
overcov.gkf.means <-
    lapply(1:length(S0), function(i){
        overcov.per.model(results.gkf[((i+1)*iterations + 1):((i+2)*iterations)],
                      S, S0[[i]], "t-GKF")
    }) %>% do.call(rbind, .)
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
    foreach(data = data.all[1:(3*iterations)],
            .packages = c("dplyr", "ConfSetRF")) %dopar% {

        lapply(samplesize.list, function(N) {
            ConfSet(data[, 1:N], S, partitions, alpha, level,
                    testcase = 1, pmethod = "mboot", mb.iter = 500)
        })

    }
stopCluster(cluster)
toc()

names(results.mb) <- names(data.all[1:(3*iterations)])

# Coverage je Modell ABC
cov.mb.ABC <-
    cov.per.model(results.mb, S0[[1]], "Multiplier Bootstrap")
summarise.cov.mb.ABC <- cov.mb.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")
# Overcoverage je Modell ABC
overcov.mb.ABC <-
    overcov.per.model(results.mb, S, S0[[1]], "Multiplier Bootstrap")
summarise.overcov.mb.ABC <- overcov.mb.ABC %>%
    group_by(method, model, samplesize) %>%
    summarise(across(.fns = list(mean = mean, sd = sd)), .groups = "drop")

rm(results.mb) # Speicherplatz frei machen
gc()


# Auswertung ----------------------------------------------------------------
pdf(file = "Auswertungen_twosided.pdf", width = 6, height = 4)

## Coverage -----------------------------------------------------------------
cov.plot.ABC <- rbind(summarise.cov.gkf.ABC, summarise.cov.mb.ABC) %>%
    arrange() %>%
    mutate(method = factor(method,level = c("t-GKF", "Multiplier Bootstrap")))

pcov <- ggplot(cov.plot.ABC, aes(x = samplesize)) +
    labs(x = "Anzahl an Samples", y = "Überdeckungsrate", color = "Modell",
         linetype = "Methode", shape = "Modell") +
    scale_linetype(labels = c("t-GKF", "MB")) +
    theme(legend.position = "right") +
    geom_hline(yintercept = 1 - alpha, col = "red3", size = 0.75,
               linetype = "dashed") +
    geom_hline(yintercept = 1 - alpha + 1.96*sqrt((1-alpha)*alpha/iterations),
               col = "black", size = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 1 - alpha - 1.96*sqrt((1-alpha)*alpha/iterations),
               col = "black", size = 0.5, linetype = "dashed")

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
               aes(y = value_mean, color = model, shape = model)) +
    ggtitle("t-GKF") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
pcovmean

## Overcoverage -----------------------------------------------
overcov.plot.ABC <-
    rbind(summarise.overcov.gkf.ABC, summarise.overcov.mb.ABC) %>%
    arrange() %>%
    mutate(method = factor(method,level = c("t-GKF", "Multiplier Bootstrap")))

povercov <- ggplot(overcov.plot.ABC, aes(x = samplesize)) +
    labs(x = "Anzahl an Samples",
         y = "Falsch-Negativ-Rate",
         color = "Modell", linetype = "Methode", shape = "Modell") +
    scale_linetype(labels = c("t-GKF", "MB")) +
    theme(legend.position = "right") +
    geom_hline(yintercept = alpha, col = "red3", size = 0.75,
               linetype = "dashed")
# Overcoverage je Modell ABC Durchschnittswerte
povercovABC <- povercov +
    geom_line(aes(y = value_mean, color = model, linetype = method)) +
    geom_point(aes(y = value_mean, color = model, shape = model))
povercovABC

# Overcoverage je Modell ABC GKF
povercovABC.gkf <- povercov +
    geom_boxplot(data = overcov.gkf.ABC,
                 aes(x = reorder(as.character(samplesize), samplesize),
                     y = value, color = model)) +
    ggtitle("t-GKF") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovABC.gkf

# Overcoverage je Modell ABC MB
povercovABC.mb <- povercov +
    geom_boxplot(data = overcov.mb.ABC,
                 aes(x = reorder(as.character(samplesize), samplesize),
                     y = value, color = model)) +
    ggtitle("Multiplier Bootstrap") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovABC.mb

# Overcoverage je Erwartungswertfunktion Durchschnittswerte
povercovmean <- povercov +
    geom_line(data = summarise.overcov.gkf.means,
              aes(y = value_mean, color = model)) +
    geom_point(data = summarise.overcov.gkf.means,
               aes(y = value_mean, color = model, shape = model))+
    ggtitle("t-GKF") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovmean

povercovmean.gkf <- povercov +
    geom_boxplot(data = overcov.gkf.means,
                 aes(x = reorder(as.character(samplesize), samplesize),
                     y = value, color = model)) +
    ggtitle("t-GKF") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
povercovmean.gkf

dev.off()

pdf(file = "resultsABC_twosided.pdf", width = 9, height = 4)
pcovABC+ theme(legend.position = "none") + povercovABC
dev.off()

pdf(file = "resultsA1-6_twosided.pdf", width = 9, height = 4)
pcovmean+ theme(legend.position = "none") + povercovmean
dev.off()

# save.image(file = "results_twosided_FNR.RData")

