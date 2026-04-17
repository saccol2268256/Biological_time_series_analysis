setwd("C:/Users/sacco/Desktop/lavoro/25/nascite")
library(readxl)
library(tidyverse)
library(forecast)
library(MTS)
library(dplyr)
library(tidyr)
library(corrplot)
library(ggplot2)
library(glmmTMB)
library(emmeans)

library(performance)
library(MuMIn)

df <- read_excel("nnFinal.xlsx", col_types = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
summary(df)
head(df)
df$born_volume[df$born_volume==0] <- 0.0000000000000001
df$dove=as.factor(df$dove)
df$box=as.factor(df$box)
boxplot(df$born_volume~df$dove)
boxplot(df$born_volume~df$box)
hist(log(df$born_volume*100),breaks = 15)
par(mfrow = c(1, 2)) 
livelli_dove <- unique(df$dove)
for (val in livelli_dove) {
  #dati formato largo
  df_subset <- df %>%
  filter(dove == val) %>%
  select(box, day, born_volume)
  df_wide <- df_subset %>%
  pivot_wider(names_from = day, 
              values_from = born_volume, 
              names_prefix = "Day_") %>%
  select(-box) 
  M <- cor(df_wide, use = "pairwise.complete.obs")
  labels <- colnames(M)
  new_labels <- sapply(seq_along(labels), function(i) {
    if ((i - 1) %% 5 == 0) return(labels[i]) else return("")
  })
  rownames(M) <- new_labels
  colnames(M) <- new_labels
  corrplot(M, 
           method = "ellipse", 
           type = "upper", 
           order = "original", 
           tl.col = "black", 
           tl.cex = 0.8,      
           tl.srt = 45,       
           diag = FALSE, 
           title = paste("Struttura Temporale | Dove =", val),
           mar = c(0,0,2,0))  # Margini per il titolo
}
##emerge dal grafico correlazione di tipo autoregressivo in entrambi i casi, 
#sottolineamo la correlazione positiva, chiara e facilmente comprensibile
#nel casi dove=1 di località breganze vediamo emergere pattern più complesso
#indicante nascite asincrone nell 'arco temporale tale da non identifivare un vero e proprio picco
#come accade nella popolazione di susegana a prescindere dalla densitaà di nascite
#come emergera dopo vediamo un piccolo rimbalzo nelle crescite intorno al giorno 60
par(mfrow = c(1, 1))
# Correlation Coefficients per Group
cor_by_group <- df %>%
  group_by(dove) %>%
  summarise(Correlation = cor(born_volume, day, use = "complete.obs"))
print(cor_by_group) 



#ANALISI TRAMITE SERIE STORICHE 
# Qui avviene il cambio logico: ignoriamo i 'box' e sommiamo per 'dove'
df_grouped <- df %>%
  group_by(day, dove) %>%
  summarise(total_volume = sum(born_volume, na.rm = TRUE), .groups = "drop") %>%
  arrange(day)
df_wide <- df_grouped %>%
  pivot_wider(names_from = dove, values_from = total_volume)
head(df_wide)
nomi_luoghi <- setdiff(names(df_wide), "day")
# Creazione delle 2 Serie Temporali (TS)
ts_luogo_1 <- ts(df_wide[[nomi_luoghi[1]]])
ts_luogo_2 <- ts(df_wide[[nomi_luoghi[2]]])
day <- df_wide$day
#ANALISI ARIMA STOCASTIC PROCESS
ts_diff_1 <- diff(ts_luogo_1, differences = 1)
fit_arima <- auto.arima(ts_diff_1, seasonal = F, stepwise = TRUE)
summary(fit_arima)
par(mfrow=c(2,1))
acf(ts_diff_1, main=paste("ACF -", nomi_luoghi[1]))
pacf(ts_diff_1, main=paste("PACF -", nomi_luoghi[1]))
ts_diff_2 <- diff(ts_luogo_2, differences = 1)
fit_arima2 <- auto.arima(ts_diff_2, seasonal = F, stepwise = TRUE)
summary(fit_arima2)
acf(ts_diff_2, main=paste("ACF -", nomi_luoghi[2]))
pacf(ts_diff_2, main=paste("PACF -", nomi_luoghi[2]))
par(mfrow=c(1,1))
#gli errori sono un processo stocastico bianco se visto sotto effetto di medie mobili in entrambe le zone
#notiamo come il processo stocastico alla base delle nascite sia lo stesso in entrambe le zone

# APPROCCIO CLASSICO (Decomposizione su 2 serie)
n <- length(ts_luogo_1)
tempo <- seq(1, n, 1)
tempo2 <- tempo^2 
tempo3 <- tempo^3 
filtro <- c(rep(1/10, 2),rep(2/10, 3),rep(2/10, 2))
analizza_luogo <- function(ts_data) {
  T_hat <- stats::filter(ts_data, filter=filtro, sides=2)
  # Isolamento Stagionalità
  CS <- ts_data - T_hat
  CS_vec <- as.numeric(CS)
  periodo <- 7
  pad <- (ceiling(length(CS_vec)/periodo)*periodo) - length(CS_vec)
  CS_mat <- matrix(c(CS_vec, rep(NA, pad)), ncol=periodo, byrow=TRUE)
  CG <- apply(CS_mat, 2, mean, na.rm=TRUE)
  # Centratura degli indici
  CI <- CG - mean(CG)
  # Destagionalizzazione
  deseasonalized <- ts_data - rep(CI, length.out=n)
  # Fit Finale Polinomiale (Cubico)
  fit_final <- lm(deseasonalized ~ -1 + tempo + tempo2+tempo3)
  print(summary(fit_final))
  return(list(trend_mm = T_hat, forecast = fitted(fit_final), seasonality_indices = CI))
}
# Eseguo l'analisi sulle 2 serie
res_L1 <- analizza_luogo(ts_luogo_1)
res_L2 <- analizza_luogo(ts_luogo_2)

build_df <- function(obs, res, nome) {
  data.frame(
    date = day,
    Category = nome,
    Observed = as.numeric(obs),
    Forecast = as.numeric(res$forecast),
    Trend = as.numeric(res$trend_mm)
  )
}
df_res_final <- bind_rows(
  build_df(ts_luogo_1, res_L1, nomi_luoghi[1]),
  build_df(ts_luogo_2, res_L2, nomi_luoghi[2])
) %>%
  pivot_longer(cols = c("Observed", "Forecast", "Trend"), 
               names_to = "Type", values_to = "Value")
ggplot(df_res_final, aes(x = date, y = Value, color = Type, linetype = Type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Category, scales = "free_y") +
  labs(title = "Analisi Serie Temporali Aggregate per Luogo",
       x = "Data", y = "Somma Born Volume",
       color = "Tipo Serie", linetype = "Tipo Serie") +
  theme_minimal() +
  scale_color_manual(values = c("Observed"="black", "Forecast"="red", "Trend"="blue"))
#"L'analisi grafica delle serie storiche evidenzia una morfologia a campana ("single-humped") per entrambe le locazioni. 
#La componente di Trend (blu) rivela che il modello previsionale classico (rosso) tende a sottostimare l'intensità del picco (bias negativo al vertice) e a sovrastimare la durata della coda (bias positivo in chiusura). 
#Ciò suggerisce che le nascite, pur essendo spalmate, si concentrano con una virulenza maggiore del previsto nella fase centrale (giorni 15-25).



#UTILIZZO DI MODELLI CON EFFETTI CASUALI
##identificazioni di ic e confronto tra i coeficienti 
day2=df$day^2
day3=df$day^3
df$gruppo_unico <- interaction(factor(df$box), factor(df$dove))
fit_ar1 <- glmmTMB(born ~ day*dove + day2 + day3 +volume +offset(log(volume)) + ar1(factor(day) | box),
                   family = poisson,
                   data = df)
summary(fit_ar1)#try without volume

fit_arrr <- glmmTMB(born ~ day*dove + day2 + day3 +volume +offset(log(volume)) + ar1(factor(day) | gruppo_unico),
                   family = poisson,
                   data = df)
summary(fit_arrr)
##1. Struttura di Correlazione ed Effetti CasualiIl modello conferma statisticamente le ipotesi sulla dinamica temporale ("nascite spalmate") discusse in precedenza attraverso due componenti fondamentali:Correlazione Autoregressiva ar1(factor(day) | box):Il parametro di correlazione stimato è 0.82.Significato: Si tratta di una correlazione positiva molto forte. Questo valore quantifica matematicamente l'effetto "trascinamento" o "momentum" descritto prima.Interpretazione Biologica: Esiste una fortissima "memoria" a breve termine: il numero di nascite/insetti al giorno $t$ è fortemente predittivo del giorno $t+1$. 
#   Questo smentisce l'ipotesi di eventi casuali o impulsivi (picchi isolati) e conferma che l'infestazione ha un'inerzia biologica notevole (una volta partita, si auto-sostiene per giorni).
#Variabilità tra Gruppi (box):La varianza dell'intercetta casuale è 0.3514 (Std.Dev 0.59).Significato: C'è eterogeneità tra le diverse "scatole" (box) di osservazione. Anche a parità di luogo e giorno, alcune scatole tendono ad avere sistematicamente più insetti di altre. Questo suggerisce l'esistenza di micro-variabili locali non misurate (es. microclima specifico della scatola) che influenzano il livello base dell'infestazione. 
#Analisi del Coefficiente dove1Il parametro fisso dove1 (effetto del luogo) è l'elemento cruciale per distinguere l'intensità del fenomeno.Stima: Il coefficiente è -2.267.Significatività: Il p-value è < 2e-16 (estremamente significativo).Interpretazione Quantitativa:Essendo un modello log-lineare, l'effetto sul conteggio atteso si calcola esponendo il coefficiente: $e^{-2.267} \approx 0.103$.Questo ci dice che, a parità di altre condizioni, nel luogo Dove=1 ci aspettiamo circa il 10% del volume di insetti rispetto al luogo di riferimento (Dove=0). 
#Altri Elementi RilevantiPolinomiali Temporali (day, day2, day3):Tutti i termini sono altamente significativi ($p < 2e-16$). La significatività del termine quadratico (day2, negativo) e cubico (day3) conferma che la relazione temporale non è lineare ma curvilinea (la famosa forma a "campana" o curva asimmetrica vista nei grafici). Il modello sta catturando efficacemente la fase di crescita, il picco e la decrescita.Interazione day:dove1:Il termine è significativo e positivo (0.023).Ciò indica che la curva temporale nel luogo 1 non è solo "più bassa", ma ha anche una forma leggermente diversa (pendenza differente) rispetto al luogo 0. La dinamica di evoluzione dell'infestazione varia leggermente a seconda dell'ambiente.
# Sempre dal pacchetto performance


# Calcola R2
#r2_vals <- r.squaredGLMM(fit_ar1)
#print(r2_vals)

icc(fit_ar1)
#Un ICC > 0.5 indica che la struttura "Box/Dove" è predominante rispetto alle variabili che hai misurato (volume, day).
# Calcola R2
r2_nakagawa <- r2_nakagawa(fit_ar1)
r2_nakagawa
#Se $R^2_m = 0.30$ e $R^2_c = 0.35$: Gli effetti casuali incidono poco (5%).
#Se $R^2_m = 0.30$ e $R^2_c = 0.80$: Gli effetti casuali sono fondamentali (spiegano il 50% della variabilità extra).

# 1. Varianza degli Effetti Fissi (Fixed Effects)
# Calcoliamo la varianza dei predetti (sul link scale)
var_f <- var(predict(fit_ar1, type = "link", re.form = NA))

# 2. Varianza degli Effetti Random (Random Effects)
# Estraiamo la somma delle varianze delle componenti random
vc <- VarCorr(fit_ar1)

get_variance <- function(matrix_structure) {
  # Se è una matrice (come nel caso AR1), prendiamo la media della diagonale.
  # Essendo stazionaria (valori uguali 0.3514), la media è uguale al singolo valore.
  # Usiamo diag() per estrarre la diagonale e mean() per averne uno solo.
  if (is.matrix(matrix_structure)) {
    return(mean(diag(matrix_structure))) 
  } else {
    # Caso fallback (se fosse uno scalare o struttura diversa)
    return(attr(matrix_structure, "stddev")^2)
  }
}
var_r <- sum(sapply(vc$cond, function(x) {
  attr(x, "stddev")^2
}))
# Applichiamo la funzione a ogni componente random (cond) e sommiamo TRA I GRUPPI
# Sommiamo Box + Dove + Interazione. NON sommiamo i giorni.
var_r <- sum(sapply(vc$cond, get_variance))
# Nota: "attr(x, 'stddev')^2" prende la varianza associata alla struttura ar1/box

# 3. Varianza della Distribuzione (Residual Variance)
# Per Poisson con link log, usiamo l'approssimazione log-normale (Nakagawa et al. 2017)
# Si calcola come log(1 + 1/mu). Usiamo la media dei valori predetti come mu.
mu_mean <- predict(fit_ar1, type = "response")
var_d <- mean(log(1 + 1/mu_mean))

# --- Calcolo R2 ---
var_f
var_r
# R2 Marginale (Spiegato solo dai fissi: day, dove, volume...)
r2_marg <- var_f / (var_f + var_r + var_d)

# R2 Condizionale (Spiegato da fissi + struttura random ar1)
r2_cond <- (var_f + var_r) / (var_f + var_r + var_d)

# Output
cat("R2 Marginale:", round(r2_marg, 4), "\n")
cat("R2 Condizionale:", round(r2_cond, 4), "\n")


library(DHARMa)
res <- simulateResiduals(fit_ar1)
plot(res)




df$predicted_values <- predict(fit_ar1, type = "response")
# grafico adattamento
p2 <- ggplot(df, aes(x = day, group = box)) +
  geom_point(aes(y = born, color = dove), alpha = 0.3, size = 1) +
  geom_line(aes(y = predicted_values, color = dove), alpha = 0.6) +
  facet_wrap(~ dove) + # Separa i grafici per Luogo
  labs(
    title = "Fit del Modello nel Tempo",
    subtitle = "Linee: Predizioni del modello | Punti: Dati reali",
    x = "Giorno",
    y = "Nati (born)"
  ) +
  theme_bw()
print(p2)

#inizia a creare un report da questa analisi 1.i dati sono una serie temporale in diverse aree dove= 0 localita susegana, dove=1 di località breganze, la variabile born segue una distribuzione di Poisson per ragioni naturali della variabile e la distribuzione discreta di questa 2.valutazione della matrice di correlazione sottolinea emerge dal grafico correlazione di tipo autoregressivo in entrambi i casi, sottolineamo la correlazione positiva ttra i vari soggetti appartenenti alla stessa categoria, chiara e facilmente comprensibile nel casi dove=1 di località breganze vediamo emergere pattern più complesso indicante nascite asincrone nell 'arco temporale tale da non identifivare un vero e proprio picco uniforme in tutte le osservazioni raccolte come accade nella popolazione di susegana a prescindere dalla densità di nascite come emergera dopo vediamo un piccolo rimbalzo nelle crescite intorno al giorno 60 3.possiamo vedere il fenomeno sotto un approccio classico di serie temporali che sottostimano il picco epidemiologico massimo con un regressione polinomiale cubica e commenta la emminente stima forviante per la non correttezza del modello parametrico alla base dell'inferenza. le nascite, pur essendo spalmate, si concentrano con una virulenza maggiore del previsto nella fase centrale (giorni 15-25) 4.approfondisci l'interessante fatto che anche se su fenomeni di infestazione su scale differenti nei 2 appezzamenti il processo stocastico basato sulla media delle nascite puo essere considerato congruo essendo che la serie differenziata per raggiungere stazionarietà e quindi le deviazioni rispetto al white noise sono riassumibili tramite medie mobili spiegane le implicazione ed interessante vedere la non componente autoregressiva temporale degli errori 5. approfondisci a livello teorico la struttura del modello con effetti casuali utilizzato qui fit_ar1 <- glmmTMB(born ~ day*dove + day2 + day3 +volume +offset(log(volume)) + ar1(factor(day) | factor(df$box) * factor(df$dove)),
# Struttura di Correlazione ed Effetti CasualiIl modello conferma statisticamente le ipotesi sulla dinamica temporale ("nascite spalmate") discusse in precedenza attraverso due componenti fondamentali:Correlazione Autoregressiva ar1(factor(day) | box):Il parametro di correlazione stimato è 0.82.Significato: Si tratta di una correlazione positiva molto forte. Questo valore quantifica matematicamente l'effetto "trascinamento" o "momentum" descritto prima.Interpretazione Biologica: Esiste una fortissima "memoria" a breve termine: il numero di nascite/insetti al giorno $t$ è fortemente predittivo del giorno $t+1$. # Questo smentisce l'ipotesi di eventi casuali o impulsivi (picchi isolati) e conferma che l'infestazione ha un'inerzia biologica notevole (una volta partita, si auto-sostiene per giorni).#Variabilità tra Gruppi (box):La varianza dell'intercetta casuale è 0.3514 (Std.Dev 0.59).Significato: C'è eterogeneità tra le diverse "scatole" (box) di osservazione. Anche a parità di luogo e giorno, alcune scatole tendono ad avere sistematicamente più insetti di altre. Questo suggerisce l'esistenza di micro-variabili locali non misurate (es. microclima specifico della scatola) che influenzano il livello base dell'infestazione. #Analisi del Coefficiente dove1Il parametro fisso dove1 (effetto del luogo) è l'elemento cruciale per distinguere l'intensità del fenomeno.Stima: Il coefficiente è -2.267.Significatività: Il p-value è < 2e-16 (estremamente significativo).Interpretazione Quantitativa:Essendo un modello log-lineare, l'effetto sul conteggio atteso si calcola esponendo il coefficiente: $e^{-2.267} \approx 0.103$.Questo ci dice che, a parità di altre condizioni, nel luogo Dove=1 ci aspettiamo circa il 10% del volume di insetti rispetto al luogo di riferimento (Dove=0). #Altri Elementi RilevantiPolinomiali Temporali (day, day2, day3):Tutti i termini sono altamente significativi ($p < 2e-16$). La significatività del termine quadratico (day2, negativo) e cubico (day3) conferma che la relazione temporale non è lineare ma curvilinea (la famosa forma a "campana" o curva asimmetrica vista nei grafici). Il modello sta catturando efficacemente la fase di crescita, il picco e la decrescita.Interazione day:dove1:Il termine è significativo e positivo (0.023).Ciò indica che la curva temporale nel luogo 1 non è solo "più bassa", ma ha anche una forma leggermente diversa (pendenza differente) rispetto al luogo 0. La dinamica di evoluzione dell'infestazione varia leggermente a seconda dell'ambiente (correggi i coeficienti di questa bozza)
