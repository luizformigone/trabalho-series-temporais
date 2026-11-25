library(quantmod) # puxar dados do Yahoo
library(tseries) #para séries temporais
library(forecast) #para auto.arima
library(urca) # para raíz unitaria e teste de Dickey-Fuller Aumentado (ADF)

# Pasta para guardar os gráficos (cria só se ainda não existir)
dir.create("graficos", showWarnings = FALSE)

# Salva o gráfico que está na aba Plots como PNG dentro da pasta "graficos"
salvar_grafico <- function(nome, largura = 1600, altura = 900) {
  dev.copy(png, filename = file.path("graficos", paste0(nome, ".png")),
           width = largura, height = altura, res = 150)
  dev.off()
}

# Dados fixos: baixa do Yahoo só na primeira vez e salva; nas próximas, lê do arquivo salvo
# (o Yahoo recalcula o Adjusted a cada novo provento, então baixar de novo poderia mudar os resultados)
if (file.exists("BBAS3_dados_brutos.rds")) {
  BBAS3.SA <- readRDS("BBAS3_dados_brutos.rds")
} else {
  getSymbols("BBAS3.SA", src = "yahoo", from = "1900-01-01", to = "2026-11-24", auto.assign = TRUE) # 03/01/2000 até 23/11/2026 (to é exclusivo)
  saveRDS(BBAS3.SA, "BBAS3_dados_brutos.rds")
}
range(index(BBAS3.SA)) # conferência: deve ser 2000-01-03 até 2026-11-23
anyNA(BBAS3.SA) # TRUE. Tem valores faltantes
sum(is.na(BBAS3.SA)) # Quantos? 6
index(BBAS3.SA)[is.na(BBAS3.SA)] # Todos em uma só data. 02/11/2004. Vamos retirar esta linha
BBAS3.SA <- na.omit(BBAS3.SA)
sum(Vo(BBAS3.SA)==0) # Quantos pregões tiveram Volume 0. Yahoo finance só copia do dia anterior. Pode dar ruim com log-retorno. Tiveram 122
index(BBAS3.SA)[Vo(BBAS3.SA)==0] # Datas que o volume foram 0. A maioria feriados. Vamos retirar também
BBAS3.SA <- BBAS3.SA[Vo(BBAS3.SA)>0] #Retirado
preço_ajustado <- BBAS3.SA$BBAS3.SA.Adjusted # Selecionando apenas a coluna de preço ajustado.
datas <- index(preço_ajustado) # Extraindo as datas dos pregões para colocarmos como eixo horizontal nos gráficos
ts_preço_ajustado <- ts(as.numeric(preço_ajustado)) # Transformando numa série temporal. as.numeric() Claude recomendou para transformar num vetor
n <- length(ts_preço_ajustado) # 6627 pregões (observações)
ts_preço_ajustado_treino <- ts(as.numeric(ts_preço_ajustado[1:(n - 30)])) # Tirando 30 observações que o enunciado pediu para formar a série de treino
datas_treino <- datas[1:(n - 30)]
ts_preço_ajustado_teste <- ts(as.numeric(ts_preço_ajustado[(n - 29):n])) # Considerando apenas as ultimas 30 observações. n - 29 pq é inclusivo
datas_teste <- datas[(n - 29):n]

## Etapa 1: Identificação (Metodologia Box-Jenkins)
plot(datas_treino,ts_preço_ajustado_treino,type = "l",main = "BBAS3", xlab = "Data", ylab = "Valor preço ajustado") #plotar ts treino para visualizarmos
salvar_grafico("01_preco_treino")
modelo_auto_arima <- auto.arima(ts_preço_ajustado_treino)
modelo_auto_arima # ARIMA(0,1,0). Random-walk sem drift. Esperado pela hipótese de eficiência dos mercados

# FAC e FACP da série de treino e de sua 1a diferença
par(mfrow = c(2,3)) # R gerar figura com plots em 2 linhas e 3 colunas.
plot(datas_treino, ts_preço_ajustado_treino, type = "l", main = "Nível", xlab = "", ylab = "R$")
acf(ts_preço_ajustado_treino, main = "FAC — nível"); pacf(ts_preço_ajustado_treino, main = "FACP — nível") # Plotar FAC e FACP da série de treino. Pela FAC nitidamente é não estacionária. Não decai exponencialmente
ts_preço_ajustado_treino_diff1 <- diff(ts_preço_ajustado_treino) # Função diff() para fazermos a primeira diferença da série temporal
datas_treino_diff1 <- datas_treino[-1] # Tiramos a primeira data para bater com a série temporal
plot(datas_treino_diff1, ts_preço_ajustado_treino_diff1, type = "l", main = "1ª diferença", xlab = "", ylab = "") # média constante porém variância vai aumentando
acf(ts_preço_ajustado_treino_diff1, main = "FAC — diferença"); pacf(ts_preço_ajustado_treino_diff1, main = "FACP — diferença") # FAC = MA(0) e FACP = AR(0): Random-Walk
salvar_grafico("02_fac_facp_nivel_diferenca", largura = 2000, altura = 1200) # salvar depois dos 6 gráficos e antes de voltar ao normal
par(mfrow = c(1, 1)) # Voltar para a visualização normal

# Testes de raíz unitária para formalizar d = 1
# Se |t| > |valor crítico|, então rejeita-se 𝐻0, ou seja, a série é estacionária. Caso contrário, há evidência de raiz unitária e, assim, deve-se tomar a diferença. (Slide 14 aula Séries Temporais Não Estacionárias)
# (O urca mostra os críticos negativos: equivale a rejeitar H0 se t < valor crítico)
# type = "trend": inclui constante e tendência na regressão do teste, pois o preço sobe ao longo do tempo
# lags = 20: no máximo 20 defasagens da diferença para eliminar autocorrelação dos resíduos ("aumentado")
# selectlags = "BIC": escolhe quantas dessas defasagens usar pelo critério BIC (escolheu 1)
adf_ts_preço_ajustado_treino <- ur.df(ts_preço_ajustado_treino, type = "trend", lags = 20, selectlags = "BIC")
summary(adf_ts_preço_ajustado_treino) # tau3 = -2,58 > -3,41 (5%) -> não rejeita H0 -> raiz unitária; phi3 = 3,74 < 6,25 (5%)  -> sem tendência determinística; phi2 = 3,15 < 4,68 (5%)  -> sem drift: compatível com random walk puro (confirma auto.arima)
# 1ª diferença: só com constante
adf_ts_preço_ajustado_treino_diff1 <- ur.df(ts_preço_ajustado_treino_diff1, type = "drift", lags = 20, selectlags = "BIC")
summary(adf_ts_preço_ajustado_treino_diff1) # tau2 = -56.49 < -2.86 (5%) -> rejeita H0 -> 1ª diferença é estacionária
# KPSS complementar ao ADF. H0: yt ~ I(0); H1: yt ~ I(1)
summary(ur.kpss(ts_preço_ajustado_treino, type = "tau", lags = "long")) # t = 2.6 > 0.15 (5%). Rejeita H0. Precisa integrar em 1. Claude fez função
summary(ur.kpss(ts_preço_ajustado_treino_diff1, type = "mu", lags = "long")) # t = 0.09 < 0.46 (5%). Não rejeita H0. Já está integrado

## Etapa 2: Estimação
# Testar 32 modelos ARIMA(p,1,q): p e q de 0 a 3, com e sem drift. p e q de 0 a 3 pq, como sabemos que é random-walk com AR(0) e MA(0), irmos até 3 é o suficiente.
resultados_modelos <- data.frame() # tabela vazia para guardar os resultados
for (p in 0:3) for (q in 0:3) for (drift in c(FALSE, TRUE)) {
  modelo_candidato <- try(Arima(ts_preço_ajustado_treino, order = c(p,1,q), include.drift = drift), silent = TRUE) # Arima com A maiúsculo (forecast): o arima minúsculo não aceita include.drift. try() para não parar se der erro. silent = TRUE é para não retornar as mensagens de erro
  if (!inherits(modelo_candidato, "try-error")) # se o modelo não(!) der erro ("try-error")
    resultados_modelos <- rbind(resultados_modelos, data.frame(p, q, drift = drift, AIC = AIC(modelo_candidato), BIC = BIC(modelo_candidato))) # rbind add ao data frame uma linha contendo: p, q, se drift, AIC e BIC
}
resultados_modelos[order(resultados_modelos$BIC), ][1:10, ] # 10 melhores pelo BIC (menor = melhor). Random-walk sem drift ganhou
# Testar um SARIMA com sazonalidade semanal
modelo_sazonal_semanal <- Arima(ts(ts_preço_ajustado_treino, frequency = 5), order = c(0,1,0), seasonal = list(order = c(1,0,1), period = 5)) # 5 pq ñ considera fds. seasonal = list(order = c(1,0,1), period = 5): acrescenta um AR e um MA sazonais
modelo_sazonal_semanal # Coeficientes sazonais não são significativos (t ≈ 0,02). Coeficientes (sar1 e sma1) + ou - desvio-padrão engloba 0. AIC e BIC pioram -> sem sazonalidade semanal
# Finalistas para comparação detalhada
modelo_random_walk       <- Arima(ts_preço_ajustado_treino, order = c(0,1,0))                        # melhor pelo BIC
modelo_random_walk_drift <- Arima(ts_preço_ajustado_treino, order = c(0,1,0), include.drift = TRUE)  # testar o drift
modelo_arima_110         <- Arima(ts_preço_ajustado_treino, order = c(1,1,0))                        # melhor ARMA simples pelo BIC
# Estatística t: |t| > 2 -> significativo
coef(modelo_random_walk_drift) / sqrt(diag(modelo_random_walk_drift$var.coef)) # 1.39. Não significativo
coef(modelo_arima_110) / sqrt(diag(modelo_arima_110$var.coef)) # -0.211. Não significativo

## Etapa 3: Diagnóstico
# Ljung-Box no random walk em nível (H0: resíduos sem autocorrelação; quer p > 0,05)
checkresiduals(modelo_random_walk, lag = 10) # p = 0,013. Não passa
checkresiduals(modelo_random_walk, lag = 20) # p = 0,023. Não passa
salvar_grafico("03_residuos_random_walk_nivel")
# Os resíduos em R$ têm variância crescente (o preço foi de ~R$ 0,40 para ~R$ 25): heterocedasticidade.
# As autocorrelações restantes são muito pequenas (|ρ| < 0,04) e a rejeição é amplificada pelo tamanho da amostra (6.597 obs.)
# Solução: modelar o log do preço (lambda = 0). Os resíduos viram variações percentuais (log-retornos)
modelo_random_walk_log <- Arima(ts_preço_ajustado_treino, order = c(0,1,0), lambda = 0)
modelo_random_walk_log # sigma^2 = 0,000682 -> desvio-padrão ≈ 2,6% ao dia
checkresiduals(modelo_random_walk_log, lag = 10) # p = 0,083. Passa
checkresiduals(modelo_random_walk_log, lag = 20) # p = 0,090. Passa
salvar_grafico("04_residuos_random_walk_log")
# Modelo final da parte (i): random walk no log do preço -> log(P_t) = log(P_t-1) + e_t
# Normalidade dos resíduos (H0: normal)
jarque.bera.test(residuals(modelo_random_walk_log)) # rejeita: caudas pesadas (leptocúrtica). IC da previsão deve ser lido com cautela
# Resíduos ao quadrado: autocorrelação indica volatilidade variável no tempo -> motiva o GARCH da parte (ii)
acf(residuals(modelo_random_walk_log)^2, main = "FAC dos resíduos² — random walk (log)")
salvar_grafico("05_fac_residuos_quadrado")

## Etapa 4: Previsão (backtest de 30 pregões)
# Função que calcula os erros da previsão de 30 passos contra a série de teste (em R$)
avaliar_previsao <- function(modelo) {
  previsao <- forecast(modelo, h = 30)
  erro <- as.numeric(ts_preço_ajustado_teste) - as.numeric(previsao$mean)
  c(RMSE = sqrt(mean(erro^2)),                                          # raiz do erro quadrático médio (R$)
    MAE  = mean(abs(erro)),                                             # erro absoluto médio (R$)
    MAPE = 100 * mean(abs(erro / as.numeric(ts_preço_ajustado_teste)))) # erro percentual médio (%)
}
rbind(random_walk       = avaliar_previsao(modelo_random_walk),
      random_walk_drift = avaliar_previsao(modelo_random_walk_drift),
      random_walk_log   = avaliar_previsao(modelo_random_walk_log)) # menor = melhor
# RW e RW log: RMSE 6,55 / MAE 5,25 / MAPE 33,0%. RW com drift: RMSE 6,61 / MAE 5,30 / MAPE 33,3% (drift inclina para cima e o preço caiu)
# Obs.: a previsão do random walk é o último preço observado (previsão ingênua). A diferença entre nível e log
# aparece no intervalo de confiança, que no log é assimétrico e proporcional ao preço
# Erro grande por causa de 2 choques de notícia no período de teste:
# 26/10/2026: -29,3% no 1º pregão após o 2º turno da eleição presidencial (circuit breaker na B3)
# 12/11/2026: -8,7% após o balanço do 3T26 (lucro -66,5% e corte do payout de 30% para 25%)

# Gráfico: últimos 120 pregões do treino + 30 do teste + previsão com IC 95%
previsao_random_walk_log <- forecast(modelo_random_walk_log, h = 30, level = 95)
plot(c(tail(datas_treino, 120), datas_teste),
     c(tail(as.numeric(ts_preço_ajustado_treino), 120), as.numeric(ts_preço_ajustado_teste)),
     type = "l", main = "BBAS3 — previsão de 30 pregões (random walk no log)", xlab = "Data", ylab = "Preço ajustado (R$)",
     ylim = range(c(tail(as.numeric(ts_preço_ajustado_treino), 120), as.numeric(ts_preço_ajustado_teste),
                    previsao_random_walk_log$lower, previsao_random_walk_log$upper)))
lines(datas_teste, previsao_random_walk_log$mean, col = "red", lwd = 2)  # previsão pontual
lines(datas_teste, previsao_random_walk_log$lower, col = "red", lty = 2) # limite inferior IC 95%
lines(datas_teste, previsao_random_walk_log$upper, col = "red", lty = 2) # limite superior IC 95%
abline(v = datas_teste[1], col = "gray", lty = 3)                        # início do teste
legend("topleft", legend = c("Observado", "Previsão", "IC 95%"), col = c("black", "red", "red"), lty = c(1, 1, 2), cex = 0.8)
salvar_grafico("06_previsao_30_pregoes")
# Quantos dos 30 preços reais ficaram dentro do IC 95%? (esperado ~28 a 30)
sum(as.numeric(ts_preço_ajustado_teste) >= as.numeric(previsao_random_walk_log$lower) &
      as.numeric(ts_preço_ajustado_teste) <= as.numeric(previsao_random_walk_log$upper)) # 11 de 30. IC supõe normalidade, mas os retornos têm caudas pesadas (Jarque-Bera)

# Erro antes e depois do choque da eleição (pregões 1 a 11: até 23/10; 12 a 30: a partir de 26/10)
erro_random_walk_log <- as.numeric(ts_preço_ajustado_teste) - as.numeric(previsao_random_walk_log$mean)
c(MAPE_antes_eleicao  = 100 * mean(abs(erro_random_walk_log[1:11]  / as.numeric(ts_preço_ajustado_teste)[1:11])),
  MAPE_depois_eleicao = 100 * mean(abs(erro_random_walk_log[12:30] / as.numeric(ts_preço_ajustado_teste)[12:30])))