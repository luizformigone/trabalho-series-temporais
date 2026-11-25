############################################################
# TRABALHO FINAL — PARTE (i): PREÇO DA BBAS3 (Box-Jenkins)
############################################################
library(tseries)
library(forecast)
library(car)

# Leitura: preço de fechamento diário da BBAS3, em R$ (03/01/2022 a 13/11/2026, 1216 pregões).
bd <- read.csv("bbas3.csv") # pode pular se já importou pelo RStudio
dados <- ts(bd$Fechamento)

# Reservar os últimos 30 pregões para avaliar as previsões (enunciado).
treino <- head(dados, length(dados)-30) # 1186 pregões, até 30/09/2026
treino <- ts(as.numeric(treino))
teste <- tail(as.numeric(dados),30) # 01/10/2026 a 13/11/2026
datas <- as.Date(bd$Data) # datas dos pregões, para o eixo dos gráficos
datas_treino <- head(datas, length(datas)-30)
datas_teste <- tail(datas, 30)

#----------------------------------------------------------#
# a) Identificação e estimação dos modelos
#----------------------------------------------------------#

# 1. Comparar o nível e a primeira diferença.
par(mfrow=c(2,3))
plot(treino, main="BBAS3", ylab="R$")
acf(treino, main="FAC - nível")
pacf(treino, main="FACP - nível")
plot(diff(treino), main="Primeira diferença", ylab="Diferença (R$)")
acf(diff(treino), main="FAC - diferença")
pacf(diff(treino), main="FACP - diferença")

# Vai precisar integrar pq FAC - nível decai muito devagar
# Na diferença, FAC e FACP ficam dentro das bandas: nem AR nem MA (cara de random walk)

# 2. Complementar a inspeção gráfica com testes.
# ADF: H0 = raiz unitária. KPSS: H0 = estacionariedade.
adf.test(treino)                     # p = 0,49: não rejeita raiz unitária
adf.test(diff(treino))               # p < 0,01: rejeita, a diferença é estacionária
kpss.test(treino, null="Level")      # p < 0,01: rejeita estacionariedade
kpss.test(diff(treino), null="Level") # p > 0,10: não rejeita, a diferença é estacionária

# p-valor sem integrar deu ruim nos dois testes
# p-valor dps de integrar só uma vez deu o suficiente

# Hipótese com base nos testes: d=1.
d <- 1

# 4. Candidatos parcimoniosos: justificar pelas FAC/FACP.
# FACP cortando sugere AR; FAC cortando sugere MA; ambas decaindo, ARMA.
# Aqui nenhuma das duas tem defasagem relevante na diferença -> random walk é o principal candidato
modelo1 <- arima(treino, order=c(0,d,0), method="ML") # random walk
modelo2 <- arima(treino, order=c(1,d,0), method="ML")
modelo3 <- arima(treino, order=c(0,d,1), method="ML")
modelo4 <- auto.arima(treino) # ARIMA(0,1,0): também escolhe o random walk
modelo5 <- arima(ts(treino, frequency=5), order=c(0,d,0), seasonal = list(order = c(1,0,1), period = 5), method="ML") # sazonalidade semanal (5 pregões)
modelo6 <- arima(treino, order=c(2,d,2), method="ML") # melhor pelo AIC numa grade p,q de 0 a 3
modelo1
modelo2
modelo3
modelo4
modelo5
modelo6 # aviso de convergência; AR e MA quase iguais com sinal trocado: as raízes se cancelam (superparametrizado)


# Candidato que será examinado nos itens b/c.
modelo <- modelo1 # rodar mesmo bloco em baixo e trocar

#----------------------------------------------------------#
# b) Significância dos coeficientes
#----------------------------------------------------------#
estimativas <- coef(modelo)
erros_padrao <- sqrt(diag(modelo$var.coef))
z <- estimativas/erros_padrao
data.frame(estimativas, erros_padrao, p_valor=2*pnorm(-abs(z)))
# Repetir para os outros candidatos.
# modelo1: random walk não tem coeficiente para testar
# modelo2: ar1 = 0,044 (p = 0,13). Não significativo
# modelo3: ma1 = 0,045 (p = 0,13). Não significativo
# modelo5: sar1 e sma1 ≈ 0 (p ≈ 1,00). Sem sazonalidade semanal
# modelo6: p ≈ 0 em tudo, mas é ilusão: as raízes AR e MA se cancelam

#----------------------------------------------------------#
# c) Diagnóstico e capacidade de previsão
#----------------------------------------------------------#
tsdiag(modelo, gof.lag=24)
residuos <- as.numeric(residuals(modelo))
padronizados <- residuos/sqrt(modelo$sigma2)

par(mfrow=c(2,2))
plot(padronizados, type="l", main="Resíduos padronizados")
abline(h=0)
acf(residuos, main="FAC dos resíduos")
acf(residuos^2, main="FAC dos resíduos²")
qqPlot(padronizados, distribution="norm", id=FALSE)

# Resíduos sem autocorrelação (FAC dentro das bandas)
# QQ-plot foge da reta nas pontas: caudas pesadas. Jarque-Bera confirma (p < 0,001)
jarque.bera.test(residuos)

# p+q: número de parâmetros AR/MA do candidato examinado.
numero_parametros <- 0 # random walk: nenhum (modelo2 e modelo3: 1; modelo6: 4)
Box.test(residuos, lag=12, type="Ljung-Box", fitdf=numero_parametros) # p = 0,79
Box.test(residuos, lag=24, type="Ljung-Box", fitdf=numero_parametros) # p = 0,73
# H0: ausência de autocorrelação. Não rejeita: o random walk passa (modelos 2, 3 e 6 também passam)

previsao <- forecast(modelo, h=30) # usando o modelo 1
par(mfrow=c(1,1))
# Zoom: o modelo continua estimado sem os 30 últimos; só o gráfico mostra menos pregões do treino
pregoes_zoom <- 60 # quantos pregões do treino aparecem antes da previsão (diminua para dar mais zoom)
plot(tail(datas_treino, pregoes_zoom), tail(as.numeric(treino), pregoes_zoom), type="l",
     xlim=range(c(tail(datas_treino, pregoes_zoom), datas_teste)), # eixo x em datas: últimos pregões do treino + 30 do teste
     ylim=range(c(tail(treino, pregoes_zoom), teste, previsao$lower, previsao$upper)), # eixo y: cabe o teste e o IC
     main="Previsão de 30 pregões x observado", ylab="R$", xlab="", xaxt="n")
axis.Date(1, at=seq(min(tail(datas_treino, pregoes_zoom)), max(datas_teste), by="2 weeks"), format="%d/%m") # datas no formato dia/mês
polygon(c(datas_teste, rev(datas_teste)), c(previsao$lower[,2], rev(previsao$upper[,2])), col="gray85", border=NA) # IC 95%
polygon(c(datas_teste, rev(datas_teste)), c(previsao$lower[,1], rev(previsao$upper[,1])), col="#B0B5D0", border=NA) # IC 80%
lines(datas_teste, previsao$mean, col="#2297E6", lwd=2) # previsão
lines(datas_teste, teste, col="red", lwd=2) # observado no teste
abline(v=datas_teste[1], lty=3) # início do teste
legend("bottomleft", c("Treino","Previsão","Observado (teste)"), col=c("black","#2297E6","red"), lwd=2, bty="n")
# A queda de ~29% em 26/10/2026 sai do intervalo de confiança: nenhum modelo prevê esse choque

# Comparar com a previsão ingênua: último valor observado.
erro <- teste-as.numeric(previsao$mean)
erro_ingenuo <- teste-tail(as.numeric(treino),1)
data.frame(modelo=c("ARIMA examinado","Ingênuo"),
           RMSE=c(sqrt(mean(erro^2)),sqrt(mean(erro_ingenuo^2))),
           MAE=c(mean(abs(erro)),mean(abs(erro_ingenuo))))
# Iguais (RMSE = 4,55; MAE = 3,63): a previsão do random walk É a previsão ingênua
# Repetir para candidatos definidos no treino.
# modelo2: RMSE 4,57 | modelo3: RMSE 4,57 | modelo6: RMSE 4,55. Nenhum supera o ingênuo

#----------------------------------------------------------#
# d) Comparação por AIC e BIC
#----------------------------------------------------------#
data.frame(modelo=c("ARIMA(0,d,0)","ARIMA(1,d,0)","ARIMA(0,d,1)","SARIMA(0,d,0)(1,0,1)[5]","ARIMA(2,d,2)"),
           AIC=c(AIC(modelo1),AIC(modelo2),AIC(modelo3),AIC(modelo5),AIC(modelo6)),
           BIC=c(BIC(modelo1),BIC(modelo2),BIC(modelo3),BIC(modelo5),BIC(modelo6)))
# Menor é preferível, sem ignorar o diagnóstico.
# AIC: todos empatados (diferenças < 1). BIC: random walk ganha com folga (1082,7)

# ARIMA(0,d,0) - random walk - melhor modelo: mais simples, passa no diagnóstico e prevê igual aos outros
# Compatível com a hipótese de eficiência de mercado (forma fraca)

