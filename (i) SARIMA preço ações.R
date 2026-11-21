library(quantmod) #importar dados BB
library(tseries) #para séries temporais
library(forecast) #para auto.arima
library(urca) # para raíz unitaria e teste de Dickey-Fuller Aumentado (ADF)

# baseado código exc 4 da lista 3

####### 1. Importação #########
bb <- getSymbols("BBAS3.SA", src = "yahoo",
                 from = "2005-01-01", to = "2026-10-21",
                 auto.assign = FALSE)

# Salvar as datas (que são o índice dos dados) Usado para plotar no eixo 
# horizontal o tempo e não o número dos pregões
datas <- index(bb)

# Série temporal
dados <- ts(as.numeric(bb$BBAS3.SA.Adjusted))

# Reservar os últimos 30 dias para avaliar as previsões.
n <- length(dados) # Série com 5423 observações
treino <- ts(as.numeric(dados[1:(n - 30)])) # série de treino, tirando as ultimas 30 observações
teste <- ts(as.numeric(dados[(n - 29):n])) # Série teste, começando pelas ultimas 30 observações. Começo é inclusivo (n-29)


# E as cortar das datas para montar o gráfico de treino e teste
datas_treino <- datas[1:(n - 30)]
datas_teste <- datas[(n - 29):n]

####### 2. Visualização #############
par(mfrow=c(1,2)) # Muda como o R mostra os gráficos. Neste caso, 1 linha e 2 colunas

#plot(dados eixo x, dados eixo y, tipo do gráfico em linhas, main é o título, xlab é o nome do eixo x, y lab eixo y)
plot(datas_treino, treino, type = "l",
     main = "BBAS3", xlab = "Data", ylab = "Valor preço ajustado")
# Preços sobem com queda na pandemia
# Certamente há efeito de inflação, mas não vamos deflacionar pq o enunciado não pede... E o Claude falou que o padrão é assim

# Fac do modelo
acf(treino, main="FAC — nível")
# FAC praticamente constante, série NÃO ESTACIONÁRIA.

########## 3. auto.arima ####################
modelo <- auto.arima(treino)
modelo # Retornou ARIMA(0,1,0). É um random-walk puro, sem drift

####### 4. Justificar com robustez pq é um random-walk #########

# 4.1 Mostrar que a série não é estacionária e que d = 1.
#autoplot(modelo)
adf <- ur.df(dados) #  Dickey-Fuller Aumentado -> Verificar random-walk
summary(adf) # Valor crítco tau1 em 5pct é -1.95. Valor de 0.7493 não é maior
# Não rejeita H0. Série é Random-Walk

kpss <- ur.kpss(dados) # KPSS -> Verificar se I(1)
summary(kpss) # Valor crítico 5% é 0.463, valor teste 39.9176.
# Rejeita H0. Diferencia 1 vez. d = 1






