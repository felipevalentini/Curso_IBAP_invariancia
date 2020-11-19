#configurar o diretório de trabalho
setwd("C:/Users/X-Casa/OneDrive/congressos_cursos/IBAP/IBAP_jornada2020")

#carregar os pacotes que serão utilizados na análise
library(lavaan)
library(semTools)

# exemplo Multigrupos simples
bancoPerson<-read.csv2("banco5Person.csv")
# 5 facetas da abertura 
# imaginação, estética, curiosidade, preferência pelo novo e ideias
# todas já padronizadas em média 50 e DP 10

#################
#modelar o básico
#################
modelo_Person_Basic <-"
f1 =~ i1 + i2 + i3 + i4 + i5
"

#rodar o modelo
fit_modelo_Person_Basic <- cfa(modelo_Person_Basic, bancoPerson)
summary(fit_modelo_Person_Basic, fit.measures=TRUE,
        standardized=TRUE)


#modelo configural
fit_configural <- cfa(model=modelo_Person_Basic, 
                      data=bancoPerson,
                      group="grupo")
summary(fit_configural)

#modelo metrico
fit_metrico <- cfa(model=modelo_Person_Basic, 
                      data=bancoPerson,
                      group="grupo",
                      group.equal = c("loadings"))
summary(fit_metrico)


#modelo escalar
fit_escalar <- cfa(model=modelo_Person_Basic, 
                   data=bancoPerson,
                   group="grupo",
                   group.equal = c("loadings", "intercepts"))
                
summary(fit_escalar)                   

#comparar os modelos de invariancia
compar_AFMG_Persona <- compareFit(fit_configural, fit_metrico,
                                  fit_escalar)
summary(compar_AFMG_Persona, 
        fit.measures=c("cfi","rmsea","srmr","mfi"))



#Semtools oferece um modo automático de testar a inariância

fit_Person_autom <- measurementInvariance(model=modelo_Person_Basic, 
                                          data=bancoPerson,
                                          group="grupo")





##################################
# para dados categóricos/ordinais
#################################

# importar o banco em txt/dat
banco <- read.table("banco_1.dat") 
names(banco)<-c("it1","it2","it3","it4","it5","it6","it7","it8","it9","it10",
                "idade", "sexo")

modelo2<- "
f1 =~ it1 + it2 + it3 +  it4 +  it5
f2 =~ it6 + it7 + it8 +  it9 +  it10"


#configural
configural_sexo <-measEq.syntax(configural.model = modelo2,
                                data=banco,
                                ordered = c("it1","it2","it3","it4","it5","it6","it7","it8","it9","it10"), 
                                estimator="WLSMV",
                                group="sexo", return.fit=TRUE,
                                ID.cat="Wu.Estabrook.2016",
                                ID.fac="std.lv", parameterization="delta")
summary(configural_sexo, fit.measures=T, standardized=T)


# thresholds (proposicao 4)
threshold_sexo <-measEq.syntax(configural.model = modelo2,
                                data=banco,
                                ordered = c("it1","it2","it3","it4","it5","it6","it7","it8","it9","it10"), 
                                estimator="WLSMV",
                                group="sexo", return.fit=TRUE,
                                ID.cat="Wu.Estabrook.2016",
                                ID.fac="std.lv", parameterization="delta",
                                group.equal="thresholds")
summary(threshold_sexo, fit.measures=T, standardized=T)


# cargas fatoriais (proposicao 7)
cargas_sexo <-measEq.syntax(configural.model = modelo2,
                               data=banco,
                               ordered = c("it1","it2","it3","it4","it5","it6","it7","it8","it9","it10"), 
                               estimator="WLSMV",
                               group="sexo", return.fit=TRUE,
                               ID.cat="Wu.Estabrook.2016",
                               ID.fac="std.lv", parameterization="delta",
                               group.equal=c("thresholds", "loadings"))
summary(cargas_sexo, fit.measures=T, standardized=T)


# escalar (proposicao 7 + interceptos)
escalar_sexo <-measEq.syntax(configural.model = modelo2,
                            data=banco,
                            ordered = c("it1","it2","it3","it4","it5","it6","it7","it8","it9","it10"), 
                            estimator="WLSMV",
                            group="sexo", return.fit=TRUE,
                            ID.cat="Wu.Estabrook.2016",
                            ID.fac="std.lv", parameterization="delta",
                            group.equal=c("thresholds","loadings", 
                                          "intercepts"))
summary(escalar_sexo, fit.measures=T, standardized=T)


#comparar modelos
invariancia_sexo_WU <- compareFit(configural_sexo, threshold_sexo, cargas_sexo, 
                               escalar_sexo)
summary(invariancia_sexo_WU , 
        fit.measures=c("cfi.scaled","rmsea.scaled","srmr","mfi"))




######################
# invariancia parcial
######################
# vimos que não podemos concluir pela invariância das cargas e dos interceptos
# vamos tentar descobrir quais parâmetros específicos não são invariântes
# e tentar um modelo de invariância parcial

# solicitar indicadores de modificação (cargas)
mi_cargas <- modindices(cargas_sexo, free.remove = FALSE)
mi_cargas[mi_cargas$op == "=~",] %>% arrange(desc(mi))
# indicadores de modificação apontam para as cargas dos itens 3 e 10


# cargas fatoriais (parcial - itens 3 e 10 livres)
cargas_sexo_parcial <-measEq.syntax(configural.model = modelo2,
                            data=banco,
                            ordered = c("it1","it2","it3","it4","it5","it6","it7","it8","it9","it10"), 
                            estimator="WLSMV",
                            group="sexo", return.fit=TRUE,
                            ID.cat="Wu.Estabrook.2016",
                            ID.fac="std.lv", parameterization="delta",
                            group.equal=c("thresholds", "loadings"),
                            group.partial = c("f1  =~ it3", "f2  =~ it10"))
summary(cargas_sexo_parcial, fit.measures=T, standardized=T)



# solicitar indicadores de modificação (interceptos)
mi_interceptos <- modindices(escalar_sexo, free.remove = FALSE)
mi_interceptos[mi_interceptos$op == "~1",] %>% arrange(desc(mi))
# indicadores de modificação apontam para o intercepto do item 3


# interceptos (parcial - intercepto do item 4 livre)
escalar_sexo_parcial <-measEq.syntax(configural.model = modelo2,
                                    data=banco,
                                    ordered = c("it1","it2","it3","it4","it5","it6","it7","it8","it9","it10"), 
                                    estimator="WLSMV",
                                    group="sexo", return.fit=TRUE,
                                    ID.cat="Wu.Estabrook.2016",
                                    ID.fac="std.lv", parameterization="delta",
                                    group.equal=c("thresholds", "loadings", "intercepts"),
                                    group.partial = c("f1  =~ it3", "f2  =~ it10",
                                                       "it4 ~1"))
summary(escalar_sexo_parcial, fit.measures=T, standardized=T)

#vamos comparar novamente
invariancia_sexo_WU_parcial <- compareFit(configural_sexo, threshold_sexo, cargas_sexo_parcial, 
                                  escalar_sexo_parcial)
summary(invariancia_sexo_WU_parcial , 
        fit.measures=c("cfi.scaled","rmsea.scaled","srmr","mfi"))




