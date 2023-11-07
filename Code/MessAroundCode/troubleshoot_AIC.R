library(MASS)
library(tidyverse)

# load in data:

load('Processed_Data/troubleshoot_AIC.Rdata')

# create list of correlation dfs:

l.cor.MLR<-df.cor %>%
  split(., df.cor$CQ_Parameter)%>%
  lapply(., \(i) i%>% 
           arrange(Spearman_Correlation)%>%  
           # slice_head(n = 7)%>%
           # bind_rows(i%>%arrange(desc(Spearman_Correlation))%>%slice_head(n = 7))%>%
           # distinct(term, .keep_all = T)%>%
           # filter(!between(Spearman_Correlation, -0.25,.25))%>%
           mutate(sig_0.05 = factor(sig_0.05, levels = c('not', 'sig')))%>%
           # mutate(term = factor(term, levels = unique(term[order(Spearman_Correlation)])))%>%
           # filter(sig_0.05 == 'sig')%>%
           as.data.frame())


x<-l.cor.MLR[[1]]

# create a list of full model dataframes from df.OLS.Sens and subseting the attributes using the names in each l.cor[[i]]$term:

l.cor.MLR.full<-lapply(1:4, \(i) df.OLS_Sens%>%
                         dplyr::select(i+1, l.cor.MLR[[i]]$term)%>%
                         as.data.frame()%>%
                         rename(term = 1))

x<-l.cor.MLR.full[[1]]

# loop through data frames and make lm and stepAIC objects (wont work on lapply for some reason):

# set up list for aic objects to append into:

l.aic<-list()

# loop: takes 10 seconds

# i<-1

for (i in 1:4){
  
  # find the best predictor name for the simple model:
  
  n<-l.cor.MLR[[i]]$term[which.max(abs(l.cor.MLR[[i]]$Spearman_Correlation))]
  
  # createformula with this predictor:
  
  n<-paste('term ~', n)
  
  # create simple model:
  
  m.simple<-lm(n,l.cor.MLR.full[[i]])
  
  # create full model:
  
  m.full<-lm(term~., l.cor.MLR.full[[i]])
  
  # create AIC object:
  
  aic<-stepAIC(m.simple, scope = list(upper=m.full, lower =~1), direction = 'both', trace = FALSE)
  
  # save aic object to list:
  
  l.aic[[i]]<-aic
  
}

# look at summary of aic objects:

lapply(l.aic, summary)

# set names of list elements:

names(l.aic)<-names(df.OLS_Sens)[2:5]

# extract best predictors for each parameter:

l<-lapply(l.aic, \(i) names(i$model)[-1])

# look at univariate plots of these. to do this:

# make a list of ggplot objects for each parameter:

l.MLR.plots<-lapply(1:4, \(i) df.OLS_Sens%>%
                      select(i+1, l[[i]])%>%
                      as.data.frame()%>%
                      pivot_longer(cols = 2:last_col(), names_to = 'Attribute', values_to = 'value')%>%
                      mutate(Attribute = factor(Attribute, levels = l[[i]]))%>%
                      ggplot(., aes(y = !!sym(names(df.OLS_Sens)[i+1]), x = value))+
                      facet_wrap('Attribute', scales = 'free')+
                      geom_point()+
                      geom_smooth(method = 'lm')+
                      ggtitle((names(df.OLS_Sens)[i+1]))
                    
                    
)

#plot:

l.MLR.plots

#