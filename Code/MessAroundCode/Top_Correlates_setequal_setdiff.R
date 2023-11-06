# determine differences in top correlates between OLS and SEns:

# intercepts, neg

OLS<-l.cor[[1]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation<0)
Sen<-l.cor[[3]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation<0)
intersect(OLS$term, Sen$term)
setdiff(OLS$term, Sen$term)
setdiff(Sen$term,OLS$term)

# intercepts, pos

OLS<-l.cor[[1]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation>0)
Sen<-l.cor[[3]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation>0)
intersect(OLS$term, Sen$term)
setdiff(OLS$term, Sen$term)
setdiff(Sen$term,OLS$term)



# slope:

# neg:

OLS<-l.cor[[2]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation<0)
Sen<-l.cor[[4]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation<0)
intersect(OLS$term, Sen$term)
setdiff(OLS$term, Sen$term)
setdiff(Sen$term,OLS$term)

# pos:

OLS<-l.cor[[2]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation>0)
Sen<-l.cor[[4]]%>%arrange(desc(Spearman_Correlation))%>%
  filter(Spearman_Correlation>0)
intersect(OLS$term, Sen$term)
setdiff(OLS$term, Sen$term)
setdiff(Sen$term,OLS$term)








