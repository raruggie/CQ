
# Picks up from NWIS.R:

df.sf.NWIS.keep.2<-left_join(df.sf.NWIS.keep, distinct(df_Seg.2, site, .keep_all = T)%>%select(.,c(site, Type, n_sample_rank)), by = c('Name'='site'))%>%
  select(Name, Type, n_sample_rank)%>%
  arrange(n_sample_rank)%>%
  mutate(n_sample_rank=1:nrow(.))%>%
  left_join(.,df.NWIS.TP%>%group_by(site_no)%>%summarise(n=n()), by = c('Name'='site_no'))%>%
  mutate(NEW = case_when(n < 50 ~ .05,
                         n >=50 & n < 100 ~ .25,
                         n >=100 & n < 500 ~ .5,
                         n >=500 ~ .9))



# merge land use from df.NWIS.CDL to df.sf.NWIS.keep.2:

df.sf.NWIS.keep.2<-left_join(df.sf.NWIS.keep.2, l.NWIS.NLCD$`2001`, by = 'Name')

# set NA to zero

df.sf.NWIS.keep.2[is.na(df.sf.NWIS.keep.2)]<-0

# combine Ag and Pasture into a single landuse for Ag
# combine developed into single for developed 

df.sf.NWIS.keep.2<-mutate(df.sf.NWIS.keep.2, Ag = `Cultivated Crops` + `Pasture/Hay`)%>%
  rowwise()%>%
  mutate(Developed = sum(c_across(starts_with("Developed")), na.rm = T))

# create the land use class column based on USGS critiera:

df.NWIS.USGS.LU<-df.sf.NWIS.keep.2%>%
  mutate(USGS.LU = 'Mixed')%>%
  mutate(USGS.LU = case_when(.default = 'Mixed',
                             Ag > .50 & Developed <= .05 ~ 'Agriculture',
                             Developed > .25 & Ag <= .25 ~ 'Urban',
                             Developed <= .05 & Ag <= .25 ~ 'Undeveloped')
  )

# merge the OLS and Sens slopes and intercepts with this df:

df.NWIS.USGS.LU<-left_join(df.NWIS.USGS.LU, df.OLS_Sens[,1:5], by = 'Name')

# pivot longer for geom_box + facet:

data<-df.NWIS.USGS.LU%>%
  pivot_longer(cols = 25:28, names_to = 'CQ_parameter', values_to = 'Value')%>%
  mutate(USGS.LU=factor(USGS.LU))

# make ggplot:

my_xlab <- paste(levels(factor(df.NWIS.USGS.LU$USGS.LU)),"\n(N=",table(factor(df.NWIS.USGS.LU$USGS.LU)),")",sep="")

ggplot(data, aes(x=USGS.LU, y=Value))+
  geom_boxplot(varwidth = TRUE, alpha=0.2)+
  scale_x_discrete(labels=my_xlab)+
  facet_wrap('CQ_parameter', scales = 'free')+
  stat_compare_means(method = "anova", label.y = max(data$Value))+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5") +
  ggtitle('Orginal USGS Thresholds using Aggregated 2001 NLCD')

