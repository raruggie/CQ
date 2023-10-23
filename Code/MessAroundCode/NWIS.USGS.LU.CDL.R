
df.sf.NWIS.keep.2<-left_join(df.sf.NWIS.keep, distinct(df_Seg.2, site, .keep_all = T)%>%select(.,c(site, Type, n_sample_rank)), by = c('Name'='site'))%>%
  select(Name, Type, n_sample_rank)%>%
  arrange(n_sample_rank)%>%
  mutate(n_sample_rank=1:nrow(.))%>%
  left_join(.,df.NWIS.TP%>%group_by(site_no)%>%summarise(n=n()), by = c('Name'='site_no'))%>%
  mutate(NEW = case_when(n < 50 ~ .05,
                         n >=50 & n < 100 ~ .25,
                         n >=100 & n < 500 ~ .5,
                         n >=500 ~ .9))
# map:

# mapview(df.sf.NWIS.keep.2, zcol = 'Type', alpha.regions = 'NEW')

#### Categorizing land use ####

# USGS criteria:
# Agricultural sites have >50% agricultural land and ≤5% urban land;
# urban sites have >25% urban and ≤25% agricultural land; 
# undeveloped sites have ≤ 5% urban and ≤ 25% agricultural land; 
# all other combinations of urban, agricultural, and undeveloped lands are classified as mixed

# IN a first pass using these thresholds the number of ag and urban sites wasvery low
# I will play with these numbers to see what happens

# merge land use from df.NWIS.CDL to df.sf.NWIS.keep.2:

df.sf.NWIS.keep.2<-left_join(df.sf.NWIS.keep.2, df.NWIS.CDL, by = 'Name')

# combine Ag and Pasture into a single landuse for Ag:

df.sf.NWIS.keep.2<-mutate(df.sf.NWIS.keep.2, Ag = Ag+Pasture)

# set NA to zero

df.sf.NWIS.keep.2[is.na(df.sf.NWIS.keep.2)]<-0

# create the land use class column based on USGS critiera:

df.NWIS.USGS.LU<-df.sf.NWIS.keep.2%>%
  mutate(USGS.LU = 'Mixed')%>%
  mutate(USGS.LU = case_when(.default = 'Mixed',
                             Ag > .30 & Developed <= .1 ~ 'Agriculture',
                             Developed > .1 & Ag <= .3 ~ 'Urban',
                             Developed <= .1 & Ag <= .1 ~ 'Undeveloped'
  )
  )

# merge the OLS and Sens slopes and intercepts with this df:

df.NWIS.USGS.LU<-left_join(df.NWIS.USGS.LU, df.OLS_Sens[,1:5], by = 'Name')

# pivot longer for geom_box + facet:

data<-df.NWIS.USGS.LU%>%
  pivot_longer(cols = 15:18, names_to = 'CQ_parameter', values_to = 'Value')%>%
  mutate(USGS.LU=factor(USGS.LU))

# make ggplot:

my_xlab <- paste(levels(factor(df.NWIS.USGS.LU$USGS.LU)),"\n(N=",table(factor(df.NWIS.USGS.LU$USGS.LU)),")",sep="")

ggplot(data, aes(x=USGS.LU, y=Value))+
  geom_boxplot(varwidth = TRUE, alpha=0.2)+
  scale_x_discrete(labels=my_xlab)+
  facet_wrap('CQ_parameter', scales = 'free')+
  stat_compare_means(method = "anova", label.y = 2)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5") +
  ggtitle('Adjusted USGS Thresholds using Aggregated 2020 CDL')

