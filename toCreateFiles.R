library("dplyr")
library("ggplot2")
#library("plyr")
library("igraph")
library("tnet")
library("Matrix")
library("cowplot")
library("tidyverse")

remove(list = ls())

setwd("~/Projects/ScientificDataFiles/")

################ Global functions for creation of directed - weighted network #########################
#Function to create directed weighted networks. Need a file with 3 columns"
#Source, Target and Weight
networkCreationDW<-function(x){
  el1 <-as.data.frame(x) #Verify that it's a data.frame
  ell3<-as.matrix(el1)
  ell3[,1]=as.character(ell3[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
  ell3[,2]=as.character(ell3[,2])
  ggg3=graph.edgelist(ell3[,1:2],directed = TRUE) #Create the directed network
  E(ggg3)$weight=as.numeric(ell3[,3]) #weights 
  ggg3<-simplify(ggg3) #Verify that there are not multiple connection or remaining loops
  compsG<-clusters(ggg3) #Check if the network is broken
  g<-induced.subgraph(ggg3,compsG$membership==which.max(compsG$csize)) #Only select the largest CC
  #If the network is not broken, keep the networks as is.
}

# These are for individual characterization directed-weighted network (all but clustering)
nodeStatsDW<-function(x){
  V(x)$Eig<-eigen_centrality(x,directed = T,weights = E(x)$weight)$vector
  V(x)$InDegree<-degree(x,mode = "in")
  V(x)$OutDegree<-degree(x,mode = "out")
  V(x)$InStrength<-strength(x,mode="in",weights = E(x)$weight)
  V(x)$OutStrength<-strength(x,mode="out",weights = E(x)$weight)
  V(x)$Betweenness<-betweenness(x, directed=T, weights=E(x)$weight,nobigint = TRUE, normalized = FALSE)
  V(x)$InCoreness<-coreness(x,mode="in")
  V(x)$OutCoreness<-coreness(x,mode="out")
  centrality <- data.frame(Node = as.numeric(V(x)),
                           OutDegree    = V(x)$OutDegree,
                           InDegree    = V(x)$InDegree,
                           OutStrength = V(x)$OutStrength,
                           InStrength   = V(x)$InStrength,
                           OutCoreness = V(x)$OutCoreness,
                           InCoreness = V(x)$InCoreness,
                           Betweenness  = V(x)$Betweenness,
                           EigenCent = V(x)$Eig)
  centrality
}

########################################################################################################

#Reading file with raw shipments data
rawDataShipments<-read.csv("shipments.csv",header = T)
head(rawDataShipments)

edgeListCompleteNetwork<-read.csv("edgelist.csv",header=T)
head(edgeListCompleteNetwork)

shipmentsNetwork<-networkCreationDW(edgeListCompleteNetwork)

grados<-degree(shipmentsNetwork,V(shipmentsNetwork),mode="in")
degree.df <- data.frame(table(grados=factor(grados, levels=seq_len(max(grados)))))
degree.df$grados <- as.numeric(as.character(degree.df$grados))
degree.df$Mode<-"In-degree"

grados<-degree(shipmentsNetwork,V(shipmentsNetwork),mode="out")
degree.df1 <- data.frame(table(grados=factor(grados, levels=seq_len(max(grados)))))
degree.df1$grados <- as.numeric(as.character(degree.df1$grados))
degree.df1$Mode<-"Out-degree"

degree.both<-rbind(degree.df,degree.df1)
str(degree.both)

cbPaletteDos <- c("#0072B2","#D55E00")
degDist<-ggplot(degree.both, aes(x=grados, y=Freq,color=factor(Mode))) +
  theme_bw()+
  geom_point(size=3)+
  scale_x_log10()+ scale_y_log10()+
  theme(legend.title = element_blank(),legend.position = c(0.8,0.8),
        text=element_text(size=20))+
  xlab("Degree")+ylab("Frequency")+
  scale_color_manual(values = cbPaletteDos)

library("sp")
#Districts
TUR2<-readRDS("gadm36_TUR_2_sp.rds")

TUR_for<-fortify(TUR2)

#All source epiunits
tt<-table(rawDataShipments$Source) 
names(tt[tt==max(tt)]) #Epiunit with most shipments sent
#Select all destinations for epiunit with most shipments sent
mostSends<-rawDataShipments %>% 
  subset(Source!=Destination) %>% 
  subset(Source==names(tt[tt==max(tt)])) %>%
  unique()
head(mostSends)
#All destination epiunits
tt1<-table(rawDataShipments$Destination) 
names(tt1[tt1==max(tt1)]) #Epiunit with most shipments received
#Select of shipments that had epiunit with most shipments received
mostRecived<-rawDataShipments %>% 
  subset(Source!=Destination) %>% 
  subset(Destination==names(tt1[tt1==max(tt1)])) %>%
  unique()
head(mostRecived)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

outD<-ggplot(TUR_for) +
  theme_void()+
  geom_polygon( aes(x = long, y = lat, group = group),
                color = "black", fill = "white") +
  geom_curve(data=mostSends,
             aes(x=Longitude_source,y=Latitude_source, xend = Longitude_destination, yend = Latitude_destination),
             curvature = 0.1, color="green4",
             alpha = 0.2)+
  geom_point(data=mostSends,
             aes(x=Longitude_destination,y=Latitude_destination),color="#E69F00",size=1)+
  geom_point(data=mostSends,
             aes(x=Longitude_source,y=Latitude_source),color="black",size=3)+
  theme(panel.background = element_rect(fill = NA))

inD<-ggplot(TUR_for) +
  theme_void()+
  geom_polygon( aes(x = long, y = lat, group = group),
                color = "black", fill = "white") +
  geom_curve(data=mostRecived,
             aes(x=Longitude_source,y=Latitude_source, xend = Longitude_destination, yend = Latitude_destination),
             curvature = 0.1, color="green4",
             alpha = 0.2)+
  geom_point(data=mostRecived,
             aes(x=Longitude_source,y=Latitude_source),color="#E69F00",size=1)+
  geom_point(data=mostRecived,
             aes(x=Longitude_destination,y=Latitude_destination),color="black",size=3)+
  theme(panel.background = element_rect(fill = NA))


mapasRed<-plot_grid(outD,inD,ncol=1,labels = "auto",label_size = 23) 

fig1<-plot_grid(mapasRed,degDist,ncol=2,rel_widths = c(1.8,1),labels = c("","c"),
                label_size = 23)

ggsave(fig1,file="Figure1.png",width = 16,height=8)

### Outbreak data

outbreaks<-read.csv("outbreaks.csv",header = T)

outbreaksFreq<-outbreaks %>% select(-Latitude,-Longitude) %>% 
  group_by(epiunitID,serotype) %>%
  dplyr::summarise(Outbreaks=n())

outbreaksFreq$serotype<-ifelse(outbreaksFreq$serotype=="1","Asia-1",outbreaksFreq$serotype)
outbreaksFreq$serotype<-factor(outbreaksFreq$serotype,levels = c("Asia-1","A","O"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

juntos<-outbreaks %>% select(-serotype,-Latitude,-Longitude) %>%
  group_by(epiunitID) %>% 
  dplyr::summarise(Outbreaks=n()) %>% 
  ggplot(aes(x=as.factor(Outbreaks),y=..count..))+
  theme_bw()+
  geom_histogram(stat = "count",color="#D55E00",fill="#D55E00") +
  theme(text=element_text(size=18))+
  xlab("Number of outbreaks") +
  ylab("Frequency")

separados<-ggplot(outbreaksFreq,aes(x=as.factor(Outbreaks),fill=factor(serotype)))+
  geom_histogram(stat = "count")+
  facet_wrap(~serotype) +
  theme_bw() +
  scale_fill_manual(values=cbPalette)+
  scale_color_manual(values = cbPalette)+
  labs(fill="Serotype",color="")+
  theme(text=element_text(size=18),legend.position = c(0.1,0.8))+
  xlab("Number of outbreaks") +
  ylab("Frequency")

uno<-plot_grid(juntos,separados,rel_widths = c(1,3),labels = c("a","b"),
               label_size = 30)

mapas<-ggplot(TUR_for) +
  theme_void()+
  geom_polygon( aes(x = long, y = lat, group = group),
                color = "gray", fill = "white") +
  geom_polygon( data=TUR_for0,aes(x = long, y = lat, group = group),
                color = "black", fill = NA) +
  geom_point(data=outbreaks,aes(x=Longitude,y=Latitude),color="#D55E00") #+

dos<-plot_grid(uno,mapas,rel_heights = c(1,1.3),ncol = 1,labels = c("","c"),label_size = 30)

ggsave(dos,file="Figure2.png",width = 16,height = 11)

### Correlation network - outbreaks
freqOfOutbreaksAll<-outbreaks %>% select(-serotype,-Latitude,-Longitude) %>%
  group_by(epiunitID) %>% 
  dplyr::summarise(Outbreaks=n()) %>% arrange(desc(Outbreaks)) %>% 
  subset(.,Outbreaks>=5)

freqOfOutbreaksA<-outbreaks %>% 
  subset(serotype=="A") %>%
  select(-serotype,-Latitude,-Longitude) %>%
  group_by(epiunitID) %>% 
  dplyr::summarise(Outbreaks=n()) %>% arrange(desc(Outbreaks)) %>% 
  subset(.,Outbreaks>=3)

freqOfOutbreaksO<-outbreaks %>% 
  subset(serotype=="O") %>%
  select(-serotype,-Latitude,-Longitude) %>%
  group_by(epiunitID) %>% 
  dplyr::summarise(Outbreaks=n()) %>% arrange(desc(Outbreaks)) %>% 
  subset(.,Outbreaks>=3)

freqOfOutbreaks1<-outbreaks %>% 
  subset(serotype=="1") %>%
  select(-serotype,-Latitude,-Longitude) %>%
  group_by(epiunitID) %>% 
  dplyr::summarise(Outbreaks=n()) %>% arrange(desc(Outbreaks)) %>% 
  subset(.,Outbreaks>=2)
#Network statistics for each epiunit.
#I can be created using igraph and correlate with outbreak file.
#I uploaded this file in github.
statsOfNodesAllStatsAndStrain<-read.csv("NodeStatsLocationAndStrainLabel.csv",header = T)
head(statsOfNodesAllStatsAndStrain)

mostTimesInfected<-
  merge(freqOfOutbreaksAll,subset(statsOfNodesAllStatsAndStrain,Strain=="All"),
        by.x="epiunitID",by.y="Node")

mostTimesInfectedA<-
  merge(freqOfOutbreaksA,subset(statsOfNodesAllStatsAndStrain,Strain=="Outbreak A"),
        by.x="epiunitID",by.y="Node") 

mostTimesInfectedO<-
  merge(freqOfOutbreaksO,subset(statsOfNodesAllStatsAndStrain,Strain=="Outbreak O"),
        by.x="epiunitID",by.y="Node") 

mostTimesInfectedAsia1<-
  merge(freqOfOutbreaks1,subset(statsOfNodesAllStatsAndStrain,Strain=="Outbreak 1"),
        by.x="epiunitID",by.y="Node")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

total1<-ggplot(mostTimesInfected,aes(x=OutDegree/mean(statsOfNodesAllStatsAndStrain$OutDegree),
                                    y=OutCoreness/mean(statsOfNodesAllStatsAndStrain$OutCoreness),
                                    size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.5,20))+ labs(size="Outbreaks")+
  theme(legend.position = c(0.8,0.3),text=element_text(size=14))+
  ylim(0.5,2.5)+
  xlab("Out-degree")+ylab("Out-coreness")

total2<-ggplot(mostTimesInfected,aes(x=EigenNormal/mean(statsOfNodesAllStatsAndStrain$EigenNormal),
                                    y=RelativeBet/mean(statsOfNodesAllStatsAndStrain$RelativeBet),
                                    size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.1,100))+ scale_y_log10(limits=c(0.5,130))+labs(size="Outbreaks")+
  theme(legend.position = "none",text=element_text(size=14))+
  xlab("Eigenvector cent.")+ylab("Betweenness cent.")

title <- ggdraw() + draw_label("All serotypes (a)", fontface='bold')

total3<-plot_grid(total1,total2)

total3.1<-plot_grid(title, total3, ncol=1, rel_heights=c(0.1, 1))

A1<-ggplot(mostTimesInfectedA,aes(x=OutDegree/mean(statsOfNodesAllStatsAndStrain$OutDegree),
                                 y=OutCoreness/mean(statsOfNodesAllStatsAndStrain$OutCoreness),
                                 size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.5,20))+ labs(size="Outbreaks")+
  theme(legend.position = c(0.8,0.3),text=element_text(size=14))+
  ylim(0.5,2.5)+
  xlab("Out-degree")+ylab("Out-coreness")

A2<-ggplot(mostTimesInfectedA,aes(x=EigenNormal/mean(statsOfNodesAllStatsAndStrain$EigenNormal),
                                 y=RelativeBet/mean(statsOfNodesAllStatsAndStrain$RelativeBet),
                                 size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.1,100))+ scale_y_log10(limits=c(0.5,130))+labs(size="Outbreaks")+
  theme(legend.position = "none",text=element_text(size=14))+
  xlab("Eigenvector cent.")+ylab("Betweenness cent.")

titleA <- ggdraw() + draw_label("Serotype A (b)", fontface='bold')

A3<-plot_grid(total1,total2)

A3.1<-plot_grid(titleA, A3, ncol=1, rel_heights=c(0.1, 1))

O1<-ggplot(mostTimesInfectedO,aes(x=OutDegree/mean(statsOfNodesAllStatsAndStrain$OutDegree),
                                 y=OutCoreness/mean(statsOfNodesAllStatsAndStrain$OutCoreness),
                                 size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.5,20))+ labs(size="Outbreaks")+
  theme(legend.position = c(0.8,0.3),text=element_text(size=14))+
  ylim(0.5,2.5)+
  xlab("Out-degree")+ylab("Out-coreness")

O2<-ggplot(mostTimesInfectedO,aes(x=EigenNormal/mean(statsOfNodesAllStatsAndStrain$EigenNormal),
                                 y=RelativeBet/mean(statsOfNodesAllStatsAndStrain$RelativeBet),
                                 size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.1,100))+ scale_y_log10(limits=c(0.5,130))+labs(size="Outbreaks")+
  theme(legend.position = "none",text=element_text(size=14))+
  xlab("Eigenvector cent.")+ylab("Betweenness cent.")

titleO <- ggdraw() + draw_label("Serotype O (c)", fontface='bold')

O3<-plot_grid(O1,O2)

O3.1<-plot_grid(titleO, O3, ncol=1, rel_heights=c(0.1, 1))

Asia1<-ggplot(mostTimesInfectedAsia1,aes(x=OutDegree/mean(statsOfNodesAllStatsAndStrain$OutDegree),
                                       y=OutCoreness/mean(statsOfNodesAllStatsAndStrain$OutCoreness),
                                       size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.5,20))+ labs(size="Outbreaks")+
  theme(legend.position = c(0.8,0.3),text=element_text(size=14))+
  ylim(0.5,2.5)+
  xlab("Out-degree")+ylab("Out-coreness")

Asia2<-ggplot(mostTimesInfectedAsia1,aes(x=EigenNormal/mean(statsOfNodesAllStatsAndStrain$EigenNormal),
                                       y=RelativeBet/mean(statsOfNodesAllStatsAndStrain$RelativeBet),
                                       size=factor(Outbreaks)))+
  theme_bw()+  geom_point(color="#0072B2")+
  geom_hline(yintercept = 1,color="#D55E00",size=1)+
  geom_vline(xintercept = 1,color="#D55E00",size=1)+
  scale_x_log10(limits=c(0.1,100))+ scale_y_log10(limits=c(0.5,130))+labs(size="Outbreaks")+
  theme(legend.position = "none",text=element_text(size=14))+
  xlab("Eigenvector cent.")+ylab("Betweenness cent.")

titleAsia <- ggdraw() + draw_label("Serotype Asia-1 (d)", fontface='bold')

Asia3<-plot_grid(Asia1,Asia2)

Asia3.1<-plot_grid(titleAsia, Asia3, ncol=1, rel_heights=c(0.1, 1))

plot_grid(total3.1,A3.1,O3.1,Asia3.1,ncol=2)

ggsave(last_plot(),file="Figure3.png",width = 18,height = 7)
