# #################################################
# R script to generate sc expression data figures in Bergheim et al, 2024 BioRxiv 
# Author: Alison G. Cole
# #################################################

# load libraries 
library(Seurat,lib='Seurat4')
library(easypackages)
libraries("Matrix", "readxl","RColorBrewer",
          'patchwork','dplyr','viridisLite','ggplot2','pals')

{
  ## load genes 
  load(file='Genes.Nv2.RData') #gene annotations and colour palettes
  load(file = 'data1.subsets.publish.Robj') # from Cole et al 2024
  
  ### general sub.cluster names ----
  for (j in 1:length(names (data1.subsets)))
  {
    data1=data1.subsets[[j]]
    clName=NULL
    for (cl in 1:length(levels(data1)))
      clName[cl] = paste(names(data1.subsets)[j],levels(data1)[cl],sep = '.')
    levels(data1@active.ident) = clName
    data1.subsets[[j]]=data1
  }
  
  ## Add labels to Alldata ----
  load(file='Alldata.Nv2.publish.Robj') # from Cole et al 2024
  {
    cl.order=NULL
    Alldata.Nv2$ID.separate = as.character(Alldata.Nv2$IDs)
    for (i in 1:length(names (data1.subsets)))
    {
      coi=NULL
      coi=colnames(data1.subsets[[i]])
      Alldata.Nv2$ID.separate[coi] = as.character(data1.subsets[[i]]@active.ident[coi])
      cl.order = c(cl.order,levels(data1.subsets[[i]]))
    }
    cl.ind = match(unique(cl.order),levels(as.factor(Alldata.Nv2$ID.separate)))
    Alldata.Nv2$ID.separate = as.factor(Alldata.Nv2$ID.separate)
    Alldata.Nv2$ID.separate = factor(Alldata.Nv2$ID.separate,levels(Alldata.Nv2$ID.separate)[cl.ind])
    DimPlot(Alldata.Nv2,group.by = 'ID.separate',cols=c(clust.cp.graded,clust.cp.separate))&NoLegend()
    levels(Alldata.Nv2$ID.separate)
  }
  
  clust.cp.all.IDs = NULL
  for(i in 1:length(data1.subsets))
    clust.cp.all.IDs = c(clust.cp.all.IDs,rep(clust.cp.separate[i],length(levels(data1.subsets[[i]]))))
  names(clust.cp.all.IDs)=cl.order
  clust.cp.all.IDs=clust.cp.all.IDs[unique(names(clust.cp.all.IDs))]
}
## gene lists ----
Matrisome=readxl::read_excel('ECM.TableS1.xlsx')
Matrisome=Matrisome[,1:6]
Matrisome=Matrisome[rowSums(is.na(Matrisome)) == 0,]

toxins=readxl::read_excel('toxins.xlsx')
toxins=toxins[-19,]

#update the gene names:
genes$gene_short_name[match(Matrisome$gene_name,genes$geneID)]=Matrisome$short_Name
genes$name=paste(genes$gene_short_name,genes$geneID,sep=' : ')

goi.all=genes$name[match(Matrisome$gene_name,genes$geneID)]

goi.ECMa=genes$name[match(Matrisome$gene_name[Matrisome$matrisome_division=="Matrisome associated"],genes$geneID)]
goi.other=genes$name[match(Matrisome$gene_name[Matrisome$matrisome_division=="Other"],genes$geneID)]
goi.core=genes$name[match(Matrisome$gene_name[Matrisome$matrisome_division=="Core matrisome"],genes$geneID)]
goi.collagen=genes$name[match(Matrisome$gene_name[Matrisome$matrisome_category=="Collagens"],genes$geneID)]
goi.polydome=genes$name[match(Matrisome$gene_name[grep('Polydom',Matrisome$short_Name)],genes$geneID)]

goi.all.sorted=c(sort(goi.core),sort(goi.ECMa),sort(goi.other))

#Prep data ----
goi=goi.all
temp=Alldata.Nv2
#try first update the gene names according to current annotations
rownames(temp@assays$RNA@counts)=genes$name
rownames(temp@assays$RNA@data)=genes$name
rownames(temp@assays$RNA@meta.features)=genes$name # otherwise fails with vs5

temp<-ScaleData(temp,features = goi,split.by = 'orig.ident')

#split the dataset into the three partitions:
levels(temp$orig.ident)
temp<-SetIdent(temp,value = 'orig.ident')
coi1=WhichCells(temp,idents = levels(temp)[1:8])
coi2=WhichCells(temp,idents = levels(temp)[8:13])
temp$lifehistory[coi1]='larvae'
temp$lifehistory[coi2]='p.polyp'
temp<-SetIdent(temp,value = 'ID.separate')

goi=genes$name[match(Matrisome$gene_name,genes$geneID)]
goi=unique(goi)

# generate working dataset:
ECM.data=subset(temp,features = goi)
ECM.data=DietSeurat(ECM.data,scale.data = T,dimreducs = 'umap')
data1=SetIdent(ECM.data,value = 'ID.separate')
# generate a named Matrisome colour pal
MAT.cp=c(cols25(25)[14:25],'black')
names(MAT.cp)=unique(Matrisome$matrisome_category)
pal.bands(MAT.cp,clust.cp.all.IDs)

# combine cnidocytes, putative stem/germ cells, and epithelial layers
ECM.data=SetIdent(ECM.data,value = 'IDs')
ECM.data$reduced=ECM.data@active.ident
levels(ECM.data$reduced)=c("pSC|PGCs", "pSC|PGCs", "neuronal", "gland.mucous", "unchar.immune", "gland.S1/S2/dig", "cnidocyte", "cnidocyte", "ectoderm", 
                           "ectoderm", "ectoderm", "retractor muscle", "mesendoderm", 
                           "mesendoderm")
ECM.data$reduced=droplevels(ECM.data$reduced)
ECM.data$reduced = factor(ECM.data$reduced,
                          levels(ECM.data$reduced)[c(1,2,3,5,4,6,9,8,7)])
# re-arange the cnidocytes to put all mature states together:
ECM.data$ID.separate = factor(ECM.data$ID.separate,levels(ECM.data$ID.separate)[c(1:58,59,60,64,61,65:67,62,68,69,70,63,71,72:110)])

# Figures.All ----

### FigS4 ----
  data1=SetIdent(ECM.data,value = 'ID.separate')
    levels(data1)
    cl.order=c(102:110,1:101)
    data1@active.ident=factor(data1@active.ident,levels(data1@active.ident)[cl.order])
    
    goi= goi.all
      FigS4=DotPlot(data1,features = goi,split.by = 'lifehistory',cols=c('darkred','darkorange','steelblue4'))
    write.csv(FigS4$data,file='ECM_matrisome.collab/pub/MatrisomeExpression.raw.cvs')

##Fig2B | Module scores ----
    {
      goi=list(goi.core,goi.ECMa,goi.other)
      names(goi)=c('CoreMatrisome','ECMassociated','Other')
      temp <- AddModuleScore(
        object = temp,
        features = goi,
        name = c('CoreMatrisome','ECMassociated','Other'),
      )
      Fig2b=FeaturePlot(temp,c('CoreMatrisome1','ECMassociated2','Other3'),order=T,cols=gene.cp,ncol=3,raster=T)&NoAxes()&NoLegend()
      Fig2b
      source.2b= as.data.frame(temp@meta.data[,10:12])
      
    }    
  
##Fig2C | DEG ----
      goi=unique(goi.core)
      all.markers <- FindAllMarkers(
        data1,
        features = goi,
        return.thresh = 0.001,
        min.pct = 0.2,
        only.pos = TRUE)
      
      write.csv(all.markers,file='DEG.celltypes.CORE.csv')
      
      #generate a collated list of unique DE genes
      list = NULL
      for (i in 1:length(unique(all.markers$cluster)))
      {
        x = all.markers[as.numeric(all.markers$cluster) == i, ][1:min(5, length(which(as.numeric(all.markers$cluster) == i))), 7]
        
        list = c(list, x)
      }
      list=unique(list)
      list=list[2:length(list)] #drop 'NA' from clusters with no DEGs

Fig2B<-DotPlot(SetIdent(ECM.data,value='reduced'),features = unique(c(list)),cols=c('darkred','darkorange','steelblue4'),split.by = 'lifehistory',dot.min = 0.01,col.min = 0, scale.by = 'size')& theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6),axis.text.y = element_text(colour ='black'))&labs(tag='B',title = 'Differentially Expressed Genes | CORE Matrisome',subtitle = 'Tissue Types, Top 5 up-regulated genes with p.val <=0.001, in at least 20% of any cell state ') &RotatedAxis()&FontSize(9,7)&coord_flip()

#generate the colour scales:
L_R=FeaturePlot(ECM.data,goi[1],cols=c('lightgrey','darkred'))
L_O=FeaturePlot(ECM.data,goi[1],cols=c('lightgrey','darkorange'))
L_B=FeaturePlot(ECM.data,goi[1],cols=c('lightgrey','steelblue4'))
legR <- ggpubr::get_legend(L_R+theme(legend.position = 'top',legend.text.align = 1))
legO <- ggpubr::get_legend(L_O+theme(legend.position = 'top',legend.text.align = 1))
legB <- ggpubr::get_legend(L_B+theme(legend.position = 'top',legend.text.align = 1))
# Convert to a ggplot and print
leg1=ggpubr::as_ggplot(legR)
leg2=ggpubr::as_ggplot(legO)
leg3=ggpubr::as_ggplot(legB)
colourscales=leg1+leg2+leg3

pdf('Fig2B.raw.pdf')    
Fig2B+colourscales+plot_layout(ncol = 1)
dev.off()
    ###FigS5 ----
FigS5=DotPlot(data1,features = rev(goi.core),cols=c('darkred','darkorange','steelblue4'),split.by = 'lifehistory',dot.min = 0.01,col.min = 0, scale.by = 'size')& theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6),axis.text.y = element_text(colour ='black'))&labs(title = 'CORE Matrisome',subtitle = 'Alldata | cell states') &RotatedAxis()&FontSize(9,7)
    pdf('FigS5.raw.pdf',width = 54,height=36)
    FigS5
    dev.off()

##FigS6 Collagens ----

    FigS6=DotPlot(data1,features = sort(goi.collagen),cols=c('darkred','darkorange','steelblue4'),split.by = 'lifehistory',dot.min = 0.01,col.min = 0, scale.by = 'size')& theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6),axis.text.y = element_text(colour ='black'))&labs(title = 'Collagens',subtitle = 'Alldata | cell states') &RotatedAxis()&FontSize(9,7)
    pdf('FigS6.raw.pdf',width = 12,height=36)
    FigS6
    dev.off()  
    
##Fig 3 Cnido----
    #set up gene lists:
{

     goi=goi.all
    #use all genes expressed in any cnidocyte cluster 
    expression.subsets=DotPlot(data1,'RNA',goi,group.by = 'ID.separate')
    expression.subsets.data=expression.subsets$data
    write.csv(expression.subsets.data,file='ECMGeneExpressionAllCellTypes.csv')
    
    # subset the table according to cnidocytes...
    expression.cnido <- expression.subsets.data[expression.subsets.data$id %in% (levels(expression.subsets.data$id)[59:71]),]
    # keep the rest of the table
    expression.NOT.cnido <- expression.subsets.data[!expression.subsets.data$id %in% (levels(expression.subsets.data$id)[59:71]),]
    
    # filter for at least 5% of the cells in any cnidocyte cluster
    cnido.ecm=expression.cnido[(expression.cnido$pct.exp > 5 & expression.cnido$avg.exp.scaled >0),]
    write.csv(cnido.ecm,file='AllCnidocyteECMGenes.csv')
    all.cnido.genes=unique(as.character(cnido.ecm$features.plot))
    expression=DotPlot(data1,'RNA',all.cnido.genes,group.by = 'IDs')
    
    # keep the rest of the table as coarse clustering
    expression.NOT.cnido <- expression$data[!expression$data$id %in% (levels(expression$data$id)[7:8]),]
    # generate a specificity index: 
    specificity.index=  expression.NOT.cnido %>% group_by(features.plot) %>% count(avg.exp.scaled <=0) %>% arrange(desc(n))  #for every gene count how many clusters have less than average expression
    
    exclusive=as.character(specificity.index$features.plot[specificity.index$n>=12 & specificity.index$`avg.exp.scaled <= 0` == T])

    exclusive=sort(exclusive)
    plot=DotPlot(ECM.data,'RNA',features = exclusive,group.by = 'IDs',  cols=c('grey90','darkred'),dot.min = 0.01,col.min = 0, scale.by = 'size')& theme(panel.grid = element_line(size=0.002,colour = 'grey'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&labs(title = 'Cnidocyte Expressed Genes | exclusive',subtitle = 'reads in >10 percent of any cnidocyte population') &RotatedAxis()&FontSize(9,7)&coord_flip()
    plot
    x.ex=plot$data[plot$data$avg.exp.scaled >= 2 & plot$data$id=='cnidocyte',]
    x2.ex=plot$data[plot$data$avg.exp.scaled >= 2 & plot$data$id=='mat.cnido',]
    specification=as.character(x.ex$features.plot)
    maturation=as.character(x2.ex$features.plot)
    exclusive=(unique(c(maturation,specification)))
    
    # ubiquity index:
    ubiquity.index=  expression.subsets.data[expression.subsets.data$features.plot %in% setdiff(all.cnido.genes,exclusive),] %>% group_by(features.plot) %>%filter(pct.exp >=5) %>% count(pct.exp >=5)%>%arrange(desc(n))
    ubiquitous=as.character(ubiquity.index$features.plot[ubiquity.index$n>54])
    ubiquitous=sort(ubiquitous)
    
    #Split these by cnido vs mature:
    plot=DotPlot(ECM.data,'RNA',features = ubiquitous,group.by = 'IDs',  cols=c('grey90','darkred'),dot.min = 0.01,col.min = 0, scale.by = 'size')& theme(panel.grid = element_line(size=0.002,colour = 'grey'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&labs(title = 'Cnidocyte Expressed Genes | exclusive',subtitle = 'reads in >10 percent of any cnidocyte population') &RotatedAxis()&FontSize(9,7)&coord_flip()
    x=plot$data[plot$data$avg.exp.scaled >=1 & plot$data$id=='cnidocyte',]
    x2=plot$data[plot$data$avg.exp.scaled >= 1 & plot$data$id=='mat.cnido',]
    x3=plot$data[plot$data$avg.exp.scaled <= 0 & plot$data$id=='mat.cnido' | plot$data$id=='cnidocyte',]
    ubiquitous=(as.character(unique(c(x$features.plot,x2$features.plot,x3$features.plot))))
    
    # consider the rest as shared:
    
    shared.1=setdiff(all.cnido.genes,c(exclusive,ubiquitous))
    
    #Split these by cnido vs mature:
    plot=DotPlot(ECM.data,features = shared.1,group.by = 'IDs',  cols=c('grey90','darkred'),
                 col.min = 0, scale.by = 'size')& theme(panel.grid = element_line(size=0.002,colour = 'grey'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&labs(title = 'Cnidocyte Expressed Genes | exclusive',subtitle = 'reads in >10 percent of any cnidocyte population') &RotatedAxis()&FontSize(9,7)&coord_flip()
    x=plot$data[plot$data$avg.exp.scaled >= 1 & plot$data$id=='cnidocyte',]
    x2=plot$data[plot$data$avg.exp.scaled >= 1 & plot$data$id=='mat.cnido',]
    
    shared=as.character(unique(c(sort(x$features.plot),sort(x2$features.plot))))
    shared=union(shared.1,shared)
    cnido.genes.sorted=unique(c(exclusive,shared,ubiquitous))

}
    
    data1=SetIdent(ECM.data,value = 'ID.separate')
    levels(data1)
    cl.order=c(102:110,1:101)
    data1@active.ident=factor(data1@active.ident,levels(data1@active.ident)[cl.order])
 ### Fig. S7 ----   
    FigS7=DotPlot(data1,features = cnido.genes.sorted,  cols=c('grey90','darkred'),dot.min = 0.05,col.min = 0, scale.by = 'size')& theme(panel.grid = element_line(size=0.002,colour = 'grey'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&labs(title = 'Cnidocyte Expressed Genes | Sorted',subtitle = 'reads in >5 percent of any cnidocyte population') &RotatedAxis()&FontSize(9,7)&coord_flip()
    pdf('FigS7.raw.pdf',width = 22,height = 26)
    FigS7
    dev.off()
###Fig3A ----
    pie(c(dist,532),col = clust.cp.separate)
dist=matrix(c(41,156,23,78),4)#filter at 10
    colnames(dist)='gene distribution'
rownames(dist)=c('ubiquitous','shared','mature-specific','specification-specific')
data_percentage <- apply((dist), 2, function(x){x*100/sum(x,na.rm=T)})
row.names(data_percentage)=c('ubiquitous','shared','maturation-specific','specification-specific')
pdf('Fig3A.raw.pdf',width = 3,height = 8)
barplot(data_percentage,col=clust.cp.separate,legend.text = row.names(data_percentage),args.legend = list(x = "topright"),width = 3,plot = T)
dev.off()

S7.expression.table=AverageExpression(data1,features = cnido.genes.sorted)

write.csv(S7.expression.table$RNA,file='AllCnidocyteECMGenesSorted.csv')

###Fig3B ----
goi=list(ubiquitous,shared,maturation,specification)
names(goi)=c('ubiquitous','shared','mature.specific','specification.specific')
ECM.data <- AddModuleScore(
  object = ECM.data,
  features = goi,ctrl=25,
  name = c('ubiquitous','shared','mature.specific','specification.specific')
)
pdf(file='Fig3B.raw.pdf')
FeaturePlot(ECM.data,rev(c('shared2','mature.specific3','specification.specific4','ubiquitous1')),order=T,cols=gene.cp,ncol=2,raster=T,min.cutoff = 0,pt.size = 2,keep.scale = 'all')&NoAxes()&NoLegend()
dev.off()

###Fig3C ----
cnido=subset(data1,idents=levels(data1)[68:80])
cnido@active.ident=factor(cnido@active.ident,levels(cnido)[c(1,2,4,8,5,12,3,6,7,9,11,10,13)])  

all.markers.c <- FindAllMarkers(
  cnido,
  logfc.threshold = 1,
  # features = goi.Matrisome,
  return.thresh = 0.001,
  min.pct = 0.2,
  only.pos = TRUE,
)
#generate a collated list of unique DE genes
list.c = NULL
for (i in 1:length(levels(cnido@active.ident)))
{
  x = all.markers.c[as.numeric(all.markers.c$cluster) == i, ][1:min(5, length(which(as.numeric(all.markers.c$cluster) == i))), 7]
  list.c = c(list.c, x)
}
list.c=unique(list.c)

pdf (file='Fig3C.raw.pdf',width = 12)
DotPlot(cnido,features = list.c,
        cols=c('grey90','darkred'),dot.min = 0.01,col.min = 0, scale.by = 'size')+ theme(panel.grid = element_line(size=0.002,colour = 'grey'),legend.title = element_text(size = 8),legend.text = element_text(size=6),legend.position = 'bottom')&labs(title = 'Differentially Expressed Genes',subtitle = 'Cnidocyte cell-states') &RotatedAxis()&FontSize(9,12)#
dev.off()

  ### Fig.S10 ----
pdf(file='polydom.raw.pdf',height = 12)
  DotPlot(data1,'RNA',sort(goi.polydome[1:8]),cols=c('grey90','darkred','darkorange','steelblue4'),dot.min = 0.01,col.min = 0, scale.by = 'size')& RotatedAxis()&FontSize(9,7) &theme(legend.title = element_text(size = 8),legend.text = element_text(size=6),legend.position = 'bottom')&labs(title = 'Polydom',subtitle = 'AllData: Cell types')
dev.off()
