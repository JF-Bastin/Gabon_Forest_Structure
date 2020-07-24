###########################
### Create a global tree db
#
# columns : Country;ID-plot;Tree-Tag;Species;Genus;Family;DBH;measured_H;E;WSG_dryad;Pheno;FCycle
# script constructs db for each country, rbind in the end
# WSG is extracted from DRYAD (Chave,Zanne 2009)
# Functional traits are extracted from CIRAD Cofortrait database (Benedet et al. 2015)
# All the floristic data (db, dryad, cofortrait) are standardized using the TNRS service


rm(list = ls())

###########################
# Build DRYAD

work_dir <- 'C:/Users/FAO_bast/Desktop/Gabon_Bastin_script_data/dryad'
setwd(work_dir)
db_DRYAD=read.csv("DRYAD.csv", sep=";")
head(db_DRYAD)
dico_dryad=read.csv("tnrs_results_Dryad.csv", sep=";")
db_corr=merge(db_DRYAD,dico_dryad,by.x=c("Binomial"),by.y=c("Name_submitted"),all.x=TRUE)
head(db_corr)
colnames(dico_dryad)

db_corr=db_corr[db_corr$Region=="Africa (tropical)",]
db_corr2_sp=tapply(db_corr$Wood.density..g.cm.3...oven.dry.mass.fresh.volume,db_corr$Name_matched,function(x)mean(x))
db_corr2_gen=tapply(db_corr$Wood.density..g.cm.3...oven.dry.mass.fresh.volume,db_corr$Genus_matched,function(x)mean(x))
db_corr2_fam=tapply(db_corr$Wood.density..g.cm.3...oven.dry.mass.fresh.volume,db_corr$Name_matched_accepted_family,function(x)mean(x))

###########################
# Data GABON
# load required data
# DBH and H data are separated in Gabon, requires to be merge
# Merge H data to global data
# Associate TNRS species/genus/family to each tree
# Merge Dryad WSG to each species/genus/family

setwd("C:/Users/FAO_bast/Desktop/Gabon_Bastin_script_data/gabon");
db_tree_gabon=read.csv("Gabon_tree_dbh_05082017.csv",sep=";");
db_treeH_gabon=read.csv("Gabon_tree_Hprinc_14012017.csv",sep=";");
db_treeHaux_gabon=read.csv("Gabon_tree_Haux_14012017.csv",sep=";");
db_gabon_type=read.csv("Gabon_plot_meta_14012017b.csv",sep=";")

dico_tnrs_sp=read.csv("dico_tnrs_sp.csv", sep=";")
dico_tnrs_gen=read.csv("dico_tnrs_gen.csv", sep=";")
dico_tnrs_fam=read.csv("dico_tnrs_fam.csv", sep=";")

Liane=tapply(db_tree_gabon$Espece%in%c("Liane Indet","LIANE Indet"),db_tree_gabon$Code,function(x)sum(x,na.rm=TRUE))
Palm=tapply(db_tree_gabon$Espece%in%c("Elais guineensis","Raphia Indet","Coffea indet"),db_tree_gabon$Code,function(x)sum(x,na.rm=TRUE))

##############################
### NEED TO EXCLUDE LIANNA ###
##############################
db_tree_gabon=db_tree_gabon[!(db_tree_gabon$Espece=="Liane Indet"|db_tree_gabon$Espece=="LIANE Indet"|
                              db_tree_gabon$Espece=="Elais guineensis"|db_tree_gabon$Espece=="Raphia Indet"),]

#####################################################
### NEED TO EXCLUDE SAVANNA and MAJOR DISTURBANCE ###
#####################################################
list_type=db_gabon_type[!db_gabon_type$VA%in%"Savane"&!db_gabon_type$PA%in%"Majeur","Code"]


#####################################################
### Calculate plot level tree height              ###
#####################################################

# Merge H data to global data
# Note: cannot merge with aux as Tag are missing for Auxilliary plots
for (n in c(1:length(db_treeH_gabon[,1]))){
  db_treeH_gabon[n,"Hmed"]=median(as.numeric(db_treeH_gabon[n,c("H1","H2","H3")]), na.rm=TRUE)
  }

db_H=merge(db_tree_gabon,db_treeH_gabon[,c("Code","Tag","Hmed")],by.x=c("Code","Tag"),by.y=c("Code","Tag"),all.x=TRUE)
summary(db_H) # 97063 NA's equals the number of trees without ground H measurment
summary(tapply(db_H$Hmed,db_H$Code,function(x)mean(x,na.rm=TRUE))) # 27 plots without H measurement

# Associate TNRS species/genus/family to each tree
db_TNRS=merge(db_H,dico_tnrs_sp,by.x=c("Espece"),by.y=c("sp_Name_submitted"),all.x=TRUE)
db_TNRS=merge(db_TNRS,dico_tnrs_gen,by.x=c("Genre"),by.y=c("Gen_Name_submitted"),all.x=TRUE)
db_TNRS=merge(db_TNRS,dico_tnrs_fam,by.x=c("Famille"),by.y=c("Fam_sub"),all.x=TRUE)

# Merge Dryad WSG to each species/genus/family
db_d=0
db_d=merge(db_TNRS,db_corr2_sp,by.x=c("sp_Name_matched"),by.y="row.names",all.x=TRUE)
db_d=merge(db_d,db_corr2_gen,by.x=c("Gen_Name_matched"),by.y="row.names",all.x=TRUE)
db_d=merge(db_d,db_corr2_fam,by.x=c("Fan_matched"),by.y="row.names",all.x=TRUE)

# Calculate WSG_new with hierachy sp/gen/fam
# Then attribute to the rest the average of the plot? (not necessary if we want the average per plot)
# except if we want the WSGweighted BA we need for all the trees

new_db=db_d[,c("Code","Tag","PlacetteAux","sp_Name_matched","Gen_Name_matched","Fan_matched","D","Hmed","y.x","y.y","y")]
colnames(new_db)<-c("Code","Tag","PlacetteAux","species","genus","family","D","Hmed","WD_sp","WD_gen","WD_fam")

#remove savannah plot and degraded plot
new_db<-new_db[new_db$Code%in%list_type,]

# Associated WSG if WSGsp is missing
# 1. Gen
# 2. Fam
# 3. Plot average

new_db$WSG<-new_db$WD_sp
new_db[is.na(new_db$WD_sp),"WSG"]<-new_db$WD_gen[is.na(new_db$WD_sp)]
new_db[is.na(new_db$WSG),"WSG"]<-new_db$WD_fam[is.na(new_db$WSG)]

WSG_plot=tapply(new_db$WSG,new_db$Code,function(x)mean(x,na.rm=TRUE))
list_plot=unique(new_db$Code)

for (i in c(1:length(list_plot))){
  value<-as.numeric(WSG_plot[names(WSG_plot)==list_plot[i]])
  filter<-is.na(new_db$WSG)&new_db$Code==list_plot[i]
  new_db[filter,"WSG"]<-value
  }

head(new_db)

####################################################################
# Height model Gabon
####################################################################

list_plot=unique(db_treeH_gabon$Code)
ID=matrix(NA,nrow=length(list_plot),1)
a_MM=matrix(NA,nrow=length(list_plot),1)
b_MM=matrix(NA,nrow=length(list_plot),1)
cor_MM=matrix(NA,length(list_plot),1)
err_MM=matrix(NA,length(list_plot),1)

all=paste(db_tree_gabon[,c("Tag")],db_tree_gabon[,c("Code")])
selected=paste(db_treeH_gabon[,c("Tag")],db_treeH_gabon[,c("Code")])
db_gabon_sub=db_tree_gabon[all%in%selected,]

db_gabon_sub2=merge(db_gabon_sub,db_treeH_gabon[,c("Hmed","Code","Tag")],c("Code","Tag"),c("Code","Tag"),all.y=TRUE)
summary(db_gabon_sub2)
list_plot=unique(db_gabon_sub2$Code)

### NLS model chosen is a Michaelis-Menten - see Molto et al. 2014
for(i in c(1:length(list_plot))){
  sub=db_gabon_sub2[db_gabon_sub2$Code==list_plot[i],]
  
  x_temp=as.numeric(as.character(sub$D))
  x=x_temp[!is.na(sub$Hmed)&!is.na(x_temp)]
  y=sub$Hmed[!is.na(sub$Hmed)&!is.na(x_temp)]
  
    m<-tryCatch(nls(y~a*x/(b+x),start=list(a=30,b=2)),error=function(e) NA )
    if (is.na(m)){
      ID[i]=as.character(list_plot[i])
      a_MM[i]<-NA
      b_MM[i]<-NA
      cor_MM[i]<-NA
      err_MM[i]<-NA
      }
    else{
      ID[i]=as.character(list_plot[i])
      a_MM[i]=coefficients(m)[1]
      b_MM[i]=coefficients(m)[2]
      cor_MM[i]=cor(y,predict(m))
      err_MM[i]=summary(m)$sigma/mean(y,na.rm=TRUE)
      }
    }


### Associate coefficient to each plot in meta table
#db_plot_gabon vient de inventories/gabon/db_tree_coord

db_plot_gabon=read.csv("Gabon_tree_coord_14012017.csv",sep=";")
head(db_plot_gabon)
gabon_princ=db_plot_gabon[db_plot_gabon$Placette=="Principale",]
coord_x=tapply(as.numeric(as.character(gabon_princ$Longitude_final)),gabon_princ$Code,function(x)mean(x))
coord_y=tapply(as.numeric(as.character(gabon_princ$Latitude_final)),gabon_princ$Code,function(x)mean(x))

df=cbind(coord_x,coord_y)
df2=cbind(ID,a_MM,b_MM,cor_MM,err_MM)
colnames(df2)<-c("ID","a_MM","b_MM","cor_MM","err_MM")
me2=merge(df,df2,by.x="row.names",by.y="ID",all=TRUE)

# Some points are missing H, for these, we compute H using the E variable of Chave et al. 2014
setwd("C:/Users/FAO_bast/Desktop/Gabon_Bastin_script_data")
library(BIOMASS)
summary(db_plot_gabon$Longitude)
temp=data.frame(cbind(coord_x,coord_y))
E=computeE(temp)

temp$Row_names=names(coord_x)
temp<-temp[complete.cases(temp),]
temp$E=E

db_H_allom=merge(me2,temp[,c("Row_names","E")],by.x="Row.names",by.y="Row_names",all.x=TRUE)
db_H_allom=db_H_allom[complete.cases(db_H_allom),]

new_db=merge(new_db,db_H_allom,by.x="Code",by.y="Row.names",all.x=TRUE)
# Calculate Hlocal from a and b - y~a*x/(b+x)
new_db$a_MM=as.numeric(as.character(new_db$a_MM))
new_db$b_MM=as.numeric(as.character(new_db$b_MM))
new_db$D=as.numeric(as.character(new_db$D))
# new_db$Hlocal[!is.na(new_db$a_MM)]=new_db$a_MM[!is.na(new_db$a_MM)]*new_db$D[!is.na(new_db$a_MM)]/(new_db$b_MM[!is.na(new_db$a_MM)]+new_db$D[!is.na(new_db$a_MM)])
new_db$Hlocal=new_db$a_MM*new_db$D/(new_db$b_MM+new_db$D)

# calculate HChave from E for each tree (study what does explain the deviation)
ln_HChave=0.893-new_db$E+(0.76*log(new_db$D))-(0.034*(log(new_db$D))^2)
new_db$HChave=exp(ln_HChave)

######################################################
# Keep only principal plots - for height and better standardization of structural values
######################################################

new_db<- new_db[new_db$PlacetteAux%in%"Principale",]

#####################################################
### ASSOCIATING TRAIT DATABASE                    ###
#####################################################

# Add Cofortrait
setwd("C:/Users/FAO_bast/Desktop/Gabon_Bastin_script_data/cofortraits")
db_trait_new=read.csv("db_trait_tnrs.csv")

temp=table(db_trait_new$Genus_matched,db_trait_new$Pheno)
temp=temp[,-1]
temp2=colnames(temp)[apply(temp,1,which.max)]
temp3=cbind(row.names(temp),temp2)
colnames(temp3)<-c("Genus","Pheno_Genus")

temp4=table(db_trait_new$Name_matched_accepted_family,db_trait_new$Pheno)
temp4=temp4[,-1]
temp5=colnames(temp4)[apply(temp4,1,which.max)]
temp6=cbind(row.names(temp4),temp5)
colnames(temp6)<-c("Family","Pheno_Family")

temp7=table(db_trait_new$Genus_matched,db_trait_new$Stol)
temp7=temp7[,-1]
temp8=colnames(temp7)[apply(temp7,1,which.max)]
temp9=cbind(row.names(temp7),temp8)
colnames(temp9)<-c("Genus","Stol_Genus")

temp10=table(db_trait_new$Name_matched_accepted_family,db_trait_new$Pheno)
temp10=temp10[,-1]
temp11=colnames(temp10)[apply(temp10,1,which.max)]
temp12=cbind(row.names(temp10),temp11)
colnames(temp12)<-c("Family","Pheno_Family")

temp13=table(db_trait_new$Genus_matched,db_trait_new$Diaspmod)
temp13=temp13[,-1]
temp14=colnames(temp13)[apply(temp13,1,which.max)]
temp15=cbind(row.names(temp13),temp14)
colnames(temp15)<-c("Genus","Stol_Genus")

temp16=table(db_trait_new$Name_matched_accepted_family,db_trait_new$Diaspmod)
temp16=temp16[,-1]
temp17=colnames(temp16)[apply(temp16,1,which.max)]
temp18=cbind(row.names(temp16),temp17)
colnames(temp18)<-c("Family","Pheno_Family")

# Merging Pheno and Stol for species
db_trait_temp=db_trait_new[,c("Name_matched","Diaspmod","Pheno","Stol")]
colnames(db_trait_temp)<-c("sp","diasp_sp","pheno_sp","stol_sp")

### 
new_db_tem=merge(new_db,db_trait_temp[,c("sp","diasp_sp","pheno_sp","stol_sp")],by.x="species",by.y="sp",all.x=TRUE)

## Associate missing trait 1.gen / 2.fam
colnames(temp3)<-c("genus","pheno_gen")
new_db_tem=merge(new_db_tem,temp3,by.x="genus",by.y="genus",all.x=TRUE)
colnames(temp6)<-c("family","pheno_fam")
new_db_tem=merge(new_db_tem,temp6,by.x="family",by.y="family",all.x=TRUE)
colnames(temp9)<-c("genus","stol_gen")
new_db_tem=merge(new_db_tem,temp9,by.x="genus",by.y="genus",all.x=TRUE)
colnames(temp12)<-c("family","stol_fam")
new_db_tem=merge(new_db_tem,temp12,by.x="family",by.y="family",all.x=TRUE)
colnames(temp15)<-c("genus","diasp_gen")
new_db_tem=merge(new_db_tem,temp15,by.x="genus",by.y="genus",all.x=TRUE)
colnames(temp18)<-c("family","diasp_fam")
new_db_tem=merge(new_db_tem,temp18,by.x="family",by.y="family",all.x=TRUE)

# Complete for genus and family
new_db_tem$pheno<-new_db_tem$pheno_sp
new_db_tem[is.na(new_db_tem$pheno_sp),"pheno"]<-new_db_tem$pheno_gen[is.na(new_db_tem$pheno_sp)]
new_db_tem[is.na(new_db_tem$pheno),"pheno"]<-new_db_tem$pheno_fam[is.na(new_db_tem$pheno)]

new_db_tem$stol<-new_db_tem$stol_sp
new_db_tem[is.na(new_db_tem$stol_sp),"stol"]<-new_db_tem$stol_gen[is.na(new_db_tem$stol_sp)]
new_db_tem[is.na(new_db_tem$stol),"stol"]<-new_db_tem$stol_fam[is.na(new_db_tem$stol)] 

new_db_tem$diasp<-new_db_tem$diasp_sp
new_db_tem[is.na(new_db_tem$diasp_sp),"diasp"]<-new_db_tem$diasp_gen[is.na(new_db_tem$diasp_sp)]
new_db_tem[is.na(new_db_tem$diasp),"diasp"]<-new_db_tem$diasp_fam[is.na(new_db_tem$diasp)]

db_Gabon=new_db_tem
head(db_Gabon)

#############################################################3
#### Aggregate structure parameters @ plot level
#############################################################3

setwd("C:/Users/FAO_bast/Desktop/Gabon_Bastin_script_data/gabon");
db_plot_gabon=read.csv("Gabon_tree_coord_14012017.csv",sep=";")
head(db_plot_gabon)
gabon_princ=db_plot_gabon[db_plot_gabon$Placette=="Principale",]
coord_x=tapply(as.numeric(as.character(gabon_princ$Longitude_final)),gabon_princ$Code,function(x)mean(x,na.rm=TRUE))
coord_y=tapply(as.numeric(as.character(gabon_princ$Latitude_final)),gabon_princ$Code,function(x)mean(x,na.rm=TRUE))
db_meta=cbind(coord_x,coord_y)
db_meta=data.frame(db_meta)

# Per tree
db_Gabon$Code<-droplevels.factor(db_Gabon$Code)
db_Gabon$BA=pi*(db_Gabon$D/2)^2
db_Gabon$WD_BA=db_Gabon$WSG*(db_Gabon$BA/10000)
db_Gabon$AGB=0.0559*((db_Gabon$WSG)*(db_Gabon$Hlocal)*(db_Gabon$D)^2)
db_Gabon$AGB_Chave=0.0559*((db_Gabon$WSG)*(db_Gabon$HChave)*(db_Gabon$D)^2)
# Where Hlocal is null, use HChave
db_Gabon$AGB[is.na(db_Gabon$AGB)]<-db_Gabon$AGB_Chave[is.na(db_Gabon$AGB)]

# Per plot
list_plot=rownames(db_meta)%in%db_Gabon$Code
db_meta=db_meta[list_plot,]
list_plot2=!unique(db_Gabon$Code)%in%rownames(db_meta)
ID=unique(db_Gabon$Code)[list_plot2]
db_Gabon=db_Gabon[!db_Gabon$Code%in%ID,]
db_Gabon$Code<-droplevels.factor(db_Gabon$Code)

db_meta$N<-tapply(db_Gabon[,"WSG"],db_Gabon[,"Code"],function(x)length(x))
db_meta$BA=tapply(db_Gabon$D,db_Gabon$Code,function(x)sum(pi*(x/2)^2, na.rm=TRUE))/10000
db_meta$Dg=tapply(db_Gabon$D,db_Gabon$Code,function(x)sqrt(sum(x^2, na.rm=TRUE)/length(x)))
db_meta$AGB=tapply(db_Gabon$AGB,db_Gabon$Code,function(x)sum(x, na.rm=TRUE))/1000
db_meta$AGB[db_meta$AGB==0]<-NA

db_meta$H95=tapply(db_Gabon$Hlocal,db_Gabon$Code,function(x)quantile(x,probs=0.95, na.rm=TRUE))
db_meta$Hmed_85=tapply(db_Gabon$Hmed,db_Gabon$Code,function(x)quantile(x,probs=0.85, na.rm=TRUE))

# For height, we have less observation with the model (Hlocal) vs. the actual measure of H (Hmed)
# We check the correlation between the model and the measure
summary(db_meta$H95)
summary(db_meta$Hmed_85)
plot(db_meta$H95,db_meta$Hmed_85)
# The two are well correlated so we keep the one with the most observations, i.e. Hmed

db_meta$WD=tapply(db_Gabon$WSG,db_Gabon$Code,function(x)mean(x, na.rm=TRUE))
db_meta$WD_BA=tapply(db_Gabon$WD_BA,db_Gabon$Code,function(x)sum(x, na.rm=TRUE))/db_meta$BA

### proportion of trait
db_meta$prop_decid=tapply(db_Gabon$pheno[db_Gabon$pheno%in%"décidue"],db_Gabon$Code[db_Gabon$pheno%in%"décidue"],function(x)length(x))/db_meta$N
db_meta$prop_old=tapply(db_Gabon$stol[db_Gabon$stol%in%"tolérante à l'ombrage"],db_Gabon$Code[db_Gabon$stol%in%"tolérante à l'ombrage"],function(x)length(x))/db_meta$N
db_meta$prop_pionn=tapply(db_Gabon$stol[db_Gabon$stol%in%"pionnière"],db_Gabon$Code[db_Gabon$stol%in%"pionnière"],function(x)length(x))/db_meta$N
db_meta$prop_wind=tapply(db_Gabon$pheno[db_Gabon$diasp%in%"anémochore"],db_Gabon$Code[db_Gabon$diasp%in%"anémochore"],function(x)length(x))/db_meta$N
db_meta$prop_zoo=tapply(db_Gabon$pheno[db_Gabon$diasp%in%"zoochore"],db_Gabon$Code[db_Gabon$diasp%in%"zoochore"],function(x)length(x))/db_meta$N

### FINAL PREPARATION OF DB
# Remove plot with pb
db_meta=db_meta[rownames(db_meta)!="MIK-407",] # oulier, problem with inventory
db_meta=db_meta[rownames(db_meta)!="OZRI.52",] # MANGROVES

setwd("C:/Users/FAO_bast/Desktop/Gabon_Bastin_script_data")
write.csv(db_meta,"db_Gabon_24072020.csv")

###
db_meta=read.csv("db_Gabon_24072020.csv")
rownames(db_meta)<-db_meta$X
db_meta=db_meta[,-1]
head(db_meta)
