###------------------###
### 0. install packages ###
###------------------###
install.packages("visNetwork")
install.packages("htmlwidgets")
install.packages("igraph")
install.packages("reshape")
install.packages("Matrix")
install.packages("RSiena")
install.packages("networkD3")
install.packages("curl")
install.packages("devtools")
library(devtools)
devtools::install_github("PABalland/EconGeo", force = T)

###------------------###
### 0. Load packages ###
###------------------###

# http://econ.geo.uu.nl/peeg/peeg1709.pdf
library (EconGeo)
library (igraph)

###-------------------------------------###
### 1. The economy as a matrix (2-mode) ###
###-------------------------------------###

# data from ?RCA
?RCA

# ingresar el dataframe original
file.choose()
m = read.csv("/Users/luz/Documents/GitHub/Temas selectos3/LYRA_LAB45 ECONGEO EN R /coop_2019.csv")
class(m)
names(m)


# obtener la matriz original # convertir un dataframe(lista) en matriz
mat = get.matrix(m)
head(mat)
mat #matriz de incidencias (original)

# practice
# transpose of the matrix
# from matrix to data frame
# back to matrix format

# solution
t(mat) #transpose of the matrix
data = get.list (mat) #from matrix to data frame
data
mat = get.matrix (data) #back to matrix format
mat

###-------------------------------------###
### 2. Spatial concentration indicators ###
###-------------------------------------###

?RCA
RCA
RCANOBIN <- as.data.frame(RCA (mat)) #location.quotient
RCABIN <- as.data.frame(RCA (mat, binary = TRUE)) #location.quotient

BALASSA_PROM <- as.data.frame(location.quotient.avg (mat))  #concentraci?n, RCA PROM EN LOS ESTADOS
LOCAT_GINI <- as.data.frame(locational.Gini (mat))  #concentraci?n, GINI EN LOS SECTORES
locational.Gini.curve (mat) #concentraci?n/desigualdad



# in-put incidence matrix
HACHMAN <- as.data.frame(Hachman (mat))  #especializaci?n por estado
HERFINDAHL <- as.data.frame(Herfindahl (mat))  #especializaci?n por estado
INV_UBICUITY <-  as.data.frame(inv.norm.ubiquity (mat))  #ubicuidad-como medida de diversidad por sectores cooperativos
KRUGMAN <- as.data.frame(Krugman.index (mat))   #especializaci?n por estado
ENTROPIA <- as.data.frame(entropy(mat))  #diversidad por estado
DIVERSITY <- as.data.frame(EconGeo::diversity (mat, RCA = TRUE))   #diversidad por estados
DIVERSIDAD <- as.data.frame( MORc(mat, RCA = T, steps = 0))
UBICUIDAD_PROMEDIO <-  as.data.frame(MORc(mat, RCA = T, steps = 1))

GINI <- as.data.frame(Gini (mat))  #concentraci?n, GINI EN LOS SECTORES
Gini (rowSums(mat))
## run the function for industry #1 only 
Gini (mat[,1])
## run the function for industry #2 only 
Gini (mat[,2])
## run the function for industry #3 only 
Gini (mat[,3])
## run the function for industry #4 only 
Gini (mat[,4])

#exportar en formato de excel
write.csv (RCABIN, file="RCABIN.csv")
write.csv (BALASSA_PROM, file="BALASSA_PROM.csv")
write.csv (LOCAT_GINI, file="LOCAT_GINI.csv")
write.csv (HACHMAN, file="HACHMAN.csv")
write.csv (HERFINDAHL, file="HERFINDAHL.csv")
write.csv (INV_UBICUITY, file="INV_UBICUITY.csv")
write.csv (KRUGMAN, file="KRUGMAN.csv")
write.csv (ENTROPIA, file="ENTROPIA.csv")
write.csv (DIVERSITY, file="DIVERSITY.csv")
write.csv (UBICUIDAD_PROMEDIO, file="UBICUIDAD_PROMEDIO.csv")
write.csv (GINI, file="GINI.csv")

#Ver la ruta en donde se guardo
getwd()

###############################################
# SACAR HEATMAP DE LA MATRIZ RCA BINARIA #####
##############################################

install.packages("pheatmap")
library(pheatmap)
# viridis, magma, plasma, cividis, inferno
install.packages("viridisLite")
library(viridis)

file.choose()

heatmap_1 <- as.matrix(
  read.csv("/Users/luz/Documents/Rstudio/ECONGEO_FULL/RCABIN.csv",
           sep = "," ,
           header = T,
           row.names = 1))


#Plotting with pheatmap!
pheatmap(heatmap_1)

colorz <- c('white', 'navyblue')

library(viridis)
pheatmap(heatmap_1, frontsize = 6, cluster_rows = T, cluster_cols = F, treeheight_row = 1, treeheight_col = 0, 
         main = "", fontsize = 12, annotation_legend = FALSE, display_numbers = FALSE, 
         fontsize_number = 6, col = colorz)



###--------------------------###
### 3. Measuring relatedness ###
###--------------------------###

# Counting & Normalizing Co-Occurrences
# remember
RCABIN <- as.data.frame(RCA (mat, binary = TRUE)) #location.quotient

#CONTINUA
co.occurrence (RCABIN)
c = co.occurrence (t(RCABIN))

#matriz de adyacencias 1 (relacionamiento con co-ocurrencias)
r1 <- relatedness(c)
R1 <- as.data.frame(r1)  


###-------------------------------###
### 4. Plotting the product space ###
###-------------------------------###

# g1, g2 y g5

# graficar con r1********************
g1 = graph_from_adjacency_matrix(r1)
plot(g1)

#MATRIZ DE ADYACENCIA NO BINANRIA CON R1
write.csv (r1, file="relatedness_nobin_R1.csv")
write.graph(g1, file="g1.gml", format="gml")
#Ver la ruta en donde se guardo
getwd()

#opcional (s?lo para trabajar con datos binarios)
# RELACIONAMIENTO BINARIO CON CO OCURRENCIAS
r2 = relatedness(c)
 
r2[r2<1] = 0
r2[r2>1] = 1

#matriz de adyacencis binarias con relacionamiento de co-ocurrencias
R2 <- as.data.frame(r2)

# PRE VIEW DE LA RED ### VER LA ESTRUCTURA DE LA RED DENTRO DE R
g2 = graph_from_adjacency_matrix (r2, mode = "undirected")
plot(g2)

#MATRIZ DE RELACIONAMIENTO BINARIA
write.csv (r2, file="relatedness_bin_R2.csv")
#aqui se guarda la red para exportala en un segundo momento en CYTOSCAPE******
write.graph(g2, file="g2.gml", format="gml")


# 3. PROXIMIDAD CON RELACIONAMIENTO DE COSINE
r3 = relatedness (c, method = "cosine")
write.csv (r3, file="relatedness_R3.csv")
# PRE VIEW DE LA RED ### VER LA ESTRUCTURA DE LA RED DENTRO DE R
g3 = graph_from_adjacency_matrix (r3, mode = "undirected")
plot(g3)
getwd()

####
# 4. PROXIMIDAD CON RELACIONAMIENTO DE JACARD
r4 = relatedness (c, method = "Jaccard")
write.csv (r4, file="relatedness_R4.csv")
# PRE VIEW DE LA RED ### VER LA ESTRUCTURA DE LA RED DENTRO DE R
g4 = graph_from_adjacency_matrix (r4, mode = "undirected")
plot(g4)
getwd()


#####
# 5. PROXIMIDAD CON RELACIONAMIENTO DE ASSOCIATION **************
r5 = relatedness (c, method = "association")
write.csv (r5, file="relatedness_R5.csv")
# PRE VIEW DE LA RED ### VER LA ESTRUCTURA DE LA RED DENTRO DE R
g5 = graph_from_adjacency_matrix (r5, mode = "undirected")
plot(g5)
getwd()

library(EconGeo)
gephi= get.list(r5)
write.csv(gephi, file="relatedness_R5lista.csv")

getwd()



# For complex networks
# Maximum Spanning Tree
# Backbone of the network
# Rule 1: keep n-1 links maximum
# Rule 2: remove the links with the lowest weight
# Rule 3: do not create isolates

M <- matrix(runif(200*200, min=0, max=100), ncol=200)
diag(M) <- 0
head (M[,1:6])

g <- graph.adjacency(M, mode="undirected", weighted=TRUE)
plot (g)

M <- - M # very important step
head (M[,1:6])
g <- graph.adjacency(M, mode="undirected", weighted=TRUE)
MST <- minimum.spanning.tree(g)
E(MST)$weight
E(MST)$weight = E(MST)$weight * (-1)
E(MST)$weight
plot (MST, vertex.shape="none", vertex.label.cex=.7)
A <- get.adjacency(MST, sparse = F)


###------------------------###
### 5. Relatedness Density ###
###------------------------###

mat = RCA(mat, binary = T) #compute binary RCA before computing relatedness density
rd = relatedness.density(mat, r)
rd
rd = get.list (rd)
rd

###---------------------###
### 6. Predicting entry ###
###---------------------###

?entry.list # grab data
set.seed(31)
mat1 <- matrix(sample(0:1,20,replace=T), ncol = 4)
rownames(mat1) <- c ("R1", "R2", "R3", "R4", "R5")
colnames(mat1) <- c ("I1", "I2", "I3", "I4")
mat2 <- mat1
mat2[3,1] <- 1

d = entry.list (mat1, mat2)
colnames (d) = c("Region", "Industry", "Entry", "Period")
d = merge (d, rd, by = c("Region", "Industry"))
summary (lm(d$Entry ~ d$Count))


###----------------------------------------------------------------------------------------###
### 6.1 growth.ind Generate a matrix of industrial growth by industries from two regions
# - industries matrices (same matrix composition from two different periods)               ###
###----------------------------------------------------------------------------------------###

## generate a first region - industry matrix with full count (period 1)
set.seed(31)
mat1 <- matrix(sample(0:10,20,replace=T), ncol = 4)
rownames(mat1) <- c ("R1", "R2", "R3", "R4", "R5")
colnames(mat1) <- c ("I1", "I2", "I3", "I4")
## generate a second region - industry matrix with full count (period 2)
mat2 <- mat1
mat2[3,1] <- 8
## run the function
growth.ind (mat1, mat2)

## run the function Generate a matrix of industrial growth in regions from two regions -industries matrices (same matrix composition from two different periods)

growth.mat (mat1, mat2)

## run the function Generate a matrix of industrial growth by regions from two regions - industries matrices (same matrix composition from two different periods)
growth.reg (mat1, mat2)


###---------------------------###
### 7. Ubiquity and diversity ###
###---------------------------###

detach("package:igraph", unload=TRUE) #unload igraph - same name function
mat
ubiquity (mat)
diversity (mat)

###----------------###
### 8. ECI and PCI ###
###----------------###

KCI(mat)
TCI(mat)

TCI.norm<- as.data.frame(TCI(mat, RCA=T))
rownames(TCI.norm) <- colnames(mat)
colnames(TCI.norm)[1] <- "TCI"
#normalise TCI so it ranges between 0 and 1
TCI.norm = (TCI.norm-min(TCI.norm))
  
MORc(mat)
MORt(mat)

mat
###-------------------------------------------------------------------------------###
### 9. Expy Compute the expy index of regions from regions - industries matrices  ###
###-------------------------------------------------------------------------------###

## generate a region - industry matrix
set.seed(31)
mat <- matrix(sample(0:100,20,replace=T), ncol = 4)
rownames(mat) <- c ("R1", "R2", "R3", "R4", "R5")
colnames(mat) <- c ("I1", "I2", "I3", "I4")
## a vector of GDP of regions
vec <- c (5, 10, 15, 25, 50)
## run the function
expy (mat, vec)

###-------------------------------------------------------------------------------###
# 10. Compute the prody index of industries from regions - industries matrices
###-------------------------------------------------------------------------------###

## generate a region - industry matrix
set.seed(31)
mat <- matrix(sample(0:100,20,replace=T), ncol = 4)
rownames(mat) <- c ("R1", "R2", "R3", "R4", "R5")
colnames(mat) <- c ("I1", "I2", "I3", "I4")
## a vector of GDP of regions
vec <- c (5, 10, 15, 25, 50)
## run the function
prody (mat, vec)








###---------------------###
### 10. Using real data ###
###---------------------###

# load from URL
pat = "https://raw.githubusercontent.com/PABalland/ON/master/pat.csv"
df = read.csv (pat)
df = subset (df, dec == 2000)
df = df[, c("Cbsa.Name", "NBER.Sub.Cat.Name", "pat.count")]
m = get.matrix (df)
r = relatedness (co.occurrence(t(m)))
rd = relatedness.density(RCA(m, binary = T), r)
m[m<10] = 0 
m = RCA (m)
m = m[rowSums(m) > 15,]
t = data.frame (colnames (m), MORt (RCA(m, binary = T)))
k = data.frame (rownames (m), MORc (RCA(m, binary = T)))

### Key references
#Balland, P.A., Boschma, R., Crespo, J. and Rigby, D. (2018)  Smart Specialization policy in the EU: Relatedness, Knowledge Complexity and Regional Diversification, Regional Studies, forthcoming
#Balland, P.A. and Rigby, D. (2017) The Geography of Complex Knowledge, Economic Geography, 93 (1): 1-23
#Boschma, R., Balland, P.A. and Kogler, D. (2015) Relatedness and Technological Change in Cities: The rise and fall of technological knowledge in U.S. metropolitan areas from 1981 to 2010, Industrial and Corporate Change, 24 (1): 223-250 
#Balland, P.A., Jara-Figueroa, C., Petralia, S., Steijn, M., Rigby, D., and Hidalgo, C. (2018) Complex Economic Activities Concentrate in Large Cities, Papers in Evolutionary Economic Geography, 18 (29): 1-10
#Balland, P.A. (2017) Economic Geography in R: Introduction to the EconGeo Package, Papers in Evolutionary Economic Geography, 17 (09): 1-75
#Hidalgo, C., Balland, P.A., Boschma, R., Delgado, M., Feldman, M., Frenken, K., Glaeser, E., He, C., Kogler, D., Morrison, A.,  Neffke, F., Rigby, D., Stern, S., Zheng, S., and Zhu, S. (2018)  The Principle of Relatedness, Proceedings of the 20th International Conference on Complex Systems, forthcoming. 
#Hidalgo, C. A., Klinger, B., Barab?si, A. L., & Hausmann, R. (2007). The product space conditions the development of nations. Science, 317(5837), 482-487.
#Hidalgo, C. A., & Hausmann, R. (2009). The building blocks of economic complexity. Proceedings of the national academy of sciences, 106(26), 10570-10575.