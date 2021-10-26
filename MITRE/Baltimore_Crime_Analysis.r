library(shiny)
library(vroom)
library(tidyverse)
library(ggplot2)
library(gmodels)
library("dplyr")
library(scatterplot3d)
library(rgl)
library(GPArotation)
library(psych)
library(dplyr)
library(ggdendro)
library(tidyr)

LENGTH <- 17 #number of crime categories
sessionInfo() #to see versions

#=====================================================
#Basis crime data and selection by post mechanism
#=====================================================
crimes <- vroom::vroom("data/BPD.csv")
#use these if doing a line graph in ggplot
#crimes$Weapon[is.na(crimes$Weapon)] <- as.integer(0)
#crimes$Weapon[crimes$Weapon=="HANDS"] <- as.integer(1) 
#crimes$Weapon[crimes$Weapon=="OTHER"] <- as.integer(2) 
#crimes$Weapon[crimes$Weapon=="KNIFE"] <- as.integer(3) 
#crimes$Weapon[crimes$Weapon=="FIRE"] <- as.integer(4) 
#crimes$Weapon[crimes$Weapon=="FIREARM"] <- as.integer(5) 
#crimes$Weapon <- as.integer(crimes$Weapon)

crimes<-as_tibble(crimes)

selected <- crimes %>% filter(Post == 111)
nrow(selected)

selected %>% count(Description, sort = TRUE)
selected %>% count(Weapon,  sort = TRUE)
selected %>% count(Premise,  sort = TRUE)
selected %>% count(Premise,Description,  sort = TRUE)
summary <- selected %>% 
  count(Premise, Description, Weapon, wt = Count) 
summary

#THIS DOESN'T WORK BECAUSE BOTH P[REMISE AND DESC ARE CHAR, 
#BOTH THE X AND Y COORS IN GGPLOT NEED TO BE NUMERIC
#WICKHAM CONVERTED AGE LEVEL TO NUMERIC
#summary %>% 
#  ggplot(aes(Weapon, n, colour = Description)) + 
#  geom_line(na.rm = TRUE) + 
#  labs(y = "Injuries ")

#===============================================
#PCA
#===============================================

#the way to do a cross tab within R is shwon after teh tabset section
crime <- read.csv("data/BPDp.csv")
crime[is.na(crime)]<-0
rownames(crime)<-crime$Post
crime<-select(crime, -Post)

str(crime)
#the curse of dimesionality
pairs(crime)
crime_pc <- princomp(crime, cor = TRUE)

#try this so can set num factors
library(psych)
crime_pc<-principal(crime, nfactors=7, score=TRUE)
#just show object to find loadings - 1 PC explains 68% of variance
crime_pc<-principal(crime, nfactors=3)
crime_pc 

summary(crime_pc, loadings = TRUE)
#18 is a magic cookie - the number of crime categories
xlim <- range(crime_pc$scores[, 1])
plot(crime_pc$scores[, 1], crime_pc$scores[, 2],
     xlab = "First PC score", ylab = "Second PC score", xlim = xlim, ylim = xlim, type = "n")
#text(crime_pc$scores[, 1], crime_pc$scores[, 2], row.names(crime), cex=1.0) 
text(crime_pc$scores[, 1], crime_pc$scores[, 2],rownames(crime), cex=1.0)

#get rid of outlier
gmrNames<-rownames(crime)
crime<-filter(crime, rownames(crime)!=141)
#crime <- read.csv("data/BPDp_no_outlier.csv")
gmrName<-gmrNames[!gmrNames %in% c(141)]
rownames(crime)<-gmrName

str(crime)
head(crime)
#convert empty cells to 0
crime[is.na(crime)]<-0
crime_pc <- princomp(crime, cor = TRUE)

#to get proportions
summary(crime_pc)


#plot again
#note 1st dimension of crime_pc$scores matches to crime$Row.labels both 179 in length
xlim <- range(crime_pc$scores[, 1])
plot(crime_pc$scores[, 1], crime_pc$scores[, 2],
     xlab = "First PC score", ylab = "Second PC score", xlim = xlim, ylim = xlim, type = "n")
#text(crime_pc$scores[, 1], crime_pc$scores[, 2], crime$Row.Labels, cex=1.0) 
text(crime_pc$scores[, 1], crime_pc$scores[, 2], rownames(crime), cex=1.0)
#crime_pc$scores[,1][crime$Row.Labels==111] to get the 1st PC score for Post 111
#least dangerous neighborhood is 
#min(crime_pc$scores[,1])
#most dangerous (excluding outlier)
#max(crime_pc$scores[,1])
#your precinct is
#crime_pc$scores[,1][crime$Row.Labels==111]
#p1<-crime_pc$scores[,1][crime$Row.Labels==111]
#p2<-crime_pc$scores[,1][crime_pc$scores[,1] > x]
#p<-length(p2)/length(crime_pc$scores[,1])*100


#3d plot
attach(crime_pc)
scatterplot3d(scores[,1],scores[,2],scores[,3], main="Neighborhood Scoring on 3 PCAs")

scatterplot3d(scores[,1],scores[,2],scores[,3], 
    pch=16,
    highlight.3d=TRUE,
    type="h",
    xlab = "PCA 1: Overall Dangerousness", ylab = "PCA 2: Relative Level of Crime on Persons",
    main="Neighborhood Scoring on 3 PCAs")
gmrs3d<-scatterplot3d(scores[,1],scores[,2],scores[,3], 
    pch=16,
    highlight.3d=TRUE,
    type="h",
    main="Neighborhood Scoring on 3 PCAs")
fit<-lm(scores[,1]~scores[,2]+scores[,3])
gmrs3d$plane3d(fit)

lambda <- crime_pc$sdev^2
Astar <- crime_pc$loadings[1:LENGTH, 1:LENGTH]
R <- Astar %*% diag(lambda) %*% t(Astar)
R

#===============================================
#EFA
#===============================================
crime_fa1 <- factanal(crime, factors = 1, rotation = "none", lower=.01)
crime_fa2 <- factanal(crime, factors = 2, rotation = "none", lower=.01)
crime_fa3 <- factanal(crime, factors = 3, rotation = "none", lower=.01)

crime_fa1
crime_fa2
crime_fa3

crime_fa3None <- fa(crime, nfactors = 3, fm = "ml", rotate = "none")
print.psych(crime_fa3None, digits = 3)

crime_fa3varimax <- fa(crime, nfactors = 3, fm = "ml", rotate = "varimax")
print.psych(crime_fa3varimax, digits = 3, sort = TRUE)

crime_fa3 <- factanal(crime, factors = 3, method = "mle", rotation = "varimax", 
                      scores = "regression", lower=.01)

dimnames(crime_fa3$scores)[[2]]
propen<-c("Propensity to Property Crime","Propensity to Crime on People","Propensity to Vehicle Crime")
dimnames(crime_fa3$scores)[[2]]<-propen

pairs(crime_fa3$scores, 
      #panel = function(x,y) text(x, y, labels = crime$Row.Labels, cex=0.8))
      panel = function(x,y) text(x, y, labels = rownames(crime), cex=0.8))
#===============================================
#Cluster
#===============================================
#rlabs <- row.names(crime)
rlabs<-crime$Row.Labels

# standardize by range
# 2 is column
#below gets range max in a column - min in column
rge <- apply(crime, 2, max) - apply(crime, 2, min)
#sweep creates an array from another array by sweeping out a summary stat
#divies values by column range
#try other ways of stadizing such a z scores

crime_std <- sweep(crime, 2, rge, FUN = "/")

# variances of the std data:
apply(crime_std, 2, var)

# plot of wgss against number of clusters:
n <- length(crime_std[, 1])
wss1 <- (n-1) * sum(apply(crime_std, 2, var))
wss <- numeric(0)
for(i in 2:6) {
    W <- sum(kmeans(crime_std, i)$withinss)
    wss <- c(wss, W)
}
wss <- c(wss1, wss)
plot(1:6, wss, type = "l", 
     xlab = "Number of groups", ylab = "Within groups sum of squares", lwd=2)


# get two-group solution from k-means and group means and membership
crime_kmean2 <- kmeans(crime_std, 2)
lapply(1:2, function(nc) apply(crime[crime_kmean2$cluster == nc, ], 2, mean))
lapply(1:2, function(nc) rlabs[crime_kmean2$cluster == nc])

crime_pc <- princomp(crime_std, cor = TRUE)
xlim <- range(crime_pc$scores[, 1])

plot(crime_pc$scores[, 1:2], type = "n", xlim = xlim, ylim = xlim)
text(crime_pc$scores[, 1:2], labels = crime_kmean2$cluster, cex=0.8)

#crime$Row.Labels
#toggle between showing two clusters and the post numbers
#WSS says two natural clsuters - 1 is dnagerous precincts, 2 is safer

plot(crime_pc$scores[, 1:2], type = "n", xlim = xlim, ylim = xlim)
#text(crime_pc$scores[, 1:2], labels = crime$Row.Labels, cex=0.8)
text(crime_pc$scores[, 1:2], labels = rownames(crime), cex=0.8)

#===============================================
#Shiny
#MUST RUN AS ADMIN BECAUSE IT OPENS PNG FILES IN RESTRICTED AREAS
#===============================================
#least dangerous neighborhood is 
#min(crime_pc$scores[,1])
#most dangerous (excluding outlier)
#max(crime_pc$scores[,1])
#your precinct is
#crime_pc$scores[,1][crime$Row.Labels==111]

#===========================================================
#Tabsets
#===========================================================
ui <- fluidPage(
  titlePanel("Principal Component Explorer"),
  sidebarLayout(
    sidebarPanel(
      #selectInput("code", "Precinct", choices = sort(unique(crimes$Post))),
      #selectInput("code2", "Precinct", choices = sort(unique(crimes$Post))),
      selectInput("code", "Precinct", choices = sort(unique(rownames(crime)))),
      selectInput("code2", "Precinct", choices = sort(unique(rownames(crime)))),
      selectInput("Map", "Pick a section", choices = c("BPD_North","BPD_South", "BPD_South_Tip")),

    ),
    mainPanel({
      tabsetPanel(
        id = "tabset",
        tabPanel("Crimes in Precinct", "Crime Index", 
           textOutput("postPCA"),
           textOutput("percPCA"),
           "---------------------------------",
           textOutput("maxPCA"),
           textOutput("minPCA"),
           "---------------------------------",
           textOutput("totCrime"),
           textOutput("vioCrime"),
           column(4,textOutput("gmrHead"),
           tableOutput("Desc")
           ),
           column(4,tableOutput("DescPct")
           )

        ),
        tabPanel("Histogram", "Crime Bar Chart", plotOutput("weapon_desc"),),
        tabPanel("3 Latent Factors", "Three Neighborhood Propensities, % Balance between pairs of Principal Components", plotOutput("EFA")),
        tabPanel("Two Principle Components", "Crime Index by Crime Shape", plotOutput("PCA"), column(4,tableOutput("Desc2PC")), column(4,tableOutput("Desc2_2PC"))),
        tabPanel("Three Principle Components", "Crime Index by Crime Shape", plotOutput("PCA3")),
        tabPanel("Baltimore PD Map", "Southern Posts",imageOutput("photo2"))

      )
    })
  )
)

# 

#shinyApp(ui, server)

server <- function(input, output, session) {
  selected <- reactive(crimes %>% filter(Post == input$code))
  selected2 <- reactive(crimes %>% filter(Post == input$code2))

  gmr <- reactive(selected() %>% count(Description))
  gmr2 <- reactive(selected2() %>% count(Description))

  output$gmrHead <- renderText({ 
         x<-"Crimes committed in this Post"
         x
  })

  output$postPCA <- renderText({
         #x<-crime_pc$scores[,1][crime$Row.Labels==input$code]
         x<-crime_pc$scores[,1][rownames(crime)==input$code]
         y<-c("This Post's Crime Index is ", round(x,2))
         y
      }
  )

  output$percPCA <- renderText({
         #x<-crime_pc$scores[,1][crime$Row.Labels==input$code]
         #y<-crime_pc$scores[,1][crime_pc$scores[,1] > x]
         x<-crime_pc$scores[,1][rownames(crime)==input$code]
         y<-crime_pc$scores[,1][crime_pc$scores[,1] > x]
         p<-length(y)/length(crime_pc$scores[,1])*100
         ret<-c(round(p,2),"% of posts have higher crime index")
         ret
      }
  )

  output$maxPCA <- renderText({
         x<-max(crime_pc$scores[,1])
         y<-c("The Highest Crime Post has index  ", round(x,2))
         y
      }
  )

  output$head1<- renderText("Crimes committed in this Post")

  output$minPCA <- renderText({
         x<-min(crime_pc$scores[,1])
         y<-c("The Lowest Crime Post has index  ", round(x,2))
         y
      }
  )

  output$totCrime <- renderText({
         x<-selected() %>% count()
         y<-c("The Total Crime Count at this Post is  ", x[[1]])
         y
      }
  )

  output$vioCrime <- renderText({
         x<-selected() %>% count(Description, sort = TRUE)
         vioN<-0
         for(i in 1:length(x$n)) {
           #if(x$Description[i]=="I") next
           #if(x$Description[i]=="O") next
           if(is.na(x$Description[i])) message(x$Description[i])
           else{
              if(str_detect(x$Description[i], "ASSAULT", negate = FALSE)){vioN<-vioN + x$n[i]}
              if(str_detect(x$Description[i], "SHOOTING", negate = FALSE)){vioN<-vioN + x$n[i]}
              if(str_detect(x$Description[i], "RAPE", negate = FALSE)){vioN<-vioN + x$n[i]}
              if(str_detect(x$Description[i], "HOMICIDE", negate = FALSE)){vioN<-vioN + x$n[i]}
           }
         }
         y<-c("The Violent Crime Count at this Post is  ", vioN)
         y
      }
  )

  output$Desc <- renderTable(
    selected() %>% count(Description, sort = TRUE)
  )

  output$Desc2PC <- renderTable(
    selected() %>% count(Description, sort = TRUE)
  )

  output$Desc2_2PC <- renderTable(
    selected2() %>% count(Description, sort = TRUE)
  )


  output$DescPct <- renderTable({
    x<-selected() %>% group_by( Description ) %>% summarise( percent = 100 * n() / nrow( selected() ) )
    y<-arrange(x, desc(percent))
    y
    #selected() %>% count(Description, sort = TRUE)

  })

  output$ggPCA <- renderPlot({
    score_df<-as.data.frame(crime_pc$score)
    ggplot(score_df, aes(Comp.1, Comp.2), geom_point())  
  }, res = 96)


  output$PCA <- renderPlot({
    plot(crime_pc$scores[, 1], crime_pc$scores[, 2],
        xlab = "PCA 1: Overall Crime Index", ylab = "PCA 2: Crime on Property to Crime on Persons", xlim = xlim, ylim = xlim, type = "n")
        #text(crime_pc$scores[, 1], crime_pc$scores[, 2], crime$Row.Labels, cex=1.0)  
        text(crime_pc$scores[, 1], crime_pc$scores[, 2], rownames(crime), cex=1.0)  
  }, res = 96)

  output$EFA <- renderPlot({
     pairs(crime_fa3$scores, 
        #panel = function(x,y) text(x, y, labels = crime$Row.Labels, cex=0.8)) 
        panel = function(x,y) text(x, y, labels = rownames(crime), cex=0.8)) 

  }, res = 96, height = 600, width = 600)

  output$data <- renderTable({
    brushedPoints(crime_pc$scores, input$plot_brush)
  })

  output$PCA3 <- renderPlot({
    attach(crime_pc)
    #scatterplot3d(scores[,1],scores[,2],scores[,3], 
    #xlab = "PCA 1: Overall Crime Index", ylab = "PCA 2: Crime on Property to Crime on Persons",
    #main="Neighborhood Scoring on 3 PCAs")
    zz <- scatterplot3d(x = scores[,1], y = scores[,2], z = scores[,3], 
                    xlab = "crime index", ylab = "Property v Person", zlab = "Auto (mobile) Crime", 
                    pch = 1, color = "blue", grid = TRUE)

    zz.coords <- zz$xyz.convert(scores[,1],scores[,2],scores[,3]) 

    text(zz.coords$x, 
     zz.coords$y,             
     labels = rownames(crime),               
     cex = .5, 
     pos = 4)  

  }, res = 96,  height = 600, width = 600)

  output$weapon_desc <- renderPlot({
    barplot(gmr()$n, names.arg=gmr()$Description)
  }, res = 96)

  output$photo <- renderImage({
    list(
      src = file.path("./", paste0("BPD_North", ".jpg")),
      contentType = "image/jpeg",
      width = 1200,
      height = 800
    )
  }, deleteFile = FALSE)

  output$photo2 <- renderImage({
    list(
      src = file.path("./", paste0(input$Map, ".jpg")),
      contentType = "image/jpeg",
      width = 1200,
      height = 800
    )
  }, deleteFile = FALSE)

}


shinyApp(ui, server)

crimes %>% group_by(Description) %>% summarize(Number = sum(Count))
crimes %>% group_by(Description, Post) %>% summarize(Number = sum(Count))
gmr<-crimes %>% group_by(Description, Post) %>% summarize(Number = sum(Count)) %>% arrange(Description, Number)
View(gmr)

#how to do cross tab table
crime<-crimes %>% group_by(Post, Description) %>% summarise(Cnt = sum(Count)) %>% spread(Description, Cnt, sep='')
crime<-as.data.frame(crime)
crime[is.na(crime)]<-0
rownames(crime)<-crime$Post
crime<-select(crime, -Post)
#one post was NA that need removal
crime<-crime[-180,]

#how to use count so don't need a oneCount column
crime<-crimes %>% group_by(Post, Description) %>% summarise(Cnt = n()) %>% spread(Description, Cnt, sep='')
crime<-as.data.frame(crime)
crime[is.na(crime)]<-0
rownames(crime)<-crime$Post
crime<-select(crime, -Post)
#one post was NA that need removal
crime<-crime[-180,]

#Graph
crime<-crimes %>% group_by(Post, Description) %>% summarise(Cnt = n()) 
ggplot(crime) + geom_bar(mapping = aes(x=Description))
