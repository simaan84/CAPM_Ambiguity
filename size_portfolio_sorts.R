library(data.table)
library(lubridate)
library(dplyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(quantmod)
library(stargazer)

rm(list = ls())
gc()

###############################################################################

#####################################
######### RECESSION PERIODS #########
#####################################

{
  bar.col <- "gray"
  ggplot_recession0 <- geom_rect(fill = bar.col,col = bar.col,
                                 aes(xmin=date("1960-04-30"),
                                     xmax=date("1961-02-28"), ymin=-Inf, ymax=Inf))
  
  ggplot_recession1 <- geom_rect(fill = bar.col,col = bar.col,
                                 aes(xmin=date("1969-12-31"),
                                     xmax=date("1970-11-30"), ymin=-Inf, ymax=Inf))
  
  ggplot_recession2 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("1973-11-30"), xmax=date("1975-03-31"), ymin=-Inf, ymax=Inf))
  ggplot_recession3 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("1980-01-31"), xmax=date("1980-07-31"), ymin=-Inf, ymax=Inf))
  ggplot_recession4 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("1981-07-31"), xmax=date("1982-11-30"), ymin=-Inf, ymax=Inf))
  ggplot_recession5 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("1990-07-31"), xmax=date("1991-03-31"), ymin=-Inf, ymax=Inf))
  ggplot_recession6 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("2001-03-31"), xmax=date("2001-11-30"), ymin=-Inf, ymax=Inf))
  ggplot_recession7 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("2007-12-31"), xmax=date("2009-06-30"), ymin=-Inf, ymax=Inf))
  ggplot_recession8 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("2020-02-01"), xmax=date("2020-04-30"), ymin=-Inf, ymax=Inf))
  
  recession_list <- list(ggplot_recession0,ggplot_recession1,
                         ggplot_recession2,ggplot_recession3,
                         ggplot_recession4,ggplot_recession5,
                         ggplot_recession6,ggplot_recession7,
                         ggplot_recession8)
  
  
  
}


###############################################################################
#########################################
###### DATA PROCESSING AND CLEANING #####
#########################################


######################
##### FAMA-FRENCH ####
######################

FF_file <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_CSV.zip"
temp <- tempfile()
download.file(FF_file,temp)
unz_files <- unzip(temp)
ds <- read.csv(unz_files,skip = 2,stringsAsFactors = FALSE)
ds <- data.frame(apply(ds,2,as.numeric))
names(ds)[1] <- "date"
ds$date <- ds$date*100  + 1
ds$date <- ymd(ds$date)
ds$date <- ceiling_date(ds$date,"m") - 1
ds[,-1] <- ds[,-1]/100
ds <- ds[,c("date","SMB")]

########################################
####### LIQUIDITY RISK FACTOR ##########

# this data comes from https://www.openassetpricing.com/data/
CZ_ds <- read.csv("Open_Source_Asset_Pricing_2024.csv")
CZ_ds <- CZ_ds[,c("date","Illiquidity")]
CZ_ds[,-1] <- CZ_ds[,-1]/100
CZ_ds$date <- date(CZ_ds$date)
CZ_ds$date <- ceiling_date(CZ_ds$date,"m") - 1
CZ_ds$Illiquidity <- CZ_ds$Illiquidity*12 # make it annual

###########################
#### LOAD VIX ###########

VIX <- getSymbols("^VIX",auto.assign = FALSE,from = "1900-01-01")
VIX <- na.omit(VIX$VIX.Adjusted)
VIX <- data.frame(date = date(VIX), VIX = as.numeric(VIX))

###########################
#### CRSP AND GU ET AL ####
###########################

# this data comes from gu et al RFS paper: https://dachxiu.chicagobooth.edu/download/datashare.zip
data_file <- "GKX94.csv"
select_var <- c("permno","date","ret","mvel1","mve_ia")
DS <- fread(data_file,select = select_var)

# there are no NAs, and all values are unique at the permno-date levels
# DS[,N := .N, by = c("permno","date")]

DS$date <- ymd(DS$date)

# see how the returns matches with one from CRSP
CRSP_file <- "CRSP_1960_2019.csv"
DT <- fread(CRSP_file,select = c("PERMNO","date","RET","PRC","SHROUT"))
DT$permno <- DT$PERMNO
DT$ret_crsp <- DT$RET
DT$PERMNO <- DT$RET <- NULL
DT$ret_crsp <- as.numeric(DT$ret_crsp)
DT$date <- ymd(DT$date)
DT$MKTCAP <- abs(as.numeric(DT$PRC)*as.numeric(DT$SHROUT))
DT$PRC <- DT$SHROUT <- NULL

# add MKTCAP to 
DS <- merge(DS,DT[,list(date,permno,MKTCAP)],by = c("date","permno"))
rm(DT);gc();


###########################################################################

##################################
######## PORTFOLIO SORTS #########
##################################

DS <- DS[order(DS$permno,DS$date),]
DS[,ret_next := shift(ret,-1), by = "permno"]

# we can group stocks based on each
G <- 10
is_equal <- FALSE
{
  DS_port <- DS
  DS_port[,size_group := ntile(mvel1,G),by = "date"]
  if(is_equal)
    DS_port <- DS_port[,list(port_ret = mean(ret_next,na.rm = TRUE)),by = c("size_group","date")]
  if(!is_equal)
    DS_port <- DS_port[,list(port_ret = sum(ret_next*MKTCAP/(sum(MKTCAP)),na.rm = TRUE)),by = c("size_group","date")]
  
  # CAST TABLE TO TIME SERIES OF PORTFOLIO RETURNS
  DS_port <- dcast.data.table(DS_port, date ~ size_group,value.var = "port_ret")
  # adjust the date
  DS_port$date <- ceiling_date(floor_date(DS_port$date,"m")  + months(1),"m") - 1
  DS_port <- na.omit(DS_port)
  names(DS_port)[-1] <- paste("portfolio_",1:G,sep = "")
  DS_port$HML <- (-DS_port[,G+1,with = FALSE] + DS_port[,2,with = FALSE])*12
  
  # add Fama-French SMB and scale both by 12 to reflect annual units
  DS_port <- merge(DS_port,ds,by = "date")
  DS_port$SMB <- DS_port$SMB*12
}

#### produce using ggplot
plot_fun <- function(DS_port,font_size = 10){
  
  ds.plot1 <- DS_port[,list(date,return = cumsum(HML),Type = "Ours" )]
  ds.plot2 <- DS_port[,list(date,return = cumsum(SMB),Type = "Fama_French" )]
  ds.plot <- rbind(ds.plot1,ds.plot2)
  ds.plot$Type <- as.factor(ds.plot$Type)
  
  p <- ggplot(ds.plot, aes(date,return,color = Type,shape = Type))
  p <- p +  ggplot_recession0 + ggplot_recession1 + ggplot_recession2 +
    ggplot_recession3 + ggplot_recession4 + ggplot_recession5 + 
    ggplot_recession6 + ggplot_recession7 + ggplot_recession8
  p <- p + geom_line(alpha = 1)
  p <- p + xlab("Date") + ylab("SMB Cumulative Return")
  p <- p + geom_abline(intercept = 0, slope = 0, color="black",  linetype="dashed", size=0.2)
  p <- p + theme(legend.title = element_text(size=font_size + 5))
  p <- p + theme(legend.text = element_text(size=font_size))
  p <- p + theme(axis.title.x = element_text(size=font_size))
  p <- p + theme(axis.title.y = element_text(size=font_size))
  p <- p + theme(axis.text.x = element_text(face="bold", size=font_size, angle=90))
  p <- p + theme(axis.text.y = element_text(face="bold", size=font_size, angle=90))
  
  p <- p + scale_x_date(limits = c(min(ds.plot$date),max(ds.plot$date)))  # or c(min_date, max_date)
  
  return(p)
}



DS_port <- merge(DS_port,CZ_ds)
{
  min_date <- "1900-01-01"
  # we can plot based on any given subset of data
  DS_port_sub <- DS_port[DS_port$date > min_date,]
  p <- plot_fun(DS_port_sub)
  p
  
  {
    # identify all dates in the recession periods
    min_date <- sapply(recession_list, function(x) as.character(x$computed_mapping$xmin[2]) )
    min_date <- strsplit(min_date,"\"")
    min_date <- sapply(min_date,function(x) x[2])
    
    max_date <- sapply(recession_list, function(x) as.character(x$computed_mapping$xmax[2]) )
    max_date <- strsplit(max_date,"\"")
    max_date <- sapply(max_date,function(x) x[2])
    recession_ds <- data.frame(min_date = date(min_date),
                               max_date = date(max_date))
    
    recession_ds$min_date <- floor_date(recession_ds$min_date,"m")
    recession_ds$max_date <- floor_date(recession_ds$max_date,"m")
    
    recession_dates <- lapply(1:nrow(recession_ds),function(i) 
      seq.Date(from = recession_ds[i,1], to = recession_ds[i,2], by = "month"))
    recession_dates <- lapply(recession_dates,as.character)
    recession_dates <- unlist(recession_dates)
    recession_dates <- date(recession_dates)
    recession_dates <- ceiling_date(recession_dates,"m") - 1
    
    DS_port_sub$Recession <- (DS_port_sub$date %in% recession_dates)*1
  }
  
  ### SUMMARIZE STATISTICAL SIGNIFICANCE ####
  lm1 <- lm(HML ~ 1, data = DS_port_sub)
  lm2 <- lm(SMB ~ 1, data = DS_port_sub)
  # add recession index
  lm3 <- lm(HML ~ Recession, data = DS_port_sub)
  lm4 <- lm(SMB ~ Recession, data = DS_port_sub)
  lm_list <- list(lm1,lm2,lm3,lm4)
  
  lm_list_nw <- lapply(lm_list, function(temp.lm) coeftest(temp.lm, vcov = NeweyWest(temp.lm, lag = 3)) )
  lm_list_nw1 <- lm_list_nw
  lm_list1 <- lm_list
  
  ### ADD ILLIQUIDITY ####
  lm1 <- lm(HML ~ Illiquidity, data = DS_port_sub)
  lm2 <- lm(SMB ~ Illiquidity, data = DS_port_sub)
  # add recession index
  lm3 <- lm(HML ~ Illiquidity + Recession, data = DS_port_sub)
  lm4 <- lm(SMB ~ Illiquidity + Recession, data = DS_port_sub)
  lm_list <- list(lm1,lm2,lm3,lm4)
  
  lm_list_nw <- lapply(lm_list, function(temp.lm) coeftest(temp.lm, vcov = NeweyWest(temp.lm, lag = 3)) )
  lm_list_nw1_liq <- lm_list_nw
  lm_list1_liq <- lm_list
  
  # REPEAT WITH VIX SUB PERIOD
  DS_port_sub <- merge(DS_port_sub,VIX)
  DS_port_sub$Uncertainty <- (DS_port_sub$VIX >= quantile(DS_port_sub$VIX,0.90))*1
  # we can plot based on any given subset of data
  p <- plot_fun(DS_port_sub)
  
  ### SUMMARIZE STATISTICAL SIGNIFICANCE ####
  lm1 <- lm(HML ~ 1, data = DS_port_sub)
  lm2 <- lm(SMB ~ 1, data = DS_port_sub)
  # add recession index
  lm3 <- lm(HML ~ Recession, data = DS_port_sub)
  lm4 <- lm(SMB ~ Recession, data = DS_port_sub)
  lm_list <- list(lm1,lm2,lm3,lm4)
  
  lm5 <- lm(HML ~ Uncertainty, data = DS_port_sub)
  lm6 <- lm(SMB ~ Uncertainty, data = DS_port_sub)
  lm_list <- c(lm_list,list(lm5,lm6))
  
  lm_list_nw2 <- lapply(lm_list, function(temp.lm) coeftest(temp.lm, vcov = NeweyWest(temp.lm, lag = 3)) )
  lm_list2 <- lm_list
  
  # ADD Illiquidity
  lm1 <- lm(HML ~ Illiquidity, data = DS_port_sub)
  lm2 <- lm(SMB ~ Illiquidity, data = DS_port_sub)
  # add recession index
  lm3 <- lm(HML ~ Illiquidity + Recession, data = DS_port_sub)
  lm4 <- lm(SMB ~ Illiquidity + Recession, data = DS_port_sub)
  lm_list <- list(lm1,lm2,lm3,lm4)
  
  lm5 <- lm(HML ~ Illiquidity + Uncertainty, data = DS_port_sub)
  lm6 <- lm(SMB ~ Illiquidity + Uncertainty, data = DS_port_sub)
  lm_list <- c(lm_list,list(lm5,lm6))
  
  lm_list_nw2_liq <- lapply(lm_list, function(temp.lm) coeftest(temp.lm, vcov = NeweyWest(temp.lm, lag = 3)) )
  lm_list2_liq <- lm_list
  
  
  # BANZ PERIOD -1980
  DS_port_sub <- DS_port[DS_port$date < "1981-01-01",]
  DS_port_sub$Recession <- (DS_port_sub$date %in% recession_dates)*1
  
  ### SUMMARIZE STATISTICAL SIGNIFICANCE ####
  lm1 <- lm(HML ~ 1, data = DS_port_sub)
  lm2 <- lm(SMB ~ 1, data = DS_port_sub)
  # add recession index
  lm3 <- lm(HML ~ Recession, data = DS_port_sub)
  lm4 <- lm(SMB ~ Recession, data = DS_port_sub)
  lm_list <- list(lm1,lm2,lm3,lm4)
  
  lm_list_nw <- lapply(lm_list, function(temp.lm) coeftest(temp.lm, vcov = NeweyWest(temp.lm, lag = 3)) )
  lm_list_banz <- lm_list_nw
  lm_list_banz <- lm_list
  
  ## ADD Illiquidity
  lm1 <- lm(HML ~ Illiquidity, data = DS_port_sub)
  lm2 <- lm(SMB ~ Illiquidity, data = DS_port_sub)
  # add recession index
  lm3 <- lm(HML ~ Illiquidity + Recession, data = DS_port_sub)
  lm4 <- lm(SMB ~ Illiquidity + Recession, data = DS_port_sub)
  lm_list <- list(lm1,lm2,lm3,lm4)
  
  lm_list_nw <- lapply(lm_list, function(temp.lm) coeftest(temp.lm, vcov = NeweyWest(temp.lm, lag = 3)) )
  lm_list_banz_liq <- lm_list_nw
  lm_list_banz_liq <- lm_list
  
  
}

stargazer(lm_list_nw1,intercept.bottom = FALSE)
stargazer(lm_list1,intercept.bottom = FALSE)

stargazer(lm_list_nw2,intercept.bottom = FALSE)
stargazer(lm_list2,intercept.bottom = FALSE)

stargazer(lm_list_banz,intercept.bottom = FALSE)
stargazer(lm_list_banz,intercept.bottom = FALSE)

### SUMMARIZE ILLIQUIDITY

stargazer(lm_list_nw1_liq,intercept.bottom = FALSE)
stargazer(lm_list1_liq,intercept.bottom = FALSE)

stargazer(lm_list_banz_liq,intercept.bottom = FALSE)
stargazer(lm_list_banz_liq,intercept.bottom = FALSE)

stargazer(lm_list_nw2_liq,intercept.bottom = FALSE)
stargazer(lm_list2_liq,intercept.bottom = FALSE)


