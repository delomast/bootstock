############################################
# START OF PROGRAM #  Bootstock 3.2
# 14 Febuary 2020

# This program was written primarily to analyze AD-clip only sport fisheries however it can be used for
# unclipped fisheries (i.e Tribal). I recommend that the Adipose clipped and adipose unclipped
# fisheries be run separately.


# Version 3.1 edits and changes ###############

# 1-10-18 fix bug changed code::: names(new.prop)<-1:b
# to:: names(new.prop)<-1:bootstraps

#1/19/18--made edits to notes. Add start and end time calculation.

# 1/23/18--fix bug in cases where only 0 or 1 GSI groups (rowSums function needs at least 2 columns)
# edit Loop 3.1 and add Loops 3.1.a and 4.1 to address this case.

# 2/6/18.  add Basin summary for GSI only and PBT only ################

# New features
# (1) add Stock and Basin to input data file (assignment.table)
# (2) Fix point estimates and Bootstrap iterations to prevent negative proportions
# (3) add adjustment loop in re-allocation loop to subtract all of the expanded PBT 
# amount from GSI groups that were sampled if the PBT amount is < rowSum of all GSI groups 
# in the group.freq (freq.table for the bootstrap iteration).  If expanded PBT >= rowSum 
# of all GSI groups--then set all GSI groups to zero and use expanded total to calculate
# harvest proportions. This is done in Loops 2.5 (for point estimates) and 3.7 (for 
# bootstrap iterations).
# (4) Summary of results writen to csv files for groups, stocks, rear (H = PBT; U = GSI),
# and basin.  Summary files have the actual sample size, expanded PBT/adjusted GSI sample
# size, harvest proportion, lower and upper CI (default = 90% but user can change), and 
# the mean harvest proportion, standard deviation, and CV of the bootstrap iterations.
############################################################################################

# Version 3.2 released 2/14/20 with change to reallocation routine

# 2/14/20. changed how reallocation table is dervied. Now only the GSI groups that were
# sampled (and the proportion assigning to PBT stocks) are used in re-allocation loop.

# changed "hatchery_fish_release_location" to "reallocate_stock" for reallocation loop
# and in the csv assignment and reallocation input files

###########################################################################################

# Remove current objects that may be in workspace

rm(list=ls())      
start.time<-Sys.time()

library(data.table)   # used in program. created in older R version.
 library(dplyr)       # used in program. created in older R version.

###################################### 
# SET WORKING DIRECTORY 
######################################

# Set the working file directory. The 3 input files listed below should
# be placed in this directory.  Output files (specified at end of program) 
# will be written to this directory.

setwd("C:/Users/........................../")

getwd() # Make sure its a valid directory.

################################################################## 
# USER-DEFINED INPUT FILES- enter filename after <- read.csv(??.
# enclose filename with ???????..xxx? where xxx = file extension
# recommend using csv format available in Excel. 3 files required.
# example files provided
###################################################################

# 1) input data file of group assignments [PBT & GSI], rear, stock, basin
# 2) PBT hatchery/brood year tagging rates,
# 3) reallocation table

assignment.table  <- read.csv("LCR Sport All 2019.csv",header=TRUE)  

tag.rate.table    <- read.csv("Tag_Rate_6-4-20.csv", header=TRUE)

reallocation.table<- read.csv("Reallocation Table_3-12-20.csv", header=TRUE)

## replace missing values with zero in reallocation.table

reallocation.table[is.na(reallocation.table)] = 0

## for naming the output files of Results create text string next

 output<-"LCR Sport 2019 All.csv"

#########################################################################
# Set number of bootstrap resamples and CI.  Default is 10000 iterations
# and 90% CI.  Enter CI as a proportion (i.e 90% = 0.9)
##########################################################################

bootstraps <- 10000
ci         <- .90
 
######################################################################
# GENERATE HARVEST PROPORTION ESTIMATES PRIOR TO CI BOOTSTRAP ROUTINE #
#######################################################################

# Tabulate actual sample count of PBT & GSI group

group.freq <- as.data.frame(matrix(table(assignment.table$group_assignment),1,length(unique(assignment.table$group_assignment))))

# Add headers (group names) to the data frame

names(group.freq) <- names(table(assignment.table$group_assignment))
group.freq

# Create a list of all 'unique' group names in the group_assignment column of # the assignment.table

group.names <- unique(assignment.table$group_assignment) 
group.names

# Generate lists of hatchery and 'wild' groups using the 'rear' column in
# the assignment.table

hatchery.list <- as.character(unique(assignment.table$group_assignment[which(assignment.table$rear=="H")]))

wild.list <- as.character(unique(assignment.table$group_assignment[which(assignment.table$rear=="U")]))

### choose GSI groups that were sampled for reallocation table

reallocation.table <- reallocation.table[, colnames(reallocation.table) %in% c("reallocate_stock", wild.list)]
for(i in 1:nrow(reallocation.table))
{
  if(sum(reallocation.table[i,2:ncol(reallocation.table)]) == 0) reallocation.table[i,2:ncol(reallocation.table)] <- 1
  reallocation.table[i,2:ncol(reallocation.table)] <- reallocation.table[i,2:ncol(reallocation.table)] / 
                                                        sum(reallocation.table[i,2:ncol(reallocation.table)])
}

## create a list of distinct hatchery stocks

distinct.hatchery.stock.names <- as.character(unique(assignment.table$stock_assignment[which(assignment.table$rear=="H")]))

# Determine sample size in the assignment.table. The sample size is used to 
# set the number of individuals that will be re-sampled during
# each iteration in the bootstrapping routine below.

sample.size <- length(assignment.table$individual)  ## all samples
sample.size

wild.sample.size<-length(assignment.table$individual[which((assignment.table$group_assignment) %in% wild.list)])
wild.sample.size    ## GSI samples only


#### BEGIN PBT EXPANSION using Tag rates  Loop 1  ########################

exp.point.ests <- group.freq
for (hatchery.by in hatchery.list) 
{ #1.1
  if (hatchery.by %in% tag.rate.table[,1])
  {  #1.2
exp.point.ests[1,hatchery.by] <- group.freq[hatchery.by] / 
tag.rate.table[as.character(tag.rate.table[,1])==hatchery.by,2]
  }  #1.2  
}  #1.1

### END PBT EXPANSION LOOP 1 ############################################

# Calculate the number of fish that were added to each PBT stock (rear = H).  This is 
# the amount that we must subtract from the GSI stocks (rear = U) 

exp.point.ests.diff <- exp.point.ests[,hatchery.list] - group.freq[,hatchery.list]

## if the expanded total - sample.size > GSI sample size, set all GSI groups to 0
## proportions of hatchery groups will be expanded / sum(expanded). All GSI proportions = 0.

if(rowSums(exp.point.ests.diff) >= wild.sample.size)  # set all GSI groups to 0
  {
    exp.point.ests[,wild.list]<- 0
  }

### BEGIN GSI -> PBT REALLOCATION LOOP 2 ONLY if exp.point.ests < wild.sample.size ####

# subtract the expanded - actual PBT amount from appropriate GSI groups
# using reallocation table that allocates the Expand-Actual PBT difference
# to the GSI groups.  This amount is the amount that would have been assigned 
# to the GSI group if there was no PBT assignment.  If the GSI group is < 0 
# after the subtraction, then set GSI group to 0 to avoid negative  
# proportions. If there were no actual assignments to the GSI group it will 
# be 0%.  exp.point.ests matrix will have the counts for each PBT/GSI group
# after the reallocation,  These counts will be used to calculate the 
# proportions of all groups.  Only GSI groups that were sampled will appear
# after the re-allocation.


if(rowSums(exp.point.ests.diff) < wild.sample.size)  ## do re-allocation loop
{  #2.1

for (hatchery.by in hatchery.list)
{ #2.2

hatchery.release <- as.character(unique(assignment.table$reallocate_stock
[which(assignment.table$group_assignment==hatchery.by)]))

for (rg in wild.list)
{ #2.3
 
  if(any(names(reallocation.table)==rg))
  
 { #2.4
exp.point.ests[,rg]  <-ifelse(exp.point.ests[,rg] - 
(exp.point.ests.diff[1,hatchery.by]*reallocation.table[which(reallocation.table$reallocate_stock==hatchery.release),rg])<=0,0,exp.point.ests[,rg] - 
(exp.point.ests.diff[1,hatchery.by]*reallocation.table[which(reallocation.table$reallocate_stock==hatchery.release),rg]))
       } #2.4
     }    #2.3
   }       #2.2


 while (rowSums(exp.point.ests) > sample.size)  
  { #2.5

for (rg in wild.list)
{  #2.6

aa = (rowSums(exp.point.ests) - rowSums(group.freq))/length(exp.point.ests[which(!exp.point.ests==0 & colnames(exp.point.ests) %in% wild.list)])

exp.point.ests[,rg]<-ifelse((exp.point.ests[,rg] - aa)<=0,0,exp.point.ests[,rg] - aa)
 } #2.6
  }  # 2.5
    } # 2.1

# Here are our point estimates (rounded to 4 decimal places)

point.estimates <- exp.point.ests/rowSums(exp.point.ests)
round(point.estimates,4)                                       

#######    Bootstrap   LOOP 3        ##############################################

## BEGIN BOOTSTRAPPING LOOP TO CALCULATE CIs AND CVS FOR EACH group ######

# Create empty data frame to store frequencies from each bootstrap iteration 
# of group assignment resampling

freq.table <- as.data.frame(matrix(rep(NA,(bootstraps*length(group.names))),bootstraps,length(group.names)))

names(freq.table) <- group.names

#### An empty table s groups wide by b bootstraps long ###

head(freq.table)

#### BEGIN BOOTSTRAP LOOP -- this may take a few minutes depending on sample size and the number of bootstraps

for (b in 1:bootstraps) 
{
resample <- sample(x = as.character(assignment.table$group_assignment),size = sample.size, replace =TRUE)

freq.table[b,] <- lapply(as.character(group.names), function(x) length(resample[resample==x]))
}

# Now the freq.table is filled out with counts from each bootstrap

head(freq.table) 

# Now the PBT group counts within each bootstrap iteration must be expanded 
# by its PBT tagging rates

# Copy the freq.table to a new table called expansion.table

expansion.table <- freq.table
head(expansion.table)

### BEGIN EXPANSION LOOP ACROSS ALL ITERATIONS IN the expansion.table ####

for (hatchery.by in hatchery.list)
{
  if (hatchery.by %in% tag.rate.table[,1])
    {
expansion.table[,hatchery.by] <- as.numeric(freq.table[,hatchery.by]) / tag.rate.table[as.character(tag.rate.table[,1])==hatchery.by,2]
} 
else 
 {
print(paste(hatchery.by, "has no tagging rate defined in tag_rates.csv"))
  }
    }
#### END EXPANSION LOOP ###

### ADJUST the GSI group COUNTS####

# These are the numbers that we added to each hatchery group in each 
# bootstrap iteration

exp.diff.table <- expansion.table[,hatchery.list] - freq.table[,hatchery.list]
# count.present<-sum(freq.table[,wild.list]==0)

## if the expanded total - sample.size > GSI sample size, set all GSI groups to 0
## proportions of hatchery groups will be expanded / sum(expanded). All GSI proportions = 0.

if (length(wild.list) > 1)  #Need at least 2 GSI groups present for rowSums() function
   { 3.1                    #and to reallocate to multiple groups in loops 3.x 

for (b in 1:bootstraps) # check rowSums of each iteration. 
 {  #3.1a                 

 
if(rowSums(exp.diff.table[b,]) >= rowSums(freq.table[b,wild.list])) # set all GSI groups to 0
  {  #3.2
    expansion.table[b,wild.list]<- 0
  }  #3.2
    

#### REALLOCATION LOOP for the iteration in expansion.table #####
####    if the rowSums(exp.diff.table) < wild.sample.size) ##### 

# subtract the expanded - actual PBT amount from appropriate GSI groups
# using reallocation table to allocate the Expand-Actual PBT difference
# from the GSI groups.  This amount is the amount that would have been 
# assigned to each GSI group if there was no PBT assignment.  If the GSI
# group is < 0 after the subtraction, then set GSI group to 0 to avoid 
# negative proportions.  GSI groups that were not in the actual sample
# will not appear in the expansion.table.
###############################################################  


if(rowSums(exp.diff.table[b,]) < rowSums(freq.table[b,wild.list]))
 {   #3.3

for (hatchery.by in hatchery.list)
{ #3.4
hatchery.release <- as.character(unique(assignment.table$reallocate_stock[which(assignment.table$group_assignment==hatchery.by)]))
 
for (rg in wild.list)
  {  #3.5
  
  if(any(names(reallocation.table)==rg))
    {  #3.6

expansion.table[b,rg] <-ifelse(expansion.table[b,rg] -
(exp.diff.table[b,hatchery.by]*reallocation.table[which(reallocation.table$reallocate_stock==hatchery.release),rg])<=0,
0,expansion.table[b,rg] - (exp.diff.table[,hatchery.by]*reallocation.table[which(reallocation.table$reallocate_stock==hatchery.release),rg]))

}  #3.6
  }  #3.5
   }  #3.4
 

####  Adjustment loop for reallocation of GSI groups not sampled.  Loop 3a ********** 

## Since an itereation many not have had all GSI groups that had a "re-allocation value" to subtract the rowSums of this iteration
## will be > the sample size.  This loop subtracts this amount from the GSI groups present.  Once completed the sum of all PBT and 
## GSI groups in the expansion.table will = the sample size for the iteration.  This loop is only excuted when wild.sample.size 
## < the (PBT expanded - PBT actual count). If (PBT expanded - PBT actual count) >= wild.sample.size this loop (and Loop 3.3) 
## is not run because  all GSI groups were set to zero in Loop 3.2.


 while (rowSums(expansion.table[b,]) > sample.size)
  {  # 3.7

bb=(rowSums(expansion.table[b,])-sample.size)/length(expansion.table[which(!expansion.table[b,]== 0 & colnames(expansion.table[b,]) %in% wild.list)])

for (rg in wild.list)
{ #3.8
expansion.table[b,rg] <- ifelse(expansion.table[b,rg]-bb <=0,0,expansion.table[b,rg]-bb)
} #3.8
  }  #3.7
    } #3.3 
   }  # 3.1.a   and go to the next iteration
 } #3.1  # End adjustment loop 3


## Now deal with datasets that have 1 GSI groups in the sample.  If no
## GSI groups were present the reallocation loop 3 will not execute 
## and harvest proportions will use expanded total not sample size in
## the denominator.

 if (length(wild.list) == 1)  # subtract all from this group
 { #4.1

exptot <-data.frame(rowSums(exp.diff.table))
  expansion.table[,wild.list]<-ifelse(expansion.table[,wild.list]-exptot[,1]<=0,0,expansion.table[,wild.list]-exptot[,1])

 } # 4.1


# Create a table of group proportions for each iteration. Every value in each
# column in Row i is divided by the rowSums of Row i.  Each rowSums = sample 
# size if (PBT expansion - PBT actual)  < GSI actual.  if (PBT expansion count - PBT actual)
#  >= GSI actual then all GSI groups = 0 and PBT group proportions found by diviving
# PBT expand count by expanded rowSums (which will be > sample.size)

prop.table <- expansion.table / rowSums(expansion.table)

head(prop.table) 

####Generate summary statistics from the 'prop.table'################

# Calculate the Lower CI, Upper CI, Mean of iterations, standard deviation of 
# iterations, and coefficient of variation

results.boot    <- data.frame(group_assignment=names(prop.table), lci=apply(prop.table, 2, quantile, (1-ci)/2), uci=apply(prop.table, 2, quantile, ((1-ci)/2)+ci), mean=apply(prop.table, 2, mean), st.dev=apply(prop.table, 2, sd))
results.boot$cv <-(results.boot$st.dev/results.boot$mean)*100

### making a dataframe from the transposed point estimates, actual count, 
### and adjusted count used for point estimate

results.point <- data.frame(Proportion=t(point.estimates),Actual=t(group.freq),Adjusted=t(exp.point.ests))      

# adding a reporting group column  
results.point$group_assignment <- row.names(results.point)           

## all stats by group_assignment
results.group<-merge(results.point, results.boot, by="group_assignment")  

### Summarize prop.table by hatchery stock, rear, basin by summing release
# groups (group_assignments)

all.groups.stocks<-subset(assignment.table,select = -c(individual,reallocate_stock))

distinct.groups<-unique(all.groups.stocks, by = "group_assignment")
 rownames(distinct.groups)<-distinct.groups$group_assignment

## bootstrap iterations by group_assignment
new.prop<-data.frame(t(prop.table))  
 names(new.prop)<-1:bootstraps

## new column with group_assignment names
new.prop$group_assignment<-rownames(new.prop) 

stock.rear.groups<-merge(x=new.prop, y = distinct.groups, by = "group_assignment", all.x = TRUE)

 stock.prop<-subset(stock.rear.groups,select = -c(group_assignment,rear,basin))
   
rear.prop<-subset(stock.rear.groups,select = -c(group_assignment,stock_assignment,basin))
    basin.prop<-subset(stock.rear.groups,select = -c(group_assignment,stock_assignment,rear))


##sum all stocks in each iteration- to get CIs by stock  ###

stock.ci=stock.prop %>% group_by(stock_assignment) %>% summarise_all(funs(sum))

### get lci and uci, mean, std of bootstrap iterations for each stock

## delete column with stock names, keep everything else the same as stock.ci

new<-subset(stock.ci, select = -c(stock_assignment))  

stock.lci.uci<-data.frame(stock=stock.ci$stock_assignment,lci=apply(new,1,quantile,(1-ci)/2),uci=apply(new,1,quantile, ((1-ci)/2)+ci),mean=apply(new,1,mean),st.dev=apply(new,1,sd))

stock.lci.uci$cv<-(stock.lci.uci$st.dev/stock.lci.uci$mean)*100


## by rear type - find CIs###

rear.ci=rear.prop %>% group_by(rear) %>% summarise_all(funs(sum))

## delete column with rear type, keep everything else the same as rear.ci
 new2<-subset(rear.ci, select = -c(rear))  

rear.lci.uci<-data.frame(rear=rear.ci$rear,lci=apply(new2,1,quantile,(1-ci)/2),uci=apply(new2,1,quantile, ((1-ci)/2)+ci),mean=apply(new2,1,mean),st.dev=apply(new2,1,sd))

rear.lci.uci$cv<-(rear.lci.uci$st.dev/rear.lci.uci$mean)*100


## by basin - find Cis using all samples###

basin.ci=basin.prop %>% group_by(basin) %>% summarise_all(funs(sum))

# delete column with basin, keep everything else the same as basin.ci
 newbasin<-subset(basin.ci, select = -c(basin)) 

basin.lci.uci<-data.frame(basin=basin.ci$basin,lci=apply(newbasin,1,quantile,(1-ci)/2),uci=apply(newbasin,1,quantile, ((1-ci)/2)+ci),mean=apply(newbasin,1,mean),st.dev=apply(newbasin,1,sd))

basin.lci.uci$cv<-(basin.lci.uci$st.dev/basin.lci.uci$mean)*100



## Get point estimates by stock, rear type, and basin (all samples)

new.point<-data.frame(t(point.estimates))    

# new column with group_assignment names
 new.point$group_assignment<-rownames(new.point)    

point.stock.rear.groups<-merge(x=new.point, y = distinct.groups, by = "group_assignment", all.x = TRUE)

point.stock<-subset(point.stock.rear.groups,select = -c(group_assignment,rear,basin))

point.rear<-subset(point.stock.rear.groups,select = -c(stock_assignment,group_assignment,basin))

point.basin<-subset(point.stock.rear.groups,select = -c(stock_assignment,group_assignment,rear))


stock.point.estimate=point.stock %>% group_by(stock_assignment) %>% summarise_all(funs(sum))

rear.point.estimate=point.rear %>% group_by(rear) %>% summarise_all(funs(sum))

basin.point.estimate=point.basin %>% group_by(basin) %>% summarise_all(funs(sum))


### get actual and adjusted sample size by stock #####

new.group.freq<-data.frame(t(group.freq))
  new.group.freq$group_assignment<-rownames(new.group.freq)
   new.group<-merge(x=new.group.freq, y = distinct.groups, by = "group_assignment", all.x = TRUE)
    new3<-subset(new.group,select = -c(group_assignment,rear,basin))

stock.freq=new3 %>% group_by(stock_assignment) %>% summarise_all(funs(sum))


new.exp.point.ests<-data.frame(t(exp.point.ests))
  new.exp.point.ests$group_assignment<-rownames(new.exp.point.ests)
   new.ests<-merge(x=new.exp.point.ests, y = distinct.groups, by = "group_assignment", all.x = TRUE)
     new4<-subset(new.ests,select = -c(group_assignment,rear,basin))

stock.adjust=new4 %>% group_by(stock_assignment) %>% summarise_all(funs(sum))


### get actual and adjusted sample size by rear and basin #####
## by rear actual##
new5<-subset(new.group,select = -c(group_assignment,stock_assignment,basin))  

## by basin actual##
 new6<-subset(new.group,select = -c(group_assignment,stock_assignment,rear))  

rear.freq=new5 %>% group_by(rear) %>% summarise_all(funs(sum))

basin.freq=new6 %>% group_by(basin) %>% summarise_all(funs(sum))

## by rear adjusted##
new7<-subset(new.ests,select = -c(group_assignment,stock_assignment,basin))   
 
## by basin adjusted##
new8<-subset(new.ests,select = -c(group_assignment,stock_assignment,rear))  

rear.adjust=new7 %>% group_by(rear) %>% summarise_all(funs(sum))

basin.adjust=new8 %>% group_by(basin) %>% summarise_all(funs(sum))


######  write results to data frame ####

## by stock ###

res1<-merge(stock.point.estimate,stock.freq)
 res2<-merge(res1,stock.adjust)
   setnames(res2, old = c(1,2,3,4), new = c('stock','proportion','actual','adjusted'))

 ## all results by stock (actual,adjust,proportion,lci,uci,mean,std,cv)   
 
all.results.stock<- merge(res2,stock.lci.uci)  ## all results by stock (actual,adjust,proportion,lci,uci,mean,std,cv)

## by rear ###

res3<-merge(rear.point.estimate,rear.freq)
 res4<-merge(res3,rear.adjust)
   setnames(res4, old = c(1,2,3,4), new = c('rear','proportion','actual','adjusted'))
    
all.results.rear<- merge(res4,rear.lci.uci)  ## all results by rear (actual,adjust,proportion,lci,uci,mean,std,cv)


## by basin using all samples###

res5<-merge(basin.point.estimate,basin.freq)
 res6<-merge(res5,basin.adjust)
   setnames(res6, old = c(1,2,3,4), new = c('basin','proportion','actual','adjusted'))
    
all.results.basin<- merge(res6,basin.lci.uci) 


## calculate CI's, bootstrap mean, st dev, cv for stocks identified with
## PBT (rear="H") by basin. These basin results exclude GSI assigned fish.

pbt.only<-subset(stock.rear.groups,rear=="H",select = -c(group_assignment,stock_assignment,rear))

pbt.only.basin = pbt.only %>% group_by(basin) %>% summarise_all(funs(sum))

pbt0<-subset(pbt.only.basin, select = -c(basin))

pbt.basin.lci.uci<- data.frame(basin=pbt.only.basin$basin,lci=apply(pbt0,1,quantile,(1-ci)/2),uci=apply(pbt0,1,quantile,((1-ci)/2)+ci), mean=apply(pbt0,1,mean),st.dev=apply(pbt0,1,sd))

pbt.basin.lci.uci$cv<-(pbt.basin.lci.uci$st.dev/pbt.basin.lci.uci$mean)*100

## point estimates for PBT only by basin
pbt1<-subset(point.stock.rear.groups,rear=="H",select = -c(group_assignment,stock_assignment,rear))
 pbt.only.point = pbt1 %>% group_by(basin) %>% summarise_all(funs(sum))

## actual sample size for PBT only by basin
pbt2<-subset(new.group,rear=="H",select= -c(group_assignment,stock_assignment,rear))
 pbt.only.freq = pbt2 %>% group_by(basin) %>% summarise_all(funs(sum))  

## adjusted sample size for PBT only by basin
pbt3<-subset(new.ests,rear=="H",select = -c(group_assignment,stock_assignment,rear))
 pbt.only.adjust = pbt3 %>% group_by(basin) %>% summarise_all(funs(sum))

# Write PBT only by basin results to output file

pbt.res1<-merge(pbt.only.point,pbt.only.freq)
 pbt.res2<-merge(pbt.res1,pbt.only.adjust)

 setnames(pbt.res2, old= c(1,2,3,4), new = c('basin','proportion','actual','adjusted'))

 all.results.pbt.only.by.basin<-merge(pbt.res2,pbt.basin.lci.uci)


## calculate CI's, bootstrap mean, st dev, cv for stocks identified with
# GSI (rear="U") by basin. These basin results exclude PBT assigned fish 

gsi.only<-subset(stock.rear.groups,rear=="U",select = -c(group_assignment,stock_assignment,rear))

gsi.only.basin = gsi.only %>% group_by(basin) %>% summarise_all(funs(sum))

gsi0<-subset(gsi.only.basin, select = -c(basin))

gsi.basin.lci.uci<- data.frame(basin=gsi.only.basin$basin,lci=apply(gsi0,1,quantile,(1-ci)/2),uci=apply(gsi0,1,quantile,((1-ci)/2)+ci), mean=apply(gsi0,1,mean),st.dev=apply(gsi0,1,sd))

gsi.basin.lci.uci$cv<-(gsi.basin.lci.uci$st.dev/gsi.basin.lci.uci$mean)*100


## point estimates for GSI only by basin
gsi1<-subset(point.stock.rear.groups,rear=="U",select = -c(group_assignment,stock_assignment,rear))
 gsi.only.point = gsi1 %>% group_by(basin) %>% summarise_all(funs(sum))

## actual sample size for GSI only by basin
gsi2<-subset(new.group,rear=="U",select= -c(group_assignment,stock_assignment,rear))
 gsi.only.freq = gsi2 %>% group_by(basin) %>% summarise_all(funs(sum))  

## adjusted sample size for GSI only by basin
gsi3<-subset(new.ests,rear=="U",select = -c(group_assignment,stock_assignment,rear))
 gsi.only.adjust = gsi3 %>% group_by(basin) %>% summarise_all(funs(sum))

# Write GSI only by basin results to output file

gsi.res1<-merge(gsi.only.point,gsi.only.freq)
 gsi.res2<-merge(gsi.res1,gsi.only.adjust)

 setnames(gsi.res2, old= c(1,2,3,4), new = c('basin','proportion','actual','adjusted'))

 all.results.gsi.only.by.basin<-merge(gsi.res2,gsi.basin.lci.uci)



# Write the table of proportions of each bootstrap iteration to a .csv file #

# if further analysis/summaries are needed. The prop.table was used to get
# Cis for group_assignments. Stock, rear, and basin CIs were found by summing # the hatchery release groups (for example all Dworshak groups
#  for stock) in each iteration.


######  write results to csv files ###

### summary stats with point estimate, actual group sample size, adjusted 
# group sample size used for the point estimates, Lower CI, Upper CI, Mean,
# CV, Std Dev. To name files program uses the fishery type (output) that was
# #set at beginning of program
###########################################################################

write.csv(prop.table, file = paste("Resampled group proportions_",output,sep="_"))

write.csv(results.group, file = paste("Results by Group",output,sep="_"))

write.csv(all.results.stock, file = paste("Results by Stock_",output,sep="_"))

write.csv(all.results.rear, file = paste("Results by Rear_",output,sep="_"))

write.csv(all.results.basin, file = paste("Results by Basin All Samples_",output,sep="_"))

write.csv(all.results.pbt.only.by.basin, file = paste("Results by Basin PBT only",output,sep="_"))

write.csv(all.results.gsi.only.by.basin, file = paste("Results by Basin GSI only",output,sep="_"))


end.time<-Sys.time()
 runtime<-end.time-start.time
   runtime

######### END OF PROGRAM  #####################