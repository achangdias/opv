######## setup ####
rm(list=ls())

setwd("~/Documents/04 OPV/U5MR")

library(foreign)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(ggthemes)
options("scipen" = 10)
library(ggrepel)
library(RColorBrewer)
library(viridis)


#### parameters #### 
pop <- 1000000

m_year1 <- 50/1000
m_year2 <- 7/1000
m_year3 <- 7/1000

year <- c("year1","year2","year3")
mort <- c(m_year1, m_year2, m_year3)

vx_eff_marginal <- 1-0.92
vx_eff_marginal_l <- vx_eff_marginal * (1+0.5)
vx_eff_marginal_u <- 0.08 * (1-0.5) 

vx_eff <- 1-0.90
vx_eff_l <- 1-0.94
vx_eff_u <- 1-0.49

campaign_num <- 3 

#  https://www.usinflationcalculator.com/inflation/consumer-price-index-and-annual-percent-changes-from-1913-to-2008/
cost_opv <- 0.15 # * 	258.811 / 255.657 # 2019 to 2020 USD, 0.1518 
cost_opv_l <- cost_opv * (1-0.25)
cost_opv_u <- cost_opv * (1+0.25)

cost_delivery <- 1.03 # (1.56+0.72+0.76+ 0.77+1.43+0.96)/6 # average of six
  # 1.14 * 258.811 / 188.9 # 1.56 2004 Ghana Zambia Nepal 
  # 0.56 * 258.811 / 201.6 # 0.72 Ethiopia vitamin A 
  # 3187/5000 * 258.811 / 218.056 # Zambia vitA 2010 
  # 0.65 * 258.11 / 218.056 # Portnoy 2015 2010 
  # 0.95 * 258.811/255.657 THOMPSON 2019  

  #  (1/2) * 2.09 *  258.811 / 188.9 # same two dose, CPI from 2004 Cavallier
cost_delivery_u <- 0.72
cost_delivery_l <- 1.56


# 2.46 * 258.811 / 255.657 = 2.49 # CPI 2019 Zimmerman NOT RELEVANT 
# previously (1/2) * 9.7 * 258.811 / 214.537 # OCV but this is for 2 doses so 1/2 ; Schaetti 2012 CPI from 2009 to 2020 

wastage <- 0.1
wastage_l <- 0.01
wastage_u <- 0.25 # OCV 1%  


#### effectiveness ####
data_effect <- data.table(year, mort)
data_effect[year == "year1", death_novx := pop * mort]
data_effect[year == "year2", death_novx := (pop - data_effect[year == "year1",.(death_novx)]) * mort]
data_effect[year == "year3", death_novx := (pop - data_effect[year == "year1",.(death_novx)] - data_effect[year == "year2",.(death_novx)]) * mort]
data_effect[year == "year1", death_vx := pop * mort * (1-vx_eff)]
data_effect[year == "year2", death_vx := (pop - data_effect[year == "year1",.(death_vx)]) * mort * (1-vx_eff_marginal)]
data_effect[year == "year3", death_vx := (pop - data_effect[year == "year1",.(death_vx)]  - data_effect[year == "year2",.(death_vx)]) * mort * (1-vx_eff_marginal)]

# slightly different  
# data_effect[year == "year1", death_novx * (1-vx_eff)] 
# data_effect[year == "year2", death_novx * (1-vx_eff_marginal)] 
# data_effect[year == "year3", death_novx * (1-vx_eff_marginal)] 

data_effect[, death_avt := death_novx - death_vx]
data_effect[, sum(death_avt)] # 5993
data_effect[, sum(death_avt)] / data_effect[, sum(death_novx)] # 9.4% 

nnv <- pop / data_effect[, sum(death_avt)] # 167

data_effect[, sum(death_novx)] *1000 / pop # U3MR = 63 per 1000 live births 

# lower bound 
data_effect_l <- data.table(year, mort)
data_effect_l[year == "year1", death_novx := pop * mort]
data_effect_l[year == "year2", death_novx := (pop - data_effect_l[year == "year1",.(death_novx)]) * mort]
data_effect_l[year == "year3", death_novx := (pop - data_effect_l[year == "year1",.(death_novx)] - data_effect_l[year == "year2",.(death_novx)]) * mort]
data_effect_l[year == "year1", death_vx := pop * mort * (1-vx_eff_l)]
data_effect_l[year == "year2", death_vx := (pop - data_effect_l[year == "year1",.(death_vx)]) * mort * (1-vx_eff_marginal_l)]
data_effect_l[year == "year3", death_vx := (pop - data_effect_l[year == "year1",.(death_vx)]  - data_effect_l[year == "year2",.(death_vx)]) * mort * (1-vx_eff_marginal_l)]
data_effect_l[, death_avt := death_novx - death_vx]
data_effect_l[, sum(death_avt)] 
data_effect_l[, sum(death_avt)] / data_effect_l[, sum(death_novx)] 
nnv_l <- pop / data_effect_l[, sum(death_avt)] 

# upper bound 
data_effect_u <- data.table(year, mort)
data_effect_u[year == "year1", death_novx := pop * mort]
data_effect_u[year == "year2", death_novx := (pop - data_effect_u[year == "year1",.(death_novx)]) * mort]
data_effect_u[year == "year3", death_novx := (pop - data_effect_u[year == "year1",.(death_novx)] - data_effect_u[year == "year2",.(death_novx)]) * mort]
data_effect_u[year == "year1", death_vx := pop * mort * (1-vx_eff_u)]
data_effect_u[year == "year2", death_vx := (pop - data_effect_u[year == "year1",.(death_vx)]) * mort * (1-vx_eff_marginal_u)]
data_effect_u[year == "year3", death_vx := (pop - data_effect_u[year == "year1",.(death_vx)]  - data_effect_u[year == "year2",.(death_vx)]) * mort * (1-vx_eff_marginal_u)]
data_effect_u[, death_avt := death_novx - death_vx]
data_effect_u[, sum(death_avt)] # 25687
data_effect_u[, sum(death_avt)] / data_effect_u[, sum(death_novx)] # 40.1% 
nnv_u <- pop / data_effect_u[, sum(death_avt)] # 39

# format 
results <- cbind(data_effect,data_effect_l[,.(death_vx, death_avt)],data_effect_u[,.( death_vx, death_avt)])
names(results) <- c("year","mort","death_novx","death_vx","death_avt","death_vx_l","death_avt_l","death_vx_u","death_avt_u")
results[year == "year1", col1 := paste(round(death_vx,0), " [",round(death_vx_u,0),",",round(death_vx_l,0),"]", sep="")]
results[year != "year1", col1 := paste(round(death_vx,0), " [",round(death_vx_l,0),",",round(death_vx_u,0),"]", sep="")]

results[year == "year1", col2 := paste(round(death_avt,0), " [",round(death_avt_l,0),",",round(death_avt_u,0),"]", sep="")]
results[year != "year1", col2 := paste(round(death_avt,0), " [",round(death_avt_u,0),",",round(death_avt_l,0),"]", sep="")]
results
write.csv(results, "results_table.csv")

# sum of averted deaths 
c(results[, round(sum(death_avt),0)], results[, round(sum(death_avt_l),0)],results[, round(sum(death_avt_u),0)])

# estimated campaign effectiveness 
c(results[, round(sum(death_avt),0)] / results[, round(sum(death_novx),0)], 
  results[, round(sum(death_avt_l),0)] / results[, round(sum(death_novx),0)],
  results[, round(sum(death_avt_u),0)] / results[, round(sum(death_novx),0)])

# nnv 
c(pop / data_effect[, sum(death_avt)], pop / data_effect_u[, sum(death_avt)] , pop / data_effect_l[, sum(death_avt)] )

# DATI 
c(data_effect[, sum(death_avt)] *1000/ (pop*3), data_effect_l[, sum(death_avt)] *1000/ (pop*3), data_effect_u[, sum(death_avt)] *1000/ (pop*3))

#### cost #### 
cost <- pop * (cost_opv + cost_delivery) * (1+wastage) * campaign_num 
cost_l <- pop * (cost_opv_l + cost_delivery_l) * (1+wastage_l) * campaign_num 
cost_u <- pop * (cost_opv_u + cost_delivery_u) * (1+wastage_u) * campaign_num 
c(cost, cost_l, cost_u)

#### VSL #### 
gni <- 2230 # 2019 PPP 
gni_us <- 57900 
vsl_floor <- gni * 20 
  
vsl_1 <- 9400000 * (gni/gni_us) ^1.5  # 71050 
  # gni * 160 * (gni/gni_us) ^1.5 # 2211, lower than floor 
vsl_2 <- 100 * gni # 223000
vsl_3 <- 160 * gni # 357000

LE_ratio <- 57.82 / 32.39 # LE at age 0/35  # 1.79
vsl_ratio = vsl_1 * LE_ratio 

vsl_half = vsl_1 * 0.5 

vsl_main <- vsl_1 

#### CEA and BCA ####
# avt death 
c(results[, round(sum(death_avt),0)], results[, round(sum(death_avt_l),0)],results[, round(sum(death_avt_u),0)])

# estimated campaign effectiveness 
c(round(results[, sum(death_avt)] / results[, sum(death_novx)],3) *100, 
  round(results[, sum(death_avt_l)] / results[, sum(death_novx)],3) *100, 
  round(results[, sum(death_avt_u)] / results[, sum(death_novx)],3) *100)

# deaths averted per 1000 vaccinated 
round(c(results[, round(sum(death_avt),0)]*1000/ pop, results[, round(sum(death_avt_l),0)] *1000 / pop,results[, round(sum(death_avt_u),0)]*1000 / pop),1)

# nnv 
round(c(pop / data_effect[, sum(death_avt)], pop / data_effect_u[, sum(death_avt)] , pop / data_effect_l[, sum(death_avt)] ),1)
# nnv with three campaigns 
round(c(pop * campaign_num / data_effect[, sum(death_avt)], pop * campaign_num / data_effect_u[, sum(death_avt)] , pop * campaign_num  / data_effect_l[, sum(death_avt)] ),1)



# total cost 
c(cost, cost_l, cost_u)
c(cost, cost_l, cost_u)/1000000

# C/E 
round(c(cost/results[, sum(death_avt)],
  cost_l/results[, sum(death_avt_u)],
  cost_u/results[, sum(death_avt_l)]),0) # mixed lower and upper ranges 

# B/C
# vsl_main
bca1 <- vsl_main * data_effect_l[, sum(death_avt)] / cost_l
bca2 <- vsl_main * data_effect_u[, sum(death_avt)] / cost_u
bca3 <- vsl_main * data_effect_u[, sum(death_avt)] / cost_l # mixed l and u
bca4 <- vsl_main * data_effect_l[, sum(death_avt)] / cost_u # mixed l and u
round(c(vsl_main * data_effect[, sum(death_avt)] / cost, min(bca1,bca2,bca3,bca4), max(bca1,bca2,bca3,bca4)),1)

# vsl_ratio 
bca1 <- vsl_ratio * data_effect_l[, sum(death_avt)] / cost_l
bca2 <- vsl_ratio * data_effect_u[, sum(death_avt)] / cost_u
bca3 <- vsl_ratio * data_effect_u[, sum(death_avt)] / cost_l # mixed l and u
bca4 <- vsl_ratio * data_effect_l[, sum(death_avt)] / cost_u # mixed l and u
round(c(vsl_ratio * data_effect[, sum(death_avt)] / cost, min(bca1,bca2,bca3,bca4), max(bca1,bca2,bca3,bca4)),1)

# vsl_1
bca1 <- vsl_1 * data_effect_l[, sum(death_avt)] / cost_l
bca2 <- vsl_1 * data_effect_u[, sum(death_avt)] / cost_u
bca3 <- vsl_1 * data_effect_u[, sum(death_avt)] / cost_l # mixed l and u
bca4 <- vsl_1 * data_effect_l[, sum(death_avt)] / cost_u # mixed l and u
round(c(vsl_1 * data_effect[, sum(death_avt)] / cost, min(bca1,bca2,bca3,bca4), max(bca1,bca2,bca3,bca4)),1)

# vsl_2
bca1 <- vsl_2 * data_effect_l[, sum(death_avt)] / cost_l
bca2 <- vsl_2 * data_effect_u[, sum(death_avt)] / cost_u
bca3 <- vsl_2 * data_effect_u[, sum(death_avt)] / cost_l # mixed l and u
bca4 <- vsl_2 * data_effect_l[, sum(death_avt)] / cost_u # mixed l and u
round(c(vsl_2 * data_effect[, sum(death_avt)] / cost, min(bca1,bca2,bca3,bca4), max(bca1,bca2,bca3,bca4)),1)

# vsl_3
bca1 <- vsl_3 * data_effect_l[, sum(death_avt)] / cost_l
bca2 <- vsl_3 * data_effect_u[, sum(death_avt)] / cost_u
bca3 <- vsl_3 * data_effect_u[, sum(death_avt)] / cost_l # mixed l and u
bca4 <- vsl_3 * data_effect_l[, sum(death_avt)] / cost_u # mixed l and u
round(c(vsl_3 * data_effect[, sum(death_avt)] / cost, min(bca1,bca2,bca3,bca4), max(bca1,bca2,bca3,bca4)),1)

#### DATI plots #### 
vx_eff_v <- seq(0.00, 0.5, 0.0001)
mort_v <- seq(10/1000, 120/1000, 1/5000) 
data2 <- data.table(expand.grid(vx_eff_v, mort_v))
names(data2) <- c("vx_eff", "mort")
data2[, nnv := campaign_num / (mort * vx_eff)]
data2[, nnv := round(nnv,0)]
data2[, NNV := as.factor(nnv)]
data2[, death_avt_k := mort * (1-1+vx_eff) * 1000 ]
data2[, death_avt_k_campaign := mort * (1-1+vx_eff) * 1000 / 3 ]

data2[, death_avt_k := round(death_avt_k, 2)]
data2[, death_avt_k_f := as.factor(death_avt_k)]
data2[, death_avt_k_campaign := round(death_avt_k_campaign, 2)]

data2[, death_avt_k_campaign_f := as.factor(death_avt_k_campaign)]
data2[, table(round(death_avt_k_campaign,2))]
select_avt <- c(0.1, 0.5 ,1, 1, 2, 3, 5, 7)
data_plot <- data2[death_avt_k_campaign%in% select_avt, ]
data_plot[ , `:=`(IDX = 1:.N ) , by = death_avt_k_campaign_f ]
data_plot <- data_plot[order(death_avt_k_campaign_f)]
for (i in select_avt) {
  dat <- data_plot[death_avt_k_campaign == i]
  count <- dat[,max(IDX)]
  data_plot[death_avt_k_campaign == i & IDX == count, label := death_avt_k_campaign_f]
}
data_plot[is.na(label), label:=""]

data_plot[mort == 0.0100 & vx_eff == 0.0735,]
data_plot[mort > 0.063 & mort < 0.065 & vx_eff > 0.095 & vx_eff < 0.097]

# 0.078 * (1-1+0.096) * 1000 / 3 = 2.496 
# 3 / (0.078*0.096) = 400  NNV 

jpeg("fig_dati_2021jun.jpg", width = 7, height = 7, units = 'in', res = 300)
p1 <- ggplot(data_plot, aes(x=mort * 1.2 * 1000, y=vx_eff  * 100, color=death_avt_k_campaign_f)) + # 
  # geom_point() + 
  geom_line(size=2) + 
  # scale_colour_brewer(palette = "Paired") +
  scale_colour_manual(values = viridis_pal(option = "D")(length(data_plot[, unique(death_avt_k_campaign_f)]))) + 
  geom_label_repel(aes(label = label),
                   nudge_x = 1, size = 4.5,
                   na.rm = TRUE, max.overlaps = Inf, label.size= NA) + 
  labs(colour = "death_avt_k_campaign_f", title = "A: Deaths averted per thousand campaign doses",
       x = "Underlying under-5 mortality rate (per 1000 live births)", 
       y = "Intervention effectiveness (%)",
       caption = "") +
  theme_classic() + theme(legend.position="none")  +
  # geom_text(aes(label=death_avt_k_campaign_f),hjust=0,vjust=0) +
  xlim(0.01 * 1000, 0.150  * 1000) +
  ylim(0.000001  * 100, 0.20 * 100) +
  geom_vline(xintercept = 0.079  * 1000, color = "gray", size = 1, linetype = "dashed") +
  geom_text(x=0.103 * 1000, y= 0.000001  * 100, label="Under-5 mortality = 80",size = 4.5,
            color = "gray", fontface="plain") +
  geom_hline(yintercept = 0.096 * 100, color = "gray", size = 1, linetype = "dashed") +
  geom_text(x=0.021 * 1000, y=0.089  * 100, label="Current estimate", size = 4.5,
            color = "gray", fontface="plain") +
  geom_text(x=0.092 * 1000, y= 0.089  * 100, label="DATI = 2.0", size = 4.5,
            color = "black", fontface="plain")
p1
dev.off()



data_plot[, table(round(death_avt_k,2))]
select_avt <- c(0.75, 1.5, 3, 6, 7.5, 9, 12, 15, 18, 21, 24)
data_plot2 <- data_plot[death_avt_k %in% select_avt, ]
data_plot2[ , `:=`(IDX = 1:.N ) , by = death_avt_k_f ]
data_plot2 <- data_plot2[order(death_avt_k_f)]
for (i in select_avt) {
  dat <- data_plot2[death_avt_k == i]
  count <- dat[,max(IDX)]
  data_plot2[death_avt_k == i & IDX == count, label := death_avt_k_f]
}
data_plot2[is.na(label), label:=""]


#### CEA and BCA plots ####  
vx_eff_v <- seq(0.00, 0.5, 0.001)
mort_v <- seq(10/1000, 120/1000, 1/5000) 
data3<- data.table(expand.grid(vx_eff_v, mort_v))
names(data3) <- c("vx_eff", "mort")
data3[, death_avt := vx_eff * mort]
data3[, cea := cost / (death_avt * pop)]
data3[, bca := death_avt * pop * vsl_main / cost] # using vsl_8 as main result (mid point)

data3[vx_eff == 0.1 & mort == 0.065] # check 

# CEA 
data3[, table(round(cea,0))]
data3[, cea_round := round(cea,0) ]
data3[cea_round > 3999 & cea_round < 4500, table(cea_round)]
data3[cea_round > 500 & cea_round < 750, table(cea_round)]
data3[cea_round > 8000 & cea_round < 12000, table(cea_round)]
data3[cea_round > 10000 & cea_round < 15000, table(cea_round)]

select_cea <- c(60,75,100,150,200,300,500,1000) # 2000, 4001
data_plot <- data3[cea_round %in% select_cea, ]
data_plot[ , `:=`(IDX = 1:.N ) , by = cea_round ]
data_plot <- data_plot[order(cea_round)]
for (i in select_cea) {
  dat <- data_plot[cea_round == i]
  count <- dat[,max(IDX)]
  data_plot[cea_round == i & IDX == count, label := round(cea_round,-1)]
}
data_plot[is.na(label), label:=""]
data_plot[, cea_round := as.factor(cea_round)]


jpeg("fig_cea_2021jun.jpg", width = 7, height = 7, units = 'in', res = 300)
p2 <- ggplot(data_plot, aes(x=mort * 1.2 * 1000, y=vx_eff  * 100, color=cea_round)) + 
  geom_line(size=2) +  
  # scale_colour_brewer(palette = "Paired") +
  scale_colour_manual(values = viridis_pal(option = "D")(length(data_plot[, unique(cea_round)]))) +
  geom_label_repel(aes(label = label),
                   nudge_x = 1, size = 4.5,
                   na.rm = TRUE, max.overlaps = Inf, label.size= NA) +
  labs(colour = "cea_round", title = "B: Cost-effectiveness ratio",
       x = "Underlying under-5 mortality rate (per 1000 live births)", 
       y = "Intervention effectiveness (%)",
       caption = "") +
  theme_classic() + theme(legend.position="none")  +
  # geom_text(aes(label=death_avt_k_campaign_f),hjust=0,vjust=0) +
  xlim(0.01 * 1000, 0.150  * 1000) +
  # ylim(0, 0.30 * 100) +
  geom_vline(xintercept = 0.079  * 1000, color = "gray", size = 1, linetype = "dashed") +
  geom_text(x=0.103 * 1000, y= 0.034  * 100, label="Under-5 mortality = 80", size = 4.5, 
            color = "darkgray", fontface="plain") +
  geom_hline(yintercept = 0.096 * 100, color = "gray", size = 1, linetype = "dashed") +
  geom_text(x=0.021 * 1000, y=0.08  * 100, label="Current estimate", size = 4.5, 
            color = "gray", fontface="plain") + 
  geom_text(x=0.092 * 1000, y= 0.08  * 100, label="C/E = $650", size = 4.5, 
            color = "black", fontface="plain")
p2
dev.off()
 

# BCA 
data3[, summary(round(bca,0))]
# data3[bca>600, table(round(bca,1))]
data3[, table(round(bca,0))]

data3[, bca_round := round(bca,1) ]

select_bca <- c(10,50,100,200,300,400,500)
data_plot <- data3[bca_round %in% select_bca, ]
data_plot[ , `:=`(IDX = 1:.N ) , by = bca_round ]
data_plot <- data_plot[order(bca_round)]
for (i in select_bca) {
  dat <- data_plot[bca_round == i]
  count <- dat[,max(IDX)]
  data_plot[bca_round == i & IDX == count, label := round(bca_round,-1)]
}
data_plot[is.na(label), label:=""]
data_plot[, bca_round := as.factor(bca_round)]

jpeg("fig_bca_2021jun.jpg", width = 7, height = 7, units = 'in', res = 300)
p3 <- ggplot(data_plot, aes(x=mort * 1.2 * 1000, y=vx_eff * 100, color=bca_round)) + 
  geom_line(size=2) +  
  # scale_colour_brewer(palette = "Paired") +
  scale_colour_manual(values = viridis_pal(option = "D")(length(data_plot[, unique(bca_round)]))) +
  geom_label_repel(aes(label = label),
                   nudge_x = 1, size = 4.5,
                   na.rm = TRUE, max.overlaps = Inf, label.size= NA) +
  labs(colour = "bca_round", title = "C: Benefit-cost ratio",
       x = "Underlying under-5 mortality rate (per 1000 live births)", 
       y = "Intervention effectiveness (%)",
       caption = "") +
  theme_classic() + theme(legend.position="none")  +
  # geom_text(aes(label=death_avt_k_campaign_f),hjust=0,vjust=0) +
  xlim(0.01 * 1000, 0.150  * 1000) +
  ylim(0.00001  * 100, 0.30 * 100) +
  geom_vline(xintercept = 0.079  * 1000, color = "gray", size = 1, linetype = "dashed") +
  geom_text(x=0.103 * 1000, y= 0.000001  * 100, label="Under-5 mortality = 80",size = 4.5, 
            color = "darkgray", fontface="plain") +
  geom_hline(yintercept = 0.096 * 100, color = "gray", size = 1, linetype = "dashed") +
  geom_text(x=0.021 * 1000, y=0.085  * 100, label="Current estimate",size = 4.5, 
            color = "gray", fontface="plain") + 
  geom_text(x=0.092 * 1000, y= 0.085  * 100, label="B/C = 110",size = 4.5, 
            color = "black", fontface="plain")
p3
dev.off()


jpeg("fig1_2021jun.jpg", width = 21, height = 7, units = 'in', res = 300)
(p1|p2|p3)
dev.off()


