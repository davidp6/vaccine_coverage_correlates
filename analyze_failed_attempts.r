# check how often a fully-vaccinated child experienced any failed attempts
# could indicate that the "reasons for non-vaccination" questions from other surveys are valid representations of "turn-away rate"
rm(list=ls())
library(data.table)

load("J:/Project/Evaluation/GAVI/hhs/zmb/data_production/release/household_external.rdata")
load("J:/Project/Evaluation/GAVI/hhs/zmb/data_production/release/child_external.rdata")

data = merge(householdData, childData, 'customid')

# make fully vaccinated
data[, fully_vaccinated:=either_penta3==1 & either_polio3==1 & either_measles==1]

# check frequency (there are many fully-vaccinated children who had at least one failed attempt)
table(data[child_age<2, c('attempt_trig.x','fully_vaccinated')], useNA='always')
data[child_age<2, sum(nosuc_num>0, na.rm=T), by='fully_vaccinated']

# check against the "never vaccinated why" variable
table(data[,c('attempt_trig.x', 'never_vac_why_t1_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t2_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t3_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t4_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t5_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t6_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t7_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t8_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t9_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t10_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t11_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t12_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t13_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t14_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t15_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t16_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t17_')])
table(data[,c('attempt_trig.x', 'never_vac_why_t18_')])
