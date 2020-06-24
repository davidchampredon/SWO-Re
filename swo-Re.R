###
###  Effective Reproduction Number ("Re") Estimation
###  for Selected Localities in SouthWestern Ontario.
###  Use the `EpiEsptim` package.
###  
###  Contact       : David Champredon
###  Creation date : 2020-06-20
###

source('utils.R')

# ---- Data ----

# Retrieve and filter data:

date.type = 'specimen'
df = get_data() %>%
    digest_data(date.type = date.type)
phus = unique(df$Reporting_PHU)
n    = length(phus)
pdf('plot-all-data.pdf')
plot_data(df)
dev.off()

# ---- Parameters ----

# First anchor date across all regions:
first.date = lubridate::ymd('2020-03-30')

# Sliding window (looking backward from the calculation date)
# to calculate Re over time. 
# Re will be assumed constant insidde this window.
window.size = 7 # in days

# Serial interval specification:
si.type = 'uncertain'  # parametric  or  uncertain
si_mean = 4    # Mean in days
si_stdv = 4.75 # Standard deviation in days
# Sources for serial interval:
# https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
# https://www.jwatch.org/na51171/2020/03/27/serial-interval-covid-19
# https://pubmed.ncbi.nlm.nih.gov/32145466/



# ---- Re Estimation ----

# -- For each PHU:

figname = paste0('plot-R-w',window.size,'-si',si_mean,'.pdf')
pdf(file = figname)

for(i in 1:n){  # i=3
    
    # Public Health Unit selection:
    phu = phus[i]
    print(paste(i,'/',n,':',phu))
    
    # Statistical estimation:
    Rest = calc_R_unit(phu, df, 
                       window.size,
                       first.date,
                       si_mean, 
                       si_stdv,
                       si.type = 'uncertain')
    
    # Plots:
    g.R = plot_R(Rest, first.date, title = phu )
    g   = plot_cases(phu, df, first.date, date.type)
    gridExtra::grid.arrange(g, g.R, ncol=1)
}
dev.off()

# -- Grouped PHUs:

figname = paste0('plot-group-R-w',window.size,'-si',si_mean,'.pdf')
pdf(file = figname)
for(grp in 1:2){ # grp=1
    
    message(paste('Group',grp, ':', get_group_name(grp)))
    
    Rest.group = calc_R_group(group.number = grp, 
                              df, window.size, 
                              first.date, 
                              si_mean, si_stdv, si.type = 'uncertain')
    
    g.R.grp = plot_R(Rest = Rest.group, 
                     first.date = first.date, 
                     title = paste("Group :", get_group_name(grp)) )
    g.grp   = plot_cases_group(group.number = grp, df, first.date, date.type)
    gridExtra::grid.arrange(g.grp, g.R.grp, ncol=1)
    }
dev.off()