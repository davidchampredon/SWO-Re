source('utils.R')

# ---- Data ----

# Retrieve and filter data:

date.type = c('specimen', 'episode', 'case', 'test')

rawdat = get_data()

dfl = list()

for(i in seq_along(date.type)){
    dfl[[i]] = digest_data(rawdat, date.type = date.type[i])
    dfl[[i]]$dateType = date.type[[i]]
}
df = do.call('rbind', dfl)

phus = unique(df$Reporting_PHU)
n    = length(phus)

first.date = lubridate::ymd('2020-03-30')

window.size = 7 # in days
si.type = 'uncertain'  # parametric  or  uncertain
si_mean = 4    # Mean in days
si_stdv = 4.75 # Standard deviation in days
reporting.lag = 6  # in days

pdf('sensi-date-type.pdf')
for(i in 1:n){  # i=1
    
    # Public Health Unit selection:
    phu = phus[i]
    print(paste(i,'/',n,':',phu))
    
    # Statistical estimation:
    udt = unique(df$dateType)
    Rest = list()
    for(k in seq_along(udt)){ # k=1
        print(udt[k])
        tmp = filter(df, dateType == udt[k])
        
        Rest[[k]] = calc_R_unit(phu, df = tmp, 
                                window.size,
                                first.date,
                                si_mean, 
                                si_stdv,
                                si.type = 'uncertain')
        
        Rest[[k]]$df.R$date.type = udt[k]
    }
    a = lapply(Rest, '[[','df.R')
    z = do.call('rbind', a)
    
    g.m = z %>% ggplot() + 
        geom_line(aes(x=date_end, y=m, colour = date.type),
                  alpha = 0.8,
                  size=2) + 
        xlab('date')  + ylab('R eff. (mean)')+
        ggtitle(phu)
    
    g.ci = z %>% 
        mutate(ci.width = qhi-qlo) %>%
        ggplot(aes(x=date_end, fill = date.type)) + 
        geom_boxplot(aes(x=date.type, y=ci.width)) +
        # geom_density(aes(ci.width)) +
        # geom_ribbon(aes(ymin = qlo, ymax = qhi), alpha=0.2) + 
        scale_color_brewer(palette = 'Paired') + 
        xlab('')  + ylab('CI width') #+ggtitle(phu)
    
    gridExtra::grid.arrange(g.m, g.ci, ncol=1)
}
dev.off()
