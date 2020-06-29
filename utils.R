library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(EpiEstim)

#' Retrieve the data from Ontario's MoH publicly available line list.
#' Warning: The URL may have to be updated.
#' NOTE: interactive URL:  https://data.ontario.ca/dataset/confirmed-positive-cases-of-covid-19-in-ontario/resource/455fd63b-603d-4608-8216-7d8647f43350
get_data <- function() {
    url = "https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv"
    dat = read.csv(url)
    return(dat)
}

#' Filter and reformat data.
#' @param dat Dataframe as returned by \code{get_data()}.
#' @param date.type String. Type of date selected \code{specimen} or \code{episode}.
#' @param phu.filename String. Path to the file that specifies the public health units selected.
#' @param remove.travel Logical. Remove travel related cases? Removing travel-related cases implicitly assumes that those identified cases were isolated and do not contribute to local chain of transmissions.
#' @return A dataframe.
#' 
digest_data <- function(dat, 
                        date.type = 'specimen',
                        phu.filename = 'phus.csv',
                        remove.travel = TRUE) {
    # DEBUG:
    # unique(dat$Reporting_PHU)
    
    phu.select = readLines(phu.filename)
    
    if(date.type == 'specimen') 
        dat = mutate(dat, date = lubridate::ymd(Specimen_Date))
    if(date.type == 'episode')
        dat = mutate(dat, date = lubridate::ymd(Accurate_Episode_Date))
    
    # Travel related cases:
    if(remove.travel){
        dat = filter(dat, Case_AcquisitionInfo != 'Travel')
    }
    
    # Calculate incidence for each date:    
    df = dat %>%
        filter(Reporting_PHU %in% phu.select) %>%
        group_by(date, Reporting_PHU) %>%
        summarise(n = n())
    
    df$Reporting_PHU = as.character(df$Reporting_PHU)
    
    return(df)    
}

#' Plot the entire dataset.
plot_data <- function(df) {
    g = df %>% 
        ggplot(aes(x=date, y=n))+
        geom_step()+
        xlab('')+
        facet_wrap(~Reporting_PHU, scales = 'fixed')
    plot(g)
}



#' Group the data of multiple Public Health Units.
#'
#' @param group.number Integer. Group number.
#' @param df Dataframe of data.
#' @param first.date Date. Anchor date across multiple locations.  
#'
#' @return Dataframe of grouped data.
#' 
group_phu_data <- function(group.number, df, first.date) {
    phugroup = read.csv('phu-group.csv', stringsAsFactors = F) %>%
        filter(group == group.number)
    
    # Select the PHU:
    dfs = df %>% 
        filter(Reporting_PHU %in% phugroup$phu) %>%
        filter(date > first.date) %>%
        rename(n_unit = n) %>%
        group_by(date) %>%
        summarize(n = sum(n_unit))
    
    # Fill in the missing dates:
    dj = fill_dates(dfs)
    
    return(dj)
}


get_group_name <- function(group.number) {
    
    phugroup = read.csv('phu-group.csv', stringsAsFactors = F) %>%
        filter(group == group.number)
    
    x = paste(phugroup$phu, collapse = '')
    x = str_remove_all(x,'Health')
    x = str_remove_all(x,'Unit')
    x = str_remove_all(x,'Public')
    x = str_remove_all(x,'County')
    res = x
    return(res)
}

#' Technical function for setting up serial intervals and time points.
setup_epiestim <- function(window.size, n.obs, si_mean, si_stdv) {
    
    # Sliding window to calculater R:
    t_start = seq(2, n.obs - window.size)   
    t_end   = t_start + window.size     
    
    # Serial Intervals
    si <- EpiEstim::make_config(list(mean_si = si_mean, 
                                     std_si  = si_stdv,
                                     t_start = t_start, 
                                     t_end   = t_end))
    
    usi <- make_config(list(
        t_start = t_start, 
        t_end   = t_end,
        # mean SI:
        mean_si     = si_mean, 
        std_mean_si = 2,
        min_mean_si = si_mean-2, 
        max_mean_si = si_mean+2,
        # std. dev SI:
        std_si      = si_stdv, 
        std_std_si  = 0.5,
        min_std_si  = si_stdv - 0.5, 
        max_std_si  = si_stdv + 0.5,
        # replicates:
        n1 = 100, n2 = 100))
    
    return(list(t_start = t_start, t_end=t_end,
                si=si, usi = usi, 
                window.size=window.size))
}

#' Fill-in the missing dates with 0:
fill_dates <- function(dfs) {
    dd = data.frame(date = seq(first.date, max(dfs$date), by = '1 day'))
    dj = left_join(dd,dfs, by='date')
    dj$n[is.na(dj$n)] <- 0
    return(dj)
}

#' Calculate the effective reproduction number (aka. R or Re or Rt)
#' using the R package EpiEstim.
#'
#' @param phu String. Name of the public health unit.
#' @param df Dataframe of cases data.
#' @param window.size Integer. Sliding window over which Re is calculated (and assumed constant). 
#' @param first.date Date. First anchor date (conveniently start the data and calculation after this date when analyzing several regions)
#' @param si_mean Numeric. Mean of the serial interval in days.
#' @param si_stdv Numeric. Standard deviation of the serial interval in days.
#' @param si.type String. Type of model for the serial interval. "parameteric" or "uncertain" (default).
#' @return List containing estimates for Re.
#' 
calc_R_unit <- function(phu, df, window.size,
                        first.date,
                        si_mean, 
                        si_stdv,
                        si.type = 'uncertain') {
    
    # Select the PHU:
    dfs = df %>% 
        filter(Reporting_PHU == phu) %>%
        filter(date > first.date)
    
    # Fill in the missing dates:
    dj = fill_dates(dfs)
    
    # Setup before Re estimation:
    z = setup_epiestim(window.size, 
                       n.obs   = nrow(dj), 
                       si_mean = si_mean, 
                       si_stdv = si_stdv)
    t_end   = z$t_end
    si      = z$si
    usi     = z$usi
    
    config = usi
    if(si.type=='parametric') config = si
    
    
    # Re estimation using EpiEstim:
    R.est = EpiEstim::estimate_R(incid = dj$n, 
                                 method = paste(si.type,'si',sep="_"),
                                 config = config)
    
    # Reformat results:
    df.R = R.est$R %>%
        mutate(date_end = max(dfs$date) + t_end - max(t_end)) %>%
        rename(m = `Mean(R)`, 
               qlo = `Quantile.0.025(R)`, 
               qhi = `Quantile.0.975(R)`)
    
    return(list(df.R = df.R, 
                config = config, 
                window.size = window.size))
}


#' Calculate the effective reproduction number (aka. R or Re or Rt)
#' using the R package EpiEstim.
#'
#' @param group.number Integer. Group number where several PHUs have been merged into.
#' @param df Dataframe of cases data.
#' @param window.size Integer. Sliding window over which Re is calculated (and assumed constant). 
#' @param first.date Date. First anchor date (conveniently start the data and calculation after this date when analyzing several regions)
#' @param si_mean Numeric. Mean of the serial interval in days.
#' @param si_stdv Numeric. Standard deviation of the serial interval in days.
#' @param si.type String. Type of model for the serial interval. "parameteric" or "uncertain" (default).
#' @return List containing estimates for Re.
#' 
calc_R_group <- function(group.number, 
                         df, 
                         window.size,
                         first.date,
                         si_mean, 
                         si_stdv,
                         si.type = 'uncertain') {  # group.number=1
    
    # Fill in the missing dates:
    dj = group_phu_data(group.number, df, first.date)
    
    # Setup before Re estimation:
    z = setup_epiestim(window.size, 
                       n.obs   = nrow(dj), 
                       si_mean = si_mean, 
                       si_stdv = si_stdv)
    t_end   = z$t_end
    si      = z$si
    usi     = z$usi
    
    config = usi
    if(si.type=='parametric') config = si
    
    
    # Re estimation using EpiEstim:
    R.est = EpiEstim::estimate_R(incid = dj$n, 
                                 method = paste(si.type,'si',sep="_"),
                                 config = config)
    
    # Reformat results:
    df.R = R.est$R %>%
        mutate(date_end = max(dj$date) + t_end - max(t_end)) %>%
        rename(m = `Mean(R)`, 
               qlo = `Quantile.0.025(R)`, 
               qhi = `Quantile.0.975(R)`)
    
    return(list(df.R = df.R, 
                config = config, 
                window.size = window.size))
}



#' Create a standardized date axis (for cases and Re estimates).
create_axis <- function(first.date, last.date) {
    return(scale_x_date(breaks='2 weeks', 
                        date_labels = '%b-%d', 
                        limits = c(first.date, last.date + 14), 
                        date_minor_breaks = '1 week'))
}

#' Create the gradient shaded area indicating the reporting lag.
reporting_gradient <- function(df, reporting.lag){
    td = lubridate::today()
    
    if ('n' %in% names(df)) y.max = max(df$n+1)
    if ('m' %in% names(df)) y.max = max(df$m+1)
    
    
    dd = data.frame(xmin = td - seq(1,reporting.lag,by=0.5),
                    xmax = td,
                    ymin = 0,
                    ymax = y.max)
    gr = geom_rect(data = dd, 
                   inherit.aes = FALSE,
                   aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                   fill = 'tomato2',
                   alpha = 0.06)
    return(gr)
}


#' Plot estimates of Re over time. 
plot_R <- function(Rest, first.date, 
                   reporting.lag = 7, 
                   title='') {
    
    # Unpack:
    df.R   = Rest$df.R
    config = Rest$config
    
    # Last value of Re is displayed:
    last.R    = round(df.R$m[nrow(df.R)],2)
    last.R.lo = round(df.R$qlo[nrow(df.R)],2)
    last.R.hi = round(df.R$qhi[nrow(df.R)],2)
    last.date = max(df.R$date_end)
    
    report.grad = reporting_gradient(df.R, reporting.lag)
    
    # Cosmetics: 
    last.shift = 9
    alpha.rib  = 0.3
    alpha.line = 0.9
    col.R      = 'steelblue2'
    axisdate   = create_axis(first.date, last.date = lubridate::today())
    
    g.R = ggplot(df.R, aes(x=date_end))+
        report.grad + 
        geom_hline(yintercept = 1, linetype='dashed') +
        geom_ribbon(aes(ymin = qlo,
                        ymax = qhi), 
                    alpha = alpha.rib,
                    fill = col.R) +
        geom_line(aes(y=m), size=2, colour = col.R, alpha=alpha.line) + 
        geom_point(aes(y=m), 
                   size=1.5, shape = 21,
                   colour=col.R, fill = 'white') + 
        # Last estimates:
        geom_label(data = data.frame(x=last.date, y=last.R),
                   aes(x = x + last.shift, 
                       y = y, 
                       label = paste0(last.R, ' (',
                                      last.R.lo,' ; ',
                                      last.R.hi,')')),
                   size = 2, 
                   fontface='bold',
                   colour = col.R) +
        axisdate +
        ylab('Eff. Reprod. Number') +
        xlab('') +
        ggtitle(label = title,
                subtitle = paste0('sliding window = ', 
                                  Rest$window.size, 
                                  ' days ; serial intrv: mean = ', 
                                  round(config$mean_si,2), ' days, stdv = ',
                                  round(config$std_si,2),' days'))
    # g.R
    return(g.R)
}




#' Plot cases used for a specific public health unit.
plot_cases <- function(phu, df, first.date, date.type, 
                       reporting.lag = 7) {
    dfs = df %>% 
        filter(Reporting_PHU == phu) %>%
        filter(date > first.date) %>%
        fill_dates()
    
    axisdate = create_axis(first.date, last.date = lubridate::today())
    
    report.grad = reporting_gradient(dfs, reporting.lag)
    
    g = dfs %>%
        ggplot(aes(x=date, y=n)) +
        # -- Reporting lag
        report.grad +
        # -- Data
        geom_step()+
        geom_point(size=1)+
        axisdate +
        theme(panel.grid.minor.y = element_blank()) +
        scale_y_continuous(limits=c(0,max(dfs$n+1)))+
        
        ylab('Count') +
        xlab(paste(date.type,'date')) +
        ggtitle(phu, subtitle = date.type)
    g
    return(g)
}

#' Plot cases used for a grouped public health unit.
plot_cases_group <- function(group.number, df, first.date, 
                             date.type, 
                             reporting.lag = 7) {
    
    df.grp = group_phu_data(group.number, df, first.date)
    
    axisdate = create_axis(first.date, last.date = lubridate::today())
    
    report.grad = reporting_gradient(df.grp, reporting.lag)
    
    g = df.grp %>%
        ggplot(aes(x=date, y=n))+
        report.grad + 
        geom_step() +
        geom_point(size=1)+
        axisdate +
        theme(panel.grid.minor.y = element_blank()) +
        scale_y_continuous(limits=c(0,max(df.grp$n+1)))+
        ylab('Count') +
        xlab(paste(date.type,'date')) +
        ggtitle(paste('Group: ', get_group_name(group.number)), 
                subtitle = date.type)
    return(g)
}

