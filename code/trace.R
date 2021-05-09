library(tidyverse)
library(readxl)

exp_fn = function(t, i, tau) {
  return(i * exp(-t/tau))
}

## average 4.5-sec current traces -----
avg_current_trace <- function(xlsx_paths) {
  num_files <- length(xlsx_paths)

  df <- read_excel(xlsx_paths[1])
  t <- df[1]
  c <- df[-1]
  
  for (i in 2:num_files) {
    df <- read_excel(xlsx_paths[i])
    c <- c + df[-1]
  }
  c <- c/num_files
  
  avg_trc <- bind_cols(t, c)
  return(avg_trc)
}

wt_paths <- dir('./data/Final MGAT1KO Data for Hui/JMCC paper/K Currents 14 Weeks/WT-4.5s-traces', 
                full.names = TRUE)
ko_paths <- dir('./data/Final MGAT1KO Data for Hui/JMCC paper/K Currents 14 Weeks/Mgat1KO-4.5s-traces',
                full.names = TRUE)

wt_avg <- avg_current_trace(wt_paths)
ko_paths <- avg_current_trace(ko_paths)

write_csv(wt_avg, './data/Final MGAT1KO Data for Hui/JMCC paper/K Currents 14 Weeks/4.5s-avg-wt.csv')
write_csv(ko_avg, './data/Final MGAT1KO Data for Hui/JMCC paper/K Currents 14 Weeks/4.5s-avg-ko.csv')

## data handling -----
# K+ table data
k_ko <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx')
k_wt <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx')

colnames(k_ko) <- colnames(k_wt)
k <- bind_rows(mutate(k_ko, group = 'KO'), mutate(k_wt, group = 'WT'))

agg_k <- k %>% 
  group_by(group) %>% 
  summarise(IKslow2 = mean(A1, na.rm = TRUE), IKslow1 = mean(A2, na.rm = TRUE), Ito = mean(A3, na.rm = TRUE),
            tau1 = mean(Tau1, na.rm = TRUE), tau2 = mean(Tau2, na.rm = TRUE), tau3 = mean(Tau3, na.rm = TRUE),
            Iss = mean(Iss, na.rm = TRUE), capa = mean(Cap, na.rm = TRUE))

# K+ raw traces
k_ko_trace_paths <- dir('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/KO_traces', full.names = TRUE)
k_wt_trace_paths <- dir('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/WT_traces', full.names = TRUE)

k_ko_traces <- vector('list', length = length(k_ko_trace_paths))
k_wt_traces <- vector('list', length = length(k_wt_trace_paths))

# 35 files for each KO and WT
for (i in seq_along(k_ko_trace_paths)) {
  k_ko_traces[[i]] <- read_excel(k_ko_trace_paths[i])
  k_wt_traces[[i]] <- read_excel(k_wt_trace_paths[i])
}

# average traces
trace_time <- k_ko_traces[[1]]$`Time(ms)`
k_ko_trace <- vector('numeric', length = length(k_ko_traces[[1]]$Trace))
k_wt_trace <- vector('numeric', length = length(k_wt_traces[[1]]$Trace))
for (i in seq_along(k_ko_trace_paths)) {
  k_ko_trace <- k_ko_trace + k_ko_traces[[i]]$Trace
  k_wt_trace <- k_wt_trace + k_wt_traces[[i]]$Trace
}
k_ko_trace <- k_ko_trace / length(k_ko_trace_paths)
k_wt_trace <- k_wt_trace / length(k_wt_trace_paths)

k_trace <- tibble(time = trace_time, KO = k_ko_trace, WT = k_wt_trace)

# k_trace <- k_trace %>% 
#   pivot_longer(2:3, values_to = 'trace', names_to = 'group') %>% 
#   arrange(group)
# k_trace %>% ggplot(aes(x = time, y = trace, color = group)) +
#   geom_point()


## extract K+ traces
t <- seq(0, 27499.95, by = 0.05)

# IKslow1
IKslow1_ko <- exp_fn(t, agg_k$IKslow1[1], agg_k$tau2[1]) * agg_k$capa[1]
IKslow1_wt <- exp_fn(t, agg_k$IKslow1[2], agg_k$tau2[2]) * agg_k$capa[2]

# IKslow2
IKslow2_ko <- exp_fn(t, agg_k$IKslow2[1], agg_k$tau1[1]) * agg_k$capa[1]
IKslow2_wt <- exp_fn(t, agg_k$IKslow2[2], agg_k$tau1[2]) * agg_k$capa[2]

# Ito
Ito_ko <- exp_fn(t, agg_k$Ito[1], agg_k$tau3[1]) * agg_k$capa[1]
Ito_wt <- exp_fn(t, agg_k$Ito[2], agg_k$tau3[2]) * agg_k$capa[2]

# Iss
Iss_ko <- agg_k$Iss[1] * agg_k$capa[1]
Iss_wt <- agg_k$Iss[2] * agg_k$capa[2]

# Ito trace
Ito_trace <- k_trace %>% 
  mutate(KO = KO - IKslow1_ko - IKslow2_ko - Iss_ko, WT = WT - IKslow1_wt - IKslow2_wt - Iss_wt)

# write_csv(Ito_trace, './Ito_trace.csv')

Ktrace <- read_csv('./k_trace.csv')
time <- Ktrace$time
which(time == 470)
which(time == 25*1000)
