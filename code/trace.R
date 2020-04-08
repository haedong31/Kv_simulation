library(tidyverse)
library(readxl)


## custom functions -----
exp_fn = function(t, i, tau) {
  return(i * exp(-t/tau))
}


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
            Iss = mean(Iss, na.rm = TRUE))

# K+ raw traces
k_trace <- read_excel()


## extract K+ traces
t <- seq(0, 25, by = 0.001)

# IKslow1
IKslow1_ko <- exp_fn(t, agg_k$IKslow1[1], agg_k$tau2[1])
IKslow1_wt <- exp_fn(t, agg_k$IKslow1[2], agg_k$tau2[2])

# IKslow2
IKslow2_ko <- exp_fn(t, agg_k$IKslow2[1], agg_k$tau1[1])
IKslow2_wt <- exp_fn(t, agg_k$IKslow2[2], agg_k$tau1[2])

# Ito
Ito_ko <- exp_fn(t, agg_k$Ito[1], agg_k$tau3[1])
Ito_wt <- exp_fn(t, agg_k$Ito[2], agg_k$tau3[2])

# Iss
IKslow1_ko <- agg_k$Iss[1]
IKslow1_wt <- agg_k$Iss[2]

plot(t, IKslow1)
qplot(t, IKslow1)
