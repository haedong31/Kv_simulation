library(tidyverse)
library(readxl)


k_ko <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx')
k_wt <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx')

colnames(k_ko) <- colnames(k_wt)
k <- bind_rows(mutate(k_ko, group = 'KO'), mutate(k_wt, group = 'WT'))

agg_k <- k %>% 
  group_by(group) %>% 
  summarise(IKslow2 = mean(A1, na.rm = TRUE), IKslow1 = mean(A2, na.rm = TRUE), Ito = mean(A3, na.rm = TRUE),
            tau1 = mean(Tau1, na.rm = TRUE), tau2 = mean(Tau2, na.rm = TRUE), tau3 = mean(Tau3, na.rm = TRUE))

IKslow1_trace = function(t, i, tau) {
  return(i * exp(-t/tau))
}

t <- seq(0, 25, by = 0.001)
IKslow1 <- IKslow1_trace(t, agg_k$IKslow1[1], agg_k$tau2[1])
plot(t, IKslow1)
qplot(t, IKslow1)
