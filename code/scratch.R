library(tidyverse)
library(readxl)
library(janitor)


ko <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/Ik 4.5 seconds IV protocol.xlsx', 
                 range = cell_limits(c(1, 1), c(42, 8)))
new_ko <- ko %>% 
  remove_empty(which = 'rows')
new_ko$Weeks <- new_ko$Weeks %>% 
  replace_na('')
new_ko$Control <- new_ko$Control %>% 
  str_replace_all('[[:space:]]+', '')
new_ko <- new_ko %>% 
  mutate(id = str_c(Weeks, Control)) %>% 
  select(id, `Peak A/F`, `Ito A/F`, `Ikslow A/F`, `Iss A/F`, Tau2, tau1)
names(new_ko) <- c('id', 'peak', 'peak_Ito', 'peak_IKslow', 'Iss', 'tau_Ito', 'tau_IKslow')
write_excel_csv(new_ko, './MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/Ik-4.5s-KO.csv')

wt <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/Ik 4.5 seconds IV protocol.xlsx', 
                 range = cell_limits(c(1, 10), c(45, 17)))
new_wt <- wt %>% 
  remove_empty(which = 'rows')
new_wt$Weeks <- new_wt$Weeks %>% 
  replace_na('')
new_wt$`Date Cell MGAT1KO` <- new_wt$`Date Cell MGAT1KO` %>% 
  str_replace_all('[[:space:]]+', '')
new_wt <- new_wt %>% 
  mutate(id = str_c(Weeks, `Date Cell MGAT1KO`)) %>% 
  select(id, `Peak A/F -FF`, `Ito A/F/ - FF`, `Ikslow A/F - FF`, `Iss A/F -FF`, `tau2 -FF`, `tau1 -FF`)
names(new_wt) <- c('id', 'peak', 'peak_Ito', 'peak_IKslow', 'Iss', 'tau_Ito', 'tau_IKslow')
write_excel_csv(new_wt, './MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/Ik-4.5s-WT.csv')
