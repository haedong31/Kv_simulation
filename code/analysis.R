library(tidyverse)
library(readxl)
library(ggpubr)


## custom function -----
data_processor <- function(excel_path, num_vars) {
  rst <- read_excel(excel_path, col_names = FALSE)
  for (i in 1:num_vars) {
    colnames(rst)[i] <- str_c('x', i)
  }
  colnames(rst)[num_vars + 1] <- 'SSE'
  colnames(rst)[num_vars + 2] <- 'predicted'
  
  return(rst)
}


## data pre-processing -----
# experimental data
k_ko <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx')
k_wt <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx')

colnames(k_ko) <- colnames(k_wt)
k <- bind_rows(mutate(k_ko, group = 'KO'), mutate(k_wt, group = 'WT'))

# results of GA
read_names <- dir('./results/', '.*xlsx')

ga_Iss_paths <- read_names[str_detect(read_names, 'Iss')]
ga_Iss_paths <- str_c('./results/', ga_Iss_paths)
ga_Ito_paths <- read_names[str_detect(read_names, 'Ito')]
ga_Ito_paths <- str_c('./results/', ga_Ito_paths)
ga_Ikslow1_paths <- read_names[str_detect(read_names, 'IKslow1')]
ga_Ikslow1_paths <- str_c('./results/', ga_Ikslow1_paths)

ga_Iss_ko <- data_processor(ga_Iss_paths[1], 4)
ga_Iss_wt <- data_processor(ga_Iss_paths[2], 4)
ga_Ito_ko <- data_processor(ga_Ito_paths[1], 6)
ga_Ito_wt <- data_processor(ga_Ito_paths[2], 6)
ga_Ikslow1_ko <- data_processor(ga_Ikslow1_paths[1], 8)
ga_Ikslow1_wt <- data_processor(ga_Ikslow1_paths[2], 8)

ga_Iss <- bind_rows(mutate(ga_Iss_ko, group = 'KO'), mutate(ga_Iss_wt, group = 'WT'))
ga_Ito <- bind_rows(mutate(ga_Ito_ko, group = 'KO'), mutate(ga_Ito_wt, group = 'WT'))
ga_Ikslow1 <- bind_rows(mutate(ga_Ikslow1_ko, group = 'KO'), mutate(ga_Ikslow1_wt, group = 'WT'))

# for barplots with error bars
k_4barplot <- k %>% 
  group_by(group) %>% 
  summarise(mean_Iss = mean(Iss, na.rm = TRUE), sem_Iss = sd(Iss, na.rm = TRUE)/sqrt(n()),
            mean_Ito = mean(A3, na.rm = TRUE), sem_Ito = sd(A3, na.rm = TRUE)/sqrt(n()),
            mean_Ikslow1 = mean(A2, na.rm = TRUE), sem_Ikslow1 = sd(A2, na.rm = TRUE)/sqrt(n()),
            mean_Ikslow2 = mean(A1, na.rm = TRUE), sem_Ikslow2 = sd(A1, na.rm = TRUE)/sqrt(n())) %>% 
  mutate(clf = 'Real Data')

ga_Iss_4barplot <- ga_Iss %>% 
  group_by(group) %>% 
  summarise(mean_Iss = mean(Iss_hat, na.rm = TRUE), sem_Iss = sd(Iss_hat, na.rm = TRUE)/sqrt(n())) %>% 
  mutate(clf = 'Simulated')

ga_Ito_4barplot <- ga_Ito %>% 
  group_by(group) %>% 
  summarise(mean_Ito = mean(Ito_hat, na.rm = TRUE), sem_Ito = sd(Ito_hat, na.rm = TRUE)/sqrt(n())) %>% 
  mutate(clf = 'Simulated')


## EDA -----
# box plot of fitted parameters
ga_Iss %>% 
  pivot_longer(colnames(ga_Iss)[1:4], names_to = 'param', values_to = 'value') %>% 
  ggplot(aes(x = param, y = value, color = group)) +
  geom_boxplot() +
  xlab('Tuning Parameters') +
  ylab('') +
  ggtitle('Iss Parameters - Genetic Algorithm') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom')

ga_Ito %>% 
  pivot_longer(colnames(ga_Ito)[1:6], names_to = 'param', values_to = 'value') %>% 
  ggplot(aes(x = param, y = value, color = group)) +
  geom_boxplot() +
  xlab('Tuning Parameters') +
  ylab('') +
  ggtitle('Ito Parameters - Genetic Algorithm') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom')

ga_Ikslow1 %>% 
  pivot_longer(colnames(ga_Ikslow1)[1:6], names_to = 'param', values_to = 'value') %>% 
  ggplot(aes(x = param, y = value, color = group)) +
  geom_boxplot() +
  xlab('Tuning Parameters') +
  ylab('') +
  ggtitle('Ikslow1 Parameters - Genetic Algorithm') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom')

# bar plot with error bars of mean Iss & Ito
z95 <- qnorm(0.975, mean = 0, sd = 1)
p1 <- bind_rows(select(k_4barplot, -mean_Ito, -sem_Ito), ga_Iss_4barplot) %>% 
  ggplot(aes(x = clf, y = mean_Iss, fill = group)) +
  geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_Iss - sem_Iss, ymax = mean_Iss + sem_Iss), width = 0.2,
                position = position_dodge(0.9)) +
  xlab('') +
  ylab('Density (pA/pF)')

p2 <- bind_rows(select(k_4barplot, -mean_Iss, -sem_Iss), ga_Ito_4barplot) %>% 
  ggplot(aes(x = clf, y = mean_Ito, fill = group)) +
  geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_Ito - sem_Ito, ymax = mean_Ito + sem_Ito), width = 0.2,
                position = position_dodge(0.9)) +
  xlab('') +
  ylab('Density (pA/pF)') 

ggarrange(p1, p2, labels = c("Iss", "Ito"), nrow = 1, ncol = 2)
