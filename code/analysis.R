library(tidyverse)
library(readxl)
library(ggpubr)


## data pre-processing -----
# experimental data
k_ko <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx')
k_wt <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx')

colnames(k_ko) <- colnames(k_wt)
k <- bind_rows(mutate(k_ko, group = 'KO'), mutate(k_wt, group = 'WT'))

# results of GA
read_names <- dir('./results/', '.*xlsx')
ga_rsts <- vector('list', length(read_names))
for (i in seq_along(read_names)) {
  ga_rsts[[i]] <- read_excel(str_c('./results/', read_names[i]), col_names = FALSE)
}
ga_Iss_ko <- ga_rsts[[1]]
ga_Iss_wt <- ga_rsts[[2]]
ga_Ito_ko <- ga_rsts[[3]]
ga_Ito_wt <- ga_rsts[[4]]

colnames(ga_Iss_ko) <- c('x1', 'x2', 'x3', 'x4', 'SSE', 'Iss_hat')
colnames(ga_Iss_wt) <- c('x1', 'x2', 'x3', 'x4', 'SSE', 'Iss_hat')
colnames(ga_Ito_ko) <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'SSE', 'Ito_hat')
colnames(ga_Ito_wt) <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'SSE', 'Ito_hat')

ga_Iss <- bind_rows(mutate(ga_Iss_ko, group = 'KO'), mutate(ga_Iss_wt, group = 'WT'))
ga_Ito <- bind_rows(mutate(ga_Ito_ko, group = 'KO'), mutate(ga_Ito_wt, group = 'WT'))

# for barplots with error bars
k_4barplot <- k %>% 
  group_by(group) %>% 
  summarise(mean_Iss = mean(Iss, na.rm = TRUE), mean_Ito = mean(A3, na.rm = TRUE),
            sem_Iss = sd(Iss, na.rm = TRUE)/sqrt(n()), sem_Ito = sd(A3, na.rm = TRUE)/sqrt(n())) %>% 
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
