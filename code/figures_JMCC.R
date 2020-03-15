library(tidyverse)
library(readxl)


## Fig. 3 -----
# Repolarizing K+ currents at a test potential 50 mV for 25 sec
fig3_wt <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-WT.xlsx')
fig3_ko <- read_excel('./MGAT1_Data_tidy/JMCC/K Currents 14 Weeks/potassium-KO.xlsx')

fig3_wt_agg <- fig3_wt %>% 
  select(-Weeks, -VAR2, -Seal, -Resis) %>% 
  summarise_all(~mean(., na.rm = TRUE)) %>% 
  mutate(group = 'WT')

fig3_ko_agg <- fig3_ko %>% 
  select(-Weeks, -`Date Cell FF`, -`Seal FF`, -`Resis FF`) %>% 
  summarise_all(~mean(., na.rm = TRUE)) %>% 
  mutate(group = 'KO')

colnames(fig3_ko_agg) <- colnames(fig3_wt_agg)
fig3 <- bind_rows(fig3_ko_agg, fig3_wt_agg)
fig3 <- fig3 %>% 
  select(-C)

# Fig. 3-B
fig3 %>% 
  select(2:6, group) %>% 
  pivot_longer(colnames(fig3)[2:6], names_to = 'measure') %>% 
  ggplot(aes(x = measure, y = value, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  xlab('')
  ylab('Density (pA/pF)')

# Fig. 3-C
fig3 %>% 
  select(Cap, group) %>% 
  pivot_longer(Cap, names_to = 'measure') %>% 
  ggplot(aes(x = measure, y = value, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  xlab('Capacitance') +
  ylab('pF') +
  scale_y_continuous(n.breaks = 6)


## Fig. 4 -----
Na_dv_WT <- read_excel('./MGAT1_Data_tidy/JMCC/Nav Currents/Ina GV MGAT1KO Final.xlsx',
                       range = cell_limits(c(1, 1), c(24, 20)))
Na_dv_KO <- read_excel('./MGAT1_Data_tidy/JMCC/Nav Currents/Ina GV MGAT1KO Final.xlsx',
                       range = cell_limits(c(1, 24), c(17, 43)))
Na_ssi_WT <- read_excel('./MGAT1_Data_tidy/JMCC/Nav Currents/Ina SSI MGAT1KO Final.xlsx',
                        range = cell_limits(c(1, 1), c(26, 21)))
Na_ssi_KO <- read_excel('./MGAT1_Data_tidy/JMCC/Nav Currents/Ina SSI MGAT1KO Final.xlsx',
                        range = cell_limits(c(1, 25), c(17, 45)))

# data pre-processing; density-voltage WT
Na_dv_WT_2 <- Na_dv_WT %>% 
  pivot_longer(colnames(Na_dv_WT)[7:20], names_to = 'V', values_to = 'G')
Na_dv_WT_2$V <- Na_dv_WT_2$V %>% as.numeric()
# G_max <- max(Na_dv_WT_2$G)
# Na_dv_WT_2 <- Na_dv_WT_2 %>% 
#   mutate(ss_act = G / G_max)
Na_dv_WT_3 <- Na_dv_WT_2 %>% 
  group_by(V) %>% 
  summarise(mean_G = mean(G))

# data pre-processing; density-voltage KO
Na_dv_KO_2 <- Na_dv_KO %>% 
  pivot_longer(colnames(Na_dv_KO)[7:20], names_to = 'V', values_to = 'G')
Na_dv_KO_2$V <- Na_dv_KO_2$V %>% as.numeric()
# G_max <- max(Na_dv_KO_2$G)
# Na_dv_KO_2 <- Na_dv_WT_2 %>% 
#   mutate(ss_act = G / G_max)
Na_dv_KO_3 <- Na_dv_KO_2 %>% 
  group_by(V) %>% 
  summarise(mean_G = mean(G))

# data pre-processing; SSI WT
Na_ssi_WT_2 <- Na_ssi_WT %>% 
  pivot_longer(colnames(Na_ssi_WT)[6:21], names_to = 'V', values_to = 'I')
Na_ssi_WT_2$V <- Na_ssi_WT_2$V %>% as.numeric()
# I_max <- max(Na_ssi_WT_2$I)
# Na_ssi_WT_2 <- Na_ssi_WT_2 %>% 
#   mutate(ss_inact = I / I_max)
Na_ssi_WT_3 <- Na_ssi_WT_2 %>% 
  group_by(V) %>% 
  summarise(mean_I = mean(I))

# data pre-processing; SSI KO
Na_ssi_KO_2 <- Na_ssi_KO %>% 
  pivot_longer(colnames(Na_ssi_KO)[6:21], names_to = 'V', values_to = 'I')
Na_ssi_KO_2$V <- Na_ssi_KO_2$V %>% as.numeric()
# I_max <- max(Na_ssi_KO_2$I)
# Na_ssi_KO_2 <- Na_ssi_KO_2 %>% 
#   mutate(ss_inact = I / I_max)
Na_ssi_KO_3 <- Na_ssi_KO_2 %>% 
  group_by(V) %>% 
  summarise(mean_I = mean(I))
Na_act_inact_WT <- full_join(Na_dv_WT_3, Na_ssi_WT_3, 'V')
Na_act_inact_KO <- full_join(Na_dv_KO_3, Na_ssi_KO_3, 'V')

p <- Na_act_inact_WT %>% 
  ggplot(aes(x = V, y = mean_G, color = 'WT')) +
  geom_point() +
  geom_line()
p + geom_point(data = Na_act_inact_WT, aes(x = V, y = mean_I, color = 'WT')) +
  geom_line(data = Na_act_inact_WT, aes(x = V, y = mean_I, color = 'WT')) +
  geom_point(data = Na_act_inact_KO, aes(x = V, y = mean_G, color = 'KO')) +
  geom_line(data = Na_act_inact_KO, aes(x = V, y = mean_G, color = 'KO')) +
  geom_point(data = Na_act_inact_KO, aes(x = V, y = mean_I, color = 'KO')) +
  geom_line(data = Na_act_inact_KO, aes(x = V, y = mean_I, color = 'KO')) +
  ylab('I/Imax, G/Gmax') +
  xlab('Voltage (MV)') 


## Fig. 6 -----
Ca_WT <- read_excel('./MGAT1_Data_tidy/JMCC/Ca Imaging 37 Degrees/Ca Imaging MGAT1KO Final.xlsx', 
                    range = cell_limits(c(1, 1), c(85, 11)))
Ca_KO <- read_excel('./MGAT1_Data_tidy/JMCC/Ca Imaging 37 Degrees/Ca Imaging MGAT1KO Final.xlsx', 
                    range = cell_limits(c(1, 13), c(55, 23)))

# Fig. 6 C-3
fig6_c3_WT <- Ca_WT %>% 
  select(`TimeToPeak (ms)`, 8:11) %>%
  summarise_all(mean) %>% 
  mutate(group = 'WT')
fig6_c3_WT$`Tau (s)` <- fig6_c3_WT$`Tau (s)` * 1000  # s to ms
fig6_c3_WT <- fig6_c3_WT %>% 
  rename(`Tau (ms)` = `Tau (s)`)

fig6_c3_KO <- Ca_KO %>% 
  select(`TimeToPeak (ms)`, 8:11) %>% 
  summarise_all(mean) %>% 
  mutate(group = 'KO')
fig6_c3_KO$`Tau (s)` <- fig6_c3_KO$`Tau (s)` * 1000  # s to ms
fig6_c3_KO <- fig6_c3_KO %>% 
  rename(`Tau (ms)` = `Tau (s)`)

fig6_c3 <- bind_rows(fig6_c3_WT, fig6_c3_KO)
fig6_c3 %>% 
  pivot_longer(colnames(fig6_c3)[1:5], names_to = 'measures') %>% 
  ggplot(aes(x = measures, y = value, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  xlab('Measures') +
  ylab('Time')


## current-densities (normalized) ----
path <- './MGAT1_Data/Final MGAT1KO Data for Hui/JMCC paper/Nav Currents/IV Recordings/EB docs 3-22-19/Mgat1KO-Control INa Vr - 3-22-19.xlsx'

# wild type (control)
normI_wt <- read_excel(path, sheet = 'Control I-V', range = cell_limits(c(1, 1), c(24, 22)))  # cellranger::as.cell_limits('V24')
colnames(normI_wt)[2] <- 'cell_ID'
normI_wt$`Normalized (pA/pF)`[3:5] <- 16.90
normI_wt$`Normalized (pA/pF)`[7] <- 16.90
normI_wt$`Normalized (pA/pF)`[9:11] <- 17.30
normI_wt$`Normalized (pA/pF)`[13:15] <- 17.40
normI_wt$`Normalized (pA/pF)`[17:19] <- 18.90
normI_wt$`Normalized (pA/pF)`[22:23] <- 18.70

normI_wt_2 <- normI_wt %>% 
  pivot_longer(colnames(normI_wt)[3:22], names_to = 'voltage', values_to = 'INa') %>% 
  arrange(cell_ID)
normI_wt_2$voltage <- normI_wt_2$voltage %>% as.numeric()

normI_wt_mean <- normI_wt_2 %>% 
  group_by(voltage) %>% 
  summarise(mean_INa = mean(INa))

# WT plottings
normI_wt_2 %>% 
  ggplot(aes(x = voltage, y = INa, group = cell_ID)) +
  geom_point(aes(color = cell_ID)) +
  geom_line(aes(color = cell_ID)) +
  xlab('Voltage (mV)') +
  ylab('Normalized INa (pA/pF) WT') +
  theme(legend.position = 'bottom')

normI_wt_mean %>% 
  ggplot(aes(x = voltage, y = mean_INa, group = 1)) +
  geom_line() +
  geom_point() +
  xlab('Voltage (mV)') +
  ylab('Mean Normalized INa (pA/pF) WT')

# MGAT1KO 
normI_ko <- read_excel(path, sheet = 'KO IV', range = cell_limits(c(1, 1), c(17, 21)))  # cellranger::as.cell_limits('U17')
normI_ko_2 <- normI_ko %>% 
  pivot_longer(colnames(normI_ko)[2:21], names_to = 'voltage', values_to = 'INa')
normI_ko_2$voltage <- normI_ko_2$voltage %>% as.numeric()

normI_ko_mean <- normI_ko_2 %>% 
  group_by(voltage) %>% 
  summarise(mean_INa = mean(INa))

# MGAT1KO plottings 
normI_ko_2 %>% 
  ggplot(aes(x = voltage, y = INa, group = Normalized)) +
  geom_point(aes(color = Normalized)) +
  geom_line(aes(color = Normalized)) +
  xlab('Voltage (mV)') +
  ylab('Normalized INa (pA/pF) MGAT1KO') +
  theme(legend.position = 'bottom')

normI_ko_mean %>% 
  ggplot(aes(x = voltage, y = mean_INa, group = 1)) +
  geom_point() + 
  geom_line() +
  xlab('Voltage (mV)') +
  ylab('Mean Normalized INa (pA/pF) MGAT1KO') 

normI_mean <- inner_join(normI_wt_mean, normI_ko_mean, "voltage")
colnames(normI_mean) <- c('voltage', 'WT', 'MGAT1KO')
normI_mean <- normI_mean %>% 
  pivot_longer(c('WT', 'MGAT1KO'), names_to = 'Group', values_to = 'density')
normI_mean %>% 
  ggplot(aes(x = voltage, y = density, group = Group)) +
  geom_point(aes(color = Group)) +
  geom_line(aes(color = Group)) +
  xlab('Voltage (mV)') +
  ylab('Density (pA/pF)')


## Nav Mgat1KO IV Data ----
path <- './MGAT1_Data/Final MGAT1KO Data for Hui/JMCC paper/Nav Currents/IV Recordings/Mgat1KO IV/01-31-2018 Na.xlsx'

# steady-state activation
mgat1ko_ivs <- read_excel(path, sheet = 'IVS Cell ', range = cell_limits(c(1, 1), c(21, 6)))
mgat1ko_ivs %>% 
  filter(mV <= -20) %>% 
  ggplot(aes(x = mV, y = `G/Gmax`, group = 1)) +
  geom_point() +
  geom_line()

# steady-state inactivation
mgat1ko_ssi <- read_excel(path, sheet = 'SSI Cell', range = cell_limits(c(1, 1), c(17, 3)))
mgat1ko_ssi %>% 
  ggplot(aes(x = mV, y = I2, group = 1)) +
  geom_line() +
  geom_point()
mgat1ko_ssi %>% 
  ggplot(aes(x = mV, y = `I2/I2max`, group = 1)) +
  geom_line() +
  geom_point()

# recovery
mgat1ko_rev <- read_excel(path, sheet = 'Recovery -100', range = cell_limits(c(1, 1), c(31, 4)))
mgat1ko_rev %>% 
  ggplot(aes(x = `Time mS`, y = I1, group = 1)) +
  geom_point() +
  geom_line()
mgat1ko_rev %>% 
  ggplot(aes(x = `Time mS`, y = I2, group = 1)) +
  geom_point() +
  geom_line()
mgat1ko_rev %>% 
  ggplot(aes(x = `Time mS`, y = `I2/I1`, group = 1)) +
  geom_point() +
  geom_line()


## Nav WT Data ----
path <- './MGAT1_Data/Final MGAT1KO Data for Hui/JMCC paper/Nav Currents/IV Recordings/WT IV/Nav_WT.xlsx'

ivs_wt <- read_excel(path, sheet = 'IVS')
ssi_wt <- read_excel(path, sheet = 'SSI')
recovery_wt <- read_excel(path, sheet = 'recovery-100')

# aggregate the data
ggmax_mean_wt <- ivs_wt %>% 
  group_by(voltage) %>% 
  summarise(mean_ggmax = mean(`G/Gmax`))

iimax_mean_wt <- ssi_wt %>% 
  group_by(voltage) %>% 
  summarise(mean_iimax = mean(`I2/I2max`))

# plotting
# activation & inactivation

ggmax_mean_wt %>%
  filter(voltage <= -20) %>% 
  ggplot(aes(x = voltage, y = mean_ggmax, group = 1)) +
  geom_point() +
  geom_line() +
  ylab('SS Activation')

# inactivation
iimax_mean_wt %>% 
  ggplot(aes(x = voltage, y = mean_iimax, group = 1)) +
  geom_point() +
  geom_line() +
  ylab('SS Inactivation')
