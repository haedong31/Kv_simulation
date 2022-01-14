library(tidyverse)
library(readxl)


## Fig 3 -----
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
  select(-Cap)

# Fig 3-B
# currents
fig3 %>% 
  select(2:6, group) %>% 
  pivot_longer(colnames(fig3)[2:6], names_to = 'measure') %>% 
  ggplot(aes(x = measure, y = value, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  xlab('')
  ylab('Density (pA/pF)')
# taus
fig3 %>% 
  select(8:9, group) %>% 
  pivot_longer(colnames(fig3)[8:9], names_to = 'measure') %>% 
  ggplot(aes(x = measure, y = value, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_y_continuous(breaks = c(0, 100, 1000, 1250, 1500, 1750, 2000)) +
  xlab('') +
  ylab('Tau Inactivation (ms)')

# Fig 3-C
fig3 %>% 
  select(Cap, group) %>% 
  pivot_longer(Cap, names_to = 'measure') %>% 
  ggplot(aes(x = measure, y = value, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  xlab('Capacitance') +
  ylab('pF') +
  scale_y_continuous(n.breaks = 6)
fig3 %>% 
  select(1:6, group) %>%
  mutate(Peak = Peak * Cap, Iss = Iss * Cap, A1 = A1 * Cap, A2 = A2 * Cap, A3 = A3 *Cap) %>% 
  rename(IPeak = Peak, IKslow1 = A2, IKslow2 = A1, Ito = A3) %>% 
  pivot_longer(2:6, names_to = 'measure') %>% 
  ggplot(aes(x = measure, y = value, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  xlab('')
  ylab('Density (pA/pF)')

  
## Fig 4 -----
fig4_bc <- read_excel('./MGAT1_Data_tidy/JMCC/Nav Currents/FF IV 01-11-2018 - 3-22-19 - final.xlsx',
                      range = cell_limits(c(1, 51), c(21, 59)))

# data pre-processing
fig4_bc_wt <- fig4_bc %>% 
  select(1:5) %>% 
  mutate(group = 'WT')
colnames(fig4_bc_wt) <- c('mV', 'norm_mean', 'mean', 'norm_SEM', 'SEM', 'group')
fig4_bc_ko <- fig4_bc %>% 
  select(c(1, 6:9)) %>% 
  mutate(group = 'KO')
colnames(fig4_bc_ko) <- c('mV', 'norm_mean', 'mean', 'norm_SEM', 'SEM', 'group')
fig4_bc <- bind_rows(fig4_bc_wt, fig4_bc_ko)

# Fig 4-B
fig4_bc %>% 
  filter(mV <= 0) %>% 
  ggplot(aes(x = mV, y = norm_mean, color = group)) +
  geom_point() +
  geom_line() +
  xlab('Voltage (mV)') +
  ylab('Density (pA/pF)')

# Fig 4-C
fig4_bc %>% 
  filter(mV <= 0) %>% 
  ggplot(aes(x = mV, y = mean, color = group)) +
  geom_point() +
  geom_line() +
  xlab('Voltage (mV)') +
  ylab('Amplitude (pA)')

# data pre-processing
fig4_ssa <- read_excel('./MGAT1_Data_tidy/JMCC/Nav Currents/Ina GV MGAT1KO Final.xlsx',
                       range = cell_limits(c(1, 1), c(24, 43)))
fig4_ssi <- read_excel('./MGAT1_Data_tidy/JMCC/Nav Currents/Ina SSI MGAT1KO Final.xlsx',
                       range = cell_limits(c(1, 1), c(26, 45)))

fig4_ssa_wt <- fig4_ssa %>% 
  select(1:20) %>% 
  mutate(group = 'WT')
colnames(fig4_ssa_wt) <- c('ID', 'Cap', 'G/Gmax', 'Slope', 'V1/2', 'Gmax/Cap', 
                           seq(-85, -20, by = 5), 'group')
fig4_ssa_ko <- fig4_ssa %>% 
  select(24:43) %>% 
  slice(1:16) %>% 
  mutate(group = 'KO')
colnames(fig4_ssa_ko) <- c('ID', 'Cap', 'G/Gmax', 'Slope', 'V1/2', 'Gmax/Cap', 
                           seq(-85, -20, by = 5), 'group')

fig4_ssa_2 <- bind_rows(fig4_ssa_wt, fig4_ssa_ko)
fig4_ssa_2 <- fig4_ssa_2 %>% 
  pivot_longer(colnames(fig4_ssa_2)[7:20], names_to = 'mV', values_to = 'ssa') %>% 
  arrange(mV)
fig4_ssa_2$mV <- fig4_ssa_2$mV %>% as.numeric()

fig4_ssi_wt <- fig4_ssi %>% 
  select(2:21) %>% 
  slice(1:25) %>% 
  mutate(group = 'WT')
colnames(fig4_ssi_wt) <- c('ID', 'Imax', 'Slope', 'V1/2', seq(-130, -55, by = 5), 'group')
fig4_ssi_ko <- fig4_ssi %>% 
  select(26:45) %>% 
  slice(1:16) %>% 
  mutate(group = 'KO')
colnames(fig4_ssi_ko) <- c('ID', 'Imax', 'Slope', 'V1/2', seq(-130, -55, by = 5), 'group')

fig4_ssi_2 <- bind_rows(fig4_ssi_wt, fig4_ssi_ko)
fig4_ssi_2 <- fig4_ssi_2 %>% 
  pivot_longer(colnames(fig4_ssi_2)[5:20], names_to = 'mV', values_to = 'ssi') %>% 
  arrange(mV)
fig4_ssi_2$mV <- fig4_ssi_2$mV %>% as.numeric()

# Fig 4-D
fig4_ssa_2 %>% 
  group_by(mV, group) %>% 
  summarise(mean_ssa = mean(ssa)) %>% 
  ggplot(aes(x = mV, y = mean_ssa, color = group)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))
fig4_ssi_2 %>% 
  group_by(mV, group) %>% 
  summarise(mean_ssi = mean(ssi)) %>% 
  ggplot(aes(x = mV, y = mean_ssi, color = group)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))

# Fig 4-G
test <- read_excel('./MGAT1_Data_tidy/JMCC/test.xlsx')
test %>% 
  ggplot(aes(x = `Time mS`, y = `I2/I1`)) +
  geom_point() +
  geom_line()


## Fig 6 -----
Ca_WT <- read_excel('./MGAT1_Data_tidy/JMCC/Ca Imaging 37 Degrees/Ca Imaging MGAT1KO Final.xlsx', 
                    range = cell_limits(c(1, 1), c(85, 11)))
Ca_KO <- read_excel('./MGAT1_Data_tidy/JMCC/Ca Imaging 37 Degrees/Ca Imaging MGAT1KO Final.xlsx', 
                    range = cell_limits(c(1, 13), c(55, 23)))

# Fig 6 C-3
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


## misc -----
fig4_ssa_ko$ID[]
!(fig4_ssa_ko$ID %in% fig4_ssa_wt$ID)
!(fig4_ssa_wt$ID %in% fig4_ssa_ko$ID)

fig4_ssa_ko$ID[!ssi_idx]
fig4_ssa_wt$ID[!ssi_idx]
