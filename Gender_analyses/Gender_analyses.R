### Analyzing correlations between phonology and syntax
### Lemmas within parts of speech

# Load packages
library(ggplot2)
library(data.table)
library(lme4)
library(lmerTest)
library(effects)
library(ggeffects)
library(sjPlot)
library(sjmisc)
library(generics)
library(dplyr)
library(bestNormalize)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(foreach)
library(doParallel)
library(DescTools)


## Creating a single file with all languages

# Load all files
allpairs = data.frame()
files <- list.files(path="/Users/pgr/Documents/Dissertation/Python_scripts/UDT_gender_paired/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
i = 1
for (file in files){
  pairs <- read.csv(file, header=TRUE)
  print(nrow(pairs))
  allpairs = rbind(allpairs, pairs)
  print(i)
  i = i + 1
}

# Keep only the columns we need
allpairs = allpairs[c('lang', 'lemma_x', 'gender_x', 'frequency_x', 'p_m_freq_x', 'dispersion_x', 'lemma_y', 'gender_y', 'frequency_y', 'p_m_freq_y', 'dispersion_y', 'lev_dist', 'sem_dist', 'syn_dist', 'syn_left_dist', 'syn_right_dist', 'syn_dep_dist', 'syn_head_dist', 'gender_sameness')]

# Save allpairs as csv
fwrite(allpairs, "/Users/pgr/Documents/Dissertation/Gender_analyses/Gender_pairs.csv")



### Loading combined data and subsetting

# Load file
p = fread("/Users/pgr/Documents/Dissertation/Gender_analyses/Gender_pairs.csv")

# Setting variable types
p$lang = as.factor(p$lang)
p$gender_x = as.factor(p$gender_x)
p$gender_y = as.factor(p$gender_y)
p$frequency_x = as.numeric(p$frequency_x)
p$frequency_y = as.numeric(p$frequency_y)
p$lev_dist = as.numeric(p$lev_dist)
p$gender_sameness = as.numeric(p$gender_sameness)


### Explore data

# Sample randomly for easier exploration and plotting
s = p %>% group_by(lang) %>% sample_n(1000)


## Languages

# Let's see how much data is coming from each language
langs = vector(mode = 'character')
sizes = vector(mode = 'numeric')
for (level in levels(p$lang)) {
  size = nrow(p[p$lang == level,])
  sizes = c(sizes, size)
  langs = c(langs, level)
}
langsizes = data.frame('langs' = langs, 'sizes' = sizes)
langsizes
sort(sizes)

# Let's get rid of any language that doesn't have at least 10,000 pairs
good_langs = langsizes[langsizes$sizes >= 10000,]
p = p[p$lang %in% good_langs$langs,]

# Drop levels that no longer exist in the data
p$lang = factor(p$lang)

# Save this meta info
write.csv(langsizes, "/Users/pgr/Documents/Dissertation/Gender_analyses/Langs_and_sizes.csv")


## Genders
table(s$lang, s$gender_x)
table(s$lang, s$gender_y)

# Get rid of underattested genders
p = p[!(p$gender_x == 'Com,Neut' | p$gender_y == 'Com,Neut' | p$gender_x == 'Fem,Masc' | p$gender_y == 'Fem,Masc')]
p$gender_x[(p$gender_x == 'Masc1' | p$gender_x == 'Masc2' | p$gender_x == 'Masc3')] = 'Masc'
p$gender_y[(p$gender_y == 'Masc1' | p$gender_y == 'Masc2' | p$gender_y == 'Masc3')] = 'Masc'
p$gender_x = factor(p$gender_x)
p$gender_y = factor(p$gender_y)


## Distances of 0
nrow(p[p$syn_dist <= 0,])
nrow(p[p$sem_dist <= 0,])

# We're going to get rid of distances of 0
p = p[!(p$syn_dist <= 0 | p$sem_dist <= 0)]


## Add variable indicating number of genders in the language
for (l in levels(p$lang)){
  print("Starting language")
  print(l)
  temp = p[p$lang == l,]
  print("Got temp")
  genders = length(unique(temp$gender_x))
  print("Got length:")
  print(genders)
  p$genders[p$lang == l] = genders
}


## Save allpairs as csv
fwrite(p, "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Gender_pairs.csv")



### Load file
p = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Gender_pairs.csv")

# Setting variable types
p$lang = as.factor(p$lang)
p$gender_x = as.factor(p$gender_x)
p$gender_y = as.factor(p$gender_y)
p$frequency_x = as.numeric(p$frequency_x)
p$frequency_y = as.numeric(p$frequency_y)
p$lev_dist = as.numeric(p$lev_dist)
p$gender_sameness = as.numeric(p$gender_sameness)
p$genders = as.factor(p$genders)

# Subset to just the columns we need for a particular purpose
p = subset(p, select = c(lang, lev_dist, sem_dist, syn_dist, genders, gender_sameness))



### Checking correlations for each language

# Getting size of data and correlations for each 
all_data = data.frame(matrix(ncol=8, nrow=0, dimnames=list(NULL, c("langs", "rows", "sem_syn", "lev_syn", "sem_lev", "sem_gen", "lev_gen", "syn_gen"))))
for (l in levels(p$lang)){
  sub = p[p$lang == l,]
  lang_data = data.frame(langs = l,
            rows = nrow(sub),
            sem_syn = cor(sub$sem_dist, sub$syn_dist),
            lev_syn = cor(sub$lev_dist, sub$syn_dist),
            sem_lev = cor(sub$sem_dist, sub$lev_dist),
            sem_gen = cor(sub$sem_dist, sub$gender_sameness),
            lev_gen = cor(sub$lev_dist, sub$gender_sameness),
            syn_gen = cor(sub$syn_dist, sub$gender_sameness)
  )
  all_data = rbind(all_data, lang_data)
  print(l)
}

# See the overall patterns
par(mfrow = c(3,3))
hist(all_data$sem_syn)
hist(all_data$lev_syn)
hist(all_data$sem_lev)
hist(all_data$sem_gen)
hist(all_data$lev_gen)
hist(all_data$syn_gen)
par(mfrow = c(1,1))

# Save this information
fwrite(all_data, "/Users/pgr/Documents/Dissertation/Gender_analyses/All_data.csv")



## Checking different syntax measures
for (l in levels(p$lang)){
  sub = p[p$lang == l,]
  stuff = c(cor(sub$syn_dist, sub$gender_sameness),
            cor(sub$syn_left_dist, sub$gender_sameness, use = "pairwise.complete.obs"),
            cor(sub$syn_head_dist, sub$gender_sameness, use = "pairwise.complete.obs")
  )
  print(l)
  print(stuff)
}

# Find out how many NAs there are for the special syntax measures
sum(is.na(p$syn_right_dist))



### Simulation approach to correlation testing


## Get all secondary correlations
all_data = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/All_data.csv")
hist(all_data$lev_syn)
hist(all_data$sem_syn)
hist(all_data$lev_gen)
hist(all_data$sem_gen)
data_long <- gather(all_data, variables, value, c(sem_syn, lev_syn, sem_gen, lev_gen), factor_key=TRUE)
data_long$first = c(rep("Semantic distance", 32), rep("Orthographic distance", 32), rep("Semantic distance", 32), rep("Orthographic distance", 32))
data_long$second = c(rep("Syntactic distance", 64), rep("Gender sameness", 64))
data_long$first = factor(data_long$first, levels = c("Semantic distance", "Orthographic distance"))
data_long$second = factor(data_long$second, levels = c("Syntactic distance", "Gender sameness"))

# Plot it
tiff("correlations.tiff", units="in", width=8, height=6, res=300)
ggplot(data_long, aes(x = value)) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  geom_histogram(fill = "darkred", binwidth = 0.04, boundary = 0) +
  facet_grid(rows = vars(second), cols = vars(first)) +
  scale_x_continuous(limits = c(-0.45, 0.45), breaks=c(-0.25, 0, 0.25), labels=c("\u20130.25","0.00","0.25")) +
  scale_y_continuous(limits = c(0, 13)) +
  labs(x = "Correlation", y = "Number of languages") +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = 0), linetype = 2)
dev.off()



## OPTIONAL: Restricting to comparisons of equal word lengths

p$length_x = nchar(p$lemma_x)
p$length_y = nchar(p$lemma_y)
p = p[p$length_x == p$length_y,]


## OPTIONAL: Changing syn_left_dist to syn_dist to test only leftward distributions
p$real_syn = p$syn_dist
p$syn_dist = p$syn_left_dist
p = p[!(p$syn_dist <= 0)]

  
## Parallelize the inner loop

# We're going to do this in parallel
cores = detectCores()
cl = makeCluster(cores - 1)
registerDoParallel(cl)

# Define dataframe for simulated correlations
full_df = data.frame(matrix(ncol=2, nrow=0, dimnames=list(NULL, c("langs", "sim_cors"))))

# Define function for getting simulated correlations for a single language
get_lang_sims = function(test, tolerance, sem_target, lev_target, sem_target_range, lev_target_range){
  
  # Subset to just 10,000 rows
  subtest = sample_n(test, 10000)
  
  # Randomly permute syntactic distances
  subtest$syn_dist = sample(subtest$syn_dist)
  
  # Get permuted correlations
  sem_cor = cor(subtest$sem_dist, subtest$syn_dist)
  lev_cor = cor(subtest$lev_dist, subtest$syn_dist)
  
  # Loop to indicate correlations aren't matched yet
  while ((sem_cor < sem_target_range[1]) | (sem_cor > sem_target_range[2]) | (lev_cor < lev_target_range[1]) | (lev_cor > lev_target_range[2])){
    
    # Find out difference between target and current correlations (positive means we need the correlation value to go up)
    sem_target_diff = sem_target - sem_cor
    lev_target_diff = lev_target - lev_cor
    
    # Loop to try different samples
    i = 0
    while (i == 0){
      
      # Sample two random rows
      rows = sample(1:nrow(subtest), size = 2)
      
      # Identify syntactic distance values
      syn_1 = subtest$syn_dist[rows[1]]
      syn_2 = subtest$syn_dist[rows[2]]
      
      # Find out differences in values between rows
      sem_diff = subtest$sem_dist[rows[1]] - subtest$sem_dist[rows[2]]
      lev_diff = subtest$lev_dist[rows[1]] - subtest$lev_dist[rows[2]]
      
      # Determine whether they should be switched
      if (((syn_1 > syn_2) & (sign(sem_target_diff) == sign(-sem_diff)) & (sign(lev_target_diff) == sign(-lev_diff))) |
          ((syn_1 < syn_2) & (sign(sem_target_diff) == sign(sem_diff)) & (sign(lev_target_diff) == sign(lev_diff)))){
        
        # Switch syntactic distance values
        subtest$syn_dist[rows[1]] = syn_2
        subtest$syn_dist[rows[2]] = syn_1
        
        # Get new correlations
        sem_cor = cor(subtest$sem_dist, subtest$syn_dist)
        lev_cor = cor(subtest$lev_dist, subtest$syn_dist)
        
        # Change i to exit loop
        i = 1
      }
    }
  }
  return(cor(subtest$syn_dist, subtest$gender_sameness))
}

# Loop through languages
for (l in levels(p$lang)){
  
  # Get subset of dataframe corresponding to the language
  test = p[p$lang == l,]
  
  # Define true correlation between syntax and gender sameness
  true_cor = cor(test$syn_dist, test$gender_sameness)
  
  # Get target correlations between syntax and both semantics and levenshtein
  sem_target = cor(test$sem_dist, test$syn_dist)
  lev_target = cor(test$lev_dist, test$syn_dist)
  
  # Set tolerance
  tolerance = .001
  
  # Get target ranges
  sem_target_range = c(sem_target - tolerance, sem_target + tolerance)
  lev_target_range = c(lev_target - tolerance, lev_target + tolerance)
  
  # Apply function in parallel
  x <- foreach(k = 1:10000, .combine = 'c', .packages='dplyr') %dopar% {
    get_lang_sims(test, tolerance, sem_target, lev_target, sem_target_range, lev_target_range)
  }
  
  # Progress report and combining language data
  print("Done with another language:")
  print(l)
  langs = rep(l, 10000)
  df = data.frame(langs, x)
  full_df = rbind(full_df, df)
}

# Close cluster
stopCluster(cl)

# Save allpairs as csv
fwrite(full_df, "/Users/pgr/Documents/Dissertation/Gender_analyses/Sim_cors_4.csv")


## Collect true correlations in dataframe

# Define dataframe for real correlations
all_df = data.frame(matrix(ncol=2, nrow=0, dimnames=list(NULL, c("langs", "x"))))

# Loop through to get correlations
for (l in levels(p$lang)){
  
  # Get subset of dataframe corresponding to the language
  test = p[p$lang == l,]
  
  # Define true correlation between syntax and gender sameness
  real_cor = cor(test$syn_dist, test$gender_sameness)
  
  # Combine into dataframe
  df = data.frame(langs = l, x = real_cor)
  
  # Merge with previous dataframe
  all_df = rbind(all_df, df)
  
  # Progress
  print("Done with another language")
}

all_df

# Save real correlations
fwrite(all_df, "/Users/pgr/Documents/Dissertation/Gender_analyses/Same_length_real_cors.csv")


## Getting significance values for each language

# Load real and simulated correlations
cors = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Sim_cors.csv")
cors$langs = as.factor(cors$langs)
real = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Real_cors.csv")
p_values = c()

# Loop through each language, calculating the proportion of simulated correlations that are greater than or equal to the real one
for (l in levels(sim_cors$langs)){
  sub = sim_cors[sim_cors$langs == l,]
  p_value = (sum(sub$x >= real_cors$x[real_cors$langs == l]) + 1)/10001
  p_values = c(p_values, p_value)
  print(l)
}

significance = data.frame(langs = levels(sim_cors$langs), p_values = p_values)
significance

sum(significance$p_values < 0.05)
fwrite(significance, "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/P_values.csv")


## Plotting permutation results

# Prep data for plotting
all_data = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/All_data.csv")
cors$language <- all_data$languages[match(cors$langs, all_data$langs)]
cors$family <- all_data$family[match(cors$langs, all_data$langs)]
real$language <- all_data$languages[match(real$langs, all_data$langs)]
cors = cors[order(cors$family)]
real$y = 200
p_values = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/P_values.csv")
real$p_values <- p_values$significance[match(real$langs, p_values$langs)]

# Generate and save plot
tiff("simulations.tiff", units="in", width=10, height=6, res=300)
ggplot(cors, aes(x = x)) +
  scale_fill_viridis(discrete = TRUE, name = "Family") +
  geom_histogram(aes(fill = family), bins = 50) +
  geom_point(data = real, aes(y = y), color = "red", size = 2) +
  facet_wrap(~factor(language, levels=c('Lithuanian', 'Latvian', 'Welsh', 'Irish',
                                        'Gaelic', 'Belarusian', 'Russian', 'Ukrainian',
                                        'Greek', 'Hindi', 'Urdu', 'Danish',
                                        'Icelandic', 'Norwegian', 'Swedish', 'Catalan',
                                        'Spanish', 'French', 'Italian', 'Latin',
                                        'Portuguese', 'Romanian', 'Hebrew', 'Bulgarian',
                                        'Croatian', 'Slovenian', 'Serbian', 'German',
                                        'Dutch', 'Czech', 'Polish', 'Slovak')), ncol = 4)+
  labs(x = "Correlation between syntactic distance and gender sameness", y = "Number of simulations") +
  scale_x_continuous(limits = c(-0.15, 0.15), breaks=c(-0.1, 0.0, 0.1), labels=c("\u20130.1","0.0","0.1")) +
  #scale_y_continuous(limits = c(0, 3000)) +
  geom_text(data = real, aes(label = p_values), size = 3.5, x = 0.1, y = 1500)
dev.off()



### Modeling

# Re-sample after data processing (using seed 77)
set.seed(77)
s = p %>% group_by(lang) %>% sample_n(10000)

# OPTIONAL: Run regression on entire dataset rather than sample
s = p
rm(p)

## Language and family names
all_data = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/All_data.csv")
s$language <- all_data$languages[match(s$lang, all_data$langs)]
s$family <- all_data$family[match(s$lang, all_data$langs)]
s$language = as.factor(s$language)
s$family = as.factor(s$family)


## Let's check out each variable

## syn_dist
ggplot(s, aes(x = syn_dist)) +
  geom_density(aes(group = lang, color = lang))

# Let's normalize and scale within each language
s = s %>%
  group_by(lang) %>%
  mutate(syn_dist_bc = boxcox(syn_dist, standardize = TRUE)$x.t)
summary(s$syn_dist_bc)
ggplot(s, aes(x = syn_dist_bc)) +
  geom_density(aes(group = lang, color = lang))


## lev_dist
ggplot(s, aes(x = lev_dist)) +
  geom_density(aes(group = lang, color = lang))

# Let's normalize and scale within each language
s = s %>%
  group_by(lang) %>%
  mutate(lev_dist_bc = boxcox(lev_dist, standardize = TRUE)$x.t)
summary(s$lev_dist_bc)
ggplot(s, aes(x = lev_dist_bc)) +
  geom_density(aes(group = lang, color = lang))


## sem_dist
ggplot(s, aes(x = sem_dist)) +
  geom_density(aes(group = lang, color = lang))

# Let's normalize and scale within each language
s = s %>%
  group_by(lang) %>%
  mutate(sem_dist_bc = boxcox(sem_dist, standardize = TRUE)$x.t)
summary(s$sem_dist_bc)
ggplot(s, aes(x = sem_dist_bc)) +
  geom_density(aes(group = lang, color = lang))

# Exploring correlation between intercepts and slopes
ggplot(s, aes(x = syn_dist_bc, y = gender_sameness)) +
  geom_smooth(aes(group = language), method = "lm", se = FALSE)
ggplot(s, aes(x = sem_dist_bc, y = gender_sameness)) +
  geom_smooth(aes(group = language), method = "lm", se = FALSE)
ggplot(s, aes(x = lev_dist_bc, y = gender_sameness)) +
  geom_smooth(aes(group = language), method = "lm", se = FALSE)
# These suggest to me that I should do uncorrelated slopes and intercepts, along with theoretical considerations


## Save allpairs as csv
fwrite(s, "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Gender_pairs_transformed_sampled.csv")


## Load allpairs
s = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Gender_pairs_transformed_sampled.csv")
s$language = as.factor(s$language)
s$family = as.factor(s$family)
s$genders = as.factor(s$genders)


## Models
m5 = glmer(gender_sameness ~ lev_dist_bc + sem_dist_bc + syn_dist_bc +
             (1|family/language) + 
             (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language), 
           s, family = "binomial")
m5

m6 = glmer(gender_sameness ~ sem_dist_bc + syn_dist_bc +
             (1|family/language) + 
             (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language),
           s, family = "binomial")
m6

m7 = glmer(gender_sameness ~ lev_dist_bc + syn_dist_bc +
             (1|family/language) + 
             (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language),
           s, family = "binomial")
m7

m8 = glmer(gender_sameness ~ lev_dist_bc + sem_dist_bc +
             (1|family/language) + 
             (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language),
           s, family = "binomial")
m8

# Compare models
anova(m5, m6)
anova(m5, m7)
anova(m5, m8)
summary(m5)

# Save models
saveRDS(m5, file = "/Users/pgr/Documents/Dissertation/Gender_analyses/Model_full.rds")
saveRDS(m6, file = "/Users/pgr/Documents/Dissertation/Gender_analyses/Model_wo_lev.rds")
saveRDS(m7, file = "/Users/pgr/Documents/Dissertation/Gender_analyses/Model_wo_sem.rds")
saveRDS(m8, file = "/Users/pgr/Documents/Dissertation/Gender_analyses/Model_wo_syn.rds")

# Load models
m5 = readRDS(file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Model_full.rds")
m6 = readRDS(file = "/Users/pgr/Documents/Dissertation/Gender_analyses/Model_wo_lev.rds")
m7 = readRDS(file = "/Users/pgr/Documents/Dissertation/Gender_analyses/Model_wo_sem.rds")
m8 = readRDS(file = "/Users/pgr/Documents/Dissertation/Gender_analyses/Model_wo_syn.rds")

# Check residuals
rs = DHARMa::simulateResiduals(m5)
plot(rs)
DHARMa::plotResiduals(rs, form = s$syn_dist_bc, quantreg = T)
DHARMa::plotResiduals(rs, form = s$sem_dist_bc, quantreg = T)
DHARMa::plotResiduals(rs, form = s$lev_dist_bc, quantreg = T)


## Plot results

# Quick look
plot(allEffects(m10))

# Get predictions to do it manually
syn = as.data.frame(predictorEffect("syn_dist_bc", m10, focal.levels=100))
syn$effect = "syn_dist_bc"
colnames(syn)[1] = "x"
sem = as.data.frame(predictorEffect("sem_dist_bc", m10, focal.levels=100))
sem$effect = "sem_dist_bc"
colnames(sem)[1] = "x"
lev = as.data.frame(predictorEffect("lev_dist_bc", m10, focal.levels=100))
lev$effect = "lev_dist_bc"
colnames(lev)[1] = "x"
preds = rbind(syn, sem, lev)

# Create long form of distances for rugs
s_long = gather(s, effect, value, syn_dist_bc:sem_dist_bc, factor_key = TRUE)

# Create labels for facets
facet_labels = c(syn_dist_bc = "Syntactic", sem_dist_bc = "Semantic", lev_dist_bc = "Orthographic")

# Plot and save
tiff("gender_glmer_effects.tiff", units="in", width=10, height=4, res=300)
ggplot(preds, aes(x = x)) +
  geom_smooth(aes(color = effect, y = fit), method = "lm") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = effect), alpha = 0.3) +
  facet_grid(cols = vars(effect), labeller = as_labeller(facet_labels)) +
  scale_color_manual(values = c("#06ffcf", "#00a384", "#004134")) +
  scale_fill_manual(values = c("#06ffcf", "#00a384", "#004134")) +
  labs(x = "Scaled distance", y = "Probability of gender sameness") +
  coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(0.15, 0.8)) +
  scale_x_continuous(breaks=c(-2, 0, 2), labels=c("\u20132","0","2")) +
  theme(legend.position = "none") +
  geom_rug(data = s_long, aes(x = value), alpha = 0.005)
dev.off()


## Make data from which to predict
newdata = expand.grid(lang = levels(s$lang), lev_dist_bc = seq(-5, 5, 0.2), sem_dist_bc = seq(-5, 5, 0.2), syn_dist_bc = seq(-5, 5, 0.2))
all_data = fread("/Users/pgr/Documents/Dissertation/Gender_analyses/All_data.csv")
newdata$family <- all_data$family[match(newdata$lang, all_data$langs)]
newdata$family = as.factor(newdata$family)

# Make preds
newdata$preds = predict(m1, newdata, type = 'response')

# Plot each effect for each language
ggplot(newdata, aes(x = syn_dist_bc, y = preds)) +
  geom_smooth(aes(color = lang), method = "lm", se = FALSE) +
  scale_color_viridis(discrete = TRUE, option = "turbo")
ggplot(newdata, aes(x = sem_dist_bc, y = preds)) +
  geom_smooth(aes(color = lang), method = "lm", se = FALSE) +
  scale_color_viridis(discrete = TRUE, option = "turbo")
ggplot(newdata, aes(x = lev_dist_bc, y = preds)) +
  geom_smooth(aes(color = lang), method = "lm", se = FALSE) +
  scale_color_viridis(discrete = TRUE, option = "turbo")


## New model with number of genders
m10 = glmer(gender_sameness ~ lev_dist_bc + sem_dist_bc + syn_dist_bc + genders +
             (1|family/language) + 
             (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language) + (0 + genders|family/language), 
           s, family = "binomial")
m10
anova(m5, m10)
# Save model
saveRDS(m10, file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Model_w_genders.rds")
m10 = readRDS("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Model_w_genders.rds")

## Comparison models

# No lev
m20 = glmer(gender_sameness ~ sem_dist_bc + syn_dist_bc + genders +
              (1|family/language) + 
              (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language) + (0 + genders|family/language), 
            s, family = "binomial")
m20
saveRDS(m20, file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Upd_model_wo_lev.rds")
anova(m10, m20)

# No sem
m21 = glmer(gender_sameness ~ lev_dist_bc + syn_dist_bc + genders +
              (1|family/language) + 
              (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language) + (0 + genders|family/language), 
            s, family = "binomial")
m21
saveRDS(m21, file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Upd_model_wo_sem.rds")
anova(m10, m21)

# No syn
m22 = glmer(gender_sameness ~ lev_dist_bc + sem_dist_bc + genders +
              (1|family/language) + 
              (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language) + (0 + genders|family/language), 
            s, family = "binomial")
m22
saveRDS(m22, file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Upd_model_wo_syn.rds")
anova(m10, m22)

# No genders
m23 = glmer(gender_sameness ~ lev_dist_bc + sem_dist_bc + syn_dist_bc +
              (1|family/language) + 
              (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language) + (0 + genders|family/language), 
            s, family = "binomial")
m23
saveRDS(m23, file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Upd_model_wo_genders.rds")
anova(m10, m23)


## Trying a correlation between genders and syn_dist_bc to see if its licensed
m11 = glmer(gender_sameness ~ lev_dist_bc + sem_dist_bc + syn_dist_bc + genders + genders*syn_dist_bc +
              (1|family/language) + 
              (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language) + (0 + genders|family/language) + (0 + genders*syn_dist_bc|family/language), 
            s, family = "binomial")
m11
saveRDS(m11, file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Model_w_genders_&_interaction.rds")

m12 = glmer(gender_sameness ~ lev_dist_bc + sem_dist_bc + syn_dist_bc + genders +
              (1|family/language) + 
              (0 + lev_dist_bc|family/language) + (0 + sem_dist_bc|family/language) + (0 + syn_dist_bc|family/language) + (0 + genders|family/language) + (0 + genders*syn_dist_bc|family/language), 
            s, family = "binomial")
m12
anova(m11, m12)
saveRDS(m12, file = "/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/Model_w_genders_&_interaction_comparison.rds")


## Get confidence intervals of model
confint(m10, method = "Wald")


## Checking correlation between coefficients and corpus size
coefs = coef(m10)$language
coefs$language = unlist(lapply(strsplit(row.names(coefs), ":"), '[', 1))

langs = vector(mode = 'character')
sizes = vector(mode = 'numeric')
for (level in levels(p$lang)) {
  size = nrow(p[p$lang == level,])
  sizes = c(sizes, size)
  langs = c(langs, level)
  print("Done with another one...")
}
langsizes = data.frame('langs' = langs, 'sizes' = sizes)

all_data = fread("/Users/pgr/Documents/Dissertation/syntactic_information_in_the_lexicon/Gender_analyses/All_data.csv")
langsizes$language <- all_data$languages[match(langsizes$langs, all_data$langs)]
langsizes$language = as.factor(langsizes$language)

coefs$sizes = langsizes$sizes[match(coefs$language, langsizes$language)]
coefs$sizes_bc = boxcox(coefs$sizes, standardize = TRUE)$x.t

# Plot relationship
ggplot(coefs, aes(x = sizes_bc, y = syn_dist_bc)) +
  geom_point() +
  geom_smooth(method = "lm")
ggplot(coefs, aes(x = sizes_bc, y = syn_dist_bc)) +
  geom_point() +
  geom_smooth()
# Visually, it doesn't look like much of a correlation

# Get correlation
cor(coefs$sizes_bc, coefs$syn_dist_bc)



### Checking correlations for each language

# Checking various correlations
for (l in levels(p$lang)){
  sub = p[p$lang == l,]
  stuff = c(cor(sub$sem_dist, sub$syn_dist),
            cor(sub$lev_dist, sub$syn_dist),
            cor(sub$sem_dist, sub$gender_sameness),
            cor(sub$lev_dist, sub$gender_sameness),
            cor(sub$syn_dist, sub$gender_sameness)
  )
  print(l)
  print(stuff)
}

# Checking different syntax measures
for (l in levels(p$lang)){
  sub = p[p$lang == l,]
  stuff = c(cor(sub$syn_dist, sub$gender_sameness),
            cor(sub$syn_left_dist, sub$gender_sameness, use = "pairwise.complete.obs"),
            cor(sub$syn_head_dist, sub$gender_sameness, use = "pairwise.complete.obs")
  )
  print(l)
  print(stuff)
}

# Find out how many NAs there are for the special syntax measures
sum(is.na(p$syn_right_dist))




## Checking Arabic (first language that missed the size cut-off)
arabic = fread("/Users/pgr/Documents/Dissertation/Python_scripts/UDT_gender_paired/UDT_gender_paired-ar.csv")
plot(as.factor(arabic$gender_sameness), arabic$syn_dist)
cor(arabic$lev_dist, arabic$gender_sameness)




### Generating some Spanish examples
span = p[p$lang == "es",]
span = span[(frequency_x > 24) & (frequency_y > 24),]
summary(span$syn_dist)
summary(span$sem_dist)


# Close syntactically and far semantically
span[(span$syn_dist < 0.1) & (span$sem_dist > 0.9),]


# Far syntactically and close semantically
span[(span$syn_dist > 0.5) & (span$sem_dist < 0.6),]
# oro (gold, 98) vs medalla (medal, 72) [7 0.4961874 0.5304857]
# censura (censorship) vs amenaza (threat)
# embargo vs duda (doubt)
# cabida (capacity/space) vs hueco (hole/opening)/cuenta (count)
# doña vs. esposa
# municipio (town/municipality) vs cantón (district)

# Let's find other words to pair with oro
oro = span[(span$lemma_x == "oro") | (span$lemma_y == "oro"),]
oro[(oro$syn_dist < 0.2) & (oro$sem_dist > 0.8),]
# oro (gold, 98) vs paz (peace, 132) [3, 0.8472961 0.1944767]

# Summary of syntactic distance in Spanish: 0.06656 0.23514 0.27768 0.28995 0.32916 0.81259
# Summary of semantic distance in Spanish: 0.04228 0.76553 0.82537 0.81647 0.87745 1.14353


## Plot Spanish example syntactic vectors

# Generate variables
dependencies = c("dep_nsubj", "head_det", "head_nmod", "dep_nmod", "head_case", "head_amod", "head_conj", "dep_obj", "head_numm", "dep_obl")
dependencies = c(dependencies, dependencies, dependencies)
probs = c(0.027, 0.212, 0.235, 0.02, 0.054, 0.054, 0.02, 0.17, 0.035, 0.027, 0.013, 0.062, 0.013, 0.348, 0.357, 0.04, 0.031, 0.04, 0.009, 0.005, 0.018, 0.166, 0.018, 0.285, 0.324, 0.021, 0.042, 0.057, 0, 0.021)
freqs = c(7, 55, 61, 5, 14, 14, 5, 44, 9, 7, 3, 14, 3, 79, 81, 9, 7, 9, 2, 1, 6, 56, 6, 96, 109, 7, 14, 19, 0, 7)
lemma = c(rep("medalla", 10), rep("oro", 10), rep("paz", 10))

# Create dataframe
span = data.frame(dependencies = dependencies, probs = probs, freqs = freqs, lemma = lemma)
span$dependencies = as.factor(span$dependencies)
span$lemma = as.factor(span$lemma)
fwrite(span, "/Users/pgr/Documents/Dissertation/Gender_analyses/Span_examples.csv")

# Generate and save plot
tiff("span_examples_classic.tiff", units="in", width=10, height=4, res=300)
ggplot(span, aes(y = probs, x = dependencies)) +
  geom_col(aes(fill = lemma)) +
  facet_grid(vars(rows = lemma)) +
  labs(x = "Syntactic dependency role (head or dependent) and relation", y = "Probability") +
  scale_fill_manual(values = c("#13306dff", "#90548bff", "#f68f46ff"), name = "Lemma") +
  #geom_text(aes(label = freqs), vjust = -0.2, size = 3.5) +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()



