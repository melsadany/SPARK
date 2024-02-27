################################################################################
#            COWAT analysis in spark to replicate findings from RPOE           #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
################################################################################
################################################################################
# read the demographics and remapping
demo <- read_csv("data/raw/ID_mapping_and_demographics.csv") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("NASD", "ASD")))
# read the COWAT
cowat <- read_csv("data/raw/COWAT_merged_word_level_data.qc.csv")
cowat.clean <- cowat %>%
  inner_join(demo) %>%
  distinct(spid, prompt, word, .keep_all = T)
write_rds(cowat.clean, "data/derivatives/cowat-clean.rds")
################################################################################
################################################################################
# read the IQ data
spark.iq <- read_csv("/sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v10_2023-07-17/iq-2023-07-21.csv")
spark.predicted.iq <- read_csv("/sdata/Simons/SPARK/DATA/phenotypes/SPARK_collection_v10_2023-07-17/predicted_iq_experimental-2023-07-21.csv")
iq.clean <- spark.iq %>%
  filter(subject_sp_id %in% cowat.clean$spid) %>%
  select(spid = subject_sp_id, fsiq, fsiq_score, viq_score, nviq_score) %>%
  drop_na(fsiq_score)
write_rds(iq.clean, "data/derivatives/iq-clean.rds")
# education
edu <- read_csv("data/raw/educational_attainment_phenotypes.csv")
# pgs
pgs <- read_csv("/wdata/lcasten/spark/prs/HapMap3_plus/LDPred2-inf-full/gathered_LDPred2-inf_pc_corrected_long.csv")
pgs.clean <- pgs %>%
  filter(grepl("cog|attainment|EA-|CP-|subcortical|school|conscientiousness|SurfaceArea|genlang|ADHD-PGC-2019", pgs_name)) %>%
  filter(IID %in% cowat.clean$spid) %>%
  pivot_wider(names_from = "pgs_name", values_from = pgs_pc_corrected, id_cols = IID) %>%
  drop_na() %>%
  rename(spid=IID)
write_rds(pgs.clean, "data/derivatives/pgs-clean.rds")
################################################################################
################################################################################
################################################################################
# get embeddings from text package
library(text)
emb.text <- textEmbed(unique(cowat.clean$word))
emb.text.m <- emb.text$texts$texts
emb.text.m <- left_join(cowat.clean[,2:ncol(cowat.clean)], 
                        cbind(word = unique(cowat.clean$word), 
                              emb.text.m))
save(emb.text.m, file = "data/derivatives/word-embedding-from-text-package.rda")
################################################################################
################################################################################
# read embeddings and keep samples of interest
load("data/derivatives/word-embedding-from-text-package.rda")
emb.iq <- emb.text.m %>%
  filter(spid %in% iq.clean$spid)
emb.edu <- emb.text.m %>%
  filter(spid %in% edu$spid)
emb.pgs <- emb.text.m %>%
  filter(spid %in% pgs.clean$spid)
rm(emb.text.m);gc()
# choose which category of these 3 above, you'd like to keep
emb.text.m <- emb.pgs %>%
  distinct(spid, prompt, word, .keep_all = T)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# start of analysis
#####
# look at categories of words said
#####
cowat.analyzed <- cbind(cowat.clean,
                     nrc = syuzhet::get_nrc_sentiment(cowat.clean$word),
                     sentimentr::profanity(cowat.clean$word)%>%select(profanity_count),
                     lingmatch = lingmatch::lma_meta(cowat.clean$word),
                     lingmatch = lingmatch::lma_termcat(cowat.clean$word)) %>%
  select(-c(lingmatch.words, lingmatch.unique_words, lingmatch.clauses, lingmatch.sentences, 
            lingmatch.words_per_clause,
            lingmatch.words_per_sentence, lingmatch.characters_per_word, lingmatch.syllables_per_word,
            lingmatch.type_token_ratio))
write_rds(cowat.analyzed, "data/derivatives/cowat-analyzed.rds", compress = "gz")
#### read if done before
cowat.analyzed <- read_rds("data/derivatives/cowat-analyzed.rds") %>%
  filter(spid %in% emb.text.m$spid)
#####
# get words count per participant
#####
word.count <- cowat.analyzed[,2:11] %>%
  group_by(spid) %>%
  dplyr::summarise(word_count = n())
#####
# number of characters
#####
chr.wise <- cowat.analyzed %>%
  pivot_longer(cols = c(lingmatch.characters, lingmatch.syllables, lingmatch.reading_grade), 
               names_to = "cat1") %>%
  group_by(spid, cat1) %>%
  dplyr::summarise(avg = mean(value, na.omit = T)) %>%
  pivot_wider(names_from = "cat1", values_from = "avg", id_cols = "spid")
#####
# ratio of words category being said
#####
word.wise <- cowat.analyzed %>%
  pivot_longer(cols = c(starts_with("nrc"), profanity_count, 
                        lingmatch.sixltr, lingmatch.ppron, lingmatch.ipron, lingmatch.adverb, 
                        lingmatch.conj, lingmatch.auxverb, lingmatch.prep, lingmatch.negate, lingmatch.quant), 
               names_to = "cat2") %>%
  group_by(spid, cat2) %>%
  dplyr::summarise(count = sum(value)) %>%
  left_join(word.count) %>%
  mutate(cat_ratio = count / word_count) %>%
  pivot_wider(names_from = "cat2", values_from = "cat_ratio", id_cols = "spid")
#####
# get average of waiting time between consecutive pairs of words
#####
wait.time <- cowat.analyzed %>%
  drop_na(word) %>%
  mutate(word = tolower(word)) %>%
  select(spid, word, prompt, start_time, end_time) %>%
  distinct(spid, prompt, word, .keep_all = T) %>% # only keep unique words and drop repeated by the same participant in the same task/word
  drop_na(start_time, end_time) # drop the words with no timestamps
wait.time2 <- cbind(wait.time[-nrow(wait.time),]%>%select(1:3),
                    w1_index = rownames(wait.time)[-nrow(wait.time)],
                    w2_index = rownames(wait.time)[-1],
                    wait=wait.time$start_time[-1] - wait.time$end_time[-nrow(wait.time)]) %>%
  filter(wait>0) %>%
  group_by(spid) %>%
  dplyr::summarise(avg_wait = mean(wait))
#####
# save summ for PS-VC
#####
ttt <- inner_join(inner_join(chr.wise, word.wise),
                  inner_join(word.count, wait.time2))
write_rds(ttt, file = "data/derivatives/cowat-pgs-summary-data.rds", compress = "gz")
###############
###############
###############
# get the correlation between PGS and word count/ums count/language features from PS-VC audio
###############
###############
m124 <- inner_join(pgs.clean, 
                   inner_join(inner_join(chr.wise, word.wise),
                              inner_join(word.count, wait.time2)))
corr.table(m124 %>% select(colnames(pgs.clean)[-1]),
           m124 %>% select(colnames(chr.wise), colnames(word.wise), 
                           word_count, avg_wait, -spid),
           method = "spearman") %>%
  mutate(FDR = p.adjust(pval, method = "fdr")) %>%
  filter(V1 %in% c(colnames(pgs.clean)), 
         V2 %in% c("word_count", "avg_wait",colnames(chr.wise), colnames(word.wise))) %>%
  mutate(cat1 = ifelse(grepl(paste(c("characters",  "reading_grade", "syllables"), 
                                   collapse = "|"), V2), 
                       "average", 
                       ifelse(V2 %in% c("word_count", "avg_wait"), 
                              "total", 
                              "ratio")),
         V2 = sub("lingmatch\\.", "", V2),
         V2 = sub("nrc\\.", "", V2),
         V2 = factor(V2, levels = unique(V2))) %>%
  ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(FDR < 0.05, "**", ifelse(pval<0.05, "*",""))))+
  geom_tile()+
  geom_text(size = 3, color = "white")+
  ggh4x::facet_grid2(rows = vars(cat1), 
                     scales = "free", space = "free") +
  scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1]) +
  labs(x = "", y = "language metrics from recorded COWAT",
       caption = paste0("n(samples): ", nrow(m124), "\n",
                        "**   FDR<0.05", "\n",
                        "*    pval<0.05")) +
  my.guides
ggsave(filename = paste0("figs/corr-pgs-COWAT-language-features.png"),
       width = 6, height = 8, units = "in", dpi = 320, bg = "white")
#####
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
##################### similarity between pairs of words ########################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# set up
rm(list=ls());gc();source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(text)
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
################################################################################
################################################################################
# read the demographics and remapping
demo <- read_csv("data/raw/ID_mapping_and_demographics.csv") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("NASD", "ASD")))
# read the COWAT
cowat.clean <- read_rds("data/derivatives/cowat-clean.rds")
# read the IQ data
iq.clean <- read_rds("data/derivatives/iq-clean.rds")
# education
edu <- read_csv("data/raw/educational_attainment_phenotypes.csv")
# pgs
pgs.clean <- read_rds("data/derivatives/pgs-clean.rds")
# 
m1.m2 <- pgs.clean
##########
# read embeddings and keep samples of interest
load("data/derivatives/word-embedding-from-text-package.rda")
emb.iq <- emb.text.m %>%
  filter(spid %in% iq.clean$spid)
emb.edu <- emb.text.m %>%
  filter(spid %in% edu$spid)
emb.pgs <- emb.text.m %>%
  filter(spid %in% pgs.clean$spid)
rm(emb.text.m);gc()
# choose which category of these 3 above, you'd like to keep
emb.text.m <- emb.pgs %>%
  distinct(prompt, spid, word, .keep_all = T)
################################################################################
# get the cosine similarity between each pair of words for each participant
registerDoMC(cores = 4)
pairs.sim <- foreach(i = 1:length(unique(emb.text.m$spid)), .combine = rbind) %dopar% {
  id <- unique(emb.text.m$spid)[i]
  # get the prompt that we have a transcription for
  prompts <- unique(emb.text.m$prompt)
  # loop over these prompts and extract similiraity between every possible pair
  p.sim.2 <- foreach(k = 1:length(prompts), .combine = rbind) %dopar% {
    # identify prompt words said
    w0 <- prompts[k]
    # extract participant's response
    tmp.emb <- emb.text.m %>% 
      filter(spid ==id, prompt ==w0)
    # make sure you have a response for that prompt
    if (nrow(tmp.emb) == 0) {
      return(NULL)
    }
    # build a dataframe that has every possible combination of pairs from words said
    sim.df <- data.frame(spid = id,
                         prompt = w0,
                         w_order = rep(c(1:nrow(tmp.emb)), each = nrow(tmp.emb)),
                         w1 = rep(tmp.emb$word, each = nrow(tmp.emb)),
                         w2 = rep(tmp.emb$word, nrow(tmp.emb)),
                         cos_similarity = NA)
    # loop over these pairs, and calculate the cosine similarity between their prompt embeddings
    for (j in 1:nrow(sim.df)) {
      # j=1
      w1 <- sim.df$w1[j]
      w2 <- sim.df$w2[j]
      # get embeddings of word 1
      e1 <- tmp.emb %>%
        distinct(word, .keep_all = T) %>%
        filter(word == w1) %>%
        select(-c(1:5))
      # get embeddings of word 2
      e2 <- tmp.emb %>%
        distinct(word, .keep_all = T) %>%
        filter(word == w2) %>%
        select(-c(1:5))
      # calculate the cosine similarity and save it
      sim.df$cos_similarity[j] <- text::textSimilarity(e1,e2, method = "cosine")
    }
    return(sim.df)
  }
  return(p.sim.2)
}
w2_order <- pairs.sim %>% 
  select(1,2,w2_order=w_order,w2=w1) %>% distinct()
pairs.sim <- pairs.sim %>%
  left_join(w2_order, relationship = "many-to-many")
# save the similarity between pairs
write_rds(pairs.sim, paste0("data/derivatives/pairs-sim-by-prompt-by-participant-pgs.rds"), compress = "gz")
pairs.sim <- read_rds(("data/derivatives/pairs-sim-by-prompt-by-participant-pgs.rds"))
################################################################################
# histogram of cosine similarity distribution for pairs
# get a df of pairs sim, just for consec pairs
cons.pairs <- pairs.sim %>%
  mutate(consec = ifelse(w2_order==w_order+1, T, F)) %>%
  filter(consec==T) %>%
  distinct() %>%
  filter(w1 != w2)
write_rds(cons.pairs, "data/derivatives/cons-pairs.rds", compress = "gz")
pairs.sim %>%
  mutate(consec = ifelse(w2_order==w_order+1, T, F)) %>%
  distinct() %>%
  filter(w_order != w2_order) %>%
  # filter(w1 != w2) %>%
  ggplot(aes(x=cos_similarity, fill = consec)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = redblack.col, name = "consecuetive pairs only") +
  labs(title = "faceted by task")
ggsave(bg = "white", filename = "figs/distribution-of-cos-similarity-of-all-pairs.png",
       width = 4, height = 4, units = "in", dpi = 320)
################################################################################
################################# correlations #################################
################################################################################
######
# correlation between consec. pairs similarity and PGS
######
# get a df of pairs sim, just for consec pairs
m125 <- left_join(cons.pairs, m1.m2) %>%
  drop_na() %>%
  left_join(demo[,2:5])
# write_rds(m125, "data/derivatives/lmer-inputs/cons-pairs.rds")
####
# lmer for demo
lm <- lmerTest::lmer(cos_similarity ~ age + sex + age:sex + diagnosis + (1|spid) + (1|prompt),
                     data = m125 %>%
                       select(cos_similarity, spid, prompt, age, sex, diagnosis))
# get summ
demo.lmer <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p) %>%
  mutate(var = "demo")
write_rds(demo.lmer, "data/derivatives/demo-lmer/pairs-sim.rds")
####
# lmer for PGS
library(lmerTest)
registerDoMC(cores = 6)
lm.results <- foreach(i=9:32, .combine = rbind) %dopar% {
  var <- colnames(m125)[i]
  # predict cosine similarity of the pair using the PGS. 
  # adding participant's ID as a random variable, and the prompt/mini-task
  lm <- lmerTest::lmer(cos_similarity ~ xx + age + sex + age:sex + diagnosis + (1|spid) + (1|prompt),
                       data = cbind(m125 %>% 
                                      select(cos_similarity, spid, prompt, age, sex, diagnosis),
                                    xx=m125[,i]) %>%
                         rename(xx=7))
  gc()
  # combine results in a df, and save
  df <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
    as.data.frame() %>%
    rownames_to_column("fixed") %>%
    filter(fixed != "(Intercept)") %>%
    rename(Estimate = `Est.`,
           confint_min = `2.5%`,
           confint_max = `97.5%`,
           pval = p) %>%
    mutate(var = var)
  write_rds(df, paste0("data/derivatives/pairs-lmer/", var, ".rds"))
  gc()
  return(df)
}
# combine the saved lmer results
lm.results <- foreach(i=9:ncol(m125), .combine = rbind) %dopar% {
  var <- colnames(m125)[i]
  if (file.exists(paste0("data/derivatives/pairs-lmer/", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/pairs-lmer/", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
lm.results <- lm.results %>% 
  mutate(FDR = p.adjust(pval, method = "fdr"))
# write_rds(lm.results, "data/derivatives/pairs-lmer/all-lmer-results.rds", compress = "gz")
# make plot for results
p1 <- lm.results %>%
  filter(fixed %in% c("xx")) %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05"),
         fixed = ifelse(fixed == "xx", "PGS", fixed)) %>%
  ggplot(aes(x=Estimate, y=var)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey"),
        strip.text.y.right = element_text(angle = 0)) +
  labs(x = "Estimate for predicting cosine similarity of consec. pairs", y="",
       caption = paste0("n(samples): ", length(unique(m125$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    lmer(cos_similarity ~ X + age + sex + age:sex + diagnosisASD + (1|spid) + (1|prompt))", "\n",
                        "    where X is a selected variable from the PGS list"))
# demographics plot
p2 <- demo.lmer %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=fixed)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate for predicting cosine similarity of consec. pairs", y="",
       caption = paste0("n(samples): ", length(unique(m125$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    lmer(cos_similarity ~ age + sex + age:sex + diagnosisASD + (1|spid) + (1|prompt))"))
patchwork::wrap_plots(p2,p1,ncol = 1,heights = c(1,5))
ggsave(filename = "figs/lmer-cos-sim-consec-pairs-by-pgs-random-id-and-prompt.png",
       width = 10, height = 14, units = "in", bg = "white", dpi = 360)
#############
# get correlation between consec. pairs similarity and the waiting time between these pairs
#############
c.clean <- cowat.clean %>%
  filter(prompt != word) %>% #drop the words that are exactly the same as the prompt prompt
  distinct(spid, prompt, word, .keep_all = T) %>% # only keep unique words and drop repeated by the same participant in the same prompt
  drop_na(start_time, end_time) # drop the words with no timestamps
wait.time <- cbind(c.clean[-nrow(c.clean),] %>% select(-word),
                   w1 = c.clean$word[-nrow(c.clean)],
                   w2 = c.clean$word[-1],
                   w1_index = rownames(c.clean)[-nrow(c.clean)],
                   w2_index = rownames(c.clean)[-1],
                   wait=(c.clean$start_time[-1] - c.clean$end_time[-nrow(c.clean)])) %>%
  filter(wait>0)
tmp <- inner_join(cons.pairs, wait.time)
tmp %>% 
  ggplot(aes(x=cos_similarity, y=log2(wait))) +
  geom_point(size=1)+
  geom_smooth(method = "loess", color = boxplot.colors[2]) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(color = "red") +
  labs(x="cosine similarity of consec. pair",
       y="log2(silence time between consec. pairs)",
       caption = paste0("n(samples): ", length(unique(tmp$spid))))
ggsave(filename = "figs/corr-between-consec-pairs-cos-similarity-and-thinking-time.png", bg = "white",
       width = 4, height = 4, units = "in", dpi = 360)
################################################################################
################################################################################
################################################################################
################### stats table for participants' performance ##################
################################################################################
################################################################################
################################################################################
################################################################################
# set up and clean
rm(list=ls());gc();source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
################################################################################
################################################################################
# load files before
# get the clean transcription for data
all <- read_rds("data/derivatives/cowat-analyzed.rds") %>%
  distinct(spid, prompt, word, .keep_all = T)
################################################################################
#####
# table of how many words said by participant in each task
#####
tmp <- all[,2:7] %>%
  group_by(spid, prompt) %>%
  dplyr::summarise(count = n()) %>%
  pivot_wider(names_from = prompt, values_from = count, id_cols = spid)
table <- kableExtra::kable(tmp, format="html") %>%
  kableExtra::kable_styling(full_width = T, protect_latex = T)
table
#####
# plot these stats
#####
tmp %>% 
  pivot_longer(cols = colnames(tmp)[-1], names_to = "word") %>%
  mutate(task = ifelse(nchar(word)==1, 3,1)) %>%
  group_by(word, value,task) %>%
  dplyr::summarise(count = n()) %>%
  ggplot(aes(x=value, y=count))+
  geom_histogram(stat = "identity")+
  facet_wrap(~reorder(word, desc(task))) +
  geom_vline(xintercept = 5, color = "red", linetype=2) +
  labs(x="count of words said by participant",
       y="count of participants")
ggsave(filename = "figs/word-count-stats-per-task.png", bg = "white",
       width = 6, height = 5, units = "in", dpi = 320)
#####
################################################################################
################################################################################
################################################################################
#################### calc Euclidean distance and correlate it ##################
################################################################################
################################################################################
################################################################################
# set up and clean
rm(list=ls());gc();source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
################################################################################
################################################################################
# load files before
# get the clean transcription for data
all <- read_rds("data/derivatives/cowat-analyzed.rds") %>%
  distinct(spid, word, prompt, .keep_all = T)
# read the demographics and remapping
demo <- read_csv("data/raw/ID_mapping_and_demographics.csv") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("NASD", "ASD")))
# get word embeddings from text package
load("data/derivatives/word-embedding-from-text-package.rda")
emb.text.m <- emb.text.m %>% distinct(spid,prompt,word, .keep_all = T) # keep unique words per participant for each mini-task
# pgs
pgs.clean <- read_rds("data/derivatives/pgs-clean.rds")
# 
m1.m2 <- pgs.clean
################################################################################
#####
# calculate the full euclidean distance traveled per task, and normalize it by number of words
#####
# this also considers at least 2 words said per task to get the euc distance

# get summary per participant per task
tmp <- all[,2:7] %>%
  group_by(spid, prompt) %>%
  dplyr::summarise(count = n()) %>%
  pivot_wider(names_from = prompt, values_from = count, id_cols = spid)
# loop over participants and select prompts that have enough data for you
# then, get the embeddings for these words, and calculate the distance
# after getting the full euclidean distance per that task, divide it over number of words
euc <- foreach(i=1:length(unique(tmp$spid)), .combine = rbind) %dopar% {
  id <- unique(tmp$spid)[i]
  # get counts of response per task for the selected participant
  df1 <- tmp %>% 
    filter(spid == id) %>%
    pivot_longer(cols = colnames(tmp)[-1]) %>%
    filter(value >=2)
  # identify unique tasks
  words.to.keep <- unique(df1$name)
  df2 <- emb.text.m %>%
    filter(spid==id,
           prompt %in% words.to.keep)
  # loop over words here 
  word.opt <- foreach(j = 1:length(words.to.keep), .combine = rbind) %dopar% {
    w0 <- words.to.keep[j]
    df3 <- df2 %>%
      filter(prompt==w0)
    ######
    # get the overall euclidean distance and divide it by number of words
    text.to.look <- unique(df3$word)
    full.dist <- 0
    for (g in 1:(length(text.to.look)-1)) {
      t0 <- text.to.look[g]
      t1 <- text.to.look[g+1]
      euc.dist <- dist(rbind(as.numeric(emb.text.m%>%
                                          filter(spid==id,prompt==w0,word==t0)%>%
                                          select(starts_with("Dim"))),
                             as.numeric(emb.text.m%>%
                                          filter(spid==id,prompt==w0,word==t1)%>%
                                          select(starts_with("Dim")))))
      full.dist <- full.dist + euc.dist
    }
    full.dist.n <- as.numeric(full.dist)/nrow(df3)
    ######
    df6 <- data.frame(spid = id,
                      prompt = w0,
                      full_euc_dist = as.numeric(full.dist),
                      full_euc_dist_normalized = as.numeric(full.dist.n))
    return(df6)
  }
  return(word.opt)
}
# save euc distances
write_rds(euc, "data/derivatives/euc-distance-w-minimum-of-2-words-per-task.rds")
# euc <- read_rds("data/derivatives/euc-distance-w-minimum-of-2-words-per-task.rds")
#####
# predict normalized Euclidean distance traveled 
#####
# use the PGS as a major predictor
# add the spid and the prompt as random variables
m123 <- inner_join(euc, m1.m2) %>% 
  left_join(demo[,2:5])
# write_rds(m123, "data/derivatives/lmer-inputs/euc.rds")
library(lmerTest)
####
# lmer for demo
lm <- lmerTest::lmer(full_euc_dist_normalized ~ age + sex + age:sex + diagnosis + (1|spid) + (1|prompt),
                     data = m123 %>%
                       select(full_euc_dist_normalized, spid, prompt, age, sex, diagnosis))
# get summ
demo.lmer <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p) %>%
  mutate(var = "demo")
write_rds(demo.lmer, "data/derivatives/demo-lmer/euc.rds")
####
# lmer for PGS
registerDoMC(cores = 6)
lm.results <- foreach(i=5:28, .combine = rbind) %dopar% {
  var <- colnames(m123)[i]
  # predict euclidean distance using the PGS. 
  # adding participant's ID as a random variable, and the prompt
  lm <- lmerTest::lmer(full_euc_dist_normalized ~ xx + age + sex + age:sex + diagnosis + (1|spid) + (1|prompt),
                       data = cbind(m123 %>% 
                                      select(full_euc_dist_normalized, spid, prompt, age, sex, diagnosis),
                                    xx=m123[,i]) %>%
                         rename(xx=7))
  gc()
  # combine results in a df, and save
  df <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
    as.data.frame() %>%
    rownames_to_column("fixed") %>%
    filter(fixed != "(Intercept)") %>%
    rename(Estimate = `Est.`,
           confint_min = `2.5%`,
           confint_max = `97.5%`,
           pval = p) %>%
    mutate(var = var)
  write_rds(df, paste0("data/derivatives/euc-lmer/", var, ".rds"))
  gc()
  return(df)
}
# combine the saved lmer results
lm.results <- foreach(i=5:ncol(m123), .combine = rbind) %dopar% {
  var <- colnames(m123)[i]
  if (file.exists(paste0("data/derivatives/euc-lmer/", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/euc-lmer/", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
lm.results <- lm.results %>% mutate(FDR = p.adjust(pval, method = "fdr"))
# write_rds(lm.results, "data/derivatives/euc-lmer/all-lmer-results.rds", compress = "gz")
# make plot for results
p1 <- lm.results %>%
  filter(fixed == "xx") %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=var,)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate for predicting normalized full Euclidean distance", y="",
       caption = paste0("n(samples): ", length(unique(m123$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    lmer(normalized_euc_distance ~ X + age + sex + age:sex + diagnosisASD + (1|spid) + (1|prompt))", "\n",
                        "    where X is a selected variable from the PGS list", "\n",
                        "Derivation of Euclidean distance was as follows:", "\n",
                        "    distance = traveled distance by participant in each prompt","\n",
                        "    (i.e., sum of distance between consec. pairs)","\n",
                        "    Normalized by word count said by a participant in this prompt","\n",
                        "    Only kept participants with at least 2 words in response to prompt"))
# demographics plot
p2 <- demo.lmer %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=fixed)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate for predicting normalized full Euclidean distance", y="",
       caption = paste0("n(samples): ", length(unique(m123$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    lmer(normalized_euc_distance ~ age + sex + age:sex + diagnosisASD + (1|spid) + (1|prompt))"))
patchwork::wrap_plots(p2,p1,ncol = 1,heights = c(1,5))
ggsave(filename = "figs/lmer-full-euc-distance-by-pgs-random-id-and-word.png",
       width = 10, height = 14, units = "in", bg = "white", dpi = 360)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
####################### calc divergence and correlate it #######################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# set up and clean
rm(list=ls());gc();source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
################################################################################
################################################################################
# load files before
# read the demographics and remapping
demo <- read_csv("data/raw/ID_mapping_and_demographics.csv") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("NASD", "ASD")))
# get the clean transcription for data
all <- read_rds("data/derivatives/cowat-analyzed.rds") %>%
  filter(word != prompt) %>%
  distinct(spid, word, prompt, .keep_all = T)
# pgs
pgs.clean <- read_rds("data/derivatives/pgs-clean.rds")
# 
m1.m2 <- pgs.clean
all <- all %>% filter(spid %in% m1.m2$spid) %>% distinct(spid, word, prompt, .keep_all = T)
# get word embeddings from text package
load("data/derivatives/word-embedding-from-text-package.rda")
emb.text.m <- emb.text.m %>% distinct(spid,word,prompt, .keep_all = T) # keep unique words per participant for each mini-task
# load the pairs similarity values
pairs.sim <- read_rds(paste0("data/derivatives/pairs-sim-by-prompt-by-participant-pgs.rds")) %>%
  filter(w1 != w2)
# get a df of pairs sim, just for consec pairs
cons.pairs <- read_rds("data/derivatives/cons-pairs.rds")
################################################################################
#####
# calculate the optimal and actual trajectories per participant for first 10 words per task/word
#####
# get summary per participant per task
tmp <- all[,2:7] %>%
  distinct(spid, prompt, word, .keep_all = T) %>%
  group_by(spid, prompt) %>%
  dplyr::summarise(count = n()) %>%
  pivot_wider(names_from = prompt, values_from = count, id_cols = spid)
# loop over participants and select prompts that have enough data for you
# then, get the embeddings for these words, and calculate the hamiltonian distance, and actual distance
registerDoMC(cores = 6)
divergence <- foreach(i=1:length(unique(tmp$spid)), .combine = rbind) %dopar% {
  id <- unique(tmp$spid)[i]
  # get data for selected id
  df1 <- tmp %>% 
    filter(spid == id) %>%
    pivot_longer(cols = colnames(tmp)[-1]) %>%
    filter(value >=3) # make sure the participant has said at least 3 words
  # identify words to keep
  words.to.keep <- unique(df1$name)
  df2 <- emb.text.m %>%
    filter(spid==id,
           prompt %in% words.to.keep)
  # make sure the participant has an actual response for each prompt
  if (length(words.to.keep)==0) {
    return(NULL)
  }
  # loop over words here to get traveled, and optimal distances
  word.opt <- foreach(j = 1:length(words.to.keep), .combine = rbind) %dopar% {
    w0 <- words.to.keep[j]
    
    df3 <- pairs.sim %>%
      filter(spid ==id, prompt==w0) %>%
      pivot_wider(names_from = "w2", values_from = "cos_similarity", id_cols = "w1") %>%
      column_to_rownames("w1")
    df3 <- df3[,rownames(df3)]
    ######
    # get the optimal trajectory
    library(TSP)
    distance.matrix <- as.dist(1 - df3) # convert cosine similarity matrix to a distance matrix
    tsp.dist <- TSP(distance.matrix) # travelling salesman problem
    tsp.sol <- as.integer(solve_TSP(tsp.dist, method = "repetitive_nn"))
    df4 <- data.frame(w_order = tsp.sol, text = rownames(df3)[tsp.sol])
    df5 <- data.frame(w1 = df4$text[-(nrow(df4))], # make a df with optimal path order
                      w1_order = df4$w_order[-(nrow(df4))],
                      w2 = df4$text[-1],
                      w2_order = df4$w_order[-1],
                      distance = NA)
    for (k in 1:nrow(df5)) { # get optimal path distances
      df5$distance[k] <- (1-df3[df5$w1_order[k],df5$w2_order[k]])
    }
    opt.dist <- sum(df5$distance) # this is the optimal distance for this word for this participant
    act.order <- data.frame(w1 = rownames(df3)[-(nrow(df3))], # make a df with actual path order
                            w1_order = c(1:(nrow(df3)-1)),
                            w2 = colnames(df3)[-1],
                            w2_order = c(2:nrow(df3)),
                            distance = NA)
    for (m in 1:nrow(act.order)) { # get actual path distances
      act.order$distance[m] <- (1-df3[act.order$w1_order[m],act.order$w2_order[m]])
    }
    act.dist <- sum(act.order$distance) # this is the actual distance for this word for this participant
    #######
    df6 <- data.frame(spid = id,
                      prompt = w0,
                      opt_dist = opt.dist,
                      act_dist = act.dist,
                      global_divergence = (act.dist - opt.dist),
                      global_divergence_normalized = (act.dist - opt.dist)/nrow(df3),
                      word_count = nrow(df3))
    return(df6)
  }
  return(word.opt)
}
# save the divergence data
write_rds(divergence, "data/derivatives/divergence-w-minimum-of-3-words-per-task.rds")
# divergence <- read_rds("data/derivatives/divergence-w-minimum-of-3-words-per-task.rds")
#####
# predict divergence
#####
# use the PGS as a major predictor
# add the spid and the prompt as random variables
m122 <- inner_join(divergence, m1.m2) 
# write_rds(m122, "data/derivatives/lmer-inputs/divergence.rds")
# add the average cosine similarity per prompt for each participant
avg.sim <- cons.pairs %>%
  group_by(spid, prompt) %>%
  dplyr::summarise(avg_cos_sim = mean(cos_similarity))
m123 <- left_join(m122, avg.sim) %>%
  left_join(demo[,2:5])
library(lmerTest)
####
# lmer for demo
lm <- lmerTest::lmer(global_divergence ~ age + sex + age:sex + diagnosis + (1|spid) + (1|prompt),
                     data = m123 %>%
                       select(global_divergence, spid, prompt, age, sex, diagnosis))
# get summ
demo.lmer <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p) %>%
  mutate(var = "demo")
write_rds(demo.lmer, "data/derivatives/demo-lmer/divergence.rds")
####
# lmer for PGS
registerDoMC(cores = 6)
lm.results <- foreach(i=8:31, .combine = rbind) %dopar% {
  var <- colnames(m123)[i]
  # predict divergence using the PGS. 
  # adding participant's ID as a random variable, and the prompt
  lm <- lmerTest::lmer(global_divergence ~ xx + word_count + age + sex + age:sex + diagnosis + (1|spid) + (1|prompt),
                       data = cbind(m123 %>% 
                                      select(global_divergence, spid, prompt, word_count, age, sex, diagnosis) %>%
                                      mutate(global_divergence = scale(global_divergence, scale = T, center = T)[,1]),
                                    xx=m123[,i]))
  gc()
  # combine results in a df, and save
  df <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
    as.data.frame() %>%
    rownames_to_column("fixed") %>%
    filter(fixed != "(Intercept)") %>%
    rename(Estimate = `Est.`,
           confint_min = `2.5%`,
           confint_max = `97.5%`,
           pval = p) %>%
    mutate(var = var)
  write_rds(df, paste0("data/derivatives/divergence-lmer/", var, ".rds"))
  gc()
  return(df)
}
# combine the saved lmer results
lm.results <- foreach(i=9:ncol(m123), .combine = rbind) %dopar% {
  var <- colnames(m123)[i]
  if (file.exists(paste0("data/derivatives/divergence-lmer/", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/divergence-lmer/", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
lm.results <- lm.results %>% mutate(FDR = p.adjust(pval, method = "fdr"))
# write_rds(lm.results, "data/derivatives/divergence-lmer/all-lmer-results.rds", compress = "gz")
# make plot for results
p1 <- lm.results %>%
  filter(fixed=="xx") %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=var,)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate for predicting z-standardized divergence", y="",
       caption = paste0("n(samples): ", length(unique(m123$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    lmer(z-standardized_divergence ~ X + word_count + age + sex + age:sex + diagnosisASD + (1|spid) + (1|prompt))", "\n",
                        "    where X is a selected variable from the PGS list", "\n",
                        "        and word_count is how many points/words were in the path", "\n",
                        "Derivation of divergence was as follows:", "\n",
                        "    divergence = (actual_path - optimal_path) in each prompt","\n",
                        "    optimal path is the sequence which visits each item exactly once in an order ", "\n",
                        "        that minimizes the total semantic distance traveled","\n",
                        "    actual path is by following the order the participant said","\n",
                        "    (i.e., negative divergence reflects increasignly optimal word selection)", "\n", 
                        "    Only kept participants with at least 3 words in response to prompt"))
# demographics plot
p2 <- demo.lmer %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=fixed)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate for predicting z-standardized divergence", y="",
       caption = paste0("n(samples): ", length(unique(m123$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    lmer(z-standardized_divergence ~ age + sex + age:sex + diagnosisASD + (1|spid) + (1|prompt))"))
patchwork::wrap_plots(p2,p1,ncol = 1,heights = c(1,5))
ggsave(filename = "figs/lmer-divergence-by-pgs-and-wc-random-id-and-word.png",
       width = 10, height = 14, units = "in", bg = "white", dpi = 360)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#################### calc vocabulary depth and correlate it ####################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# set up and clean
rm(list=ls());gc();source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
library(text)
library(umap)
################################################################################
################################################################################
# load files before
# read the demographics and remapping
demo <- read_csv("data/raw/ID_mapping_and_demographics.csv") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("NASD", "ASD")))
# get the clean transcription for data
all <- read_rds("data/derivatives/cowat-analyzed.rds") %>%
  distinct(spid, word, prompt, .keep_all = T)
# get word embeddings from text package
load("data/derivatives/word-embedding-from-text-package.rda")
emb.text.m <- emb.text.m %>% distinct(spid,prompt,word, .keep_all = T) # keep unique words per participant for each mini-task
# pgs
pgs.clean <- read_rds("data/derivatives/pgs-clean.rds")
# 
m1.m2 <- pgs.clean
################################################################################
#####
# calculate the vocabulary depth as the total volume between embeddings/UMAP points in space
#####
#####
# get UMAP for embeddings
tmp.emb <- emb.text.m %>% distinct(word, .keep_all = T)
umap.o <- umap(tmp.emb %>% 
                 select(starts_with("Dim")) %>% 
                 as.matrix(), 
               n_components = 3, preserve.seed = T)
umap.dim <- cbind(tmp.emb[,1:6],
                  umap.o$layout %>%
                    as.data.frame() %>%
                    rename_all(.funs = function(x) paste0("Dim", sub("V", "", x))))
write_rds(umap.dim, "data/derivatives/umap-3d-from-cowat-embeddings.rds", compress = "gz")
# umap.dim <- read_rds("data/derivatives/umap-3d-from-cowat-embeddings.rds")
# 
umap.dim <- full_join(emb.text.m[,1:6], umap.dim[,c(5,7:9)])
registerDoMC(cores = 6)
chulls <- foreach(i=1:nrow(demo), .combine = rbind) %dopar% {
  id <- demo$spid[i]
  # get data for selected id
  df <- umap.dim %>%
    filter(spid==id) %>%
    distinct(word, .keep_all=T)
  # filter(!grepl(paste("two", "four", "fourth", "five", "fives", "six", "seven", "seventeen", collapse = "|"), text))
  if (nrow(df)>=4) {
    vol.all <- cxhull::cxhull(df[,7:9]%>%as.matrix())$volume
    c.all <- nrow(df)
  } else {
    vol.all = 0;c.all = 0
  }
  # what if we did the chull for each letter independently? possibly helping identify outliers source
  df.L <- df %>% filter(prompt == "L") %>% select(starts_with("Dim"))
  df.F <- df %>% filter(prompt == "F") %>% select(starts_with("Dim"))
  df.C <- df %>% filter(prompt == "C") %>% select(starts_with("Dim"))
  df.A <- df %>% filter(prompt == "A") %>% select(starts_with("Dim"))
  df.S <- df %>% filter(prompt == "S") %>% select(starts_with("Dim"))
  if (nrow(df.L)>=4) {
    vol.L <- cxhull::cxhull(df.L%>%as.matrix())$volume
    c.L = nrow(df.L)
  } else {
    vol.L = 0;c.L = 0
  }
  if (nrow(df.F)>=4) {
    vol.F <- cxhull::cxhull(df.F%>%as.matrix())$volume
    c.F = nrow(df.F)
  } else {
    vol.F = 0;c.F = 0
  }
  if (nrow(df.C)>=4) {
    vol.C <- cxhull::cxhull(df.C%>%as.matrix())$volume
    c.C = nrow(df.C)
  } else {
    vol.C = 0;c.C = 0
  }
  if (nrow(df.A)>=4) {
    vol.A <- cxhull::cxhull(df.A%>%as.matrix())$volume
    c.A = nrow(df.A)
  } else {
    vol.A = 0;c.A = 0
  }
  if (nrow(df.S)>=4) {
    vol.S <- cxhull::cxhull(df.S%>%as.matrix())$volume
    c.S = nrow(df.S)
  } else {
    vol.S = 0;c.S = 0
  }
  return(data.frame(spid = id, vol_all = vol.all, 
                    count_all = c.all,
                    vol_L = vol.L, count_L = c.L,
                    vol_F = vol.F,count_F = c.F,
                    vol_C = vol.C,count_C = c.C,
                    vol_A = vol.A,count_A = c.A,
                    vol_S = vol.S,count_S = c.S))
}
write_rds(chulls, "data/derivatives/chulls.rds")
####
# check distribution
####
# outliers are there for the total volume from all tasks combined
chulls %>%
  pivot_longer(cols = starts_with("vol"), names_to = "vol_source", values_to = "vol_value") %>%
  pivot_longer(cols = starts_with("count"), names_to = "count_source", values_to = "count_value") %>%
  mutate(vol_source = sub("vol_", "", vol_source),
         count_source = sub("count_", "", count_source)) %>%
  filter(vol_source == count_source) %>%
  ggplot(aes(x=vol_value)) +
  geom_histogram() + facet_wrap("vol_source", scales = "free") +
  labs(title = "distribution convex hull volume of word embeddings")
ggsave(filename = "figs/distribution-of-vocab-depth.png",
       width = 7, height = 6, units = "in", bg = "white", dpi = 360)
############
# correlate with PGS
############
inner_join(m1.m2, chulls) %>%
  pivot_longer(cols = c(colnames(m1.m2), -spid), names_to = "measure") %>%
  # filter(vol_all<50) %>%
  ggplot(aes(x=vol_all, y=value)) +
  geom_point() +geom_smooth(method = "lm") + ggpubr::stat_cor(color = "red") +
  facet_wrap(~measure, scales = "free") +
  labs(caption = paste0("no correction done for anything here", "\n",
                        "only dropped outliers because of their calculated vocabulary depth"))
ggsave(filename = "figs/corr-pgs-vocab-depth-no-correction.png",
       width = 12, height = 10, units = "in", bg = "white", dpi = 360)
######################################
# probably need to correct for randomeness coming from task
######################################
m123 <- inner_join(m1.m2, 
                   chulls %>%
                     pivot_longer(cols = starts_with("vol"), names_to = "vol_source", values_to = "vol_value") %>%
                     pivot_longer(cols = starts_with("count"), names_to = "count_source", values_to = "count_value") %>%
                     mutate(vol_source = sub("vol_", "", vol_source),
                            count_source = sub("count_", "", count_source)) %>%
                     filter(vol_source == count_source)) %>%
  # filter(!spid %in% chulls$spid[which(chulls$vol_all>50)]) %>%
  # filter(vol_value<40) %>%
  filter(vol_source != "all") %>%
  left_join(demo[,2:5])
# write_rds(m123, "data/derivatives/lmer-inputs/chulls.rds")
####
# lmer for demo
lm <- glm(vol_value ~ count_value + vol_source + age + sex + age:sex + diagnosis,
          data = m123 %>%
            select(vol_value, count_value, vol_source, spid, age, sex, diagnosis))
# get summ
demo.lmer <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p) %>%
  mutate(var = "demo")
write_rds(demo.lmer, "data/derivatives/demo-lmer/chulls.rds")
####
# lmer for PGS
registerDoMC(cores = 6)
lm.results <- foreach(i=2:25, .combine = rbind) %dopar% {
  var <- colnames(m123)[i]
  lm <- glm(vol_value ~ xx + count_value + vol_source + age + sex + age:sex + diagnosis,
            data = cbind(m123 %>% 
                           select(vol_value, vol_source, spid, count_value, age, sex, diagnosis) %>%
                           mutate(vol_value = scale(vol_value, scale = T, center = T)[,1]),
                         xx=m123[,i]) %>%
              rename(xx = 8))
  gc()
  # combine results in a df, and save
  df <- jtools::summ(lm, confin = T, pval = T)$coeftable %>%
    as.data.frame() %>%
    rownames_to_column("fixed") %>%
    filter(fixed != "(Intercept)") %>%
    rename(Estimate = `Est.`,
           confint_min = `2.5%`,
           confint_max = `97.5%`,
           pval = p) %>%
    mutate(var = var)
  write_rds(df, paste0("data/derivatives/chulls-lmer/glm-all-union-", var, ".rds"))
  gc()
  return(df)
}
# combine the saved lmer results
lm.results <- foreach(i=3:29, .combine = rbind) %dopar% {
  var <- colnames(m123)[i]
  if (file.exists(paste0("data/derivatives/chulls-lmer/glm-all-union-", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/chulls-lmer/glm-all-union-", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
lm.results <- lm.results %>% mutate(FDR = p.adjust(pval, method = "fdr"))
# write_rds(lm.results, "data/derivatives/chulls-lmer/all-lmer-results.rds", compress = "gz")
# make plot for results
p1 <- lm.results %>%
  filter(fixed=="xx") %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=var,)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate for predicting z-standardized vocabulary depth", y="",
       caption = paste0("n(samples): ", length(unique(m123$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    glm(z-standardized_vocab_depth ~ X + word_count + prompt + age + sex + age:sex + diagnosisASD)", "\n",
                        "    where X is a selected variable from the PGS list", "\n",
                        "        and word_count is how many points/words in participant's response in this task", "\n",
                        "Derivation of vocabulary depth was as follows:", "\n",
                        "    vocabulary_depth = volume enclosed between 3d semantic space of words said","\n",
                        "    765 word embeddings were derived from 'text' package and then", "\n",
                        "        UMAP was utilized to reduce dimensions for 3D","\n",
                        "    Only kept participants with at least 4 words in response to task/word"))
# demographics plot
p2 <- demo.lmer %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=fixed)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate for predicting z-standardized vocabulary depth", y="",
       caption = paste0("n(samples): ", length(unique(m123$spid)), "\n",
                        "the estimates are derived from the model below:", "\n",
                        "    glm(z-standardized_vocab_depth ~ word_count + prompt + age + sex + age:sex + diagnosisASD)"))
patchwork::wrap_plots(p2,p1,ncol = 1,heights = c(1,5))
ggsave(filename = "figs/glm-vocab-depth-by-pgs-and-wc-and-word.png",
       width = 10, height = 14, units = "in", bg = "white", dpi = 360)
################################################################################
################################################################################
################################################################################
###################### combine all lmer results in 1 plot ######################
################################################################################
################################################################################
################################################################################
# set up and clean
rm(list=ls());gc();source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
################################################################################
################################################################################
# load files before
# pgs
pgs.clean <- read_rds("data/derivatives/pgs-clean.rds")
# 
m1.m2 <- pgs.clean
#####
# pairs lmer
#####
pairs <- foreach(i=2:25, .combine = rbind) %dopar% {
  var <- colnames(m1.m2)[i]
  if (file.exists(paste0("data/derivatives/pairs-lmer/", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/pairs-lmer/", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
pairs <- pairs %>% 
  mutate(FDR = p.adjust(pval, method = "fdr"), 
         source = "cos_similarity for consecutive pairs")
pairs.demo <- read_rds("data/derivatives/demo-lmer/pairs-sim.rds") %>%
  mutate(source = "cos_similarity for consecutive pairs")
#####
# euc lmer
#####
euc <- foreach(i=2:25, .combine = rbind) %dopar% {
  var <- colnames(m1.m2)[i]
  if (file.exists(paste0("data/derivatives/euc-lmer/", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/euc-lmer/", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
euc <- euc %>% 
  mutate(FDR = p.adjust(pval, method = "fdr"), 
         source = "normalized full Euclidean distance")
euc.demo <- read_rds("data/derivatives/demo-lmer/euc.rds") %>%
  mutate(source = "normalized full Euclidean distance")
#####
# vocab depth lmer
#####
v.depth <- foreach(i=2:25, .combine = rbind) %dopar% {
  var <- colnames(m1.m2)[i]
  if (file.exists(paste0("data/derivatives/chulls-lmer/glm-all-union-", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/chulls-lmer/glm-all-union-", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
v.depth <- v.depth %>% 
  mutate(FDR = p.adjust(pval, method = "fdr"), 
         source = "vocabulary depth")
v.depth.demo <- read_rds("data/derivatives/demo-lmer/chulls.rds") %>%
  mutate(source = "vocabulary depth")
#####
# divergence lmer
#####
divergence <- foreach(i=2:25, .combine = rbind) %dopar% {
  var <- colnames(m1.m2)[i]
  if (file.exists(paste0("data/derivatives/divergence-lmer/", var, ".rds"))) {
    df <- read_rds(paste0("data/derivatives/divergence-lmer/", var, ".rds"))
    return(df)
  } else {
    return(NULL)
  }
}
divergence <- divergence %>% 
  mutate(FDR = p.adjust(pval, method = "fdr"), 
         source = "divergence")
divergence.demo <- read_rds("data/derivatives/demo-lmer/divergence.rds") %>%
  mutate(source = "divergence")
################################################################################
################################################################################
# combine all in one major df
all <- rbind(euc,
             divergence,
             pairs,
             v.depth %>% mutate(`d.f.` = NA))
all.demo <- rbind(euc.demo,
                  divergence.demo,
                  pairs.demo,
                  v.depth.demo %>% mutate(`d.f.` = NA))
# plot
p1 <- all %>%
  filter(fixed=="xx") %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=var,)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  ggh4x::facet_grid2(cols = vars(source), scales = "free") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey"),
        strip.text.y.right = element_text(angle = 0)) +
  labs(x = "Estimate", y="")
p1
p2 <- all.demo %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=fixed)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  ggh4x::facet_grid2(cols = vars(source), scales = "free") +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey")) +
  labs(x = "Estimate", y="")
p2
patchwork::wrap_plots(p2,p1, ncol = 1, heights = c(1,5))
ggsave(filename = "figs/lmer-all-metrics-by-pgs.png",
       width = 16, height = 14, units = "in", bg = "white", dpi = 360)
#
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#################### use years of education instead of pgs #####################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# set up and clean
rm(list=ls());gc();source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SPARK/language"
setwd(project.dir)
################################################################################
################################################################################
# load files before
# read the demographics and remapping
demo <- read_csv("data/raw/ID_mapping_and_demographics.csv") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("NASD", "ASD")))
# get the clean transcription for data
all <- read_rds("data/derivatives/cowat-analyzed.rds") %>%
  filter(word != prompt) %>%
  distinct(spid, word, prompt, .keep_all = T)
# edu
edu <- read_csv("data/raw/educational_attainment_phenotypes.csv")
# 
m1.m2 <- edu[,c(2,4)]
all <- all %>% filter(spid %in% m1.m2$spid) %>% distinct(spid, word, prompt, .keep_all = T)
# get word embeddings from text package
load("data/derivatives/word-embedding-from-text-package.rda")
emb.text.m <- emb.text.m %>% distinct(spid,word,prompt, .keep_all = T) # keep unique words per participant for each mini-task
# load the pairs similarity values
pairs.sim <- read_rds(paste0("data/derivatives/pairs-sim-by-prompt-by-participant-pgs.rds")) %>%
  filter(w1 != w2)
# get a df of pairs sim, just for consec pairs
cons.pairs <- read_rds("data/derivatives/cons-pairs.rds")
################################################################################
################################################################################
################################################################################
# get pairs sim
p.sim <- read_rds("data/derivatives/pairs-sim-by-prompt-by-participant-pgs.rds") %>%
  left_join(edu[,c(2,4)]) %>%
  left_join(demo[,2:6]) %>%
  drop_na()
p.sim.ea <- jtools::summ(lmerTest::lmer(cos_similarity ~ ea_num + age + sex + age:sex + diagnosis + (1|spid) + (1|prompt), data = p.sim), confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p)
###
# get euc
euc <- read_rds("data/derivatives/euc-distance-w-minimum-of-2-words-per-task.rds") %>%
  left_join(edu[,c(2,4)]) %>%
  left_join(demo[,2:6]) %>%
  drop_na()
euc.ea <- jtools::summ(lmerTest::lmer(full_euc_dist_normalized ~ ea_num + age + sex + age:sex + diagnosis + (1|spid) + (1|prompt), data = euc), confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p)
###
# get divergence
div <- read_rds("data/derivatives/divergence-w-minimum-of-3-words-per-task.rds") %>%
  left_join(edu[,c(2,4)]) %>%
  left_join(demo[,2:6]) %>%
  drop_na()
div.ea <- jtools::summ(lmerTest::lmer(global_divergence ~ ea_num + word_count + age + sex + age:sex + diagnosis + (1|spid) + (1|prompt), data = div), confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p)
###
# get chulls
chulls <- read_rds("data/derivatives/chulls.rds") %>%
  left_join(edu[,c(2,4)]) %>%
  left_join(demo[,2:6]) %>%
  drop_na() %>%
  pivot_longer(cols = starts_with("vol"), names_to = "vol_source", values_to = "vol_value") %>%
  pivot_longer(cols = starts_with("count"), names_to = "count_source", values_to = "count_value") %>%
  mutate(vol_source = sub("vol_", "", vol_source),
         count_source = sub("count_", "", count_source)) %>%
  filter(vol_source == count_source) %>%
  filter(vol_source != "all")
chulls.ea <- jtools::summ(glm(vol_value ~ ea_num + count_value + age + sex + age:sex + diagnosis + count_source, data = chulls), confin = T, pval = T)$coeftable %>%
  as.data.frame() %>%
  rownames_to_column("fixed") %>%
  filter(fixed != "(Intercept)") %>%
  rename(Estimate = `Est.`,
         confint_min = `2.5%`,
         confint_max = `97.5%`,
         pval = p)
###
# combine all and plot
all <- rbind(p.sim.ea %>% mutate(source = "consec sim", n = length(unique(p.sim$spid))),
             euc.ea %>% mutate(source = "full euc", n = length(unique(euc$spid))),
             div.ea %>% mutate(source = "divergence", n = length(unique(div$spid))),
             chulls.ea %>% mutate(`d.f.`=NA, source = "vocab depth", n = length(unique(chulls$spid))))

all %>%
  mutate(sig = ifelse(pval<0.05, "pval < 0.05", "pval \u2265 0.05")) %>%
  ggplot(aes(x=Estimate, y=fixed)) +
  geom_point(aes(alpha = sig),  position = position_dodge(width = 0.6), size =2.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  scale_alpha_manual(values = c("pval < 0.05" = 1, "pval \u2265 0.05" = 0.3), name ="") +
  geom_errorbarh(aes(xmin = confint_min, xmax = confint_max, alpha = sig), 
                 linewidth = 0.4, height = 0, 
                 position = position_dodge(width = 0.6)) +
  facet_wrap(~source, scales = "free") +
  theme(panel.grid = element_line(linewidth = 0.1, colour = "grey"),
        strip.text.y.right = element_text(angle = 0)) +
  labs(x = "Estimate", y="")
ggsave(filename = "figs/lmer-all-metrics-by-education-years-num.png", bg = "white",
       width = 8, height = 8, units = "in", dpi = 360)
################################################################################
################################################################################
################################################################################
