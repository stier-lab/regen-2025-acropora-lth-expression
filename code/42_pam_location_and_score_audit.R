# Paired tip/base analysis and morphology audit for the September 21 follow-up.
# Locations are spatial subsamples, not interchangeable technical replicates.
# No morphology scores are corrected or carried forward here.
source(here::here("code", "00_setup.R"))

raw_pam <- suppressWarnings(read_csv(file.path(DATA_RAW, "pam", "PAM_data.csv"),
              col_types=cols(.default=col_character()), show_col_types = FALSE))
write_csv(problems(raw_pam),file.path(TBL_DIR,"42_pam_csv_parse_audit.csv"))
# Ragged CSV rows omit trailing notes; complete measurement pairs and exact
# agreement with the canonical cleaned means are checked below.
p <- raw_pam |> janitor::clean_names() |>
  filter(!is.na(id)) |>
  rename(fv_fm = matches("^fv_fm"), thicket = matches("^thicket")) |>
  mutate(fv_fm = coalesce(suppressWarnings(as.numeric(fv_fm)), as.numeric(y)/1000),
         location = str_to_lower(str_trim(location)),
         treatment = factor(as.numeric(treatment), c(28,31), c("28C","31C")),
         wound = factor(wound, c("no","yes")), thicket = factor(thicket),
         tank = factor(as.integer(tank)), id = factor(as.integer(id)), day=as.integer(day))
stopifnot(all(is.finite(p$fv_fm)), all(p$fv_fm > 0 & p$fv_fm < 1),
          all(p$location %in% c("top","bottom")),
          !anyDuplicated(p[c("id","day","location")]))
pairs <- p |> select(id, tank, treatment, wound, thicket, day, location, fv_fm) |>
  pivot_wider(names_from = location, values_from = fv_fm) |>
  mutate(gap = top-bottom, average = (top+bottom)/2)
stopifnot(nrow(pairs)==336, all(complete.cases(pairs)), n_distinct(pairs$id)==48)
old <- readRDS(file.path(DATA_PROC,"pam_clean.rds")) |>
  mutate(id = as.character(id))
check <- pairs |> mutate(id = as.character(id)) |>
  left_join(old |> select(id,day,fv_fm), by=c("id","day"))
stopifnot(all(abs(check$average-check$fv_fm)<1e-12))
write_csv(pairs, file.path(TBL_DIR,"42_pam_location_pairs.csv"))

# Exact tank-label tests preserve temperature's experimental unit. The primary
# contrast averages the six post-clipping visits; day-specific tests are secondary.
tank_day <- pairs |> group_by(tank,treatment,day) |>
  summarise(gap=mean(gap), .groups="drop")
exact_heat <- function(d, scope) {
  stopifnot(nrow(d)==8, sum(d$treatment=="31C")==4)
  observed <- mean(d$gap[d$treatment=="31C"])-mean(d$gap[d$treatment=="28C"])
  null <- apply(combn(8,4),2,function(j) mean(d$gap[j])-mean(d$gap[-j]))
  tt <- t.test(gap ~ treatment, data=d)
  tibble(scope=scope, estimate_heated_minus_ambient=observed,
         lower=-tt$conf.int[2],upper=-tt$conf.int[1],
         p_exact=mean(abs(null)>=abs(observed)-1e-12),n_tanks=8)
}
tank_mean <- tank_day |> filter(day>=0) |> group_by(tank,treatment) |>
  summarise(gap=mean(gap),.groups="drop")
heat_tests <- bind_rows(exact_heat(tank_mean,"Post-clipping average"),
  map_dfr(sort(unique(tank_day$day)), function(day_i)
    exact_heat(filter(tank_day,day==day_i),paste("Day",day_i)))) |>
  mutate(p_holm_secondary=c(NA,p.adjust(p_exact[-1],"holm")))
write_csv(heat_tests,file.path(TBL_DIR,"42_pam_location_heat_tests.csv"))
tank_changes <- tank_day |> filter(day %in% c(0,14)) |>
  pivot_wider(names_from=day,values_from=gap,names_prefix="day") |>
  mutate(gap=day14-day0)
write_csv(exact_heat(tank_changes,"Change in gap: Day 14 minus Day 0"),
          file.path(TBL_DIR,"42_pam_location_change_test.csv"))
# Test whether the shape of the gap trajectory differs between temperatures,
# using whole tank trajectories rather than treating visits as independent.
curves <- tank_day |> filter(day>=0) |>
  pivot_wider(names_from=day,values_from=gap,names_prefix="day")
mat <- as.matrix(select(curves,starts_with("day")))
mat <- mat-rowMeans(mat)
stat <- function(j) sum((colMeans(mat[j,,drop=FALSE])-colMeans(mat[-j,,drop=FALSE]))^2)
obs <- stat(which(curves$treatment=="31C"))
write_csv(tibble(statistic=obs,p_exact=mean(apply(combn(8,4),2,stat)>=obs-1e-12),
                 assignments=70),file.path(TBL_DIR,"42_pam_location_trajectory_test.csv"))
site_changes <- pairs |> filter(day %in% c(0,14)) |>
  group_by(tank,treatment,day) |> summarise(across(c(top,bottom,average),mean),.groups="drop") |>
  pivot_longer(c(top,bottom,average),names_to="site",values_to="gap") |>
  pivot_wider(names_from=day,values_from=gap,names_prefix="day") |>
  mutate(gap=day14-day0)
write_csv(map_dfr(c("top","bottom","average"),function(s)
  exact_heat(filter(site_changes,site==s),s)),
  file.path(TBL_DIR,"42_pam_location_heat_decline_by_site.csv"))
robust_tanks <- pairs |> filter(day>=0) |> group_by(tank,treatment) |>
  summarise(gap=median(gap),.groups="drop")
write_csv(exact_heat(robust_tanks,"Median paired gap within each tank"),
          file.path(TBL_DIR,"42_pam_location_robust_test.csv"))
group_gaps <- tank_mean |> group_by(treatment) |>
  summarise(mean_gap=mean(gap),se=sd(gap)/sqrt(n()),n_tanks=n(),
            lower=mean_gap-qt(.975,n_tanks-1)*se,
            upper=mean_gap+qt(.975,n_tanks-1)*se,.groups="drop")
write_csv(group_gaps,file.path(TBL_DIR,"42_pam_location_group_gaps.csv"))

# Repeated-measures model of within-fragment differences; day is categorical.
# The exact tank-level test above is the primary temperature comparison.
post <- pairs |> filter(day>=0) |> mutate(day_f=factor(day)) |>
  arrange(tank,id,day)
m <- lmerTest::lmer(gap ~ treatment*wound*day_f + thicket +
                      (1|tank)+(1|id),data=post,REML=TRUE)
saveRDS(m,file.path(MOD_DIR,"42_pam_location_gap_lmm.rds"))
write_csv(as.data.frame(anova(m)) |> rownames_to_column("term"),
          file.path(TBL_DIR,"42_pam_location_model_tests.csv"))
emm <- emmeans(m,~treatment|wound*day_f)
write_csv(as_tibble(summary(contrast(emm,"revpairwise"),infer=TRUE,adjust="none")) |>
            mutate(p_holm=p.adjust(p.value,"holm")),
          file.path(TBL_DIR,"42_pam_location_model_contrasts.csv"))
res <- post |> mutate(fitted=fitted(m),residual=residuals(m))
lag_pairs <- res |> group_by(id) |> arrange(day,.by_group=TRUE) |>
  mutate(previous=lag(residual)) |> ungroup() |> filter(!is.na(previous))
diagnostics <- tibble(model="paired-gap LMM",singular=lme4::isSingular(m),
  convergence=paste(m@optinfo$conv$lme4$messages,collapse="; "),
  residual_sd=sd(res$residual),
  adjacent_visit_residual_correlation=cor(lag_pairs$residual,lag_pairs$previous),
  max_absolute_standardized_residual=max(abs(res$residual/sigma(m))))
write_csv(diagnostics,file.path(TBL_DIR,"42_pam_location_diagnostics.csv"))
write_csv(res |> group_by(treatment,day) |> summarise(sd_residual=sd(residual),
            .groups="drop"),file.path(TBL_DIR,"42_pam_location_residual_spread.csv"))
loo <- map_dfr(levels(tank_mean$tank),function(t) {
  z <- filter(tank_mean,tank!=t)
  tibble(omitted_tank=t, heat_gap=mean(z$gap[z$treatment=="31C"])-
           mean(z$gap[z$treatment=="28C"]))
})
write_csv(loo,file.path(TBL_DIR,"42_pam_location_leave_one_tank_out.csv"))
# Unequal residual spread and temporal dependence sensitivity, same fixed effects.
ar <- tryCatch(nlme::lme(gap ~ treatment*wound*day_f + thicket,
  random=~1|tank/id, correlation=nlme::corCAR1(form=~day|tank/id),
  weights=nlme::varIdent(form=~1|treatment),data=post,method="REML",
  control=nlme::lmeControl(opt="optim",msMaxIter=300)),error=function(e)e)
if(inherits(ar,"error")) {
  write_csv(tibble(status="failed",message=conditionMessage(ar)),
            file.path(TBL_DIR,"42_pam_location_correlated_sensitivity.csv"))
} else {
  saveRDS(ar,file.path(MOD_DIR,"42_pam_location_correlated_lme.rds"))
  write_csv(as_tibble(summary(contrast(emmeans(ar,~treatment),"revpairwise"),infer=TRUE)),
            file.path(TBL_DIR,"42_pam_location_correlated_sensitivity.csv"))
}

plot_data <- p |> mutate(site=if_else(location=="top","Near tip","Near base"),
                         injury=if_else(wound=="yes","Clipped","Unwounded")) |>
  group_by(tank,treatment,injury,day,site) |>
  summarise(score=mean(fv_fm),.groups="drop")
means <- plot_data |> group_by(treatment,injury,day,site) |>
  summarise(score=mean(score),.groups="drop")
p1 <- ggplot(plot_data,aes(day,score,colour=site)) +
  geom_line(aes(group=interaction(tank,site)),alpha=.22,linewidth=.35) +
  geom_line(data=means,linewidth=.8)+geom_point(data=means,size=1.4)+
  facet_grid(injury~treatment)+scale_colour_manual(values=c("Near tip"="#009E73","Near base"="#7B3294"))+
  scale_x_continuous(breaks=c(-1,0,3,6,9,12,14))+
  labs(x="Day relative to clipping",y="Photosynthesis efficiency (Fv/Fm)",
       colour=NULL,title="Photosynthesis near the tip and base")+theme_pub(10)
gap_plot <- pairs |> mutate(injury=if_else(wound=="yes","Clipped","Unwounded")) |>
  group_by(tank,treatment,injury,day) |> summarise(gap=mean(gap),.groups="drop")
gap_means <- gap_plot |> group_by(treatment,injury,day) |>
  summarise(gap=mean(gap),.groups="drop")
p2 <- ggplot(gap_plot,aes(day,gap,colour=treatment))+
  geom_hline(yintercept=0,linetype="dashed",colour="grey50")+
  geom_line(aes(group=tank),alpha=.35,linewidth=.4)+
  geom_line(data=gap_means,linewidth=.9)+geom_point(data=gap_means,size=1.5)+
  facet_wrap(~injury)+scale_colour_manual(values=PAL_TEMP)+
  scale_x_continuous(breaks=c(-1,0,3,6,9,12,14))+
  labs(x="Day relative to clipping",y="Near-tip minus near-base score",colour="Temperature",
       title="Above zero: higher score near the tip")+theme_pub(10)
save_fig(p1/p2,"42_pam_location_comparison",width=190,height=215)
q1 <- ggplot(res,aes(fitted,residual,colour=treatment))+geom_point(alpha=.5)+
  geom_hline(yintercept=0)+scale_colour_manual(values=PAL_TEMP)+theme_pub(10)
q2 <- ggplot(res,aes(sample=residual))+stat_qq()+stat_qq_line()+theme_pub(10)
save_fig(q1+q2,"42_pam_location_diagnostics",width=180,height=90)

# Audit all wounded fragments, including missing calls, without imposing an
# irreversible biological sequence or treating missing observations as absence.
d <- readRDS(file.path(DATA_PROC,"physio_clean.rds")) |> filter(wound=="yes")
stopifnot(!anyDuplicated(d[c("id","day")]))
missing <- d |> group_by(treatment,day) |>
  summarise(n=n(),pigment_missing=sum(is.na(pigment_over_wound)),
            tip_missing=sum(is.na(tip_exist)),.groups="drop")
write_csv(missing,file.path(TBL_DIR,"42_morphology_missing_calls.csv"))
review <- d |> group_by(id,treatment,tank,thicket) |> arrange(day,.by_group=TRUE) |>
  summarise(first_tip=if(any(tip_exist==1,na.rm=TRUE)) min(day[which(tip_exist==1)]) else NA_real_,
    tip_positive_days=paste(day[which(tip_exist==1)],collapse=","),
    absent_after_first=paste(day[which(tip_exist==0 & day>first_tip)],collapse=","),
    terminal_absence=any(tip_exist==1,na.rm=TRUE) &&
      any(tip_exist==0 & day>max(day[which(tip_exist==1)]),na.rm=TRUE),
    pigment_missing_days=paste(day[is.na(pigment_over_wound)],collapse=","),.groups="drop")
write_csv(review,file.path(TBL_DIR,"42_morphology_fragment_review.csv"))
print(group_gaps); print(heat_tests); print(diagnostics); print(loo)
print(review |> filter(absent_after_first!="")); print(missing |> filter(pigment_missing>0))
