#' ## relative change in FID ~ period
    d[, n_c_sp :=.N, by = sp_country]
    d[, FID_z :=scale(FID), by = sp_country]
    d[, SD_z :=scale(SD), by = sp_country]

    ms_z <- lmer(scale(FID_z) ~
      scale(SD_z) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Covid) +
      (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | Country) + (scale(Covid) | IDLocality) + (1 | sp_loc),
    data = d, REML = FALSE,
    control <- lmerControl(
      optimizer = "optimx", optCtrl = list(method = "nlminb")
    )
    )

    est_ms_z <- est_out(ms_z, "All: (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (1|sp_loc)")


#+ gsfig, fig.width=4.5, fig.height = 3.5
col2_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col2__ <- col2_[3:7]
ggplot(s, aes(x = StringencyIndex, y = parks_percent_change_from_baseline, col = Country)) +
  stat_smooth(method = "lm") +
  stat_cor(method = "pearson", size = 2) +
  geom_point() +
  scale_color_manual(values = col2__) +
  labs(subtitle = "simple lm & Pearson's R")


# Stringency

  # full
  mss <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
  data = s, REML = FALSE,
  control = lmerControl(
    optimizer = "optimx", optCtrl = list(method = "nlminb")
  )
  )
  est_mss <- est_out(mss, "ALL: (scale(StringencyIndex) |genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc)")
  est_mss[, control_for_starting_distance := "yes"]

  msx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex) | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
  data = s, REML = FALSE,
  control = lmerControl(
    optimizer = "optimx", optCtrl = list(method = "nlminb")
  )
  )
  est_msx <- est_out(msx, "ALL: (scale(StringencyIndex) |genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc)")
  est_msx[, control_for_starting_distance := "no"]


  # CZ - singular fits only due to genera estimated as zero (removing it changes no results)
  css <- lmer(scale(log(FID)) ~
    # scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Czechia"], REML = FALSE
  )
  est_css <- est_out(css, "Czechia: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_css[, control_for_starting_distance := "yes"]

  csx <- lmer(scale(log(FID)) ~
    # scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Czechia"], REML = FALSE
  )
  est_csx <- est_out(csx, "Czechia: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_csx[, control_for_starting_distance := "no"]

  # FI
  fss <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  data = s[Country == "Finland"], REML = FALSE
  )
  est_fss <- est_out(fss, "Finland: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(scale(StringencyIndex)|IDLocality)+(1|sp_loc)")
  est_fss[, control_for_starting_distance := "yes"]

  fsx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  data = s[Country == "Finland"], REML = FALSE
  )
  est_fsx <- est_out(fsx, "Finland: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(scale(StringencyIndex)|IDLocality)+(1|sp_loc)")
  est_fsx[, control_for_starting_distance := "no"]

  # HU - singular fits only due to sp_loc estimated as zero (removing it changes no results)
  hss <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Hungary"], REML = FALSE
  )
  est_hss <- est_out(hss, "Hungary: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_hss[, control_for_starting_distance := "yes"]

  hsx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  # (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Hungary"], REML = FALSE
  )
  est_hsx <- est_out(hsx, "Hungary: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_hsx[, control_for_starting_distance := "no"]

  # AU - singular fits only due to Year and random slope estimated as zero (removing those changes no results)
  ass <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Australia"], REML = FALSE
  )
  est_ass <- est_out(ass, "Australia: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_ass[, control_for_starting_distance := "yes"]

  asx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
  # (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Australia"], , REML = FALSE,
  control = lmerControl(
    optimizer = "optimx", optCtrl = list(method = "nlminb")
  )
  )
  est_asx <- est_out(asx, "Australia: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_asx[, control_for_starting_distance := "no"]

  # PL
  pss <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Poland"], REML = FALSE
  )
  est_pss <- est_out(pss, "Poland: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)")
  est_pss[, control_for_starting_distance := "yes"]

  psx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (scale(StringencyIndex)  | genus) + (1 | Species) + (1 | sp_day_year),
  # (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s[Country == "Poland"], REML = FALSE
  )
  est_psx <- est_out(psx, "Poland: (scale(StringencyIndex) |genus)+(1|Species)+(1|sp_day_year)")
  est_psx[, control_for_starting_distance := "no"]

  # combine
  est_mss[, Country := "All\n(mixed model)"]
  est_msx[, Country := "All\n(mixed model)"]
  est_ass[, Country := "Australia"]
  est_asx[, Country := "Australia"]
  est_css[, Country := "Czechia"]
  est_csx[, Country := "Czechia"]
  est_hss[, Country := "Hungary"]
  est_hsx[, Country := "Hungary"]
  est_pss[, Country := "Poland"]
  est_psx[, Country := "Poland"]
  est_fss[, Country := "Finland"]
  est_fsx[, Country := "Finland"]

  os = rbind(
    est_mss, est_msx,
    est_ass, est_asx,
    est_css, est_csx,
    est_hss, est_hsx,
    est_pss, est_psx,
    est_fss, est_fsx
  )
  save(os, file = here::here("Data/dat_est_Stringency_genus-r-slope__rev.Rdata"))

  #+ est_str, fig.width=3, fig.height = 2.5
  load(here::here("Data/dat_est_Stringency_genus-r-slope__rev.Rdata"))
  os[predictor %in% c("scale(StringencyIndex)"), predictor := "Stringency Index"]
  oso <- os[predictor %in% c("Stringency Index")]
  oso[, N := as.numeric(sub(".*N = ", "", model))]
  # add meta-analytical mean
  oso_s = oso[control_for_starting_distance == "yes"]
  met = summary(meta.summaries(d = oso_s$estimate, se = oso_s$sd, method = "fixed", weights = oso_s$N))$summci
  oso_met = data.table(predictor = "Period", estimate = met[2], lwr = met[1], upr = met[3], sd = NA, model = NA, control_for_starting_distance = "yes", Country = "Combined\n(metanalytical)", N = NA)

  oso_sx = oso[control_for_starting_distance == "no"]
  metx = summary(meta.summaries(d = oso_sx$estimate, se = oso_sx$sd, method = "fixed", weights = oso_sx$N))$summci
  oso_metx = data.table(predictor = "Period", estimate = metx[2], lwr = metx[1], upr = metx[3], sd = NA, model = NA, control_for_starting_distance = "no", Country = "Combined\n(metanalytical)", N = NA)

  oso = rbind(oso, oso_met, oso_metx)

  oso[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia", "Combined\n(metanalytical)", "All\n(mixed model)")))]

  width_ <- .5 # spacing between error bars

  # col_ <- c(brewer.pal(n = 12, name = "Paired"), "grey30", "grey80")
  # Tol_bright <- c("#EE6677", "#228833", "#4477AA", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
  # Tol_muted <- c("#88CCEE", "#44AA99", "#117733", "#332288", "#DDCC77", "#999933", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
  # Tol_light <- c("#BBCC33", "#AAAA00", "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#DDDDDD")

  # From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
  # Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  # col_ = Okabe_Ito[7:1]
  # JAMA and LocusZoom modified order
  # col_ =  c("#374E55FF", "#374E55FF", "#DF8F44FF", "#79AF97FF", "#00A1D5FF", "#B24745FF",  "#80796BFF") #"#6A6599FF",
  # col_ <- c("#357EBDFF", "#9632B8FF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#D43F3AFF", "#D43F3AFF")[7:1] # "#D43F3AFF", "#B8B8B8FF"
  col_ = c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1] # "#D43F3AFF", "#B8B8B8FF"
  # show_col(col_)

  # p =
  ggplot(oso, aes(x = estimate, y = Country, col = Country, shape = control_for_starting_distance)) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
    # geom_point(position = ggstance::position_dodgev(.6)) +
    geom_point(position = position_dodge(width = width_), bg = "white", size = 1.1) +
    # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
    # geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
    # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+

    scale_shape_manual(name = "Controlled for\nstarting distance", guide = guide_legend(reverse = TRUE), values = c(21, 19)) +
    # scale_color_jama(guide = "none")+ #, palette = 'light'
    scale_color_manual(guide = "none", values = col_) + # guide_legend(reverse = TRUE)
    scale_x_continuous(breaks = round(seq(-0.3, 0.4, by = 0.1), 1)) +
    ylab("") +
    xlab("Standardized effect size of Stringency Index\n[on flight initiation distance]") +
    # coord_cartesian(xlim = c(-.15, .15)) +
    # scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      # legend.spacing.y = unit(0.1, 'cm'),
      legend.key.height = unit(0.5, "line"),
      legend.margin = margin(0, 0, 0, 0),
      # legend.position=c(0.5,1.6),
      plot.title = element_text(color = "grey", size = 7),
      plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = ax_lines, size = 0.25),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
      # axis.text.x = element_text()
      axis.ticks.length = unit(1, "pt"),
      axis.text.x = element_text(, size = 6),
      axis.text.y = element_text(colour = "black", size = 7),
      axis.title = element_text(size = 7)
    )

  ggsave("Outputs/Fig_Stringency_rev_width_CustomLocusZoom_v2_genu-r-slope.png", width = 8, height = 6, unit = "cm", dpi = 600)

# using random slope
#' ## predictions for FID ~ Google Mobiliity
# predictions
# full
mgs <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | genus) + (1 | Species) + (1 | sp_day_year) +
  (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
data = ss, REML = FALSE,
control = lmerControl(
  optimizer = "optimx", optCtrl = list(method = "nlminb")
)
)
est_mgs <- est_out(mgs, "ALL: (scale(Google)|genus) + (1|Species)  + (1|sp_day_year) + (scale(Google)|Country) + (1|IDLocality) +(1|sp_loc)")
est_mgs[, control_for_starting_distance := "yes"]

mgx <- lmer(scale(log(FID)) ~
  scale(Year) +
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(StringencyIndex) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | genus) + (1 | Species) + (1 | sp_day_year) +
  (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
data = ss, REML = FALSE,
control = lmerControl(
  optimizer = "optimx", optCtrl = list(method = "nlminb")
)
)
est_mgx <- est_out(mgx, "ALL: (scale(Google)|genus) + (1|Species) + (1|sp_day_year) + (scale(Google)|Country) + (1|IDLocality) +(1|sp_loc)")
est_mgx[, control_for_starting_distance := "no"]


# CZ - singular fits only due to genera estimated as zero (removing it changes no results)
cgs <- lmer(scale(log(FID)) ~
  # scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Czechia"], REML = FALSE
)
est_cgs <- est_out(cgs, "Czechia: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_cgs[, control_for_starting_distance := "yes"]

cgx <- lmer(scale(log(FID)) ~
  # scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Czechia"], REML = FALSE
)
est_cgx <- est_out(cgx, "Czechia: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_cgx[, control_for_starting_distance := "no"]

# FI
fgs <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
data = ss[Country == "Finland"], REML = FALSE
)
est_fgs <- est_out(fgs, "Finland: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(scale(StringencyIndex)|IDLocality)+(1|sp_loc)")
est_fgs[, control_for_starting_distance := "yes"]

fgx <- lmer(scale(log(FID)) ~
  scale(Year) +
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
data = ss[Country == "Finland"], REML = FALSE
)
est_fgx <- est_out(fgx, "Finland: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(scale(StringencyIndex)|IDLocality)+(1|sp_loc)")
est_fgx[, control_for_starting_distance := "no"]

# HU - singular fits only due to sp_loc estimated as zero (removing it changes no results)
hgs <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Hungary"], REML = FALSE
)
est_hgs <- est_out(hgs, "Hungary: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hgs[, control_for_starting_distance := "yes"]

hgx <- lmer(scale(log(FID)) ~
  scale(Year) +
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Hungary"], REML = FALSE
)
est_hgx <- est_out(hgx, "Hungary: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hgx[, control_for_starting_distance := "no"]

# AU - singular fits only due to Year and random slope estimated as zero (removing those changes no results)
ags <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Australia"], REML = FALSE
)
est_ags <- est_out(ags, "Australia: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_ags[, control_for_starting_distance := "yes"]

agx <- lmer(scale(log(FID)) ~
  scale(Year) +
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Australia"], , REML = FALSE,
control = lmerControl(
  optimizer = "optimx", optCtrl = list(method = "nlminb")
)
)
est_agx <- est_out(agx, "Australia: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_agx[, control_for_starting_distance := "no"]

# PL
pgs <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Poland"], REML = FALSE
)
est_pgs <- est_out(pgs, "Poland: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)")
est_pgs[, control_for_starting_distance := "yes"]

pgx <- lmer(scale(log(FID)) ~
  scale(Year) +
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline)  | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Poland"], REML = FALSE
)
est_pgx <- est_out(pgx, "Poland: (scale(parks_percent_change_from_baseline) |genus)+(1|Species)+(1|sp_day_year)")
est_pgx[, control_for_starting_distance := "no"]

# combine
est_mgs[, Country := "All\n(mixed model)"]
est_mgx[, Country := "All\n(mixed model)"]
est_ags[, Country := "Australia"]
est_agx[, Country := "Australia"]
est_cgs[, Country := "Czechia"]
est_cgx[, Country := "Czechia"]
est_hgs[, Country := "Hungary"]
est_hgx[, Country := "Hungary"]
est_pgs[, Country := "Poland"]
est_pgx[, Country := "Poland"]
est_fgs[, Country := "Finland"]
est_fgx[, Country := "Finland"]

og <- rbind(
  est_mgs, est_mgx,
  est_ags, est_agx,
  est_cgs, est_cgx,
  est_hgs, est_hgx,
  est_pgs, est_pgx,
  est_fgs, est_fgx
)
save(og, file = here::here("Data/dat_est_Google_rev_genus-ran.Rdata"))

#+ est_str, fig.width=3, fig.height = 2.5
load(here::here("Data/dat_est_Google_rev_genus-ran.Rdata"))
og[predictor %in% c("scale(parks_percent_change_from_baseline)"), predictor := "Google Mobility"]
ogo <- og[predictor %in% c("Google Mobility")]
ogo[, N := as.numeric(sub(".*N = ", "", model))]
# add meta-analytical mean
ogo_s <- ogo[control_for_starting_distance == "yes"]
met <- summary(meta.summaries(d = ogo_s$estimate, se = ogo_s$sd, method = "fixed", weights = ogo_s$N))$summci
ogo_met <- data.table(predictor = "Period", estimate = met[2], lwr = met[1], upr = met[3], sd = NA, model = NA, control_for_starting_distance = "yes", Country = "Combined\n(metanalytical)", N = NA)

ogo_sx <- ogo[control_for_starting_distance == "no"]
metx <- summary(meta.summaries(d = ogo_sx$estimate, se = ogo_sx$sd, method = "fixed", weights = ogo_sx$N))$summci
ogo_metx <- data.table(predictor = "Period", estimate = metx[2], lwr = metx[1], upr = metx[3], sd = NA, model = NA, control_for_starting_distance = "no", Country = "Combined\n(metanalytical)", N = NA)

ogo <- rbind(ogo, ogo_met, ogo_metx)

ogo[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia", "Combined\n(metanalytical)", "All\n(mixed model)")))]

width_ <- .5 # spacing between error bars

# col_ <- c(brewer.pal(n = 12, name = "Paired"), "grey30", "grey80")
# Tol_bright <- c("#EE6677", "#228833", "#4477AA", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
# Tol_muted <- c("#88CCEE", "#44AA99", "#117733", "#332288", "#DDCC77", "#999933", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
# Tol_light <- c("#BBCC33", "#AAAA00", "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#DDDDDD")

# From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
# Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
# col_ = Okabe_Ito[7:1]
# JAMA and LocusZoom modified order
# col_ =  c("#374E55FF", "#374E55FF", "#DF8F44FF", "#79AF97FF", "#00A1D5FF", "#B24745FF",  "#80796BFF") #"#6A6599FF",
# col_ <- c("#357EBDFF", "#9632B8FF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#D43F3AFF", "#D43F3AFF")[7:1] # "#D43F3AFF", "#B8B8B8FF"
col_ = c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1] # "#D43F3AFF", "#B8B8B8FF"
# show_col(col_)

# p =
ggplot(ogo, aes(x = estimate, y = Country, col = Country, shape = control_for_starting_distance)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
  # geom_point(position = ggstance::position_dodgev(.6)) +
  geom_point(position = position_dodge(width = width_), bg = "white", size = 1.1) +
  # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
  # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
  # geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
  # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+

  scale_shape_manual(name = "Controlled for\nstarting distance", guide = guide_legend(reverse = TRUE), values = c(21, 19)) +
  # scale_color_jama(guide = "none")+ #, palette = 'light'
  scale_color_manual(guide = "none", values = col_) + # guide_legend(reverse = TRUE)
  scale_x_continuous(breaks = round(seq(-0.3, 0.2, by = 0.1), 1)) +
  ylab("") +
  xlab("Standardized effect size of Google Mobility (human presence)\n[on flight initiation distance]") +
  # coord_cartesian(xlim = c(-.15, .15)) +
  # scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.margin = margin(0, 0, 0, 0),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text.x = element_text(, size = 6),
    axis.text.y = element_text(colour = "black", size = 7),
    axis.title = element_text(size = 7)
  )

ggsave("Outputs/Fig_Google_rev_width_CustomLocusZoom_v2__genus-ran.png", width = 8, height = 6, unit = "cm", dpi = 600)