library(Rcpp) 
library(RcppGSL)
library(HDInterval)
library(ggplot2)
library(data.table)
library(socialmixr)
library(lubridate)
library(stringr)
library(Hmisc)
library(extraDistr)
library(nloptr)
library(qs)
library(rlang)
library(readxl)
library(cowplot)
library(spatstat)

txt_theme = theme(plot.title = element_text(face = "plain", size = 7, hjust = 0.5))
path = function(x, prefix = "~/Dropbox/nCoV/Analyses/") { paste0(prefix, x); }

ccols = c("Milan" = "#623b84", "Birmingham" = "#518c9c", "Bulawayo" = "#a3d55a")
c_colours = function(a = c("colour", "fill")) scale_colour_manual(values = ccols, aesthetics = a)


# covidm options
cm_path = "~/Dropbox/nCoV/covidm/";
cm_force_rebuild = F;
cm_verbose = F
if (Sys.info()["nodename"] %like% "lshtm") {
   cm_build_dir = paste0(cm_path, "build/lshtm");
}
source(path("R/covidm.R", cm_path))

# Plot one or several contact matrices. mat must be a named list
plot_matrices = function(titl, mat, legpos = "right")
{
    data = NULL;
    mNS = mat$home + mat$work + mat$other;
    mS = mat$school;
    
    mm = rbind(
        cbind(data.table(reshape2::melt(mNS)), name = "Non-school"),
        cbind(data.table(reshape2::melt(mS)), name = "School")
    );
    
    names(mm) = c("Participant", "Contact", "N", "Matrix");

    ggplot(mm) + 
        geom_tile(aes(x = Contact, y = Participant, fill = N)) + 
        facet_wrap(~Matrix) + 
        scale_fill_viridis_c() +
        labs(title = titl, x = "Age of individual", y = "Age of contacts", fill = "Number of\ncontacts") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.background = element_rect(fill = "white", colour = "white"),
              legend.position = legpos)
}

theme_set(theme_cowplot(font_size = 5, font_family = "Helvetica", line_size = 0.25))

pm1 = plot_matrices("Wuhan (derived from Zhang et al.)", cm_matrices[["China | Wuhan"]])
pm2 = plot_matrices("Beijing (derived from Zhang et al.)", cm_matrices[["China | Beijing"]])
pm3 = plot_matrices("Shanghai (derived from Zhang et al.)", cm_matrices[["China | Shanghai"]])
pm4 = plot_matrices("Japan (from Prem et al.)", cm_matrices[["Japan"]])
pm5 = plot_matrices("South Korea (from Prem et al.)", cm_matrices[["Republic of Korea"]])
pm6 = plot_matrices("Singapore (from Prem et al.)", cm_matrices[["Singapore"]])
pm7 = plot_matrices("Canada (from Prem et al.)", cm_matrices[["Canada"]])

pm8 = plot_matrices("Milan (derived from Mossong et al.)", cm_matrices[["Italy | Milan"]])
pm9 = plot_matrices("Birmingham (derived from Mossong et al.)", cm_matrices[["UK | Birmingham"]])
pm10 = plot_matrices("Bulawayo (derived from Melegaro et al.)", cm_matrices[["Zimbabwe | Bulawayo"]])


plot_grid(pm1, pm2, pm3, pm4, pm5, pm6, pm7, pm8, pm9, pm10, nrow = 5, ncol = 2, labels = letters[1:10], label_size = 12)
ggsave("~/Dropbox/nCoV/Submission/Supp Figs/S-matrices.pdf", width = 30, height = 40, unit = "cm", useDingbats = F)
