#### CPR Data Anomaly Explorer  ####
# Adam A. Kemberling
# 9/22/2020
# App for displaying fits and metrics of seasonal splines from cpr data
# zooplankton data is from continuous plankton recorder
# Full time series a combination of NOAA and SAHFOS data products

# Run app from command line with: 
# docker run -p 3838:3838 --name cpr gmri/bashful-badger:1.0

####  Packages  ####
library(lubridate)
library(ggpmisc)
library(mgcv)
library(patchwork)
library(sf)
library(shiny)
library(tidyverse)
library(here)
library(gmRi)


####  Functions  ####

# Load CPR Helper Functions
source(here::here("R/cpr_helper_funs.R"))

# Use GMRI Stylesheet
gmRi::use_gmri_style_shiny(stylesheet = "gmri rmarkdown")

#Set ggplot theme
theme_set(theme_bw() + theme(plot.background = element_rect(fill = "transparent"), 
                             legend.box.background = element_rect(fill = "transparent"),
                             legend.title = element_blank(),
                             legend.position = "bottom",
                             axis.text = element_text(size = 12, color = "gray10"),
                             axis.title = element_text(size = 14, color = "gray10"), 
                             legend.text = element_text(size = 12, color = "gray10")))

# Coastlines from rnaturalearth
# northeast <- rnaturalearth::ne_states(country = "united states of america") %>% st_as_sfc(crs = 4326)
# canada <- rnaturalearth::ne_states(country = "canada") %>% st_as_sfc(crs = 4326)

# save them, load them from Data/
# sf::write_sf(obj = northeast, dsn = here::here("Data/rnatearth_ne_usa.geojson"))
# sf::write_sf(obj = canada, dsn = here::here("Data/rnatearth_canada.geojson"))
northeast <- read_sf(here::here("Data/rnatearth_ne_usa.geojson"))
canada    <- read_sf(here::here("Data/rnatearth_canada.geojson"))

####_________________####
####  Load and Prep CPR Data  ####




####__ Gulf of Maine CPR  __####
# Combined dataset from NOAA/SAHFOS, concentrations in common units: # / meters cubed
cpr <- read_csv(here::here("Data/zooplankton_combined.csv"), 
                guess_max = 1e6, col_types = cols())



#eliminate sample_id column that is mixed in, and format the dates
cpr <- cpr %>% 
    mutate(
        cal_date = as.Date(str_c(year, month, day, sep = "-")),
        jday = lubridate::yday(cal_date)) %>% 
    rename(
        lat = `latitude (degrees)`, 
        lon = `longitude (degrees)`) %>% 
    mutate(lon = ifelse(lon > 0, lon * -1, lon)) %>% 
    select(`Data Source`, cruise, station, year, month, day, hour, minute, lat, lon, cal_date, jday,
           `phytoplankton color index`, everything())





####  1. Make Taxa Lists  ####


#Identify the columns that represent abundances
start_name <- which(names(cpr) == "acartia danae")
taxa_cols <- names(cpr)[start_name:ncol(cpr)]
names(taxa_cols) <- taxa_cols


# Make a list with details on each taxa
taxa_list <- map(taxa_cols, function(x){
    taxa_name <- sym(x)
    taxa_subset <- cpr %>%
        select(year, jday, lat, lon, abundance = !!taxa_name)

}) 





#### 2. Remove Incomplete Time Series  ####

#Find those pesky NA taxa
na_counts <- map(taxa_list, function(x){ sum(is.na(x$abundance))}) %>% 
    bind_cols() %>%
    pivot_longer(names_to = "taxa", values_to = "total NA's", cols = everything()) %>%
    mutate(status = case_when(
        `total NA's` == 290 ~ "NOAA only",
        `total NA's` == 4799 ~ "SAHFOS only",
        `total NA's` == 0 ~ "Found in both",
        TRUE ~ "Only NA's"
    ))


#Taxa with full time series
keepers <- filter(na_counts, status == "Found in both")
fullts_taxa <- taxa_list[names(taxa_list) %in% keepers$taxa]


####  3. Selection Options  ####

# Update the list of choices to full time series only
taxa_names <- taxa_cols[which(names(taxa_list) %in% keepers$taxa)]

# List of Display options for each
display_names <- c(
    "Spatial Distribution"      = "map_plot",
    "Abundance Timeline"        = "timeline_plot",
    "Seasonal Variation"        = "seasonal_spline",
    "Seasonal Spline Residuals" = "resid_hist",
    "Anomaly Timeline"          = "anom_plot")



####  4. Calculate Detrended Abundances  ####
anomaly_list <- map(.x = fullts_taxa,
                    .f = cpr_spline_fun,
                    spline_bins = 10,       # k for smoothers
                    season_bins = 4,        # date breaks in a year (seasons/months)
                    study_area = "GOM_new")





####  5. Format Data for Plots  ####

# Add a new date column for displaying month x labels
# Format factor levels
anomaly_list <- map(anomaly_list, function(x){
    # Provide a "flat_date" to label with months for seasonal cycle
    base_date <- base_date <- as.Date("2000-01-01")
    x$cprdat_predicted <- x$cprdat_predicted %>%
        mutate(`Taxa Presence` = ifelse(abundance > 0, "Present", "Absent"),
               datebounds = factor(datebounds, levels = c("1-365", "0-92", "92-184", "184-276", "276-365")),
               `Anomaly Direction` = ifelse(anomaly > 0, "Positive Anomaly", "Negative Anomaly"),
               date = as.Date(str_c(year, "-01-01")),
               date = date + jday,
               flat_date = as.Date(jday - 1, origin = base_date))
    
    
    x$period_summs <- x$period_summs %>%
        mutate(
            datebounds = factor(datebounds, levels = c("1-365", "0-92", "92-184", "184-276", "276-365")),
            `Anomaly Direction` = ifelse(period_anom_mu > 0, "Positive Anomaly", "Negative Anomaly"))
    
    
    return(x)
})


####  6. Phytoplankton Index  ####

# Make a Phytoplankton Color index list as well
phyto <- cpr %>% 
    select(year, jday, lat, lon, abundance = 'phytoplankton color index') %>% 
    mutate(`Taxa Presence` = ifelse(abundance > 0, "Present", "Absent"),
           date = as.Date(str_c(year, "-01-01")),
           date = date + jday,
           flat_date = as.Date(jday - 1, origin = as.Date("2000-01-01")))


# Phytoplankton index plot
phyto_avg <- phyto %>% 
    group_by(date) %>% 
    summarise(mean_abund = mean(abundance, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(jday = lubridate::yday(date),
           flat_date = as.Date(jday - 1, origin = as.Date("2000-01-01")))


# Color key
phyto_color_key <- c("Individual Transect Abundance" = "gray",
                     #"Daily Mean Abundance" = "black")
                     #"Daily Mean Abundance" = "#00608A")
                     "Daily Mean Abundance" = "#00608A")

# Plot it
phyto_season <- ggplot(phyto) +
    geom_point(aes(flat_date, abundance, color = "Individual Transect Abundance"), alpha = 0.5, size = 0.5) +
    geom_point(data = phyto_avg,
               aes(flat_date, mean_abund, color = "Daily Mean Abundance"), size = 1) +
    geom_smooth(aes(flat_date, abundance),
                formula = y ~ s(x, bs = "cc", k = 10),
                method = "gam",
                #color = "lime green",
                color = gmri_cols("gmri green"),
                size = 2) +
    scale_x_date(date_labels = "%b", date_breaks = "1 month") +
    scale_color_manual(values = phyto_color_key) +
    labs(x = "", y = "Phytoplankton Abundance Index")








####_____####
####__ Plot Assembly  __####


# From there we have the:
# 1. Predictions with actual original values
# 2. Summaries over 91-day quarters
# 3. the model itself


####  1. Distribution Map  ####
anomaly_list <- map(anomaly_list, function(x){
    
    # Make Simple Features Data set
    data <- x$cprdat_predicted 
    data_sf <- data %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

    # Build Map
    x$map_plot <- ggplot() +
        geom_sf(data = filter(data_sf, `Taxa Presence` == "Present"), aes(shape = `Taxa Presence`, color = `Taxa Presence`), alpha = 0.8) +
        geom_sf(data = filter(data_sf, `Taxa Presence` == "Absent"), aes(shape = `Taxa Presence`, color = `Taxa Presence`), alpha = 0.6) +
        geom_sf(data = northeast) +
        geom_sf(data = canada) +
        scale_shape_manual(values = c("Absent" = 3, "Present" = 16)) +
        scale_color_manual(values = c("Absent" = "black", "Present" = "royalblue")) +
        guides(fill = guide_legend(title = "", label = FALSE),
               shape = guide_legend(title.position = "top", title.hjust = 0.5)) +
        coord_sf(xlim = c(-71, -66), ylim = c(41.5, 44.3)) +
        facet_wrap(~datebounds)

    return(x)


})

# # Tester
# anomaly_list$`calanus finmarchicus v-vi`$map_plot




#### 2. Abundance Timelines  ####
anomaly_list <- map(anomaly_list, function(x){

    # Grab data
    timeline_data <- x$cprdat_predicted 
    
    # Get daily_averages to plot vs. raw concentrations
    daily_avgs <- timeline_data %>% 
        group_by(date) %>% 
        summarise(mean_conc = mean(abundance, na.rm = T), .groups = "keep")

    # Color key
    color_key <- c("Individual Transect Abundance" = "gray",
                   "Daily Mean Abundance" = "#00608A")
    
    # gam period, how wiggly do we want it
    gam_period <- round(length(unique(timeline_data$year)) /5, digits = 0)
    
    
    # Build Plot
    x$timeline_plot <- timeline_data %>%
        ggplot(aes(date, abundance) ) +
        geom_point(aes(color = "Individual Transect Abundance"), size = 0.5, alpha = 0.5) +
        geom_point(data = daily_avgs, 
                   aes(date, mean_conc, color = "Daily Mean Abundance"), size = 1) +
        geom_smooth(method = "gam", 
                    formula = y ~ s(x, bs = "cs", k = gam_period),
                    color = gmri_cols("orange", as_char = TRUE)) +
        scale_y_log10(labels = scales::comma_format()) +
        scale_x_date(limits = as.Date(c("1959-01-01", "2020-01-01")), date_labels = "%Y", date_breaks = "5 year") +
        scale_color_manual(values = color_key) +
        labs(x = NULL, y = "Observed Abundance (individuals / cubic meter)") 

    return(x)


})

# # Test Plot
# anomaly_list$`calanus finmarchicus v-vi`$timeline_plot




#### 3. Seasonal Splines  ####
anomaly_list <- map(anomaly_list, function(x){
    timeline_data <- x$cprdat_predicted
        
    
    # Daily averages - actual date not day of year, build back the flat date
    daily_avg <- timeline_data %>% 
        group_by(date) %>% 
        summarise(mean_conc = mean(abundance),
                  .groups = "keep") %>% 
        ungroup() %>% 
        mutate(jday = lubridate::yday(date),
               flat_date = as.Date(jday - 1, origin = as.Date("2000-01-01")))
    
    # Color key
    color_key <- c("Individual Transect Abundance" = "gray",
                   "Daily Mean Abundance" = "#00608A")
    
    
    # Seasonal Pattern
    zoo_seasonal <- timeline_data %>% 
        ggplot(aes(flat_date, abundance)) +
        geom_point(aes(color = "Individual Transect Abundance"), size = 0.5, alpha = 0.5) +
        geom_point(data = daily_avg, aes(flat_date, mean_conc, color = "Daily Mean Abundance"), size = 1) +
        geom_smooth(method = "gam",
                    formula = y ~  s(x, bs = "cc", k = 10),
                    color = gmri_cols("orange", as_char = T),
                    size = 2) +
        
        scale_y_log10(labels = scales::comma_format()) +
        scale_x_date(date_labels = "%b", date_breaks = "1 month") +
        scale_color_manual(values = color_key) +
        labs(x = "", y = "Concentration (# / cubic meter)") +
        theme(legend.position = "none")
    
    
    x$seasonal_spline <- (zoo_seasonal / phyto_season)
    
    return(x)
})    

# Tester    
# anomaly_list$`calanus i-iv`$seasonal_spline
    

####  4. Seasonal Anomaly Timeseries  ####
anomaly_list <- map(anomaly_list, function(x){
    
    # Pull the period summaries
    timeline_data <- x$period_summs 
    
    # gam period, how wiggly do we want it
    gam_period <- round(length(unique(timeline_data$year)) /5, digits = 0)
    
    anom_facets <- function(timeline_data){
        timeline_data %>%
        ggplot(aes(year, period_anom_mu)) +
        geom_hline(yintercept = 0, linetype = 2, color = "gray60") +
        geom_smooth(method = "gam", 
                    formula = y ~ s(x, bs = "cs", k = gam_period),
                    size = 2, 
                    color = "gray50") +
        geom_point(aes(color = `Anomaly Direction`), size = 1) +
        scale_color_manual(values = c("Positive Anomaly" = as.character(gmri_cols("orange")),
                                      "Negative Anomaly" = as.character(gmri_cols("gmri blue")))) +
        facet_wrap(~datebounds, ncol = 2) +
        xlim(c(1960, 2010)) +
        labs(x = NULL, y = "Seasonal Anomaly Magnitude") +
        guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) +
        theme(legend.position = "none")}
    
    
    full_year <- timeline_data %>% 
        filter(datebounds == "1-365") %>% 
        mutate(datebounds = "Annual Average") %>% 
        anom_facets()
    four_seasons <- timeline_data %>% 
        filter(datebounds != "1-365") %>% 
        mutate(datebounds = paste0("Day ", datebounds),
               datebounds = factor(datebounds,
                levels = c("Day 0-92", "Day 92-184", "Day 184-276", "Day 276-365"))) %>% 
        anom_facets()
    
    
    
    x$anom_plot <- full_year | four_seasons
    
    return(x)
    
}) 


# Tester
anomaly_list$`calanus i-iv`$anom_plot





####  5. Histogram of residuals  ####
anomaly_list <- map(anomaly_list, function(x){

    
    residual_data <- tibble(`Spline Residuals` = x$spline_model$residuals)
    
    x$resid_hist <- ggplot(residual_data) +
        geom_histogram(aes(`Spline Residuals`), bins = 30, fill = gmri_cols("gmri blue")) +
        geom_vline(xintercept = 0, linetype = 2, color = gmri_cols("orange"), size = 2) +
        labs(x = "Residuals from Seasonal Average - As Predicted from GAM",
             y = "Count")
    
    return(x)
    
    
})


# Tester
anomaly_list$`calanus i-iv`$resid_hist


####_________________####
####  Define UI ####
ui <- fluidPage(

    ####  1. Page Details  ####
    title = "CPR Explorer",
    theme = "styles/gmri_rmarkdown.css",

    ####  2. Application Banner with Logo  ####
    fluidRow(style = "padding: 20px;",

             column(12,
                    img(src = "GMRI_logo.png")
             )
    ),

    ####  3. Title and Introduction  ####

    fluidRow(style = 'padding-bottom: 20px;',

             #Buffer column
             column(1),
             column(width = 10,
                tags$h1("Exploring Zooplankton Distribution and Abundance Patterns"),
                tags$hr(),
                
                tags$h2("Continuous Plankton Recorder Survey: Gulf of Maine"),
                tags$br(),
                
                tags$h3(strong("About the Data:")),
                tags$p("The Gulf of Maine Continuous Plankton Recorder Survey is a scientific resource
                       used to identify spatial and temporal patterns of zooplankton in near-surface waters (<50m).
                       The CPR sampler is towed behind a ship of opportunity several times a year. Transects extend
                       from Nova Scotia to either Boston MS, or more recently Portland ME. Data for the 
                       years 1961 to 2013 were received through communications with NOAA scientists, 
                       with data from 2014-2017 received through communication with SAHFOS. Once received, zooplankton 
                       concentration data from both sources were converted to the same measurement 
                       units (zooplankton density per cubic meter of seawater), 
                       and all differences and inconsistencies in taxon identification records were compared to ensure direct
                       matches between taxon-stage groups. The two CPR datasets were then joined to create a complete 
                       time-series covering the years 1961-2017."),
                
                
                tags$h3(strong("Navigating the Displays:")),
                tags$p("The CPR Web Explorer app has several display options for viewing spatial and 
                       temporal variations in abundances for individual zooplankton taxa. To explore
                       different patterns in the zooplankton data, begin by choosing a taxa for each plot,
                       and the display for the data.")
                ),
             column(1)
    ),

    ####  Selection Options  ####
    fluidRow(style = 'padding-bottom: 20px;',

             #Buffer column
             column(1),
             column(4,
                    style = "color: black",
                    align = "left",


                    # Taxa 1 Selector
                    selectInput(inputId = "taxa1",
                                label   = "1. Select a Taxa",
                                choices = taxa_names),

                    # Buoy 1 Selector
                    selectInput(inputId = "display1",
                                label   = "2. Select a Display",
                                choices = display_names),

                    actionButton(inputId = "left_button", label = "Update Left-Side", icon("refresh"))

             ),
             column(2),

             column(4,
                    style = "color: black",
                    align = "left",


                    # Taxa 2 Selector
                    selectInput(inputId = "taxa2",
                                label   = "1. Select a Taxa",
                                choices = taxa_names),

                    # Buoy 2 Selector
                    selectInput(inputId = "display2",
                                label   = "2. Select a Display",
                                choices = display_names),

                    actionButton(inputId = "right_button", label = "Update Right-Side", icon("refresh"))

             ),

             column(1)





    ),

    ####  Plotting Space  ####
    fluidRow(
        column(width = 6,
               wellPanel(
                   align = "center",
                   style = 'float: center;',
                   plotOutput("left_plot",
                              height = "600px",
                              width = "auto")
               )
        ),
        column(width = 6,
               wellPanel(
                   align = "center",
                   style = 'float: center;',
                   plotOutput("right_plot",
                              height = "600px",
                              width = "auto")
               )
        )

    ),

    ####  Additional Comments  ####
    fluidRow(
        column(1),
        column(10,
               tags$strong("NOTE:"),
               tags$p("The data and figures above are a product of CPR data collected from both NOAA
                       and the Sir Alister Hardy Foundation (SAHFOS). Pre-processing to common units, and consistent
                       group classifications occurred prior to any plotting. Spatial extent was also limited to the
                       study area of the Gulf of Maine."),
               # Add footer with github info
              gmRi::insert_gmri_footer(footer_file = "akemberling_gmri_footer.html")),
        column(1)
    )
)

####_________________####
####  Define Server Logic  ####
server <- function(input, output) {


    ####  Observe Left Side  ####
    observeEvent(input$left_button, {
        ####Left Side Plot  ####
        output$left_plot <- renderPlot({
            #ggplot()
            anomaly_list[[input$taxa1]][input$display1]
        }, bg="transparent")
    })

    ####  Observe Right Side  ####
    observeEvent(input$right_button, {
        ####Right Side Plot  ####
        output$right_plot <- renderPlot({
            #ggplot()
            anomaly_list[[input$taxa2]][input$display2]
        }, bg="transparent")
    })




}

#####  Run the application   ####
shinyApp(ui = ui, server = server)