library(shiny)
library(bslib)
library(tidyverse)
library(geomtextpath)
library(evd)
library(rhandsontable)
library(rfishbase)
library(prompter)
library(scales)

# Define UI ----
ui <- page_sidebar(
    use_prompt(),
    title = "Reconstruction of unfished body size distributions",
    nav_panel("tt"),
    sidebar = sidebar(title = "Input parameters", 
                      width = "500px",
                      radioButtons("metric", "Available size metric:",
                                   c("Asymptotic (theoretical)", 
                                     "Maximum (observed)", 
                                     "Maximum (estimated from maxima)", 
                                     "Mean (observed)")),
                      conditionalPanel(
                          condition = "input.metric == 'Asymptotic (theoretical)'",
                          sliderInput("slider_linf", "Asymptotic length:",
                                      value = 60,
                                      min = 1,
                                      max = 200,
                                      post  = "cm"
                          )
                      ),
                      conditionalPanel(
                          condition = "input.metric == 'Maximum (observed)'",
                          sliderInput("slider_maxsize", "Maximum observed length:",
                                      value = 50,
                                      min = 1,
                                      max = 200,
                                      post  = "cm"
                          )
                      ),
                      
                      conditionalPanel(
                          condition = "input.metric == 'Mean (observed)'",
                          sliderInput("slider_meansize", "Mean observed length:",
                                      value = 25,
                                      min = 1,
                                      max = 200,
                                      post  = "cm"
                          )
                      ),
                      hr(),
                      checkboxGroupInput("dist", 
                                         label = "Distribution:",
                                         choiceNames = list(
                                             span("Truncated normal", style = "color: rgb(181, 144, 19);"),
                                             span("Log-normal", style = "color: rgb(29, 84, 128);")
                                         ),
                                         choiceValues = list(
                                             "p_norm", 
                                             "p_lnorm"
                                         ), 
                                         selected = "p_norm"),
                      hr(),
                      sliderInput("xaxis",
                                  label = "X-axis limit",
                                  value = c(0, 100),
                                  min = 0,
                                  max = 400,
                                  post  = "cm"
                      ), 
                      hr(),
                      p("Overlay observational data:"),
                      downloadButton("downloadData", "Download example data"),
                      fileInput("datafile", "Upload .csv file", accept = ".csv"),
                      uiOutput("species_dropdown"),
                      conditionalPanel(
                          condition = "input.spp_selected != 'None'",
                          radioButtons("population_scale", label = "At which scale?",
                                       choices = list("Species level" = "spp",
                                                      "Population (gridcell) level" = "pop"),
                                       selected = "spp"), 
                          radioButtons(inputId = "radio_scaled", 
                                       label = "Scale the observed data",
                                       choiceNames =  c(
                                           "No scaling (relative abundance)", 
                                           "Scale based on maximum value", 
                                           "Observed as a subsample", 
                                           "Choose own scaling"
                                       ), 
                                       choiceValues = c(
                                           "rel_abun", 
                                           "scale_max", 
                                           "scale_subsample",
                                           "choose_scaling"
                                       ), 
                                       selected = "rel_abun"
                          ),
                          conditionalPanel(
                              condition = "input.radio_scaled == 'choose_scaling'",
                              sliderInput(inputId = "y_scaling",
                                          label = "",
                                          value = 1,
                                          step = 0.01,
                                          min = 0,
                                          max = 2,
                                          post  = "x"
                              )
                          )
                      ),
                      hr(),
                      checkboxInput("additional", 
                                    label = "Show advanced settings",
                                    value = FALSE),
                      conditionalPanel(
                          condition = "input.additional == true",
                          actionButton("reset", label = "Reset to default"),
                          h5(withMathJax("$$L_{max} = aL_{\\infty}$$"), align="center"),
                          sliderInput("slider_linf_lmax",
                                      label = "a",
                                      value = 95,
                                      min = 50,
                                      max = 150,
                                      post  = "%"
                          ),
                          h5(withMathJax("$$log(L_{mean}) = b + c \\cdot log(L_{\\infty})$$"), 
                             align="center"),
                          sliderInput("slider_linf_lmean_b",
                                      label = "b",
                                      value = 0.37,
                                      min = 0.01,
                                      max = 1.5
                          ),
                          sliderInput("slider_linf_lmean_c",
                                      label = "c",
                                      value = 0.78,
                                      min = 0.01,
                                      max = 2.00
                          ),
                          h5("Coefficient of Variation in Mean size"),
                          sliderInput("slider_cv",
                                      label = "CV",
                                      value = 0.34,
                                      min = 0.2,
                                      max = 0.8
                          ),
                      )
                      
    ),
    
    conditionalPanel(
        condition = "input.metric == 'Maximum (estimated from maxima)'",
        navset_card_underline(
            
            nav_panel("Extreme value theory", 
                      uiOutput("evt_explanation")),
            
            nav_panel("Estimation of max length", 
                      layout_columns(
                          card(
                              card_header(
                                  tags$span(
                                      "Input length maxima",
                                      tags$span(
                                          icon(
                                              name = "question-circle",
                                          ) 
                                      ) |>
                                          add_prompt(message = "Each row corresponds to an  indvidual fish", 
                                                     position = "bottom")
                                  )
                              ),
                              card_body(
                                  p("Add extreme length measures here:"),
                                  rHandsontableOutput("hot")
                              )),
                          card(plotOutput("evt_plot")), 
                          col_widths = c(4,8)
                      ),
                      
                      div(textOutput("max_length"), 
                          style = "text-align: center; 
                    font-weight: bold;")),
            
            nav_panel("Unfished body size distribution", 
                      plotOutput("plot"), 
                      
                      
            ),
            
        )
    ),
    conditionalPanel(
        condition = "input.metric != 'Maximum (estimated from maxima)'",
        navset_card_underline(
            
            # Panel with summary ----
            nav_panel("Unfished body size distribution", 
                      plotOutput("plot2"),
            
                      conditionalPanel(
                          condition = "input.spp_selected != 'None'",
                          plotOutput("r_relchange_plot"),
                          plotOutput("r_relchange_plot_pc")

                          
                      )
            ),
            nav_panel("Understanding",
                      uiOutput("deepdive_p1"),
                      textOutput("deepdive_p2")),
            
        )
    )
)



# Define server logic ----
server <- function(input, output, session) {
    
    
    
    
    # General input data ----
    snapper_data <- 
        tibble(length_cm = c(91.3,
                             NA,
                             102, 112,
                             NA,
                             NA,
                             107, 107, 99.2,
                             95, 82.2,
                             NA),
               weight_kg = c(NA, 11.8, NA, NA, 18.4, 16.5, rep(NA, 5), 17.2))
    
    rls_bin_breaks <- 
        c(2.5, 5.0, 7.5,  10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 
          50.0, 62.5, 75.0, 87.5, 100.0, 112.5, 125.0, 137.5, 150.0, 
          162.5, 175.0, 187.5, 200.0, 250.0, 300.0, 350.0, 400.0)
    
    rls_bin_table <-
        tibble(size_class = c(0, rls_bin_breaks, 500)) %>% 
        mutate(
            size_indx  = 0:(length(size_class)-1),
            size_min = (size_class + lag(size_class))/2,
            size_max = lead(size_min)
        ) %>% 
        filter(size_class %in% c(rls_bin_breaks, 500))
    
    available_species <- 
        c(
            "Acanthopagrus australis",
            "Arothron hispidus",
            "Arothron nigropunctatus",
            "Arothron stellatus",
            "Chirodactylus spectabilis",
            "Choerodon rubescens",
            "Choerodon schoenleinii",
            "Pagrus auratus",
            "Kyphosus bigibbus",
            "Lethrinus miniatus",
            "Lethrinus nebulosus",
            "Lutjanus carponotatus",
            "Lutjanus sebae",
            "Mendosoma lineatum",
            "Nemadactylus douglasii",
            "Notolabrus tetricus",
            "Othos dentex",
            "Plectropomus leopardus",
            "Plectropomus maculatus",
            "Pseudocaranx georgianus",
            "Rhabdosargus sarba",
            "Scarus altipinnis",
            "Scarus ghobban",
            "Scarus prasiognathos",
            "Scarus rivulatus",
            "Scarus schlegeli"
        )
    fb_spp <- 
        fb_tbl("species") %>% 
        mutate(species_name = paste(Genus, Species)) %>% 
        select(SpecCode, species_name)  
    
    fb_meansize <- 
        fb_tbl("poplf") %>%
        rename(SpecCode = Speccode) %>% 
        left_join(fb_spp, by = join_by(SpecCode)) %>% 
        select(species_name, 
               lmean = MeanLength) %>% 
        filter(species_name %in% available_species) %>% 
        drop_na() %>%
        summarise(lmean = median(lmean, na.rm = TRUE), 
                  .by = species_name) 
    
    fb_linf <- 
        fb_tbl("popgrowth") %>% 
        left_join(fb_spp, by = join_by(SpecCode))%>% 
        select(species_name, 
               linf = Loo, 
               linf_TL = TLinfinity, 
               sex = Sex) %>% 
        filter(species_name %in% available_species) %>% 
        summarise(linf = median(linf, na.rm = TRUE), 
                  .by = species_name)
    
    
    # Functions ----
    kg_2_cm <- function(w, a = 0.01, b = 3) {((w*1000)/a)^(1/b)}
    ecdf_func <- function(data) { 
        Length <- length(data) 
        sorted <- sort(data) 
        
        ecdf <- rep(0, Length) 
        for (i in 1:Length) { 
            ecdf[i] <- sum(sorted <= data[i]) / Length 
        } 
        return(ecdf) 
    } 
    
    g <- function(x, shape) {gamma(1-(x*shape))}
    var <- function(scale, shape){
        if(shape == 0) v <- (scale^2)*((pi^2)/6)
        if(shape != 0 & shape < 0.5) v <- (scale^2)*(g(2, shape)-(g(1, shape))^2)/(shape^2)
        if(shape >= 0.5) v <- Inf
        return(v)
    }
    plot_exception <-function(
        ...,
        sep=" ",
        type=c("message","warning","cat","print"),
        color="auto",
        console=FALSE ,
        size = 6){      
        type=match.arg(type)
        txt = paste(...,collapse=sep)
        if(console){
            if(type == "message") message(txt)
            if(type == "warning") warning(txt)
            if(type == "cat") cat(txt)
            if(type == "print") print(txt)
        }
        if(color =="auto") color <- if(type == "cat") "black" else "red"
        if(txt == "warning") txt <- paste("warning:",txt)
        return(ggplot2::ggplot() +
                   ggplot2::geom_text(ggplot2::aes(x=0,y=0,label=txt),color=color,size=size) + 
                   ggplot2::theme_void())
        invisible(NULL)
    }
    
    # Reactive values ----
    reactives <- reactiveValues(evt_data = snapper_data)
    
    
    # Reactive ----
    est_max <- 
        reactive({
            evt_table() %>%
                pull(l) %>%
                as.numeric() %>%
                round(2)
        })
    
    # Update the obs_data when then dropdown is changed
    obs_data <- reactive({
        
        if(input$spp_selected %in% c(available_species)){
            
            
            input_data <- 
                paste0("example_data/", input$spp_selected, ".csv") %>% 
                read_csv(show_col_types = FALSE)  %>% 
                separate(col = gridyr, into = c("grid", "year"), sep = "__")
            
            if(input$population_scale == "spp"){
                input_data_rename <- 
                    input_data %>% 
                    rename(population = species_name) 
                
            } else {
                input_data_rename <- 
                    input_data %>% 
                    rename(population = grid) 
            }
            
            input_data_rename %>% 
                summarise(n = sum(mean_abundance), 
                          .by = c(population, size_class)) %>% 
                filter(size_class %in% rls_bin_breaks) %>% 
                add_count(population, wt = n, name = "total_n") %>% 
                mutate(p_obs = n/total_n) %>% 
                add_count(population, name = "n_classes") %>%
                filter(n_classes >= 4) %>% 
                left_join(rls_bin_table, by = join_by(size_class)) %>% 
                return()
            
        } else {
            req(input$datafile)
            input$datafile$datapath %>% 
                read_csv(show_col_types = FALSE) %>% 
                add_count(population, wt = n, name = "total_n") %>% 
                mutate(p_obs = n/total_n) %>% 
                return()
        }
        
    })
    
    meansize <- reactive({
        
        req(input$metric)
        
        if(input$metric == 'Asymptotic (theoretical)'){
            linf <- input$slider_linf
        }
        if(input$metric == 'Maximum (observed)'){
            linf <- input$slider_maxsize/(input$slider_linf_lmax/100)
        } 
        if(input$metric == 'Maximum (estimated from maxima)'){
            linf <- est_max()/(input$slider_linf_lmax/100)
            
        }
        if(input$metric == 'Mean (observed)'){
            meansize <- input$slider_meansize
        } else {
            linf_b <- input$slider_linf_lmean_b
            linf_c <- input$slider_linf_lmean_c
            meansize <- exp(linf_b + (linf_c*log(linf)))
        }
        meansize
    })
    
    
    plot_data <- 
        reactive({
            
            c <- input$slider_cv
            binwidth <- 0.1
            meansize <- meansize()
            
            # if data is input - use the size classes given
            if(input$spp_selected %in% c(available_species)){
                
                dat <- obs_data()
                
                main <-
                    rls_bin_table %>% 
                    select(size_class, size_min, size_max)
                
            } else {
                main <-
                    tibble(size_class = seq(1, input$xaxis[2], by = binwidth)) %>%
                    mutate(size_min = size_class - (binwidth/2),
                           size_max = size_class + (binwidth/2))
            }
            
            # main <- 
            #   tibble(size_class = seq(1, input$xaxis[2], by = binwidth)) %>% 
            #   mutate(size_min = size_class - (binwidth/2), 
            #          size_max = size_class + (binwidth/2)) 
            
            main %>%
                mutate(mu = meansize,
                       sd = mu*c,
                       sdlog = sqrt(log((c^2)+1)),
                       meanlog = log(meansize) - ((sdlog^2)/2)) %>%
                mutate(p_norm =
                           pnorm(size_max, mean = mu, sd = sd) -
                           pnorm(size_min, mean = mu, sd = sd),
                       p_lnorm =
                           plnorm(size_max, meanlog = meanlog, sdlog = sdlog) -
                           plnorm(size_min, meanlog = meanlog, sdlog = sdlog)) %>%
                pivot_longer(cols = contains("p_")) %>%
                filter(name %in% input$dist) %>% 
                filter(size_class >= input$xaxis[1], 
                       size_class <= input$xaxis[2]) %>% 
                return()
            
        })
    
    # Label to be used in the geom_vline()
    data_source_label <- 
        reactive({
            
            case_when(input$metric == 'Mean (observed)' ~ "Mean length = ", 
                      input$metric == 'Asymptotic (theoretical)' ~ "Asymptotic length = ", 
                      input$metric == 'Maximum (observed)' ~ "Maximum length = ",
                      input$metric == 'Maximum (estimated from maxima)' ~ "Maximum length = ")
        })
    
    # The value of the input metric slider
    data_source_number <- reactive({
        
        metric <- input$metric
        
        if(metric == 'Maximum (estimated from maxima)'){
            
            est_max <- 
                evt_table() %>% 
                pull(l) %>% 
                as.numeric() %>% 
                round(2)
            
            case_when(metric == 'Mean (observed)' ~ input$slider_meansize, 
                      metric == 'Asymptotic (theoretical)' ~ input$slider_linf, 
                      metric == 'Maximum (observed)' ~ input$slider_maxsize,
                      metric == 'Maximum (estimated from maxima)' ~ est_max)
        } else {
            case_when(metric == 'Mean (observed)' ~ input$slider_meansize, 
                      metric == 'Asymptotic (theoretical)' ~ input$slider_linf, 
                      metric == 'Maximum (observed)' ~ input$slider_maxsize)
        }
        
    })
    
    
    edited_evt_table <- reactive({
        req(input$hot)
        hot_to_r(input$hot) %>% 
            as_tibble() %>% 
            mutate(length_est = case_when(!is.na(weight_kg) ~ kg_2_cm(weight_kg),
                                          TRUE ~ length_cm))
    })
    
    
    
    # EVT calculations based on evt_table
    calc_table <- 
        reactive({
            
            l_m <- 
                edited_evt_table() %>% 
                pull(length_est) %>% 
                na.omit()
            
            n_samples <- length(l_m)
            
            if(n_samples >= 10){
                gev_fit <- evd::fgev(x = l_m, 
                                     std.err = TRUE)
                gev_shape = gev_fit$estimate[3]
            } else if(n_samples > 2){
                gev_fit <- evd::fgev(x = l_m, 
                                     shape = 0, # force gumbel
                                     std.err = FALSE)
                gev_shape = 0
            } else {
                gev_fit <- NULL
            }
            
            if(!is.null(gev_fit)){
                gev_location = gev_fit$estimate[1]
                gev_scale = gev_fit$estimate[2]
                
                tibble(l = seq(0, max(l_m)*1.5), by = 0.001) %>% 
                    rowwise() %>% 
                    mutate(gev = pgev(l, gev_location, gev_scale, gev_shape)) %>% 
                    rowwise() %>% 
                    mutate(var_gev = var(gev_scale, gev_shape),
                           sd_gev = sqrt(var_gev), 
                           se_gev = sd_gev/sqrt(length(l_m))) %>% 
                    mutate(lwr = l-(1.96*se_gev), 
                           upr = l+(1.96*se_gev)) 
            } else {
                NULL
            }
            
            
            
        })
    
    raw_table <- reactive({
        
        l_m <- 
            edited_evt_table() %>% 
            pull(length_est) %>% 
            na.omit()
        
        if(length(l_m) >= 3){
            raw_data <- 
                tibble(l_m = l_m) %>% 
                mutate(rank = rank(l_m)) %>% 
                arrange(rank) %>% 
                mutate(m = rank/max(rank+1)) %>% 
                mutate(ecdf = ecdf_func(l_m))
        } else {
            NULL
        }
        
        
        
    })
    
    evt_table <- reactive({
        
        raw_data <- raw_table()
        calculation_table <- calc_table()
        
        percentile <- max(raw_data$m)
        # percentile <- 0.99
        
        if(!is.null(calculation_table)){
            calculation_table %>% 
                mutate(dist = abs(gev - percentile)) %>% 
                arrange(dist) %>% 
                head(1)
        } else {
            NULL
        }
        
    })
    
    # Unfished distribution plot
    plt1 <- reactive({
        
        req(input$spp_selected)
        req(input$xaxis)
        req(input$dist)
        
        pdat <-
            plot_data() %>%
            mutate(name = case_when(name == "p_norm" ~ "Normal",
                                    name == "p_lnorm" ~ "Log-normal",
                                    TRUE ~ name))
        
        p <-
            pdat %>%
            ggplot() +
            aes(x = size_class, y = value, color = name) +
            {if(data_source_number() >= input$xaxis[1] &
                data_source_number() <= input$xaxis[2]) {
                geom_textvline(aes(xintercept = xint,
                                   label = lab),
                               size = 5,
                               data = tibble(xint = data_source_number(),
                                             lab = paste(data_source_label(),
                                                         data_source_number(),
                                                         "cm")))
            }} +
            {if("p_norm" %in% input$dist) {geom_line(lty = 2, linewidth = 1.5,
                                                     col = rgb(181, 144, 19, maxColorValue=255),
                                                     data = filter(pdat, name == "Normal"))}} +
            {if("p_lnorm" %in% input$dist) {geom_line(lty = 2, linewidth = 1.5,
                                                      col = rgb(29, 84, 128, maxColorValue=255),
                                                      data = filter(pdat, name == "Log-normal"))}} +
            scale_x_continuous(limits = input$xaxis) +
            theme_classic(20) +
            labs(y = "Relative abundance",
                 x = "Body length (cm)") +
            guides(col=guide_legend(title = "Population"))
        
        if((!is.null(input$datafile)) | input$spp_selected != "None"){
            
            req(input$radio_scaled)
            req(input$population_scale)
            req(input$dist)
            
            obs_data_v2 <-
                obs_data() %>%
                select(size_class,
                       p_obs,
                       population) %>%
                left_join(plot_data() %>%
                              mutate(name = case_when(name == "p_norm" ~ "Normal",
                                                      name == "p_lnorm" ~ "Log-normal",
                                                      TRUE ~ name)),
                          by = join_by(size_class)) %>%
                {if("p_norm" %in% input$dist) filter(., name == "Normal") else filter(., name == "Log-normal")} %>%
                {if(input$population_scale == "pop") group_by(., population) else .} %>%
                mutate(
                    scaling = p_obs/value, # observed scaling value
                    max_diff = max(scaling),
                    value_maxdiff = max(case_when(scaling==max_diff ~ value,
                                                  TRUE ~ 0)),
                    value_atmax = max(case_when(p_obs==max(p_obs) ~ value,
                                                TRUE ~ 0)),
                    max_pobs = max(p_obs)) %>%
                mutate(
                    p_obs_scaled =  case_when(
                        input$radio_scaled == "rel_abun" ~ p_obs,
                        input$radio_scaled == "scale_max" ~ p_obs/(max(p_obs)/value_atmax),
                        input$radio_scaled == "scale_subsample" ~ p_obs/max_diff,
                        input$radio_scaled == "choose_scaling" ~ p_obs*input$y_scaling
                    )) %>%
                {if(input$population_scale == "pop") ungroup(.) else .}
            
            p +
                geom_line(aes(y = p_obs_scaled,
                              col = population),
                          data = obs_data_v2) +
                {if(input$radio_scaled == "scale_max") geom_point(
                    aes(y = p_obs_scaled),
                    data = obs_data_v2 %>% filter(p_obs==max_pobs),
                    col = "black",
                    size = 1.5)} +
                theme(legend.position = "none") +
                {if(input$radio_scaled != "rel_abun") labs(y = "Scaled abundance")}
        } else {
            p
        }
        
    })
    
    relchange_data <- reactive({
        
        req(input$spp_selected)
        req(input$xaxis)
        req(input$dist)
        
        pdat <-
            plot_data() %>%
            mutate(name = case_when(name == "p_norm" ~ "Normal",
                                    name == "p_lnorm" ~ "Log-normal",
                                    TRUE ~ name))
        
        if((!is.null(input$datafile)) | input$spp_selected != "None"){
            
            req(input$radio_scaled)
            req(input$population_scale)
            req(input$dist)
            
            obs_data_v2 <-
                plot_data() %>%
                mutate(name = case_when(name == "p_norm" ~ "Normal",
                                        name == "p_lnorm" ~ "Log-normal",
                                        TRUE ~ name)) %>% 
                left_join(obs_data() %>%
                              select(size_class,
                                     p_obs,
                                     population),
                          by = join_by(size_class)) %>%
                mutate(p_obs = replace_na(p_obs, 0)) %>% 
                {if(input$population_scale == "pop") group_by(., population) else .} %>%
                mutate(
                    scaling = p_obs/value, # observed scaling value
                    max_diff = max(scaling),
                    value_maxdiff = max(case_when(scaling==max_diff ~ value,
                                                  TRUE ~ 0)),
                    value_atmax = max(case_when(p_obs==max(p_obs) ~ value,
                                                TRUE ~ 0)),
                    max_pobs = max(p_obs)) %>%
                mutate(
                    p_obs_scaled =  case_when(
                        input$radio_scaled == "rel_abun" ~ p_obs,
                        input$radio_scaled == "scale_max" ~ p_obs/(max(p_obs)/value_atmax),
                        input$radio_scaled == "scale_subsample" ~ p_obs/max_diff,
                        input$radio_scaled == "choose_scaling" ~ p_obs*input$y_scaling
                    )) %>%
                {if(input$population_scale == "pop") ungroup(.) else .}
            
        } else {
            NULL
        }
        
    })
    
    # Relative change by length ----
    relchange_plot <- reactive({
        
        plot_data <- relchange_data()
        plot_data %>% 
            mutate(rel_change = p_obs - value) %>% 
            ggplot(aes(xmin = size_min, 
                       xmax = size_max,
                       ymax = rel_change, 
                       fill = name))  +
            geom_hline(yintercept = 0, lty = 2) +
            geom_rect(ymin = 0, col = "black") +
            facet_wrap(~name, ncol = 1) +
            theme_classic(20) +
            theme(strip.background = element_blank(),
                  strip.text.x = element_blank(), 
                  legend.title=element_blank()) +
            labs(x = "Body length (cm)",
                 y = "Relative change") +
            scale_fill_manual(values = c("Normal" = rgb(181, 144, 19, maxColorValue=255), 
                                         "Log-normal" = rgb(29, 84, 128, maxColorValue=255))) +
            scale_y_continuous(limits = c(-1, 1)*max(plot_data$p_obs)) +
            theme(legend.position = "bottom")
        
        
    })
    
    # Relative change by percentile ----
    relchange_plot_pc <- reactive({

        plot_data <- relchange_data()
        c <- input$slider_cv
        mu = meansize()
        sd = mu*c
        sdlog = sqrt(log((c^2)+1))
        meanlog = log(mu) - ((sdlog^2)/2)

        plot_data %>%
            mutate(percentile = case_when(name == "Normal" ~ pnorm(size_class, mean = mu, sd = sd),
                                    name == "Log-normal" ~ plnorm(size_class, meanlog = meanlog, sdlog = sdlog))) %>%
            mutate(rel_change = p_obs - value) %>%
            ggplot(aes(x = percentile,
                       y = rel_change,
                       col = name)) +
            geom_hline(yintercept = 0, lty = 2) +
            # geom_rect(ymin = 0, col = "black") +
            geom_path() +
            facet_wrap(~name, ncol = 1) +
            theme_classic(20) +
            theme(strip.background = element_blank(),
                  strip.text.x = element_blank(),
                  legend.title=element_blank()) +
            labs(x = "Percentile",
                 y = "Relative change") +
            scale_color_manual(values = c("Normal" = rgb(181, 144, 19, maxColorValue=255),
                                         "Log-normal" = rgb(29, 84, 128, maxColorValue=255))) +
            scale_y_continuous(limits = c(-1, 1)*max(plot_data$p_obs)) +
            scale_x_continuous(label = label_percent())+
            theme(legend.position = "none")

        
    })
    
    
    # Observe ----
    # Create the length_est col when the evt_table is modified
    observe({
        if(!is.null(input$hot)){
            temp_data <- as_tibble(hot_to_r(input$hot))
            temp_data <- temp_data %>% 
                mutate(length_est = case_when(!is.na(weight_kg) ~ kg_2_cm(weight_kg),
                                              TRUE ~ length_cm))
            # only update the evt_data if same as temp data
            # this avoid the infinite loop is when the table is quickly modified
            if (!identical(temp_data, reactives$evt_data)) {
                reactives$evt_data <- temp_data
                output$hot <- renderRHandsontable({
                    temp_data %>% 
                        rhandsontable() %>%
                        # hot_col("date", readOnly = FALSE) %>%
                        hot_col("length_est", readOnly = TRUE)
                })
            }
        }
    })
    
    
    
    # Observe Events ----
    observeEvent(input$reset, {
        updateSliderInput(session, "slider_linf_lmax", value = 95)
        updateSliderInput(session, "slider_linf_lmean_b", value = 0.37)
        updateSliderInput(session, "slider_linf_lmean_c", value = 0.78)
        updateSliderInput(session, "slider_cv", value = 0.34)
    })
    
    
    
    
    # Table outputs ----
    output$hot <- 
        renderRHandsontable({
            rhandsontable(reactives$evt_data)
        })
    
    
    # Plot outputs ----
    output$plot <- renderPlot({
        plt1()
    })
    output$plot2 <- renderPlot({
        plt1()
    })
    
    output$r_relchange_plot <- renderPlot({
        relchange_plot()
    })
    
    output$r_relchange_plot_pc <- renderPlot({
        relchange_plot_pc()
    })
    
    output$evt_plot <- renderPlot({
        
        raw_data <- raw_table()
        
        if(!is.null(calc_table())){
            calculation_table <- calc_table()
            out <- evt_table()
            
            calculation_table %>%
                ggplot(aes(x = l, y = gev)) +
                geom_segment(y = out$gev, yend = out$gev, x = -Inf, xend = out$upr, col = "grey") +
                geom_segment(y = -Inf, yend = out$gev, x = out$lwr, xend = out$lwr, col = "grey") +
                geom_segment(y = -Inf, yend = out$gev, x = out$l, xend = out$l, col = "grey") +
                geom_segment(y = -Inf, yend = out$gev, x = out$upr, xend = out$upr, col = "grey") +
                geom_line() +
                geom_line(aes(x = lwr), col = "red", lty = 2) +
                geom_line(aes(x = upr), col = "red", lty = 2) +
                geom_point(aes(x = l_m, y = m), 
                           data = raw_data) +
                labs(x = "Length (cm)", y = "Probability density") +
                theme_classic(20)
        } else {
            plot_exception("Please provide at least 3 maxima values (length or weight).")
        }
        
    })
    
    
    evt_text <- reactive({
        list(
            p("The maximum body length of a population or species (\\(L_{max}\\)) is dependant on the rate of mortality. We can use historical and contemporary estimates of the largest individuals to estimate \\(L_{max}\\), however we are unlikely to know if these estimates accurately represent the 'largest' possible body length. One approach is to use statistical methods to estimate the most likely largest size, based on the knowledge of large individuals that have been measured."),
            p("Extreme value theory (EVT) is based on the premise that if you sample from a given distribution (i.e., body length distribution), the maxima of those samples will follow a defined distribuiton. EVT is used extensively in fields such as hydrology (predicting the next flood), and stock markets (predicting the next crash), but can be applied in this case to the probability of larger individuals than those observed."), 
            p("Much like the central-limit-theorem, which says the distribution of the means from a set of samples from any dsitribution will approxmate a normal distribution, EVT shows the maxima will also follow a predictable distribution. This distribution of maxima, can be approximated by the 'Gumbel' distribution, which has two parameters. An extension of the Gumbel distribution is the generalised extreme value (GEV) distribution, which has three parameters and simplies to the Gumbel dsitribution when one of its parameters (scale) is zero."
            ),
            p("In the following tab use can apply the EVT to estimate \\(L_{max}\\) of a species given a set of length maxima. For example, these maxima could be a set of independent fishing competition largest length. "
            )
        )
    })
    
    
    # Text outputs ----
    
    output$evt_explanation <- renderUI({
        withMathJax(
            # helpText(
            p(
                evt_text()
            ))
        
    })
    
    output$max_length <- renderText({
        
        out <- evt_table()
        
        paste0("Estimated maximum length = ", round(out$l, 2), "cm (", round(out$lwr, 2), "cm-", round(out$upr, 2), "cm, 95% confidence interval)")
        
    })
    
    
    deepdive_p1 <- reactive({
        
        # Paragraph 1
        if(input$metric == 'Asymptotic (theoretical)'){
            p1 <- 
                list(p("The asymptotic body length of a species is a theorectical concept, based on the idea that as growth slows with age, fish tend to reach an asymptote in body length. This \\(L_{\\infty}\\) value is often modelled using the von Bertalanffy growth equation. Its value represents the mean body length of a sample of individuals of 'infinite age'."), 
                     p("\\(L_{a} = L_{\\infty}(1-e^{-k(a - t_0)})\\)"), 
                     p("where, \\(L_{a}\\) is length at age (\\(a\\)), \\(k\\) is the growth coefficient and \\(t_0\\) is the theoretical age when body length is zero."), 
                     p("\\(L_{\\infty}\\) is usually estimated from growth data (e.g. otolith growth rings), and fitting the von Bertalanffy growth equation. Asymptotic length (\\(L_{\\infty}\\)) differs from maximum body length (\\(L_{max}\\)), which is an observed real-world value. \\(L_{max}\\) may be greater or less than \\(L_{\\infty}\\) depending on the mortality rate of the population, as high mortality may result in few individuals reaching the asymtote.")
                )
            
        }
        if(input$metric == 'Maximum (observed)'){
            p1 <- 
                list(p("The maximum body length (\\(L_{max}\\)) of a species or population is an observed parameter, i.e., what is the largest individual ever obsered for this population or species? The largest observed individual is dependent on mortality rates, both natural and fishing morality. In situations with relatively low mortality, maximum length is expected to be greater than asymptotic length (\\(L_{\\infty}\\)), since \\(L_{\\infty}\\) is the mean value for individuals of infinite age."), 
                     p("Since mortality influences \\(L_{max}\\) we may want to estimate it statistically, using either historical values or a set of known sample maxima. See the 'Maximum (Estimated)' metric for a method of how to estimate \\(L_{max}\\) from length maxima.")
                )
        } 
        if(input$metric == 'Mean (observed)'){
            p1 <- 
                "The mean body length of a species generally unknown, and is very difficult to estimate accurately. This may be possible given appropriate (unbiased) sampling techniques with sufficient sample sizes."
        }
        
        # Paragraph 2
        if(input$metric == 'Asymptotic (theoretical)'){
            metric_simple <- "\\(L_{\\infty}\\)"
        }
        if(input$metric == 'Maximum (observed)'){
            metric_simple <- "\\(L_{max}\\)"
        } 
        if(input$metric == 'Mean (observed)'){
            metric_simple <- "\\(L_{\\mu}\\)"
        }
        
        p1 
    })
    
    deepdive_p2 <- reactive({
        
        # Paragraph 2
        if(input$metric == 'Asymptotic (theoretical)'){
            metric_simple <- "\\(L_{\\infty}\\)"
        }
        if(input$metric == 'Maximum (observed)'){
            metric_simple <- "\\(L_{max}\\)"
        } 
        if(input$metric == 'Mean (observed)'){
            metric_simple <- "\\(L_{\\mu}\\)"
        }
        
        paste0("Move the slider in the sidebar (currently set at ", data_source_number(), "cm) to change the ", metric_simple, " value and see how this impacts the unfished size distribution plot.")
        
    })
    
    
    
    
    output$deepdive_p1 <- renderUI({
        withMathJax(
            # helpText(
            p(deepdive_p1()), 
            p(deepdive_p2())
            # )
        )
        
    })
    output$deepdive_p2 <- renderText({
        
        
        
    })
    output$species_dropdown <- 
        renderUI({
            selectInput(inputId = "spp_selected", 
                        label = "Or select species (from Reef Life Survey)", 
                        choices = c("None", available_species), 
                        selected = "None")
        })
    
    # Download data
    output$downloadData <- downloadHandler(
        filename = function() paste("data-", Sys.Date(), ".csv", sep=""),
        content = function(file) write.csv(
            read.csv("example_data/example_data.csv"), 
            file)
    )
    
    output$downloadTemplate <- downloadHandler(
        filename = function() paste("length_maxima", Sys.Date(), ".csv", sep=""),
        content = function(file) write.csv(
            read.csv("example_data/length_maxima_template.csv"), 
            file)
    )
}

# Run the app ----
shinyApp(ui = ui, server = server)
