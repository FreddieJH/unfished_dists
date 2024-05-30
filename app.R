library(shiny)
library(tidyverse)


# Define UI for application
ui <- fluidPage(
    
    # Application title
    titlePanel("Reconstruction of unfished body size distributions"),
    
    # Sidebar layout with input and output definitions
    sidebarLayout(
        
        # Sidebar panel for inputs
        sidebarPanel(
            
            h3("Available body size metric"),
            selectInput("datatype", "Data type:", choices = c("Asymptotic (theoretical)", "Maximum (observed)", "Mean (observed)")),
            conditionalPanel(
                condition = "input.datatype == 'Asymptotic (theoretical)'",
                sliderInput("slider_linf", "Asymptotic length (cm)",
                            value = 60,
                            min = 1,
                            max = 200, 
                            post  = "cm"
                )
            ),
            conditionalPanel(
                condition = "input.datatype == 'Maximum (observed)'",
                sliderInput("slider_maxsize", "Maximum observed size (cm)",
                            value = 50,
                            min = 1,
                            max = 200, 
                            post  = "cm"
                )
            ),
            conditionalPanel(
                condition = "input.datatype == 'Mean (observed)'",
                sliderInput("slider_meansize", "Mean observed size (cm)",
                            value = 25,
                            min = 1,
                            max = 200, 
                            post  = "cm"
                )
            ),
            hr(),
            h3("Assumed distribution"),
            checkboxGroupInput("dist", label = "Distribution of unfished", 
                               choices = list("Truncated normal" = "p_norm",
                                              "Log-normal" = "p_lnorm"),
                               selected = "p_norm"),
            hr(),
            h3("Observed data"),
            checkboxInput("inc_obsdata", label = "Include observations", value = TRUE),
            conditionalPanel(
                condition = "input.inc_obsdata == true",
                downloadButton("downloadData", "Download example data"), 
                fileInput("datafile", "Choose CSV File", accept = ".csv"),
            ),
            
           
            
            hr(),
            checkboxInput("additional", label = "Show advanced settings", value = FALSE),
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
                h5(withMathJax("$$log(L_{mean}) = b + c \\cdot log(L_{\\infty})$$"), align="center"),
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
        
        # Main panel for displaying outputs
        mainPanel(
            # Plot output
            plotOutput("plot")
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    
    observeEvent(input$reset, {
        updateSliderInput(session, "slider_linf_lmax", value = 95)
        updateSliderInput(session, "slider_linf_lmean_b", value = 0.37)
        updateSliderInput(session, "slider_linf_lmean_c", value = 0.78)
        updateSliderInput(session, "slider_cv", value = 0.34)
    })
    
    # Reactive expression for the data
    data <- reactive({
        
        meansize <-
            case_when(input$datatype == 'Mean (observed)' ~ input$slider_meansize, 
                      input$datatype == 'Asymptotic (theoretical)' ~ exp(input$slider_linf_lmean_b + (input$slider_linf_lmean_c*log(input$slider_linf))), 
                      input$datatype == 'Maximum (observed)' ~ exp(input$slider_linf_lmean_b + (input$slider_linf_lmean_c*log(input$slider_maxsize/(input$slider_linf_lmax/100)))))
        
        # req(input$datafile)
        # read.csv(input$datafile$datapath)
        
        c <- input$slider_cv
        binwidth <- 0.1
        
        
            if(!is.null(input$datafile)){
                
                dat <- obs_data()
                main <- 
                    tibble(size_class = dat$size_class %>% unique() %>% sort()) %>% 
                    mutate(size_min = dat$size_min %>% unique() %>% sort(), 
                           size_max = dat$size_max %>% unique() %>% sort()) 
            } else {
                main <- tibble(size_class = seq(1, 100, by = binwidth)) %>% 
                    mutate(size_min = size_class - (binwidth/2), 
                           size_max = size_class + (binwidth/2)) 
            }
        
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
            filter(name %in% input$dist)
    })
    
    obs_data <- reactive({
        
        req(input$datafile)
        read.csv(input$datafile$datapath) %>% 
            add_count(population, wt = n, name = "total_n") %>% 
            mutate(p_obs = n/total_n)
        
    })
    
    
    # Render the plot
    output$plot <- renderPlot({
        p <- 
            data() %>% 
            mutate(name = case_when(name == "p_norm" ~ "Normal", 
                                    name == "p_lnorm" ~ "Log-normal",
                                    TRUE ~ name)) %>% 
            ggplot() +
            aes(x = size_class, y = value, col = name) +
            geom_line(lty = 2, linewidth = 1.5) +
            theme_classic(15) +
            labs(y = "Relative abundance", 
                 x = "Body length (cm)") +
            guides(col=guide_legend(title = "Population"))
        
        if(!is.null(input$datafile)){
            p + geom_line(aes(y = p_obs, col = population), data = obs_data())
        } else {
            p
        }
       
    })
    
    # Download data
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(read.csv("example_data.csv"), 
                      file)
        }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
