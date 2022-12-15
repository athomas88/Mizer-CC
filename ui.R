# ui.R
##Testing
shinyUI(fluidPage(
  navbarPage(("Myanmar Fishery Simulation"),
             tabPanel("Model Inputs and Outputs",
                      sidebarLayout(
                        sidebarPanel(
                          h2("Model Inputs"),
                          selectInput("siteAnalysis",label="Select site",choices=c("Rakhine", "Delta", "Tanintharyi"),selected="Rakhine"),
                          sliderInput(inputId = "fisheryAge", label = "Age of fishery for burn-in period [Years]", 1, 100, 50, step = 1),
                          sliderInput(inputId = "managementStart", label = "When should new management start?", 2019, 2039, 2019, step = 1, sep=""),
                          sliderInput(inputId = "tProjection", label = "Years projected into future", 1, 100, 50, step = 1),
                          sliderInput(inputId = "effortCreep", label = "Increase in fishing over time [% per year]", 0, 5, 1, step = 0.1),
                          sliderInput(inputId = "compliance", label = "What percentage of the fishery will obey management rules?", 0, 100, 80, step = 5),
                          checkboxGroupInput("interventionGroup", label = "Select potential management options", 
                                             choices = list("No Change" =  1, 
                                                            "Limit F relative to natural mortality" = 2,
                                                            "Single F across species" = 3,
                                                            "Single minimum size limit across species" = 4,
                                                            "Minimum size limit (for each species)" = 5,
                                                            "Size limits and F limits" = 6,
                                                            "No-take zone" = 7,
                                                            "Seasonal closure" = 8),
                                             selected = c(1)),
                          hr(),
                          uiOutput("customControlF"),
                          uiOutput("customControlFsize"),
                          hr(),
                          uiOutput("customControlNTZ"),
                          hr(),
                          uiOutput("customControlSingleF"),
                          uiOutput("customControlSingleSize"),
                          hr(),
                          uiOutput("customControlNTseason"),
                          uiOutput("customControlNTF"),
                          hr(),
                          actionButton(inputId = "runSim", label = "Run")
                        ),
                        mainPanel(
                          h2("Model Outputs"),
                          selectInput("plotSelect", label = "Select plot types", 
                                      choices = list("Total catch and population of all species" = 1, "Species-specific catch and population" = 2,"Projected size spectrum" = 3), 
                                      selected = 1),
                          conditionalPanel(
                            condition = "input.plotSelect == 1",
                            plotOutput("PlotProjectionsAggregated",height="1200px",width="1000px")
                          ),
                          conditionalPanel(
                            condition = "input.plotSelect == 2",
                            uiOutput("speciesSelect"),
                            plotOutput("PlotProjectionsSpecies",height="1200px",width="1000px")
                          ),
                          conditionalPanel(
                            condition = "input.plotSelect == 3",
                            plotOutput("PlotSizeSpectrum",height="400px",width="1000px")
                          ),
                        )
                      )
             ),
             tabPanel("Species-specific inputs",
                      rHandsontableOutput("hot"),
                      actionButton(inputId = "saveParams", label = "Save"),
                      actionButton(inputId = "revertParams", label = "Revert")
             )
  )
)
)
