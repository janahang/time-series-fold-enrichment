library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
  # Application title
  headerPanel("Fold Enrichment Analyses for a Time Series Experiment"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    helpText(p("Choose a condition you wish to analyze"),
             p("Up-regulated = 1"),
             p("Unchanged = 0"),
             p("Down-regulated = -1")),
    #sliderInput("obs", "Number of observations:", min = 0, max = 1000, value = 500)
    selectInput("condition", "Choose a condition", choices = c("1", "0", "-1"))
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("distPlot",width = "120%", height = "600px")
  )
))