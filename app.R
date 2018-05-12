
library(tidyverse)
library(depmixS4)
library(shiny)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Uploading Files"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput(
        "file1", "Select bed file",
        multiple = T,
        accept = c("text/tsv", ".bed")
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
    
      # Output: Data file ----
      tableOutput("contents")
    
    )
  )   
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  options(shiny.maxRequestSize = 500 * 1024 ^ 2)
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        
        numfiles = nrow(input$file1) 
        file_list = list()
        
        
        for (i in 1:numfiles)
        {
          
          JSON_csv = read.csv(input$file1[[i, 'datapath']], header = TRUE)
          lastrow = nrow(JSON_csv)
          shift = function(x, n){
            c(x[-(seq(n))], rep(NA, n))
          }
          JSON_csv$companyID1 = shift(JSON_csv$companyID1, 1)
          kata_csv1[[i]] = JSON_csv[-lastrow, ]
          
        }
        
        
        col_names <- c(
          "chrom", "start",
          "end", "name",
          "win_id", "strand",
          "count"
        )
        
        df <- read_tsv(input$file1$datapath, col_names)
                       
        #df <- read.csv(input$file1$datapath,
        #               header = input$header,
        #               sep = input$sep,
        #               quote = input$quote)
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    #return(df %>% head())
    return("done")
    
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)





# Horizontal line ----
#     tags$hr(),

# Input: Select separator ----
#      radioButtons("sep", "Separator",
#                   choices = c(Comma = ",",
#                               Semicolon = ";",
#                               Tab = "\t"),
#                   selected = ","),

# Input: Select quotes ----
#      radioButtons("quote", "Quote",
#                   choices = c(None = "",
#                               "Double Quote" = '"',
#                               "Single Quote" = "'"),
#                   selected = '"'),

# Horizontal line ----
#     tags$hr(),

# Input: Select number of rows to display ----
#      radioButtons("disp", "Display",
#                   choices = c(Head = "head",
#                              All = "all"),
#                   selected = "head")

#    ),
#  )
#)









