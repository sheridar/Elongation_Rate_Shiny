
library(depmixS4)
library(DT)
library(shiny)
library(tidyverse)
library(magrittr)

# Function to add "key" columns to list of dfs ----
add_key <- function(input) {
  
  tbl_list <- list()
  
  for (i in 1:length(input)) {
    x <- input[i] 
    y <- data.frame(x)
    
    new_names <- str_replace(colnames(y), str_c(names(x), "."), "")
    colnames(y) <- str_replace(new_names, "count", names(x))
    
    tbl_list <- c(tbl_list, list(y))
  }
  
  tbl_list
}

# Function to merge tables ---- 
tbl_merge <- function(input, merge_by = "name") {
  
  tbl_list <- add_key(input) 
  
  # Merged tables 
  res <- purrr::reduce(tbl_list, function(x, y) {
    left_join(x, y, by = merge_by)
  }) %>% 
    na.omit()
  
  res
}

# Function to merge tables ---- 
DRB_merge <- function(input, win_min = 1, win_max = 200) {
  
  # Selected columns and filtered 
  res <- map(input, function(x) {
    x %>% 
      dplyr::select(chrom, start, end, name, win_id, count) %>%
      group_by(name) %>%
      filter(
        win_id >= win_min,
        win_id <= win_max
      ) %>%
      filter(sum(count) > 0) %>%
      ungroup()
  })
  
  # Merged tables 
  res <- tbl_merge(res, merge_by = c("chrom", "start", "end", "name", "win_id")) %>%
    gather(key, count, -chrom, -start, -end, -name, -win_id)
  
  res
}

# Function to normalize signal ----
DRB_norm <- function(input_file, win_max) {
  
  win_cutoff <- win_max - 5
  
  res <- input_file %>% 
    # Added pseudo count 
    group_by(key, name) %>% 
    mutate(zero = ifelse(count == 0, T, F)) %>% 
    group_by(key, name, zero) %>% 
    mutate(min_count = min(count)) %>%
    group_by(key, name) %>% 
    mutate(count = ifelse(count == 0, max(min_count) / 2, count)) %>% 
    ungroup() %>% 
    dplyr::select(-zero, -min_count) %>% 
    
    # Normalized by -DRB signal 
    separate(key, sep = "_", into = c("treatment", "tm")) %>% 
    spread(tm, count) %>%
    gather(tm, count, -chrom, -start, -end, -name, -win_id, -treatment, -con) %>% 
    mutate(count = count / con) %>%
    dplyr::select(-con) %>% 
    
    # Normalized ratios by the average ratio for the last 5 bins 
    unite(key, treatment, tm, sep = "_") %>% 
    mutate(win_type = ifelse(win_id > win_cutoff, "con_wins", "data_wins")) %>% 
    group_by(key, name, win_type) %>% 
    mutate(ave_signal = mean(count)) %>% 
    ungroup() %>% 
    spread(win_type, ave_signal) %>% 
    mutate(con_wins = ifelse(is.na(con_wins), 0, con_wins)) %>% 
    group_by(name, key) %>% 
    mutate(con_wins = max(con_wins)) %>% 
    ungroup() %>% 
    mutate(count = count / con_wins) %>% 
    dplyr::select(-data_wins, -con_wins) %>% 
    
    # Digitized ratios to a range of 0 - 2.0 and step size of 0.5
    group_by(name, key) %>%
    mutate(max_count = max(count)) %>% 
    ungroup() %>% 
    mutate(
      count = (count / max_count) * 2,
      count = floor(count / 0.05) / 20
    ) %>%
    dplyr::select(-max_count)
  
  res
}

# Function to run depmix ----
find_edge <- function(input_file) {
  tot_win <- input_file %>% 
    group_by(win_id) %>% 
    group_size() %>% 
    length()
  
  tot_gene <- input_file %>% 
    group_by(name) %>% 
    group_size() %>% 
    length()
  
  table_sort <- input_file %>% arrange(name, win_id)
  name <- table_sort$name
  count <- table_sort$count
  res <- rep(NA, tot_gene)
  
  for (i in 0:(tot_gene - 1)) {
    input <- count[ (i * tot_win + 1) : (i * tot_win + tot_win) ]
    input <- data.frame(input) 
    
    mod <- depmix(response = input ~ 1, data = input, nstates = 2, trstart = runif(4))
    
    tryCatch(
      fm <- fit(mod, emc = em.control(rand = FALSE)),
      error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }
    )
    
    esttrans <- posterior(fm)$state
    
    if (esttrans %>% unique() %>% length() == 2) {
      edge <- rep(NA, tot_win)
      
      for (j in 1:length(esttrans)) {
        if (j > 4) {
          sum_state <- sum(esttrans[ (j - 4) : j ])  
          if (sum_state == 5) {
            edge[j] <- j
          }
        }
      }
    }
    
    edge %<>%
      na.omit() %>%
      tail(1) 
    
    if (length(edge) == 1) {
      edge <- edge * 500 + 500
      res[i + 1] <- edge
    }
  }
  
  res <- data.frame(unique(name), res)
  colnames(res) <- c('name', 'wave_edge')
  res
}

# Function to run find_edge() ----
run_find_edge <- function(input_file) {
  
  df <- input_file %>% 
    dplyr::select(-chrom, -start, -end) %>% 
    spread(key, count) 
  
  name_cols <- df %>% dplyr::select(name, win_id) 
  
  res <- df %>% 
    dplyr::select(name) %>% 
    unique()
  
  col_names <- "name"
  
  for (i in 3:ncol(df)) {
    input_data <- bind_cols(name_cols, df[, i])
    data_name <- names(input_data)[3]
    col_names <- c(col_names, data_name)
    input_data %<>% gather(key, count, -name, -win_id)
    
    waves <- find_edge(input_data)
    res %<>% left_join(waves, by = "name") 
    colnames(res) <- col_names
  }
  
  res
}

# Function to extract gene symbols ----
extract_gene_symbol <- function(input_file) {
  
  get_last_name <- function(gene_string) {
    res <- str_split(gene_string, "\\|")
    
    str_len <- length(res[[1]])
    
    res <- res[[1]][[str_len]]
    
    res
  }
  
  gene_names <- input_file %>%
    dplyr::select(name)
  
  other_data <- input_file %>%
    dplyr::select(-name)
  
  gene_matrix <- as.matrix(gene_names)
  
  new_names <- map(gene_matrix, get_last_name)
  
  new_names <- data.frame(name = as.matrix(new_names)) %>%
    mutate(name = as.character(name))
  
  res <- bind_cols(new_names, other_data)
}


# Define UI for data upload app ----
ui <- fluidPage(
  
  titlePanel("Elongation rate calculator"),
  
  column(12,
  
  fluidRow(
    column(4,
      fileInput(
        "file_1", "Timecourse data",
        accept = c("text/tsv", ".bed")
      )
    ),
    
    column(2, 
      numericInput(
        "time_1", "Time (min)", 
        10, min = 0, max = 999
      )
    )
  ),
  
  fluidRow(
    column(4,
      fileInput(
        "file_2", label = NULL,
        accept = c("text/tsv", ".bed")
      )
    ),
    
    column(2, 
      numericInput(
        "time_2", label = NULL, 
        20, min = 1, max = 999
      )
    )  
  ),

  fluidRow(
    column(4,
      fileInput(
        "control", "Control data",
        accept = c("text/tsv", ".bed")
      )
    ),
    
    column(2,
      numericInput(
        "win_min", "Window min",
        1, min = 1, max = 200
      )
    )
  ),
    
  fluidRow(
    column(4,
      fileInput(
        "gene_list", "Gene list",
        accept = c("text/tsv", ".txt")
      )
    ),
    
    column(2,
      numericInput(
        "win_max", "Window max",
        200, min = 1, max = 200
      )
    )
  ),
  
  fluidRow(
    column(1, actionButton("runAnalysis", "RUN")),
    column(1, downloadButton("download", "Export"))
  )
  )
  
  #column(6, fileInput("test_file", "test_file"))
  
  #column(6, 
  #  div(
  #    DT::dataTableOutput("contents")
      #style = "font-size: 75%; width: 75%; height:400px; text-overflow: ellipsis"
  #  )
  #)
  
  #mainPanel(
  #  dataTableOutput("contents")
  #)
)


# Define server logic to read selected file ----
server <- function(input, output) {
  
  options(shiny.maxRequestSize = 500 * 1024 ^ 2)
  
  tryCatch(
    {
      table_out <- eventReactive(input$runAnalysis, ignoreInit = T, {
          
        req(input$file_1)
        req(input$file_2)
        req(input$time_1)
        req(input$time_2)
        
        
        # Merged data tables 
        col_names <- c(
          "chrom", "start",
          "end", "name",
          "win_id", "strand",
          "count"
        )
    
        file_list <- list(input$file_1$datapath, input$file_2$datapath, input$control$datapath)
        df_list <- map(file_list, function(x) read_tsv(x, col_names))
        
        if (!is.null(input$gene_list)) {
          gene_list <- read_tsv(input$gene_list$datapath, "name")
          
          df_list <- map(df_list, function(x){
            x %>% semi_join(gene_list, by = "name") 
          })
        }
        
        df_list <- map(df_list, extract_gene_symbol) 
        
        name_list <- list("tm_1", "tm_2", "tm_con")
        name_list <- list("tm_con", "tm_1", "tm_2")
        names(df_list) <- name_list
        
        win_min <- input$win_min
        win_max <- input$win_max 
        
        
        # For testing win_min = 51, win_max = 195
        df_merge <- DRB_merge(df_list, win_min = win_min, win_max = win_max)
        
        # Normalized values and identified waves
        df_norm <- DRB_norm(df_merge, win_max = win_max)
        
        waves <- run_find_edge(df_norm)
        
        
        # Calculated elongation rates 
        tm <- input$time_2 - input$time_1
        
        rates <- waves %>%
          na.omit() %>% 
          filter(tm_2 > tm_1) %>%
          mutate(
            tm_1 = tm_1 / 1000,
            tm_2 = tm_2 / 1000,
            rate = (tm_2 - tm_1) / tm,
            rate = round(rate, digits = 1)
          ) 
        
        tm1_name <- str_c(input$time_1, "min")
        tm2_name <- str_c(input$time_2, "min")
        
        colnames(rates) <- c("Name", tm1_name, tm2_name, "Rate (kb/min)")
        
        rates
        
      })
        
      
      # Print table   
      output$contents <- renderDataTable({
        table_out()
      },
      options = list(
        autoWidth = TRUE,
        columnDefs = list(list(width = '75px', targets = "_all"))
      ))
      
      
      # Download table
      output$download <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".txt", sep="")
        },
        
        content = function(file) {
          write_tsv(table_out(), path = file)
        }
      )
    },
  
    
    # Return a safeError if a parsing error occurs
    error = function(e) {
      stop(safeError(e))
    }
  )
}


# Create Shiny app ----
shinyApp(ui, server)

