
library(depmixS4)
library(DT)
library(shiny)
library(shinythemes)
library(tidyverse)
library(magrittr)

# Color list 
color_list <- list(
  purple_1 = "#253494",
  orange_1 = "#fe9929",
  orange_2 = "#ec7014",
  brown_1  = "#cc4c02",
  brown_2  = "#993404",
  green_1  = "#78c679",
  green_2  = "#41ab5d",
  green_3  = "#238443",
  blue_1 = "#3288bd",
  blue_2 = "#225ea8",
  blue_3 = "#08519c",
  aqua_1 = "#41b6c4",
  grey_1 = "#969696",
  grey_2 = "#737373",
  grey_3 = "#525252",
  black  = "black",
  red_1  = "#fc4e2a",
  red_2  = "#e31a1c",
  red_3  = "#bd0026"
)

# Function to retrieve colors from color_list
get_color <- function(targets) {
  as.character(color_list[targets])
}

# Functions ----

# Function to add "key" columns to list of dfs
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

# Function to merge tables
tbl_merge <- function(input, merge_by = "name") {
  
  tbl_list <- add_key(input) 
  
  # Merged tables 
  res <- purrr::reduce(tbl_list, function(x, y) {
    left_join(x, y, by = merge_by)
  }) %>% 
    na.omit()
  
  res
}

# Function to merge tables
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
      filter(
        min(win_id) == win_min,
        max(win_id) == win_max,
        sum(count) > 0
      ) %>% 
      ungroup()
  })
  
  # Merged tables 
  res <- tbl_merge(res, merge_by = c("chrom", "start", "end", "name", "win_id")) %>%
    gather(key, count, -chrom, -start, -end, -name, -win_id)
  
  res
}

# Function to normalize signal
DRB_norm <- function(input_file, win_min = 1, win_max = 200) {
  
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
    dplyr::select(-max_count) %>%
    
    # Converted win_id to distance (kb)
    mutate(win_id = (win_id - win_min) * ((end - start) / 1000))
    
  res
}

# Function to run depmix
find_edge <- function(input_file) {
  win_tot <- input_file %>% 
    group_by(win_id) %>% 
    group_size() %>% 
    length()
  
  gene_tot <- input_file %>% 
    group_by(name) %>% 
    group_size() %>% 
    length()
  
  table_sort <- input_file %>% arrange(name, win_id)
  name <- table_sort$name
  win_id <- table_sort$win_id
  count <- table_sort$count
  res <- rep(NA, gene_tot)
  
  for (i in 0:(gene_tot - 1)) {
    count_in <- count[ (i * win_tot + 1) : (i * win_tot + win_tot) ]
    count_in <- data.frame(count_in) 
    win_in <- win_id[ (i * win_tot + 1) : (i * win_tot + win_tot) ] 
    
    mod <- depmix(response = count_in ~ 1, data = count_in, nstates = 2, trstart = runif(4))
    
    tryCatch(
      fm <- fit(mod, emc = em.control(rand = FALSE)),
      error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }
    )
    
    HMMstate <- posterior(fm)$state
    
    if (HMMstate %>% unique() %>% length() == 2) {
      edge <- rep(NA, win_tot)
      
      for (j in 1:length(HMMstate)) {
        if (j > 4) {
          sum_state <- sum(HMMstate[ (j - 4) : j ])  
          if (sum_state == 5) {
            edge[j] <- win_in[j]
          }
        }
      }
    }
    
    edge %<>%
      na.omit() %>%
      tail(1) 
    
    if (length(edge) == 1) {
      #edge <- edge * 500 + 500
      res[i + 1] <- edge
    }
  }
  
  res <- data.frame(unique(name), res)
  colnames(res) <- c("name", "wave_edge")
  res
}

# Function to run find_edge()
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

# Function to extract gene symbols
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
ui <- fluidPage(theme = "custom.css",
                
  #shinythemes::themeSelector(),
  
  titlePanel("Elongation rate calculator"),
  
  column(6,
  
    fluidRow(
      column(8,
        fileInput(
          "file_1", "Timecourse data",
          accept = c("text/tsv", ".bed")
        )
      ),
    
      column(4, 
        numericInput(
          "time_1", "Time (min)", 
          10, min = 0, max = 999
        )
      )
    ),
  
    fluidRow(
      column(8,
        fileInput(
          "file_2", label = NULL,
          accept = c("text/tsv", ".bed")
        )
      ),
    
      column(4, 
        numericInput(
          "time_2", label = NULL, 
          20, min = 1, max = 999
        )
      )  
    ),

    fluidRow(
      column(8,
        fileInput(
          "control", "Control data",
          accept = c("text/tsv", ".bed")
        )
      ),
    
      column(4,
        numericInput(
          "win_min", "Window min",
          1, min = 1, max = 200
        )
      )
    ),
    
    fluidRow(
      column(8,
        fileInput(
          "gene_list", "Gene list",
          accept = c("text/tsv", ".txt")
        )
      ),
    
      column(4,
        numericInput(
          "win_max", "Window max",
          200, min = 1, max = 200
        )
      )
    ),
    
    fluidRow(
      column(2, actionButton("runAnalysis", "RUN")),
      column(2, downloadButton("download", "Export")),
      column(2, offset = 1, actionButton("plotSelected", "Plot"))
    )
  ),
  
  column(6, 
    div(
      DT::dataTableOutput("rateTable"),
      style = "font-size: 75%; text-overflow: ellipsis"
    )
  ),
  
  
  
  fluidRow(plotOutput("genePlot"))
)


# Define server logic to read selected file ----
server <- function(input, output) {
  
  options(shiny.maxRequestSize = 500 * 1024 ^ 2)
  
  tryCatch(
    {
      tablesOut <- eventReactive(input$runAnalysis, ignoreInit = T, {
          
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
        
        #df_list <- map(df_list, extract_gene_symbol) 
        
        name_list <- list("tm_1", "tm_2", "tm_con")
        names(df_list) <- name_list
        
        win_min <- input$win_min
        win_max <- input$win_max 
        
        df_merge <- DRB_merge(df_list, win_min = win_min, win_max = win_max)
        
        
        # Normalized values and identified waves
        df_norm <- DRB_norm(df_merge, win_min = win_min, win_max = win_max)
        #wave_coords <- run_find_edge(df_norm)
        
        
        
        
        
        
        withProgress(message = "Calculating rates", {
        
          win_tot <- df_norm %>% 
            group_by(win_id) %>% 
            group_size() %>% 
            length()
          
          gene_tot <- df_norm %>% 
            group_by(name) %>% 
            group_size() %>% 
            length()
          
          data_tot <- df_norm %>%
            group_by(key) %>%
            group_size() %>%
            length()
          
          prog_tot <- data_tot * gene_tot
          prog_n <- 1
          
          input_file <- df_norm %>% 
            dplyr::select(-chrom, -start, -end) %>% 
            spread(key, count) 
          
          name_cols <- input_file %>% dplyr::select(name, win_id) 
          
          wave_coords <- input_file %>% 
            dplyr::select(name) %>% 
            unique()
          
          col_names <- "name"
          
          for (i in 3:ncol(input_file)) {
            input_data <- bind_cols(name_cols, input_file[, i])
            data_name <- names(input_data)[3]
            col_names <- c(col_names, data_name)
            input_data %<>% gather(key, count, -name, -win_id)
            
            
            #waves <- find_edge(input_data)
            
            
            table_sort <- input_data %>% arrange(name, win_id)
            name <- table_sort$name
            win_id <- table_sort$win_id
            count <- table_sort$count
            waves <- rep(NA, gene_tot)
            
            for (j in 0:(gene_tot - 1)) {
              count_in <- count[ (j * win_tot + 1) : (j * win_tot + win_tot) ]
              count_in <- data.frame(count_in) 
              win_in <- win_id[ (j * win_tot + 1) : (j * win_tot + win_tot) ] 
              
              tryCatch(
                {
                  HMMmodel <- depmix(response = count_in ~ 1, data = count_in, nstates = 2, trstart = runif(4))
                  fm <- fit(HMMmodel, emc = em.control(rand = FALSE))
                },
                error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }
              )
              
              HMMstate <- posterior(fm)$state
              
              if (HMMstate %>% unique() %>% length() == 2) {
                edge <- rep(NA, win_tot)
                
                for (n in 1:length(HMMstate)) {
                  if (n > 4) {
                    sum_state <- sum(HMMstate[ (n - 4) : n ])  
                    if (sum_state == 5) {
                      edge[n] <- win_in[n]
                    }
                  }
                }
              }
              
              edge %<>%
                na.omit() %>%
                tail(1) 
              
              if (length(edge) == 1) {
                #edge <- edge * 500 + 500
                waves[j + 1] <- edge
              }
              
              incProgress(1/prog_tot)
            }
            
            waves <- data.frame(unique(name), waves)
            colnames(waves) <- c("name", "wave_edge")
            
            wave_coords %<>% left_join(waves, by = "name") 
            colnames(wave_coords) <- col_names
          }
        })
        
        
        
        
        
        
        
        
        # Calculated elongation rates 
        tm <- input$time_2 - input$time_1
        
        rateTable <- wave_coords %>%
          na.omit() %>% 
          filter(tm_2 > tm_1) %>%
          mutate(
            rate = (tm_2 - tm_1) / tm,
            rate = round(rate, digits = 1)
          ) 
        
        colnames(rateTable) <- c(
          "Name", 
          str_c(input$time_1, "min"),
          str_c(input$time_2, "min"),
          "Rate (kb/min)"
        )
        
        list(rateTable, df_merge)
      })
        
      
      # Print table
      output$rateTable <- renderDataTable(
        datatable(
          tablesOut()[[1]],
          selection = list(mode = "single")
        )
      )
      
      rateTable_selected <- reactive({
        ids <- input$rateTable_rows_selected
        gene_name <- tablesOut()[[1]][ids, 1]
        wave_1    <- tablesOut()[[1]][ids, 2]
        wave_2    <- tablesOut()[[1]][ids, 3]
        rate      <- tablesOut()[[1]][ids, 4]
        list(gene_name, wave_1, wave_2, rate)
      })
      
      
      # Plot data for selected gene 
      plotOut <- eventReactive(input$plotSelected, ignoreInit = T, {
        
        # Gene target
        gene_text <- as.character(rateTable_selected()[[1]])
        gene_target <- as_data_frame(gene_text) %>% 
          dplyr::select(name = 1) 
        
        # Input file 
        tm1_name <- str_c(input$time_1, " min")
        tm2_name <- str_c(input$time_2, " min")
        
        input_file <- as_data_frame(tablesOut()[[2]]) %>%
          semi_join(gene_target, by = "name") %>%
          mutate(
            key = ifelse(key == "tm_1", tm1_name, key),
            key = ifelse(key == "tm_2", tm2_name, key),
            key = ifelse(key == "tm_con", "Control", key),
            key = fct_inorder(key),
            win_id = (win_id - input$win_min) * ((end - start) / 1000)
          ) %>%
          rename(Timepoint = key)
        
        max_value <- input_file %>% 
          group_by(name) %>% 
          summarize(max_value = max(count))
        
        max_value <- max_value$max_value * 0.9
        
        # Wave coordinates 
        wave_1 <- as.numeric(rateTable_selected()[[2]]) 
        wave_2 <- as.numeric(rateTable_selected()[[3]])
        wave_1_label <- str_c(wave_1, " kb")
        wave_2_label <- str_c(wave_2, " kb") 
        rate <- as.character(rateTable_selected()[[4]])
        
        # Plot colors 
        plot_colors <- get_color(c("red_1", "blue_2", "green_2"))
        
        # Plotted data for selected gene 
        input_file %>%
          ggplot(aes(win_id, count, color = Timepoint)) +
          geom_line(size = 2) +
          geom_vline(
            xintercept = c(wave_1, wave_2), 
            size = 1, linetype = 2,
            color = plot_colors[1:2]
          ) +
          scale_color_manual(values = plot_colors) +
          annotate("text", 
            x = wave_1 + 4, 
            y = max_value, 
            label = wave_1_label,
            size = 6,
            color = plot_colors[1]
          ) +
          annotate("text", 
            x = wave_2 + 4, 
            y = max_value, 
            label = wave_2_label, 
            size = 6,
            color = plot_colors[2]
          ) +
          labs(
            title = gene_text,
            subtitle = str_c(rate, " kb/min"),
            x = "Distance from TSS (kb)",
            y = ""
          ) +
          theme_classic() +
          theme(
            strip.background = element_blank(),
            plot.title = element_text(size = 35, face = "bold"),
            plot.subtitle = element_text(size = 20),
            axis.title = element_text(size = 20, face = "bold"),
            axis.line = element_line(size = 2),
            axis.ticks = element_line(size = 2),
            axis.ticks.length = unit(10, units = "point"),
            axis.text = element_text(size = 15, color = "black"),
            legend.title = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 18),
            legend.text.align = 0,
            legend.background = element_blank(),
            legend.position = c(0.8, 0.8)
          )
      })
      
      output$genePlot <- renderPlot(
          plotOut()
      )
  
      
      # Download table
      output$download <- downloadHandler(
        filename = function() {
          paste("data-", Sys.Date(), ".txt", sep="")
        },
        
        content = function(file) {
          write_tsv(tablesOut()[[1]], path = file)
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




