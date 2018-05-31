
# Install required packages ----
req_packages <- c(
  "depmixS4", "DT",
  "shiny", "tidyverse",
  "magrittr", "rlang"
)

avail_packages <- installed.packages()[, "Package"]

missing_packages <- req_packages[ !(req_packages %in% avail_packages) ]

if (length(missing_packages)) {
  install.packages(missing_packages)
}

for (i in seq_along(req_packages)) {
  library(req_packages[i], character.only = T)
}


# Define UI for data upload app ----
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css?family=Roboto:900');
      h1 {
        font-family: 'Roboto', sans-serif;
        font-weight: 500;
        font-size: 400%;
        line-height: 1.1;
        color: #cb181d;
      }
    "))
  ),
  
  headerPanel("Elongation rate calculator"),
  
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
      div(
        column(2, actionButton("runAnalysis", "RUN")),
        column(2, downloadButton("download", "Export")),
        column(2, offset = 1, actionButton("createPlot", "Plot")),
        style = "height: 75px; background-color: white;"
      )
    )
  ),
  
  column(6, 
    div(
      DT::dataTableOutput("rateTable"),
      style = "font-size: 75%; text-overflow: ellipsis"
    )
  ),
  
  fluidRow(
    #column(9, plotOutput(width = 925, "metaPlot")),
    column(9, plotOutput("metaPlot")),
    column(3, plotOutput("boxPlot"))
  )
)


# Define server logic to read selected file ----
server <- function(input, output) {
  
  options(shiny.maxRequestSize = 500 * 1024 ^ 2)
  
  tryCatch(
    {
      #####################################
      # Creates table of elongation rates #
      #####################################
      
      tablesOut <- eventReactive(input$runAnalysis, ignoreInit = T, {
          
        ################
        # Input values #
        ################
        
        file_1     <- input$file_1
        file_2     <- input$file_2
        con_path   <- input$control$datapath
        file1_path <- input$file_1$datapath
        file2_path <- input$file_2$datapath
        genes      <- input$gene_list
        genes_path <- input$gene_list$datapath
        
        time_1  <- input$time_1
        time_2  <- input$time_2
        tm1_name <- str_c(time_1, " min")
        tm2_name <- str_c(time_2, " min")
        
        win_min <- input$win_min
        win_max <- input$win_max
        
        req(file_1)
        req(file_2)
        req(time_1)
        req(time_2)
        
        
        ##################
        # Imported files #
        ##################
        
        col_names <- c(
          "chrom", "start",
          "end", "name",
          "win_id", "strand",
          "count"
        )
    
        file_list <- list(con_path, file1_path, file2_path)
        df_list <- map(file_list, function(x) read_tsv(x, col_names))
        gene_list <- read_tsv(genes_path, col_names[1:4])
        
        name_list <- list("tm_con", "tm_1", "tm_2")
        names(df_list) <- name_list
        
        
        ################
        # Merge tables #
        ################
        
        # Function to merge tables 
        DRB_merge <- function(input, gene_list, win_min, win_max, merge_by) {
          
          # Function to calculate distance from TSS 
          calc_kb <- function(input, id_col, len_col) {
            
            input_sort <- input %>% 
              ungroup() %>% 
              arrange(!!sym(id_col)) 
            
            lens <- c(input_sort[[len_col]])
            
            kb_tot <- 0
            kb_list <- vector("double", length(lens))
            
            for (i in seq_along(lens)) {
              kb_list[i] <- kb_tot
              kb_tot <- kb_tot + lens[i]
            }
            
            kb_list <- tibble(kb_dist = kb_list)
            res <- bind_cols(input_sort, kb_list)
            
            res
          }
          
          # Function to add "key" columns to list of dfs
          add_key <- function(input) {
            
            res <- list()
            
            for (i in seq_along(input)) {
              x <- input[i] 
              y <- data.frame(x)
              
              new_names <- str_replace(colnames(y), str_c(names(x), "."), "")
              colnames(y) <- str_replace(new_names, "count", names(x))
              
              res <- c(res, list(y))
            }
            
            res
          }
          
          # Function to merge tables
          tbl_merge <- function(input, ...) {
            
            tbl_list <- add_key(res) 
            
            # Merge tables 
            res <- purrr::reduce(tbl_list, function(x, y) {
              left_join(x, y, ...)
            }) %>% 
              na.omit()
            
            res
          }
          
          # Table names
          tbl_names <- names(input)
          
          # Calculate gene length
          genes <- gene_list %>%
            mutate(length = (end - start) / 1000) %>%
            dplyr::select(name, length)
          
          # Filter and calculate distance from TSS 
          res <- map(input, function(x) {
            
            # Filter by win_min and win_max
            res <- x %>% 
              left_join(genes, by = "name") %>%
              na.omit() %>% 
              dplyr::select(-strand) %>% 
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
              ungroup() %>%
              
              # Calculate distance for each window
              mutate(
                win_id = win_id - win_min,
                win_len = (end - start) / 1000
              ) %>% 
              group_by(name) %>%
              nest() %>%
              mutate(
                data = map(data, ~calc_kb(.x, id_col = "win_id", len_col = "win_len"))
              ) %>%
              unnest() %>%
              ungroup() %>%
              dplyr::select(name, length, win_id, kb_dist, count)
            
            res
          })
          
          # Merge tables 
          res <- tbl_merge(res, by = merge_by) %>%
            gather_("key", "count", tbl_names)
          
          res
        }
        
        df_merge <- DRB_merge(df_list, gene_list, win_min, win_max, merge_by = c("name", "length", "win_id", "kb_dist"))
        
        
        ####################
        # Normalize tables #
        ####################
        
        # Function to normalize signal
        DRB_norm <- function(input, win_tot = 60) {
          
          # Remove windows that are past pAS
          res <- input %>% 
            filter(kb_dist < length) %>% 
            group_by(key, name) %>% 
            mutate(win_count = n()) %>%
            filter(win_count >= win_tot) %>% 
            ungroup() %>%
            
            # Merge windows 
            mutate(
              mutate_num = round(win_count / win_tot),
              win_id = floor(win_id / mutate_num)
            ) %>% 
            group_by(key, name, win_id) %>%
            summarize(
              count = mean(count),
              kb_dist = min(kb_dist)
            ) %>%
            ungroup() %>%
            #dplyr::select(name, key, "win_id" = kb_dist, count) %>%
            dplyr::select(name, key, win_id, kb_dist, count) %>%
            
            # Add pseudo count
            group_by(key, name) %>% 
            mutate(zero = ifelse(count == 0, T, F)) %>% 
            group_by(key, name, zero) %>% 
            mutate(min_count = min(count)) %>%
            group_by(key, name) %>% 
            mutate(count = ifelse(count == 0, max(min_count) / 2, count)) %>% 
            ungroup() %>% 
            dplyr::select(-zero, -min_count) %>% 
            
            # Internally normalize signal 
            group_by(key, name) %>% 
            mutate(count = count / sum(count)) %>%
            ungroup() %>% 
            
            # Normalize by -DRB signal 
            separate(key, sep = "_", into = c("treatment", "tm")) %>% 
            spread(tm, count) %>%
            #gather(tm, count, -name, -win_id, -treatment, -con) %>% 
            gather(tm, count, -name, -win_id, -kb_dist, -treatment, -con) %>% 
            mutate(count = count / con) %>%
            dplyr::select(-con) %>% 
            unite(key, treatment, tm, sep = "_") #%>% 
            
            # Bin values using a range of 0 - 1.0 and step size of 0.025
            # group_by(name, key) %>%
            # mutate(max_count = max(count)) %>% 
            # ungroup() %>% 
            # mutate(
            #   count = count / max_count,
            #   count = floor(count / 0.025) / 40
            # ) %>%
            # dplyr::select(-max_count)
          
          res
        }
        
        df_norm <- DRB_norm(df_merge, win_tot = 60)
        
        
        #############################
        # Identify wave coordinates #
        #############################
        
        # Function to find waves using HMM
        find_HMM_waves <- function(input) {
          res <- input %>% 
            group_by(key, name) %>%
            nest() %>%
            
            mutate(
              data = map(data, function(x) {
                
                df_sort <- x %>% arrange(win_id) 
                
                wins <- df_sort$win_id
                wins_max <- wins[ (length(wins) - 5) ]
                counts <- df_sort$count
                
                trstart_vals <- c(0.7, 0.2, 0.002, 0.3) 
                HMMmod <- depmix(response = counts ~ 1, data = data.frame(counts), nstates = 2, trstart = trstart_vals)
                
                tryCatch(
                  HMMfit <- fit(HMMmod, emc = em.control(rand = FALSE)),
                  error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }
                )
                
                if (exists("HMMfit")) {
                  summary(HMMfit)
                  HMMstate <- posterior(HMMfit)$state
                  
                  if (HMMstate %>% unique() %>% length() == 2) {
                    wave_edge <- NA
                    
                    for (i in seq_along(HMMstate)) {
                      if (i > 4) {
                        sum_state <- sum(HMMstate[ (i - 4) : i ])  
                        
                        if (sum_state == 5) {
                          wave_edge <- wins[i]
                        }
                      }
                    }
                  }
                  
                  if (is.na(wave_edge)) {
                    wave_edge

                  } else if (wave_edge < wins_max) {
                    wave_edge
                    
                  } else {
                    NA
                  }
                  
                } else {
                  NA
                }
                
              })
            ) %>%

            ungroup() %>%
            mutate(type = map(data, function(x) typeof(x))) %>%
            unnest(type) %>%
            filter(type != "NULL") %>%
            dplyr::select(-type) %>%
            rename(wave_edge = data) %>% 
            unnest() 
          
          res
        }
        
        #wave_coords <- find_HMM_waves(df_norm)
        
        # Function to find waves using arbitrary cutoff
        find_waves <- function(input, sd_lim = 10) {
          
          res <- input %>%
            group_by(name, key) %>%
            mutate(
              win_max = max(win_id),
              win_min = win_max - 5,
              win_type = ifelse(win_id <= win_min, "data", "background")
            ) %>%                         
            group_by(name, key, win_type) %>% # Calculated mean count and sd for each timepoint
            mutate(
              mean_count = mean(count), 
              sd_count = sd(count),
              limit = mean_count + (sd_lim * sd_count)
            ) %>%
            group_by(name, key) %>% 
            mutate(limit = ifelse(win_id <= win_min, min(limit), limit)) %>%  
            filter(count > limit) %>% # Identified highest bin where the count is greater than the limit
            arrange(desc(win_id)) %>%
            dplyr::slice(1) %>% 
            ungroup() %>%
            filter(
              win_id > 0,
              win_id < win_max
            ) %>% 
            dplyr::select(name, key, "wave_edge" = kb_dist)
          
          res
        }
        
        wave_coords <- find_waves(df_norm, sd_lim = 10)
        
        
        ##############################
        # Calculate elongation rates #
        ##############################
        
        # Function to calculate elongation rates 
        calc_rates <- function(input, time_1, time_2, col_names, win_min = 1, win_max = 200) {
          
          # Function to extract gene symbols from dataframe
          extract_gene_symbol <- function(input) {
            
            # Function to extract gene symbol from string
            get_last_name <- function(gene_string) {
              res <- str_split(gene_string, "\\|")
              
              str_len <- length(res[[1]])
              
              res <- res[[1]][[str_len]]
              
              res
            }
            
            gene_names <- input %>%
              dplyr::select(name)
            
            other_data <- input %>%
              dplyr::select(-name)
            
            gene_matrix <- as.matrix(gene_names)
            
            new_names <- map(gene_matrix, get_last_name)
            
            new_names <- data.frame(name = as.matrix(new_names)) %>%
              mutate(name = as.character(name))
            
            res <- bind_cols(new_names, other_data)
          }
          
          # Calculate distance traveled
          tm <- time_2 - time_1
          
          # Calculate elongation rate 
          rate_table <- input %>%
            spread(key, wave_edge) %>% 
            na.omit() %>% 
            filter(tm_2 > tm_1) %>%
            mutate(
              rate = (tm_2 - tm_1) / tm,
              rate = round(rate, digits = 1),
              long_name = name
            ) %>%
            filter(rate > 0) %>% 
            dplyr::select(long_name, name, tm_1, tm_2, rate)
          
          # Extract gene symbols
          rate_table <- extract_gene_symbol(rate_table)
          
          # Update column names 
          col_names <- c(
            "Name", "Long_name", 
            tm1_name, tm2_name,
            "Rate (kb/min)"
          )
          
          colnames(rate_table) <- col_names
          
          rate_table
        }
        
        rate_table <- calc_rates(wave_coords, time_1, time_2, col_names, win_min, win_max)
        
        list(rate_table, df_merge)
      })
      
      
      # Output table 
      output$rateTable <- DT::renderDataTable(
        datatable(tablesOut()[[1]],
          options = list(
            columnDefs = list(list(visible = F, targets = c(2)))
          ),
          
          selection = list(mode = "multiple")
        )
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
      
      
      ###################
      # Create metaplot #
      ###################
      
      # Reactive to retrieve info for selected genes
      rateTable_selected <- reactive({
        ids <- input$rateTable_rows_selected
        gene_name <- tablesOut()[[1]][ids, 1]
        long_name <- tablesOut()[[1]][ids, 2]
        wave_1    <- tablesOut()[[1]][ids, 3]
        wave_2    <- tablesOut()[[1]][ids, 4]
        rate      <- tablesOut()[[1]][ids, 5]
        list(gene_name, long_name, wave_1, wave_2, rate)
      })
      
      metaplotOut <- eventReactive(input$createPlot, ignoreInit = T, {
        
        # Function to calculate mean signal 
        DRB_mean <- function(input, strand = F, relFreq = F) {
          
          if (strand == T) {
            res <- input %>%
              separate(key, sep = "_", into = c("key", "rep", "strand", "type")) %>%
              unite(key, key, rep, type, sep = "_")
          } 
          
          else res <- input
          
          if (relFreq == T) {
            res <- res %>%
              group_by(key, name) %>%
              mutate(count = count / sum(count)) %>%
              ungroup()
          }
          
          if (strand == T) {
            res <- res %>% 
              separate(key, sep = "_", into = c("key", "rep", "type")) %>%
              unite(key, key, rep, strand, type, sep = "_")
          }
          
          res <- res %>%
            group_by(key, win_id) %>%
            summarize(count = mean(count)) %>%
            ungroup()
          
          res
        }
        
        # Function to create metaplots 
        DRB_metaplot <- function(
          input, 
          plot_title = NULL, 
          sub_title = NULL, 
          y_title = NULL,
          wave_1, wave_2, rate, 
          text_pos, 
          plot_colors = c("#41ab5d", "#cb181d", "#225ea8")
          ) {
        
          wave_1_lab <- str_c(wave_1, " kb")
          wave_2_lab <- str_c(wave_2, " kb")
          
          meta_plot <- input %>%
            ggplot(aes(win_id, count, color = Timepoint)) +
            geom_line(size = 3) +
            geom_vline(
              xintercept = c(wave_1, wave_2), 
              size = 1, linetype = 2,
              color = plot_colors[2:3]
            ) +
            labs(
              subtitle = sub_title,
              x = "Distance from TSS (kb)",
              y = y_title
            ) +
            scale_color_manual(values = plot_colors) +
            annotate("text", 
              x = wave_1 + 5, 
              y = text_pos, 
              label = wave_1_lab,
              size = 6,
              color = plot_colors[2]
            ) +
            annotate("text", 
              x = wave_2 + 5, 
              y = text_pos, 
              label = wave_2_lab, 
              size = 6,
              color = plot_colors[3]
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
          
          if (!is.null(plot_title[[1]])) {
            meta_plot <- meta_plot + labs(title = plot_title)
          }
          
          meta_plot
        }
        
        # Input values 
        time_1  <- input$time_1
        time_2  <- input$time_2
        tm1_name <- str_c(time_1, " min")
        tm2_name <- str_c(time_2, " min")
        
        win_min <- input$win_min
        
        df_merge <- data.frame(tablesOut()[[2]]) %>%
          dplyr::select(name, key, "win_id" = kb_dist, count)
        
        # Create metaplot for selected genes
        if (!is.null(input$rateTable_rows_selected)) {
          
          # Gene targets
          gene_text <- as.character(rateTable_selected()[[1]])
          gene_targets <- as_data_frame(rateTable_selected()[[2]]) %>% 
            dplyr::select(name = 1) 
          
          # Wave coordinates 
          wave_1    <- round( mean( as.numeric( rateTable_selected()[[3]] )), digits = 1)
          wave_2    <- round( mean( as.numeric( rateTable_selected()[[4]] )), digits = 1)
          rate      <- as.numeric(rateTable_selected()[[5]])
          mean_rate <- round(mean(rate), digits = 1)
          med_rate  <- round(median(rate), digits = 1)
          
          # Input file 
          df_merge %<>% semi_join(gene_targets, by = "name")
          
          df_mean <- DRB_mean(df_merge)
          
          input_file <- df_mean %>% 
            mutate(
              key = ifelse(key == "tm_1", tm1_name, key),
              key = ifelse(key == "tm_2", tm2_name, key),
              key = ifelse(key == "tm_con", "Control", key),
              key = fct_relevel(key, c("Control", tm1_name, tm2_name))
            ) %>%
            rename(Timepoint = key)
          
          # Coordinates for plot labels 
          max_y <- input_file %>% 
            mutate(max_value = max(count)) %>%
            dplyr::select(max_value) %>% 
            unique()
          
          max_x <- input_file %>%
            mutate(max_value = max(win_id)) %>%
            dplyr::select(max_value) %>%
            unique()
          
          wave_text_y <- as.numeric(max_y) * 0.9
          rate_text_x <- as.numeric(max_x) * 0.745
          rate_text_y <- as.numeric(max_y) * 0.5
          
          # Created metaplots 
          if (length(gene_text) == 1) {
            DRB_metaplot(
              input_file, 
              plot_title = gene_text,
              sub_title = str_c(rate, " kb/min"),
              y_title = "",
              wave_1 = wave_1, 
              wave_2 = wave_2, 
              rate = rate, 
              text_pos = wave_text_y
            )
          
          } else {
            DRB_metaplot(
              input_file, 
              sub_title = "",
              y_title = "Average Signal",
              wave_1 = wave_1, 
              wave_2 = wave_2, 
              rate = rate, 
              text_pos = wave_text_y
            ) +
              annotate("text", 
                x = rate_text_x, 
                y = rate_text_y, 
                label = str_c("Mean: ", mean_rate, " kb/min\nMedian: ", med_rate, " kb/min"), 
                size = 6.5,
                hjust = 0
              )
          }
            
        } else {
          
          # Wave coordinates 
          gene_name <- tablesOut()[[1]][, 1]
          long_name <- tablesOut()[[1]][, 2]
          wave_1    <- as.numeric(tablesOut()[[1]][, 3])
          wave_2    <- as.numeric(tablesOut()[[1]][, 4])
          rate      <- as.numeric(tablesOut()[[1]][, 5])
          
          wave_1    <- round( mean( wave_1 ), digits = 1)
          wave_2    <- round( mean( wave_2 ), digits = 1)
          mean_rate <- round( mean( rate ), digits = 1)
          med_rate <- round( median( rate ), digits = 1)
          
          # Input file 
          df_mean <- DRB_mean(df_merge, strand = F, relFreq = F)
          
          input_file <- df_mean %>% 
            mutate(
              key = ifelse(key == "tm_1", tm1_name, key),
              key = ifelse(key == "tm_2", tm2_name, key),
              key = ifelse(key == "tm_con", "Control", key),
              key = fct_relevel(key, c("Control", tm1_name, tm2_name))
            ) %>%
            rename(Timepoint = key)
          
          # Coordinates for plot labels 
          max_y <- input_file %>% 
            mutate(max_value = max(count)) %>%
            dplyr::select(max_value) %>% 
            unique()
          
          max_x <- input_file %>%
            mutate(max_value = max(win_id)) %>%
            dplyr::select(max_value) %>%
            unique()
          
          wave_text_y <- as.numeric(max_y) * 0.9
          rate_text_x <- as.numeric(max_x) * 0.745
          rate_text_y <- as.numeric(max_y) * 0.5
          
          # Create metaplots
          DRB_metaplot(
            input_file, 
            sub_title = "",
            y_title = "Average Signal",
            wave_1 = wave_1, 
            wave_2 = wave_2, 
            rate = rate, 
            text_pos = wave_text_y
          ) +
            annotate("text", 
              x = rate_text_x, 
              y = rate_text_y, 
              label = str_c("Mean: ", mean_rate, " kb/min\nMedian: ", med_rate, " kb/min"), 
              size = 6.5,
              hjust = 0
            )
        }
      })
      
      # Output metaplot
      output$metaPlot <- renderPlot(metaplotOut())
      
      
      ##################
      # Create boxplot #
      ##################
      
      boxplotOut <- eventReactive(input$createPlot, ignoreInit = T, {
        
        # Function to create boxplot
        DRB_boxplot <- function(input_file, plot_colors = "#cb181d") {
          
          input_file %>% 
            ggplot(aes("", rate)) + 
            geom_boxplot(size = 2, color = plot_colors) +
            labs(
              title = "",
              x = "",
              y = "Elongation Rate (kb/min)"
            ) +
            theme_classic() +
            theme(
              legend.key.size   = element_blank(),
              strip.background  = element_blank(),
              axis.title        = element_text(size = 20, face = "bold"),
              axis.line         = element_line(size = 2),
              axis.ticks        = element_line(size = 2),
              axis.ticks.length = unit(10, units = "point"),
              axis.ticks.x      = element_blank(),
              axis.text         = element_text(size = 15, color = "black"),
              legend.title      = element_blank(),
              legend.text       = element_blank(),
              legend.background = element_blank()
            )
        }
        
        rate_table <- as_data_frame(tablesOut()[[1]]) %>%
          dplyr::select(
            name = 1, long_name = 2,
            wave_1 = 3, wave_2 = 4, 
            rate = 5
          ) %>%
          mutate(key = "dataset")
        
        box_len <- dim(rate_table)[[1]]
        
        if (!is.null(input$rateTable_rows_selected)) {
          new_points <- rateTable_selected()[[5]]
          
          point_len <- length(new_points)
          
          na_len <- box_len - point_len
          
          box_points <- c(new_points, rep(NA, na_len))
          
          if (point_len > 20) {
            
            rate_table <- data.frame(rate = new_points)
            
            box_output <- DRB_boxplot(rate_table)
             
          } else {
            box_output <- DRB_boxplot(rate_table) +
              geom_jitter(aes(y = box_points), color = "#225ea8", width = 0.1, height = 0, size = 4)
          }
          
        } else {
          box_output <- DRB_boxplot(rate_table)
        }
          
        box_output
      })
      
      # Output boxplot
      output$boxPlot <- renderPlot(boxplotOut())
      
    },
  
    # Return a safeError if a parsing error occurs
    error = function(e) {
      stop(safeError(e))
    }
  )
}


# Create Shiny app ----
shinyApp(ui, server)




