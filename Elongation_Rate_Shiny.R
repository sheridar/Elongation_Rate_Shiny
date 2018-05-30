
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
        gene_list <- read_tsv(genes_path, "name")
        
        # if (!is.null(genes)) {
        #   gene_list <- read_tsv(genes_path, "name")
        #   
        #   df_list <- map(df_list, function(x){
        #     x %>% 
        #       left_join(gene_list, by = "name") %>%
        #       na.omit() 
        #   })
        # }
        # 
        # name_list <- list("tm_con", "tm_1", "tm_2")
        # names(df_list) <- name_list
        
        
        #################
        # Merged tables #
        #################
        
        # # Function to merge tables
        # DRB_merge <- function(input, win_min = 1, win_max = 200) {
        #   
        #   # Function to add "key" columns to list of dfs
        #   add_key <- function(input) {
        #     
        #     tbl_list <- list()
        #     
        #     for (i in seq_along(input)) {
        #       x <- input[i] 
        #       y <- data.frame(x)
        #       
        #       new_names <- str_replace(colnames(y), str_c(names(x), "."), "")
        #       colnames(y) <- str_replace(new_names, "count", names(x))
        #       
        #       tbl_list <- c(tbl_list, list(y))
        #     }
        #     
        #     tbl_list
        #   }
        #   
        #   # Function to merge tables
        #   tbl_merge <- function(input, merge_by = "name") {
        #     
        #     tbl_list <- add_key(input) 
        #     
        #     # Merged tables 
        #     res <- purrr::reduce(tbl_list, function(x, y) {
        #       left_join(x, y, by = merge_by)
        #     }) %>% 
        #       na.omit()
        #     
        #     res
        #   }
        #   
        #   # Selected columns and filtered 
        #   res <- map(input, function(x) {
        #     x %>% 
        #       dplyr::select(chrom, start, end, name, win_id, count) %>%
        #       group_by(name) %>%
        #       filter(
        #         win_id >= win_min,
        #         win_id <= win_max
        #       ) %>%
        #       filter(
        #         min(win_id) == win_min,
        #         max(win_id) == win_max,
        #         sum(count) > 0
        #       ) %>% 
        #       ungroup()
        #   })
        #   
        #   # Merged tables 
        #   res <- tbl_merge(res, merge_by = c("chrom", "start", "end", "name", "win_id")) %>%
        #     gather(key, count, -chrom, -start, -end, -name, -win_id)
        #   
        #   res
        # }
        
        # Function to merge tables 
        DRB_merge <- function(input, gene_list, win_min = 1, win_max = 200, merge_by = c("name", "win_id")) {
          
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
            
            tbl_list <- add_key(input) 
            
            # Merge tables 
            res <- purrr::reduce(tbl_list, function(x, y) {
              left_join(x, y, ...)
            }) %>% 
              na.omit()
            
            res
          }
          
          # Table names
          tbl_names <- names(input)
          
          # Filter and calculate distance from TSS 
          res <- map(input, function(x) {
            
            # Filter by win_min and win_max
            x %>% 
              left_join(gene_list, by = "name") %>%
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
          })
          
          # Merge tables 
          res <- tbl_merge(res, by = merge_by) %>%
            gather_("key", "count", tbl_names)
          
          res
        }
        
        df_merge <- DRB_merge(df_list, gene_list, win_min = 51, win_max = 190, merge_by = c("name", "length", "win_id", "kb_dist"))
        
        # win_len <- df_merge %>% 
        #   mutate(len = (end - start) / 1000) %>% 
        #   group_by(len) %>% 
        #   summarize(n())
        
        
        #####################
        # Normalized tables #
        #####################
        
        # OLD VERSION
        # DRB_norm <- function(input_file, win_min = 1, win_max = 200, win_tot = 200) {
        #   
        #   # Function to calculate the distance from TSS 
        #   calc_kb <- function(input, id_col, len_col) {
        #     
        #     id_sym <- sym(id_col)
        #     
        #     input_sort <- input %>% 
        #       ungroup() %>% 
        #       arrange(!!id_sym) 
        #     
        #     lens <- c(input_sort[[len_col]])
        #     
        #     kb_tot <- 0
        #     kb_list <- vector("double", length(lens))
        #     
        #     for (i in seq_along(lens)) {
        #       kb_list[i] <- kb_tot
        #       kb_tot <- kb_tot + lens[i]
        #     }
        #     
        #     kb_list <- tibble(kb_dist = kb_list)
        #     res <- bind_cols(input_sort, kb_list)
        #     
        #     res
        #   }
        #   
        #   # Values used for window merging
        #   win_num <- length(seq(win_min, win_max))
        #   mutate_num <- round(win_num / win_tot)
        #   
        #   # Merged windows  
        #   res <- input_file %>% 
        #     mutate(
        #       win_id = win_id - win_min,
        #       win_id = floor(win_id / mutate_num),
        #       win_len = (end - start) / 1000
        #     ) %>% 
        #     group_by(name, key, win_id) %>% 
        #     summarize(
        #       count = mean(count), 
        #       win_len = sum(win_len)
        #     ) %>%
        #     group_by(key, name) %>%
        #     nest() %>%
        #     mutate(data = map(data, ~calc_kb(.x, id_col = "win_id", len_col = "win_len"))) %>%
        #     unnest() %>% 
        #     dplyr::select(-win_len)
        #   
        #   # Added pseudo count 
        #   win_cutoff <- max(res$win_id) - 5
        #   
        #   res %<>%  
        #     group_by(key, name) %>% 
        #     mutate(zero = ifelse(count == 0, T, F)) %>% 
        #     group_by(key, name, zero) %>% 
        #     mutate(min_count = min(count)) %>%
        #     group_by(key, name) %>% 
        #     mutate(count = ifelse(count == 0, max(min_count) / 2, count)) %>% 
        #     ungroup() %>% 
        #     dplyr::select(-zero, -min_count) %>% 
        #     
        #     # Normalized by -DRB signal 
        #     separate(key, sep = "_", into = c("treatment", "tm")) %>% 
        #     spread(tm, count) %>%
        #     gather(tm, count, -name, -win_id, -kb_dist, -treatment, -con) %>% 
        #     mutate(count = count / con) %>%
        #     dplyr::select(-con) %>% 
        #     
        #     # Normalized ratios by the average ratio for the last 5 bins 
        #     unite(key, treatment, tm, sep = "_") %>% 
        #     mutate(win_type = ifelse(win_id > win_cutoff, "con_wins", "data_wins")) %>% 
        #     group_by(key, name, win_type) %>% 
        #     mutate(ave_signal = mean(count)) %>% 
        #     ungroup() %>% 
        #     spread(win_type, ave_signal) %>% 
        #     mutate(con_wins = ifelse(is.na(con_wins), 0, con_wins)) %>% 
        #     group_by(name, key) %>% 
        #     mutate(con_wins = max(con_wins)) %>% 
        #     ungroup() %>% 
        #     mutate(count = count / con_wins) %>% 
        #     dplyr::select(-data_wins, -con_wins) %>% 
        #     
        #     # Digitized ratios to a range of 0 - 2.0 and step size of 0.5
        #     group_by(name, key) %>%
        #     mutate(max_count = max(count)) %>% 
        #     ungroup() %>% 
        #     mutate(
        #       count = (count / max_count) * 2,
        #       count = floor(count / 0.05) / 20
        #     ) %>%
        #     dplyr::select(-max_count, -win_id) %>%
        #     rename(win_id = kb_dist)
        #   
        #   res
        # }
        # 
        
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
            dplyr::select(name, key, "win_id" = kb_dist, count) %>%
            
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
            gather(tm, count, -name, -win_id, -treatment, -con) %>% 
            mutate(count = count / con) %>%
            dplyr::select(-con) %>% 
            unite(key, treatment, tm, sep = "_") %>% 
            
            # Bin values using a range of 0 - 1.0 and step size of 0.025
            group_by(name, key) %>%
            mutate(max_count = max(count)) %>% 
            ungroup() %>% 
            mutate(
              count = count / max_count,
              count = floor(count / 0.025) / 40
            ) %>%
            dplyr::select(-max_count)
          
          res
        }
        
        df_norm <- DRB_norm(df_merge, win_tot = 60)
        
        
        ###############################
        # Identified wave coordinates #
        ###############################
        
        # Function to run depmix
        find_waves <- function(input) {
          
          # Function to retrieve group size
          get_group_size <- function(input, col_name) {
            
            res <- input %>%
              group_by_(col_name) %>%
              group_size() %>% 
              length()
            
            res
          }
          
          # Function to run depmix
          run_depmix <- function(input, win_tot, gene_tot, prog_tot) {
            
            df_sort <- input %>% arrange(name, win_id)
            name <- df_sort$name
            win_id <- df_sort$win_id
            count <- df_sort$count
            res <- rep(NA, gene_tot)
            
            for (i in 0:(gene_tot - 1)) {
              count_in <- count[ (i * win_tot + 1) : (i * win_tot + win_tot) ]
              
              #count_in <- as.numeric(smooth(count_in)) # SMOOTHING
              
              win_in <- win_id[ (i * win_tot + 1) : (i * win_tot + win_tot) ] 
              
              trstart_val <- c(0.7, 0.2, 0.002, 0.3)
              HMMmod <- depmix(response = count_in ~ 1, data = data.frame(count_in), nstates = 2, trstart = trstart_val)
              
              tryCatch(
                HMMfit <- fit(HMMmod, emc = em.control(rand = FALSE)),
                error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }
              )
              
              if (exists("HMMfit")) {
                summary(HMMfit)
                
              } else {
                next
              }
              
              HMMstate <- posterior(HMMfit)$state
              
              if (HMMstate %>% unique() %>% length() == 2) {
                wave_edge <- rep(NA, win_tot)
                
                for (j in seq_along(HMMstate)) {
                  if (j > 4) {
                    sum_state <- sum(HMMstate[ (j - 4) : j ])  
                    if (sum_state == 5) {
                      wave_edge[j] <- win_in[j]
                    }
                  }
                }
              }
              
              wave_edge %<>%
                na.omit() %>%
                tail(1) 
              
              if (length(wave_edge) == 1) {
                res[i + 1] <- wave_edge
              }
              
              incProgress(1/prog_tot)
            }
            
            res <- data.frame(unique(name), res)
            colnames(res) <- c("name", "wave_edge")
            
            res
          }
          
          withProgress(message = "Calculating rates...", {
            
            # Group sizes
            win_tot <- get_group_size(input, "win_id")
            gene_tot <- get_group_size(input, "name")  
            data_tot <- get_group_size(input, "key")
            prog_tot <- data_tot * gene_tot
            
            # Spread table
            df_spread <- input %>% spread(key, count) 
            
            name_cols <- df_spread %>% dplyr::select(name, win_id) 
            
            res <- df_spread %>% 
              dplyr::select(name) %>% 
              unique()
            
            col_names <- "name"
            
            for (i in 3:ncol(df_spread)) {
              input_data <- bind_cols(name_cols, df_spread[, i])
              data_name <- names(input_data)[3]
              col_names <- c(col_names, data_name)
              input_data %<>% gather(key, count, -name, -win_id)
              
              waves <- run_depmix(input_data, win_tot, gene_tot, prog_tot)
              
              res %<>% left_join(waves, by = "name") 
              colnames(res) <- col_names
            }
          })
          
          res
        }
        
        wave_coords <- find_waves(df_norm)
        
        
        # Function to identify wave coordinates 
        # find_waves <- function(input_file) {
        #   
        #   # Function to retrieve group size
        #   get_group_size <- function(input_file, col_name) {
        #     
        #     target_col <- sym(col_name)
        #     
        #     res <- input_file %>%
        #       group_by(!!target_col) %>%
        #       group_size() %>% 
        #       length()
        #     
        #     res
        #   }
        #   
        #   # Function to run depmix
        #   run_depmix <- function(input_file, win_tot, gene_tot, prog_tot) {
        #     
        #     table_sort <- input_file %>% arrange(name, win_id)
        #     name <- table_sort$name
        #     win_id <- table_sort$win_id
        #     count <- table_sort$count
        #     res <- rep(NA, gene_tot)
        #     
        #     for (i in 0:(gene_tot - 1)) {
        #       count_in <- count[ (i * win_tot + 1) : (i * win_tot + win_tot) ]
        #       count_in <- as.numeric(smooth(count_in))
        #       win_in <- win_id[ (i * win_tot + 1) : (i * win_tot + win_tot) ] 
        #       
        #       trstart_val <- c(0.7, 0.2, 0.002, 0.3)
        #       HMMmod <- depmix(response = count_in ~ 1, data = data.frame(count_in), nstates = 2, trstart = trstart_val)
        #       #HMMmod <- depmix(response = count_in ~ 1, data = data.frame(count_in), nstates = 2, trstart = runif(4))
        #       
        #       tryCatch(
        #         HMMfit <- fit(HMMmod, emc = em.control(rand = FALSE)),
        #         error = function(e) { cat("ERROR :", conditionMessage(e), "\n") }
        #       )
        #       
        #       if (exists("HMMfit")) {
        #         summary(HMMfit)
        #         
        #       } else {
        #         next
        #       }
        #       
        #       HMMstate <- posterior(HMMfit)$state
        #       
        #       if (HMMstate %>% unique() %>% length() == 2) {
        #         edge <- rep(NA, win_tot)
        #         
        #         for (j in seq_along(HMMstate)) {
        #           if (j > 4) {
        #             sum_state <- sum(HMMstate[ (j - 4) : j ])  
        #             if (sum_state == 5) {
        #               edge[j] <- win_in[j]
        #             }
        #           }
        #         }
        #       }
        #       
        #       edge %<>%
        #         na.omit() %>%
        #         tail(1) 
        #       
        #       if (length(edge) == 1) {
        #         res[i + 1] <- edge
        #       }
        #       
        #       incProgress(1/prog_tot)
        #     }
        #     
        #     res <- data.frame(unique(name), res)
        #     colnames(res) <- c("name", "wave_edge")
        #     
        #     res
        #   }
        #   
        #   withProgress(message = "Calculating rates...", {
        #     
        #     win_tot <- get_group_size(input_file, "win_id")
        #     gene_tot <- get_group_size(input_file, "name")  
        #     data_tot <- get_group_size(input_file, "key")
        #     prog_tot <- data_tot * gene_tot
        #     
        #     shortened_df <- input_file %>% spread(key, count) 
        #     
        #     name_cols <- shortened_df %>% dplyr::select(name, win_id) 
        #     
        #     res <- shortened_df %>% 
        #       dplyr::select(name) %>% 
        #       unique()
        #     
        #     col_names <- "name"
        #     
        #     for (i in 3:ncol(shortened_df)) {
        #       input_data <- bind_cols(name_cols, shortened_df[, i])
        #       data_name <- names(input_data)[3]
        #       col_names <- c(col_names, data_name)
        #       input_data %<>% gather(key, count, -name, -win_id)
        #       
        #       waves <- run_depmix(input_data, win_tot, gene_tot, prog_tot)
        #       
        #       res %<>% left_join(waves, by = "name") 
        #       colnames(res) <- col_names
        #     }
        #   })
        #   res
        # }
        # 
        
        
        ###############################
        # Calculated elongation rates #
        ###############################
        
        col_names <- c(
          "Name", "Long_name", 
          tm1_name, tm2_name,
          "Rate (kb/min)"
        )
        
        # Function to calculate elongation rates 
        calc_rates <- function(input_file, time_1, time_2, col_names, win_min = 1, win_max = 200, win_len = 0.5) {
          
          # Function to extract gene symbols from dataframe
          extract_gene_symbol <- function(input_file) {
            
            # Function to extract gene symbol from string
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
          
          wave_max <- (win_max - win_min) * win_len - 5
          
          tm <- time_2 - time_1
          
          rate_table <- input_file %>%
            na.omit() %>% 
            filter(
              tm_2 > tm_1,
              tm_1 <= wave_max,
              tm_2 <= wave_max
            ) %>%
            mutate(
              rate = (tm_2 - tm_1) / tm,
              rate = round(rate, digits = 1),
              long_name = name
            ) %>%
            dplyr::select(long_name, name, tm_1, tm_2, rate)
          
          rate_table <- extract_gene_symbol(rate_table)
          
          colnames(rate_table) <- col_names
          
          rate_table
        }
        
        rate_table <- calc_rates(wave_coords, time_1, time_2, col_names, win_min, win_max, win_len)
        
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
        DRB_mean <- function(input_file, strand = F, relFreq = F) {
          
          if (strand == T) {
            res <- input_file %>%
              separate(key, sep = "_", into = c("key", "rep", "strand", "type")) %>%
              unite(key, key, rep, type, sep = "_")
          } 
          
          else res <- input_file
          
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
          input_file, 
          plot_title = NULL, 
          sub_title = NULL, 
          y_title = NULL,
          wave_1, wave_2, rate, 
          text_pos, 
          plot_colors = c("#41ab5d", "#cb181d", "#225ea8")
          ) {
        
          wave_1_lab <- str_c(wave_1, " kb")
          wave_2_lab <- str_c(wave_2, " kb")
          
          meta_plot <- input_file %>%
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
        
        df_merge <- as_data_frame(tablesOut()[[2]]) 
        
        win_len <- df_merge %>% 
          mutate(len = (end - start) / 1000) %>% 
          group_by(len) %>% 
          summarize(n())
        
        win_len <- win_len$len 
        
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
              win_id = (win_id - win_min) * win_len, 
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
              win_id = (win_id - win_min) * win_len,
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




