
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
        column(2, checkboxInput("HMMcheckbox", "HMM", value = TRUE)),
        column(2, checkboxInput("simpleCheckbox", "Simple", value = FALSE)),
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
  
  session$onSessionEnded(stopApp)
  options(shiny.maxRequestSize = 500 * 1024 ^ 2)
  
  tryCatch(
    {
      
      ####################################
      # Create table of elongation rates #
      ####################################
      
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
        
        win_min <- input$win_min
        win_max <- input$win_max
        
        req(file_1)
        req(file_2)
        req(time_1)
        req(time_2)
        
        
        ################
        # Import files #
        ################
        
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
          
          # NOT SURE WHY THIS DOESN'T WORK
          # Function to add "key" co lumns to list of dfs
          # add_key <- function(input) {
          # 
          #   tbl_names <- names(input)
          # 
          #   res <- map2(input, tbl_names, function(x, y) {
          #     x %>% 
          #       rename_(.dots = setNames("count", y)) %>%
          #       as_data_frame()
          #   })
          # 
          #   res
          # }
          
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
            mutate(Length = round((end - start) / 1000), digits = 1) %>%
            dplyr::select(name, Length)
          
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
              dplyr::select(name, Length, win_id, win_len, kb_dist, count)
            
            res
          })
          
          # Merge tables 
          res <- tbl_merge(res, by = merge_by) %>%
            gather_("key", "count", tbl_names)
          
          res
        }
        
        df_merge <- DRB_merge(
          df_list, gene_list, 
          win_min, win_max, 
          merge_by = c(
            "name", "Length", 
            "win_id", "win_len", 
            "kb_dist"
          ))
        
        
        ####################
        # Normalize tables #
        ####################
        
        # Function to normalize signal
        DRB_norm <- function(input, win_tot = 60) {
          
          # Remove windows that are past pAS
          res <- input %>% 
            filter(kb_dist < Length) %>% 
            group_by(key, name) %>% 
            mutate(win_count = n()) %>%
            filter(win_count >= win_tot) %>% 
            ungroup() %>%
            
            # Merge windows 
            mutate(
              mutate_num = round(win_count / win_tot),
              win_id = floor(win_id / mutate_num)
            ) %>% 
            mutate(count = count * win_len) %>% 
            group_by(key, name, win_id) %>%
            mutate(
              count = sum(count),
              win_len = sum(win_len),
              kb_dist = min(kb_dist)
            ) %>%
            unique() %>% 
            mutate(count = count / win_len) %>% 
            ungroup() %>%
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
                  
                  wave_edge <- NA
                  
                  if ( HMMstate %>% unique() %>% length() == 2 ) {
                    
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

            unnest() %>% 
            na.omit() %>% 
            dplyr::rename(wave_edge = data)
            
            # ungroup() %>%
            # mutate(type = map(data, function(x) typeof(x))) %>%
            # unnest(type) %>%
            # filter(type != "NULL") %>%
            # dplyr::select(-type) %>%
            # rename(wave_edge = data) %>% 
            # unnest() 
          
          res
        }
        
        HMM_coords <- find_HMM_waves(df_norm)
        
        # Function to find waves using arbitrary cutoff
        find_simple_waves <- function(input, sd_lim = 10) {
          
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
        
        simple_coords <- find_simple_waves(df_norm, sd_lim = 10)
        
        
        ##############################
        # Calculate elongation rates #
        ##############################
        
        # Function to calculate elongation rates 
        calc_rates <- function(input, time_1, time_2, prefix, win_min = 1, win_max = 200) {
          
          # Function to extract gene symbols from dataframe
          extract_gene_symbol <- function(input) {
            
            # Function to extract gene symbol from string
            get_last_name <- function(gene_string) {
              res <- str_split(gene_string, "\\|")
              
              str_len <- length(res[[1]])
              
              res <- res[[1]] [[str_len]]
              
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
          tm1_name <- str_c(prefix, time_1, "min", sep = " ")
          tm2_name <- str_c(prefix, time_2, "min", sep = " ")
          rate_name <- str_c(prefix, " rate (kb/min)")
          
          col_names <- c(
            "Name", "Long_name", 
            tm1_name, tm2_name,
            rate_name
          )
          
          colnames(rate_table) <- col_names
          
          rate_table
        }
        
        HMM_rates <- calc_rates(HMM_coords, time_1, time_2, prefix = "HMM", win_min, win_max)
        
        simple_rates <- calc_rates(simple_coords, time_1, time_2, prefix = "Simple", win_min, win_max)
        
        merged_rates <- left_join(HMM_rates, simple_rates, by = c("Long_name", "Name")) %>% 
          na.omit() 
        
        # Function to merge df_merge and rate tables 
        merge_tbls <- function(meta_tbl, rate_tbl) {
          meta_tbl %>%
            rename(Long_name = name) %>% 
            left_join(rate_tbl, by = c("Long_name")) %>%
            na.omit() 
        }
        
        HMM_meta_rates <- merge_tbls(df_merge, HMM_rates)
        simple_meta_rates <- merge_tbls(df_merge, simple_rates)
        merged_meta_rates <- merge_tbls(df_merge, merged_rates)
        
        list(HMM_meta_rates, simple_meta_rates, merged_meta_rates)
      })
      
      
      ####################################
      # Output table of elongation rates #
      ####################################
      
      # Function to simplify rate tables
      simplify_rate_tbls <- function(input) {
        input %>%
          dplyr::select(-win_id, -win_len, -kb_dist, -key, -count) %>%
          unique()
      }
      
      # Output table 
      output$rateTable <- DT::renderDataTable(
        if (input$HMMcheckbox == TRUE && input$simpleCheckbox == FALSE) {
          
          HMM_rates <- simplify_rate_tbls( tablesOut() [[1]] )
          
          datatable(HMM_rates,
            options = list(
              columnDefs = list(list(visible = F, targets = c(1)))
            ),
          
            selection = list(mode = "multiple")
          )
        } else if (input$HMMcheckbox == FALSE && input$simpleCheckbox == TRUE) {
          
          simple_rates <- simplify_rate_tbls( tablesOut() [[2]] )
          
          datatable(simple_rates,
            options = list(
              columnDefs = list(list(visible = F, targets = c(1)))
            ),
                    
            selection = list(mode = "multiple")
          )          
        } else if (input$HMMcheckbox == TRUE && input$simpleCheckbox == TRUE) {
          
          merged_rates <- simplify_rate_tbls( tablesOut() [[3]] )
          
          datatable(merged_rates,
            options = list(
              columnDefs = list(list(visible = F, targets = c(1)))
            ),
                    
            selection = list(mode = "multiple")
          )          
        }
      )
      
      # Download table
      output$download <- downloadHandler(
        filename = function() {
          str_c("data-", Sys.Date(), ".txt")
        },
        
        content = function(file) {
          if (input$HMMcheckbox == TRUE && input$simpleCheckbox == FALSE) {
            HMM_rates <- simplify_rate_tbls( tablesOut() [[1]] )    
            write_tsv(HMM_rates, path = file)
            
          } else if (input$HMMcheckbox == FALSE && input$simpleCheckbox == TRUE) {
            simple_rates <- simplify_rate_tbls( tablesOut() [[2]] )    
            write_tsv(simple_rates, path = file)
            
          } else if (input$HMMcheckbox == TRUE && input$simpleCheckbox == TRUE) {
            merged_rates <- simplify_rate_tbls( tablesOut() [[3]] )
            write_tsv(merged_rates, path = file)
          }
        }
      )
      
      
      ###################
      # Create metaplot #
      ###################
      
      # Reactive to create metaplots
      metaplotOut <- eventReactive(input$createPlot, ignoreInit = T, {
        
        # Input times
        time_1  <- input$time_1
        time_2  <- input$time_2
        tm1_name <- str_c(time_1, " min")
        tm2_name <- str_c(time_2, " min")
        
        # Input tables 
        if (input$HMMcheckbox == TRUE && input$simpleCheckbox == FALSE) {
          meta_rates_tbl <- tablesOut() [[1]]
          
        } else if (input$HMMcheckbox == FALSE && input$simpleCheckbox == TRUE) {
          meta_rates_tbl <- tablesOut() [[2]]
          
        } else if (input$HMMcheckbox == TRUE && input$simpleCheckbox == TRUE) {
          meta_rates_tbl <- tablesOut() [[3]]
        }
        
        # Function to simplify metaplot tables 
        simplify_meta_tbls <- function(input) {
          input %>% 
            dplyr::select(Long_name, key, kb_dist, count) %>%
            rename(win_id = kb_dist)
        }
        
        # Simplify metaplot and rate tables 
        meta_tbl <- simplify_meta_tbls(meta_rates_tbl)
        rate_tbl <- simplify_rate_tbls(meta_rates_tbl)
                
        # Function to create metaplots 
        DRB_metaplot <- function(meta_in, rate_in) {
          
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
          create_metaplots <- function(
            input, 
            plot_title = NULL, 
            sub_title = NULL, 
            y_title = NULL,
            waves, 
            line_type = 2,
            text_pos, 
            plot_colors = c("#41ab5d", "#cb181d", "#225ea8")
          ) {
            
            # Wave labels 
            wave_labels <- map(waves, function(input) {
              str_c(input, " kb")
            })
            
            meta_plot <- input %>%
              ggplot(aes(win_id, count, color = Timepoint)) +
              geom_line(size = 3) +
              scale_color_manual(values = plot_colors) +
              labs(
                subtitle = sub_title,
                x = "Distance from TSS (kb)",
                y = y_title
              ) +
              annotate("text",
                x = waves[[1]] + 5,
                y = text_pos,
                label = wave_labels[[1]],
                size = 6,
                color = plot_colors[2]
              ) +
              annotate("text",
                x = waves[[2]] + 5,
                y = text_pos,
                label = wave_labels[[2]],
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
              ) +
              geom_vline(
                xintercept = waves[1:2],
                size = 1, linetype = line_type,
                color = plot_colors[2:3]
              )
            
            if (length(waves) == 4) {
              meta_plot <- meta_plot +
                geom_vline(
                  xintercept = waves[3:4],
                  size = 1, linetype = 3,
                  color = plot_colors[2:3]
                ) +
                annotate("text",
                  x = waves[[3]] + 5,
                  y = text_pos * 0.8,
                  label = wave_labels[[3]],
                  size = 6,
                  color = plot_colors[2]
                ) +
                annotate("text",
                  x = waves[[4]] + 5,
                  y = text_pos * 0.8,
                  label = wave_labels[[4]],
                  size = 6,
                  color = plot_colors[3]
                )
            }
            
            if (!is.null(plot_title[[1]])) {
              meta_plot <- meta_plot + labs(title = plot_title)
            }
            
            meta_plot
          }
          
          # Wave coordinates 
          wave_1    <- round( mean( as.numeric( rate_in [, 4] )), digits = 1)
          wave_2    <- round( mean( as.numeric( rate_in [, 5] )), digits = 1)
          waves     <- c(wave_1, wave_2)
          rate      <- as.numeric( rate_in [, 6] )
          mean_rate <- round( mean( rate ), digits = 1)
          med_rate <- round( median( rate ), digits = 1)
          
          if (input$HMMcheckbox == TRUE && input$simpleCheckbox == TRUE) {
            sim_wave_1    <- round( mean( as.numeric( rate_in [, 7] )), digits = 1)
            sim_wave_2    <- round( mean( as.numeric( rate_in [, 8] )), digits = 1)
            sim_waves     <- c(sim_wave_1, sim_wave_2)
            sim_rate      <- as.numeric( rate_in [, 9] )
            sim_mean_rate <- round( mean( sim_rate ), digits = 1 )
            sim_med_rate  <- round( median( sim_rate ), digits = 1 )
            waves <- c(waves, sim_waves)
            mean_rate <- mean( c( mean_rate, sim_mean_rate ))
            med_rate <- mean( c( med_rate, sim_med_rate ))
          } 
          
          # Plot data 
          meta_mean <- DRB_mean(meta_in)
          
          plot_data <- meta_mean %>% 
            mutate(
              key = ifelse(key == "tm_1", tm1_name, key),
              key = ifelse(key == "tm_2", tm2_name, key),
              key = ifelse(key == "tm_con", "Control", key),
              key = fct_relevel(key, c("Control", tm1_name, tm2_name))
            ) %>%
            rename(Timepoint = key)
          
          # Coordinates for plot labels 
          max_y <- plot_data %>% 
            mutate(max_value = max(count)) %>%
            dplyr::select(max_value) %>% 
            unique()
          
          max_x <- plot_data %>%
            mutate(max_value = max(win_id)) %>%
            dplyr::select(max_value) %>%
            unique()
          
          wave_text_y <- as.numeric(max_y) * 0.9
          rate_text_x <- as.numeric(max_x) * 0.745
          rate_text_y <- as.numeric(max_y) * 0.5
          
          # Changed line type depending on wave-calling method
          if (input$HMMcheckbox == FALSE && input$simpleCheckbox == TRUE) {
            line_type <- 3
          } else {
            line_type <- 2
          }
            
          # Created metaplots 
          if (nrow(rate_in) == 1) {
            create_metaplots(
              plot_data, 
              plot_title = as.character( rate_in[, "Name"] ),
              sub_title = str_c(mean_rate, " kb/min"),
              y_title = "",
              waves = waves, 
              line_type = line_type,
              text_pos = wave_text_y
            )
            
          } else {
            create_metaplots(
              plot_data, 
              sub_title = "",
              y_title = "Average Signal",
              waves = waves,
              line_type = line_type,
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
        }
        
        # Reactive to retrieve info for selected genes
        rateTable_selected <- reactive({
          ids <- input$rateTable_rows_selected
          long_name   <- rate_tbl [ids, 1]
          gene_symbol <- rate_tbl [ids, 3]
          
          list(long_name, gene_symbol)
        })
        
        # Create metaplot for selected genes
        if (!is.null(input$rateTable_rows_selected)) {          
          
          # Selected genes 
          gene_long   <- data.frame( "Long_name" = rateTable_selected() [[1]] )
          gene_symbol <- as.character( rateTable_selected() [[2]] )
          
          meta_tbl <- meta_tbl %>% 
            semi_join(gene_long, by = "Long_name")
          
          rate_tbl <- rate_tbl %>%
            semi_join(gene_long, by = "Long_name")
        }
        
        DRB_metaplot(meta_tbl, rate_tbl)
      })
      
      # Output metaplot
      output$metaPlot <- renderPlot(metaplotOut())
      
      
      ##################
      # Create boxplot #
      ##################
      
      boxplotOut <- eventReactive(input$createPlot, ignoreInit = T, {
        
        # Function to create boxplot
        DRB_boxplot <- function(input, plot_colors) {
          
          input %>% 
            ggplot(aes(rate_key, rate)) + 
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
              axis.text         = element_text(size = 15, color = "black"),
              legend.title      = element_blank(),
              legend.text       = element_blank(),
              legend.background = element_blank()
            )
        }
        
        # Plot colors 
        plot_colors <- "#cb181d"
        
        # Input tables 
        if (input$HMMcheckbox == TRUE && input$simpleCheckbox == FALSE) {
          rate_tbl <- simplify_rate_tbls( tablesOut() [[1]] ) %>%
            gather(rate_key, rate, `HMM rate (kb/min)`)
          
        } else if (input$HMMcheckbox == FALSE && input$simpleCheckbox == TRUE) {
          rate_tbl <- simplify_rate_tbls( tablesOut() [[2]] ) %>%
            gather(rate_key, rate, `Simple rate (kb/min)`)
          
        } else if (input$HMMcheckbox == TRUE && input$simpleCheckbox == TRUE) {
          
          plot_colors <- c("#cb181d", "#225ea8")
          
          rate_tbl <- simplify_rate_tbls( tablesOut() [[3]] ) %>%
            gather(rate_key, rate, `HMM rate (kb/min)`, `Simple rate (kb/min)`)
        }
        
        rate_tbl <- rate_tbl %>% 
          dplyr::select(Long_name, rate_key, rate) %>%
          mutate(rate_key = str_replace(rate_key, "rate \\(kb/min\\)", ""))
        
        # Reactive to retrieve info for selected genes
        rateTable_selected <- reactive({
          ids <- input$rateTable_rows_selected
          long_name   <- rate_tbl [ids, 1]
          gene_symbol <- rate_tbl [ids, 3]
          
          list(long_name, gene_symbol)
        })  
        
        # Create boxplots for selected genes
        if (!is.null(input$rateTable_rows_selected)) {          
          
          # Selected genes 
          gene_long   <- data.frame( "Long_name" = rateTable_selected() [[1]] )
          gene_symbol <- as.character( rateTable_selected() [[2]] )
          select_rates <- rate_tbl %>% semi_join(gene_long, by = "Long_name")
          rate_names <- rate_tbl %>% dplyr::select(Long_name, rate_key)
          select_tbl <- rate_names %>% left_join(select_rates)
          
          if ( length(gene_symbol) < 21 ) {
            DRB_boxplot(rate_tbl, plot_colors) +
              geom_jitter(data = select_tbl, aes(rate_key, rate), color = "#41ab5d", width = 0.1, height = 0, size = 4)
              
          } else {
            DRB_boxplot(select_tbl, plot_colors)
          }
          
        } else {
          DRB_boxplot(rate_tbl, plot_colors)
        }          
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




