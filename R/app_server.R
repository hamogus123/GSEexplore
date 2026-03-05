#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @keywords internal
#' @export
app_server <- function(input, output, session) {
  # Reactive values to store data
  data_state <- shiny::reactiveValues(
    eset = NULL,
    data_loaded = FALSE,
    de_results = NULL,
    pca_result = NULL,
    batch_results = NULL,
    corrected_expr = NULL,
    batch_checked = FALSE
  )

  # ===== DATA LOADING =====
  shiny::observeEvent(input$btn_fetch, {
    shiny::showNotification("Fetching dataset...", type = "message", duration = NULL, id = "fetch_msg")

    tryCatch(
      {
        data_state$eset <- fetch_geo_dataset(input$gse_accession)
        data_state$data_loaded <- TRUE

        # Update phenotype column choices for downstream analyses
        pdata <- Biobase::pData(data_state$eset)
        phenotype_cols <- colnames(pdata)
        shiny::updateSelectInput(session, "de_phenotype_col", choices = phenotype_cols)
        shiny::updateSelectInput(session, "pca_color_by", choices = c("No coloring", phenotype_cols))
        shiny::updateSelectInput(session, "batch_column", choices = c("None", phenotype_cols))
        shiny::updateSelectInput(session, "heatmap_color_by", choices = c("None", phenotype_cols))
        shiny::updateSelectInput(session, "gene_select", choices = Biobase::featureNames(data_state$eset))

        shiny::removeNotification("fetch_msg")
        shiny::showNotification("Dataset loaded successfully!", type = "message", duration = 3)
      },
      error = function(e) {
        shiny::removeNotification("fetch_msg")
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 5)
      }
    )
  })

  # ===== DATA LOADING TAB OUTPUTS =====
  output$loading_status <- shiny::renderUI({
    if (!data_state$data_loaded) {
      shiny::p("No dataset loaded. Enter a GSE accession and click 'Fetch Dataset'.")
    } else {
      shiny::div(
        shiny::tags$span("[Loaded]", style = "color: green; font-weight: bold;"),
        shiny::br(),
        shiny::br(),
        shiny::actionButton("btn_clear", "Clear Dataset")
      )
    }
  })

  output$dataset_info <- shiny::renderUI({
    if (!data_state$data_loaded) {
      return(NULL)
    }

    summary_info <- get_dataset_summary(data_state$eset)
    shiny::div(
      shiny::h4("Dataset Summary:"),
      shiny::tags$ul(
        shiny::tags$li("Features (genes): ", summary_info$n_features),
        shiny::tags$li("Samples: ", summary_info$n_samples),
        shiny::tags$li("Expression range: ", round(summary_info$expression_range[1], 2), " to ", round(summary_info$expression_range[2], 2))
      )
    )
  })

  output$sample_metadata_table <- DT::renderDataTable({
    if (!data_state$data_loaded) {
      return(NULL)
    }
    pdata <- Biobase::pData(data_state$eset)
    DT::datatable(pdata, rownames = TRUE, options = list(pageLength = 5))
  })

  # ===== QC TAB: PCA =====
  shiny::observe({
    if (!data_state$data_loaded) {
      return(NULL)
    }

    data_list <- prepare_analysis_data(data_state$eset)
    data_state$pca_result <- compute_pca(data_list$expr_matrix)
  })

  output$pca_plot <- plotly::renderPlotly({
    if (!data_state$data_loaded || is.null(data_state$pca_result)) {
      return(NULL)
    }

    pca_scores <- data_state$pca_result$scores
    var_exp <- data_state$pca_result$var_explained

    # Extract PC1 and PC2
    pc1 <- pca_scores[, 1]
    pc2 <- pca_scores[, 2]

    # Color by phenotype if selected
    colors <- NULL
    if (input$pca_color_by != "No coloring") {
      pdata <- Biobase::pData(data_state$eset)
      colors <- as.factor(pdata[[input$pca_color_by]])
    }

    df_plot <- data.frame(
      PC1 = pc1,
      PC2 = pc2,
      Sample = names(pc1),
      Color = if (!is.null(colors)) colors else "Group"
    )

    p <- plotly::plot_ly(df_plot) %>%
      plotly::add_markers(
        x = ~PC1,
        y = ~PC2,
        text = ~Sample,
        color = ~Color,
        mode = "markers",
        marker = list(size = 8)
      ) %>%
      plotly::layout(
        title = "PCA Plot",
        xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
        yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)"))
      )

    p
  })

  output$distributions_plot <- plotly::renderPlotly({
    if (!data_state$data_loaded) {
      return(NULL)
    }

    data_list <- prepare_analysis_data(data_state$eset)
    dist_stats <- calculate_sample_distributions(data_list$expr_matrix)

    p <- plotly::plot_ly(dist_stats) %>%
      plotly::add_markers(
        x = ~sample_name,
        y = ~mean,
        error_y = list(
          type = "data",
          array = ~sd,
          visible = TRUE
        ),
        mode = "markers",
        marker = list(size = 8)
      ) %>%
      plotly::layout(
        title = "Sample Expression Distributions",
        xaxis = list(title = "Sample"),
        yaxis = list(title = "Mean Expression"),
        hovermode = "closest"
      )

    p
  })

  output$dist_stats_table <- DT::renderDataTable({
    if (!data_state$data_loaded) {
      return(NULL)
    }

    data_list <- prepare_analysis_data(data_state$eset)
    dist_stats <- calculate_sample_distributions(data_list$expr_matrix)
    DT::datatable(dist_stats, options = list(pageLength = 5))
  })

  # ===== DE TAB =====
  shiny::observeEvent(input$btn_run_de, {
    if (!data_state$data_loaded) {
      shiny::showNotification("Please load a dataset first", type = "warning")
      return()
    }

    if (input$de_phenotype_col == "" || input$de_group1 == "" || input$de_group2 == "") {
      shiny::showNotification("Please select all required fields", type = "warning")
      return()
    }

    shiny::showNotification("Running differential expression analysis...", type = "message", duration = NULL, id = "de_msg")

    tryCatch(
      {
        data_state$de_results <- run_differential_expression(
          data_state$eset,
          phenotype_col = input$de_phenotype_col,
          group1 = input$de_group1,
          group2 = input$de_group2
        )

        shiny::removeNotification("de_msg")
        shiny::showNotification("Analysis complete!", type = "message", duration = 3)
      },
      error = function(e) {
        shiny::removeNotification("de_msg")
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 5)
      }
    )
  })

  shiny::observeEvent(input$de_phenotype_col, {
    if (!data_state$data_loaded || input$de_phenotype_col == "") {
      return()
    }

    pdata <- Biobase::pData(data_state$eset)
    unique_vals <- unique(as.character(pdata[[input$de_phenotype_col]]))
    shiny::updateSelectInput(session, "de_group1", choices = unique_vals)
    shiny::updateSelectInput(session, "de_group2", choices = unique_vals)
  })

  output$de_status <- shiny::renderUI({
    if (is.null(data_state$de_results)) {
      return(NULL)
    }
    shiny::div(
      shiny::tags$span("[Complete]", style = "color: green; font-weight: bold;")
    )
  })

  output$volcano_plot <- plotly::renderPlotly({
    if (is.null(data_state$de_results)) {
      return(NULL)
    }

    de_res <- data_state$de_results
    de_res$neg_log10_pval <- -log10(de_res$P.Value)
    threshold <- -log10(0.05)

    p <- plotly::plot_ly(de_res) %>%
      plotly::add_markers(
        x = ~logFC,
        y = ~neg_log10_pval,
        text = ~gene_id,
        mode = "markers",
        marker = list(size = 6),
        hovertemplate = "<b>%{text}</b><br>logFC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>"
      ) %>%
      plotly::add_trace(
        x = c(-4, 4),
        y = c(threshold, threshold),
        mode = "lines",
        line = list(dash = "dash", color = "red"),
        name = "p=0.05",
        hovertemplate = "threshold<extra></extra>"
      ) %>%
      plotly::layout(
        title = "Volcano Plot",
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10(p-value)")
      )

    p
  })

  output$de_results_table <- DT::renderDataTable({
    if (is.null(data_state$de_results)) {
      return(NULL)
    }
    DT::datatable(data_state$de_results, options = list(pageLength = 10, order = list(list(5, "asc"))))
  })

  output$download_de_results <- shiny::downloadHandler(
    filename = function() {
      paste0("DE_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(data_state$de_results)) {
        utils::write.csv(data_state$de_results, file, row.names = FALSE)
      }
    }
  )

  # ===== EXPLORATION TAB =====
  output$gene_expr_plot <- plotly::renderPlotly({
    if (!data_state$data_loaded || input$gene_select == "") {
      return(NULL)
    }

    gene_id <- input$gene_select
    data_list <- prepare_analysis_data(data_state$eset)
    gene_idx <- which(data_list$feature_names == gene_id)[1]

    if (is.na(gene_idx)) {
      return(NULL)
    }

    expr_values <- data_list$expr_matrix[gene_idx, ]
    sample_names <- data_list$sample_names

    df_plot <- data.frame(
      Sample = sample_names,
      Expression = expr_values
    )

    p <- plotly::plot_ly(df_plot) %>%
      plotly::add_bars(
        x = ~Sample,
        y = ~Expression
      ) %>%
      plotly::layout(
        title = paste("Expression of", gene_id),
        xaxis = list(title = "Sample"),
        yaxis = list(title = "Expression Level")
      )

    p
  })

  output$gene_expr_table <- DT::renderDataTable({
    if (!data_state$data_loaded || input$gene_select == "") {
      return(NULL)
    }

    gene_id <- input$gene_select
    data_list <- prepare_analysis_data(data_state$eset)
    gene_idx <- which(data_list$feature_names == gene_id)[1]

    if (is.na(gene_idx)) {
      return(NULL)
    }

    expr_values <- data_list$expr_matrix[gene_idx, ]
    pdata <- data_list$phenotypes

    result_df <- data.frame(
      Sample = data_list$sample_names,
      Expression = expr_values,
      pdata
    )

    DT::datatable(result_df, options = list(pageLength = 5))
  })

  # ===== CLEAR DATA =====
  shiny::observeEvent(input$btn_clear, {
    data_state$eset <- NULL
    data_state$data_loaded <- FALSE
    data_state$de_results <- NULL
    data_state$pca_result <- NULL
    data_state$batch_results <- NULL
    data_state$corrected_expr <- NULL
    data_state$batch_checked <- FALSE
    shiny::updateTextInput(session, "gse_accession", value = "")
    shiny::showNotification("Dataset cleared", type = "message", duration = 2)
  })

  # ===== BATCH CORRECTION TAB =====
  shiny::observeEvent(input$btn_check_batch, {
    if (!data_state$data_loaded) {
      shiny::showNotification("Please load a dataset first", type = "warning")
      return()
    }

    if (input$batch_column == "None" || input$batch_column == "") {
      shiny::showNotification("Please select a batch variable", type = "warning")
      return()
    }

    shiny::showNotification("Checking for batch effects...", type = "message", duration = NULL, id = "batch_msg")

    tryCatch(
      {
        data_list <- prepare_analysis_data(data_state$eset)
        data_state$batch_results <- detect_batch_effects(
          data_list$expr_matrix,
          data_list$phenotypes,
          input$batch_column
        )
        data_state$batch_checked <- TRUE

        shiny::removeNotification("batch_msg")
        shiny::showNotification("Batch analysis complete!", type = "message", duration = 3)
      },
      error = function(e) {
        shiny::removeNotification("batch_msg")
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 5)
      }
    )
  })

  output$batch_status <- shiny::renderUI({
    if (!data_state$batch_checked || is.null(data_state$batch_results)) {
      return(NULL)
    }

    batch_info <- data_state$batch_results
    if (batch_info$has_batch) {
      shiny::div(
        shiny::tags$span("[Batch detected]", style = "color: orange; font-weight: bold;"),
        shiny::br(),
        shiny::p(paste("Surrogate variables:", batch_info$n_sv))
      )
    } else {
      shiny::div(
        shiny::tags$span("[No batch detected]", style = "color: green; font-weight: bold;")
      )
    }
  })

  output$has_batch_effects <- shiny::reactive({
    !is.null(data_state$batch_results) && data_state$batch_results$has_batch
  })
  shiny::outputOptions(output, "has_batch_effects", suspendWhenHidden = FALSE)

  output$batch_summary_table <- DT::renderDataTable({
    if (!data_state$batch_checked || is.null(data_state$batch_results)) {
      return(NULL)
    }

    batch_info <- data_state$batch_results
    summary_df <- data.frame(
      Metric = c("Surrogate Variables", "Batch Detected"),
      Value = c(batch_info$n_sv, ifelse(batch_info$has_batch, "Yes", "No"))
    )

    DT::datatable(summary_df, options = list(pageLength = 5, dom = "t"))
  })

  shiny::observeEvent(input$btn_correct_batch, {
    if (!data_state$data_loaded || is.null(data_state$batch_results)) {
      shiny::showNotification("Please check for batch effects first", type = "warning")
      return()
    }

    if (input$batch_column == "None" || input$batch_column == "") {
      shiny::showNotification("Please select a batch variable", type = "warning")
      return()
    }

    shiny::showNotification("Correcting batch effects...", type = "message", duration = NULL, id = "corr_msg")

    tryCatch(
      {
        data_list <- prepare_analysis_data(data_state$eset)
        batch_var <- data_list$phenotypes[[input$batch_column]]

        data_state$corrected_expr <- correct_batch_effects(
          data_list$expr_matrix,
          batch_var
        )

        shiny::removeNotification("corr_msg")
        shiny::showNotification("Batch correction complete!", type = "message", duration = 3)
      },
      error = function(e) {
        shiny::removeNotification("corr_msg")
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 5)
      }
    )
  })

  output$correction_status <- shiny::renderUI({
    if (is.null(data_state$corrected_expr)) {
      return(NULL)
    }
    shiny::div(
      shiny::tags$span("[Corrected]", style = "color: blue; font-weight: bold;")
    )
  })

  # ===== HEATMAP TAB =====
  shiny::observeEvent(input$btn_plot_heatmap, {
    if (!data_state$data_loaded) {
      shiny::showNotification("Please load a dataset first", type = "warning")
      return()
    }

    shiny::showNotification("Generating heatmap...", type = "message", duration = NULL, id = "heat_msg")

    tryCatch(
      {
        data_list <- prepare_analysis_data(data_state$eset)

        # Select which expression matrix to use
        expr_matrix <- if (!is.null(data_state$corrected_expr)) {
          data_state$corrected_expr
        } else {
          data_list$expr_matrix
        }

        # Determine gene indices
        gene_indices <- NULL
        if (input$heatmap_gene_source == "Top DE genes" && !is.null(data_state$de_results)) {
          top_n <- min(input$heatmap_n_genes, nrow(data_state$de_results))
          top_genes <- order(data_state$de_results$adj.P.Val)[1:top_n]
          gene_indices <- which(data_list$feature_names %in% data_state$de_results$gene_id[top_genes])
        } else if (input$heatmap_gene_source == "Most variable") {
          top_var <- get_top_variable_genes(expr_matrix, top_n = input$heatmap_n_genes)
          gene_indices <- top_var$row_index
        } else {
          gene_indices <- seq_len(min(input$heatmap_n_genes, nrow(expr_matrix)))
        }

        # Plot heatmap
        color_by <- if (input$heatmap_color_by != "None") input$heatmap_color_by else NULL

        p <- plot_expression_heatmap(
          expr_matrix,
          data_list$phenotypes,
          gene_indices = gene_indices,
          color_by = color_by,
          main = paste("Expression Heatmap -", input$heatmap_gene_source)
        )

        shiny::removeNotification("heat_msg")
      },
      error = function(e) {
        shiny::removeNotification("heat_msg")
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = 5)
      }
    )
  })

  output$heatmap_plot_container <- shiny::renderUI({
    if (!data_state$data_loaded) {
      return(NULL)
    }
    shiny::plotOutput("heatmap_output", height = "700px")
  })

  output$heatmap_output <- shiny::renderPlot({
    if (!data_state$data_loaded) {
      return(NULL)
    }

    data_list <- prepare_analysis_data(data_state$eset)

    # Select which expression matrix to use
    expr_matrix <- if (!is.null(data_state$corrected_expr)) {
      data_state$corrected_expr
    } else {
      data_list$expr_matrix
    }

    # Determine gene indices
    gene_indices <- NULL
    if (input$heatmap_gene_source == "Top DE genes" && !is.null(data_state$de_results)) {
      top_n <- min(input$heatmap_n_genes, nrow(data_state$de_results))
      top_genes <- order(data_state$de_results$adj.P.Val)[1:top_n]
      gene_indices <- which(data_list$feature_names %in% data_state$de_results$gene_id[top_genes])
    } else if (input$heatmap_gene_source == "Most variable") {
      top_var <- get_top_variable_genes(expr_matrix, top_n = input$heatmap_n_genes)
      gene_indices <- top_var$row_index
    } else {
      gene_indices <- seq_len(min(input$heatmap_n_genes, nrow(expr_matrix)))
    }

    # Create color annotation
    color_by <- if (input$heatmap_color_by != "None") input$heatmap_color_by else NULL
    annotation_col <- NULL
    if (!is.null(color_by) && color_by %in% colnames(data_list$phenotypes)) {
      annotation_col <- data.frame(
        Group = data_list$phenotypes[[color_by]],
        row.names = colnames(expr_matrix)
      )
    }

    # Create heatmap
    pheatmap::pheatmap(
      expr_matrix[gene_indices, ],
      scale = "row",
      annotation_col = annotation_col,
      main = paste("Expression Heatmap -", input$heatmap_gene_source),
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation"
    )
  })
}
