#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    shiny::fluidPage(
      # Header
      shiny::h1("GSEexplore: GEO Expression Data Explorer"),
      shiny::hr(),

      # Main content with tabs
      shiny::tabsetPanel(
        # Tab 1: Data Loading
        shiny::tabPanel(
          "Data Loading",
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              width = 6,
              shiny::textInput(
                inputId = "gse_accession",
                label = "GEO Accession (e.g., GSE63310):",
                placeholder = "GSE63310"
              ),
              shiny::actionButton(
                inputId = "btn_fetch",
                label = "Fetch Dataset",
                class = "btn-primary"
              ),
              shiny::br(),
              shiny::br(),
              shiny::uiOutput("loading_status")
            ),
            shiny::column(
              width = 6,
              shiny::uiOutput("dataset_info")
            )
          ),
          shiny::br(),
          shiny::h4("Sample Metadata:"),
          DT::dataTableOutput("sample_metadata_table")
        ),

        # Tab 2: QC Views
        shiny::tabPanel(
          "QC Views",
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::selectInput(
                inputId = "pca_color_by",
                label = "Color by (PCA):",
                choices = c("No coloring")
              )
            )
          ),
          shiny::fluidRow(
            shiny::column(width = 6, plotly::plotlyOutput("pca_plot")),
            shiny::column(width = 6, plotly::plotlyOutput("distributions_plot"))
          ),
          shiny::br(),
          shiny::h4("Sample Distribution Statistics:"),
          DT::dataTableOutput("dist_stats_table")
        ),

        # Tab 3: Differential Expression
        shiny::tabPanel(
          "Differential Expression",
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::selectInput(
                inputId = "de_phenotype_col",
                label = "Phenotype Column:",
                choices = c()
              )
            ),
            shiny::column(
              width = 3,
              shiny::selectInput(
                inputId = "de_group1",
                label = "Control Group:",
                choices = c()
              )
            ),
            shiny::column(
              width = 3,
              shiny::selectInput(
                inputId = "de_group2",
                label = "Treatment Group:",
                choices = c()
              )
            ),
            shiny::column(
              width = 3,
              shiny::br(),
              shiny::actionButton(
                inputId = "btn_run_de",
                label = "Run Analysis",
                class = "btn-primary"
              )
            )
          ),
          shiny::br(),
          shiny::uiOutput("de_status"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              width = 6,
              plotly::plotlyOutput("volcano_plot")
            ),
            shiny::column(
              width = 6,
              shiny::downloadButton("download_de_results", "Download Results (CSV)")
            )
          ),
          shiny::br(),
          shiny::h4("Differential Expression Results:"),
          DT::dataTableOutput("de_results_table")
        ),

        # Tab 4: Batch Correction
        shiny::tabPanel(
          "Batch Correction",
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              width = 6,
              shiny::selectInput(
                inputId = "batch_column",
                label = "Batch Variable:",
                choices = c("None")
              ),
              shiny::actionButton(
                inputId = "btn_check_batch",
                label = "Check for Batch Effects",
                class = "btn-info"
              )
            ),
            shiny::column(
              width = 6,
              shiny::uiOutput("batch_status")
            )
          ),
          shiny::br(),
          shiny::h4("Batch Effect Summary:"),
          DT::dataTableOutput("batch_summary_table"),
          shiny::br(),
          shiny::conditionalPanel(
            condition = "output.has_batch_effects",
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::actionButton(
                  inputId = "btn_correct_batch",
                  label = "Correct Batch Effects",
                  class = "btn-warning"
                )
              ),
              shiny::column(
                width = 6,
                shiny::uiOutput("correction_status")
              )
            )
          )
        ),

        # Tab 5: Heatmap Visualization
        shiny::tabPanel(
          "Heatmap",
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              width = 3,
              shiny::numericInput(
                inputId = "heatmap_n_genes",
                label = "Number of Genes:",
                value = 50,
                min = 10,
                max = 500
              )
            ),
            shiny::column(
              width = 3,
              shiny::selectInput(
                inputId = "heatmap_gene_source",
                label = "Gene Selection:",
                choices = c("Top DE genes", "Most variable", "All genes")
              )
            ),
            shiny::column(
              width = 3,
              shiny::selectInput(
                inputId = "heatmap_color_by",
                label = "Color Samples by:",
                choices = c("None")
              )
            ),
            shiny::column(
              width = 3,
              shiny::br(),
              shiny::actionButton(
                inputId = "btn_plot_heatmap",
                label = "Generate Heatmap",
                class = "btn-primary"
              )
            )
          ),
          shiny::br(),
          shiny::uiOutput("heatmap_plot_container")
        ),

        # Tab 6: Exploration (moved from position 4)
        shiny::tabPanel(
          "Exploration",
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              width = 6,
              shiny::selectInput(
                inputId = "gene_select",
                label = "Select a gene to explore:",
                choices = c(),
                selectize = TRUE
              )
            )
          ),
          shiny::br(),
          plotly::plotlyOutput("gene_expr_plot"),
          shiny::br(),
          shiny::h4("Gene Expression Values:"),
          DT::dataTableOutput("gene_expr_table")
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "GSEexplore"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
