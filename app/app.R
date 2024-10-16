library(shiny)
library(shinyjs)
library(shinycssloaders)
library(rmarkdown)
library(tidyverse)

# Define UI for application
ui <- fluidPage(
  useShinyjs(),
  
  titlePanel("RNAseq Differential Expression Analysis"),
  p("Fill out this form to run DE analysis downstream of nf-core/rnaseq pipeline."),
  hr(),
  
  # Display all inputs directly in the main panel
  fluidRow(
    column(12,
           textAreaInput("intro", "Project Introduction and Experimental Design:",
                         placeholder = "Write the project introduction and the experimental design here.",
                         width = "100%", height = "200px"),
           textAreaInput("nf_notes", "A Breif Description of the Executed Pipeline:", " - We used the nf-core/rnaseq v3.16.0 pipeline for pre-processing of raw reads. \n - We used the STAR->Salmon route for read alignment and quantification. \n - We used the GRCm39 reference genome for read mapping and Gencode vM27 for gene annotation.",
                         width = "100%", height = "200px"),
           selectInput("multiqc_path", "MultiQC report to show:",
                       choices = list.files(path = Sys.getenv("PWD"), pattern = ".html", full.names = TRUE, recursive = TRUE), width = "100%"),
           textAreaInput("multiqc_notes", "Notable Facts in the MultiQC Report:",
                         placeholder = "Write something here if any notable facts are found in the multiQC report.",
                         width = "100%", height = "200px"),
           selectInput("report_out_dir", "Output Directory for DE Analysis:",
                       choices = c(Sys.getenv("PWD"), grep("\\/\\.", list.dirs(path = Sys.getenv("PWD"), full.names = TRUE, recursive = FALSE), invert = TRUE, value = TRUE)), width = "100%"),
           selectInput("nextflow_out_dir", "Nextflow Output Directory:",
                       choices = c(Sys.getenv("PWD"), grep("\\/\\.", list.dirs(path = Sys.getenv("PWD"), full.names = TRUE, recursive = FALSE), invert = TRUE, value = TRUE)), 
                       selected = grep("\\/\\.", list.dirs(path = Sys.getenv("PWD"), full.names = TRUE, recursive = FALSE), invert = TRUE, value = TRUE)[1], width = "100%"),
           helpText("By default, the app lists all GTF files in /gpfs/data/referenceFiles/. If your GTF file is stored somewhere else on Randi, please enter the full path to the GTF file into the blank below."),
           selectizeInput("gtf_file", "GTF File Used in the nf-core/rnaseq Run:", choices = list.files(path = "/gpfs/data/referenceFiles", pattern = "\\.gtf$", full.names = TRUE, recursive = TRUE), options = list(create = TRUE), width = "100%"),
           helpText("The protein-coding gene filtering is based on the gene_type attribute at the 9th column of the GTF file. So, if you did not use Gencode references, the GTF format could be incompatible with this program and the protein-coding gene filtering won't work properly."),
           checkboxInput("protein_coding", "Retain Protein-Coding Genes Only", TRUE),
           selectInput("metadata_file", "Metadata File:",
                       choices = list.files(path = Sys.getenv("PWD"), pattern = ".txt", full.names = TRUE, recursive = TRUE), width = "100%"),
           selectizeInput("sample_remove", "Samples to Remove:", choices = NULL, multiple = TRUE, options = list(create = FALSE), width = "100%"),
           selectizeInput("color_by", "Sample Attribute for PCA Color:", choices = NULL, width = "100%"),
           selectizeInput("shape_by", "Sample Attribute for PCA Shape:", choices = NULL, selected = NULL, options = list(placeholder = "Leave empty if no secondary attribute to examine."), width = "100%"),
           textInput("top_var", "Top N Most Variable Genes for PCA:", value = NULL, placeholder = "Leave empty to use all genes", width = "100%"),
           selectInput("de_method", "DE Test Method:", choices = c("DESeq2", "limma-voom"), selected = "DESeq2", width = "100%"),
           helpText("A column named \"batch\" in the metadata table is required for batch effect correction."),
           checkboxInput("batch_correction", "Correct for Batch Effect", FALSE),
           numericInput("fdr_thres", "FDR Cutoff for DE Test Results:", value = 0.05, width = "100%"),
           numericInput("fc_thres", "Fold Change Cutoff for DE Test Results:", value = 1.5, width = "100%"),
           actionButton("add_comparison_btn", "Add Group Comparison Pairs", class = "btn btn-primary", style = "margin-bottom: 20px;"),
           div(id = "comparison_container"),  # Container for dynamically added comparison pairs
           uiOutput("dynamic_ui"),
           hr(),
           numericInput("ora_fdr_thres", "FDR Cutoff for ORA Results:", value = 0.05, width = "100%"),
           selectInput("species", "Species:", choices = c("human", "mouse"), selected = "human", width = "100%"),
           helpText("Select the \"Perform ORA with All DEGs\" option below to merge up-/down-regulated DEGs into a single list for ORA. Otherwise, they will be analyzed separately. "),
           checkboxInput("ora_all", "Perform ORA with All DEGs", TRUE),
           checkboxInput("ora_go", "Perform ORA with GO Terms", TRUE),
           checkboxInput("ora_kegg", "Perform ORA with KEGG Pathways", TRUE),
           checkboxInput("ora_reactome", "Perform ORA with Reactome Pathways", FALSE),
           checkboxInput("ora_msigdb", "Perform ORA with MSigDB Gene Sets", FALSE),
           selectInput("ora_msigdbr_category", "MSigDB Category for ORA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_cat), selected = "C2", width = "100%"),
           selectInput("ora_msigdbr_subcategory", "MSigDB Subcategory for ORA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_subcat), selected = "CP:KEGG", width = "100%"),
           checkboxInput("gsea_msigdb", "Perform GSEA with MSigDB Gene Sets:", FALSE),
           numericInput("gsea_fdr_thres", "FDR Cutoff for GSEA Results:", value = 0.05, width = "100%"),
           selectInput("gsea_msigdbr_category", "MSigDB Category for GSEA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_cat), selected = "C2", width = "100%"),
           selectInput("gsea_msigdbr_subcategory", "MSigDB Subcategory for GSEA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_subcat), selected = "CP:KEGG", width = "100%"),
           withSpinner(actionButton("submit_btn", "Submit", class = "btn btn-success btn-lg"), type = 8)  # Adding spinner
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Initialize submit button as disabled
  disable("submit_btn")
  
  # Observe the selected MSigDB category and update the subcategory options
  observeEvent(input$ora_msigdbr_category, {
    req(input$ora_msigdbr_category)  # Ensure a category is selected
    
    # Filter the subcategories based on the selected category
    filtered_subcategories <- msigdbr::msigdbr_collections() %>%
      dplyr::filter(gs_cat == input$ora_msigdbr_category) %>%
      dplyr::pull(gs_subcat) %>%
      unique()
    
    # Update the ora_msigdbr_subcategory input with the filtered subcategories
    updateSelectInput(session, "ora_msigdbr_subcategory", choices = filtered_subcategories, selected = "CP:KEGG")
  })
  
  observeEvent(input$gsea_msigdbr_category, {
    req(input$gsea_msigdbr_category)  # Ensure a category is selected
    
    # Filter the subcategories based on the selected category
    filtered_subcategories <- msigdbr::msigdbr_collections() %>%
      dplyr::filter(gs_cat == input$gsea_msigdbr_category) %>%
      dplyr::pull(gs_subcat) %>%
      unique()
    
    # Update the gsea_msigdbr_subcategory input with the filtered subcategories
    updateSelectInput(session, "gsea_msigdbr_subcategory", choices = filtered_subcategories, selected = "CP:KEGG")
  })
  
  # Reactive variable for group names
  group_names <- reactiveVal(NULL)
  
  # Reactive expression to check if metadata is valid (i.e., contains the 'group' column)
  validate_metadata <- reactive({
    req(input$metadata_file)  # Ensure a file is uploaded
    
    # Read the metadata file
    metadata <- read.delim(input$metadata_file)
    
    # Check if 'group' column exists
    if (!"group" %in% colnames(metadata)) {
      showModal(modalDialog(
        title = "Warning",
        "The metadata file does not contain the required 'group' column. Please check your file.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(FALSE)  # Metadata is invalid
    }
    
    # Assuming the first column contains the sample names and the `group` column contains group information
    sample_names <- metadata[[1]]
    group_names(metadata$group)
    
    # Update the selectizeInput for sample_remove with the sample names
    updateSelectizeInput(session, "sample_remove", choices = sample_names, server = TRUE)
    
    updateSelectizeInput(session, "color_by", choices = colnames(metadata)[-1], server = TRUE)
    updateSelectizeInput(session, "shape_by", choices = c("", colnames(metadata)[-1]), selected = "", server = TRUE)
    
    return(TRUE)  # Metadata is valid
  })
  
  # Function to check if input is a valid integer > 0 or empty
  validate_top_var <- reactive({
    if (is.null(input$top_var) || input$top_var == "") return(TRUE)  # Allow empty input
    top_var_val <- as.numeric(input$top_var)
    return(!is.na(top_var_val) && top_var_val %% 1 == 0 && top_var_val > 0)
  })
  
  # Show a warning message if the input is invalid
  observe({
    if (!validate_top_var()) showNotification("Please enter an integer greater than 0 for the 'top N most variable genes' argument or leave it empty.", type = "error", duration = 6)
  })
  
  # You can conditionally enable/disable a submit button based on validation
  observe({
    if (validate_top_var() && validate_metadata()) {
      enable("submit_btn")
    } else {
      disable("submit_btn")
    }
  })
  
  # Keep track of how many group pairs have been added
  counter <- reactiveVal(0)
  
  # Initialize an empty reactive list to store active comparison pairs
  active_comparisons <- reactiveVal(list())
  
  # Dynamically add UI for a new comparison pair
  observeEvent(input$add_comparison_btn, {
    # Increment the counter
    count <- counter() + 1
    counter(count)
    
    # Insert new input fields for the group pair
    insertUI(
      selector = "#comparison_container",
      where = "beforeEnd",
      ui = tags$div(
        id = paste0("comparison_", count),  # Unique ID for each pair
        fluidRow(
          column(5, selectizeInput(paste0("group_test_", count), "Group for Test:", choices = group_names(), multiple = FALSE)),
          column(5, selectizeInput(paste0("group_base_", count), "Group as Base:", choices = group_names(), multiple = FALSE)),
          column(2, actionButton(paste0("remove_btn_", count), "Remove", class = "btn-danger"))
        )
      )
    )
    
    # Add the newly created pair to the list of active comparisons
    current_comparisons <- active_comparisons()
    current_comparisons[[count]] <- TRUE  # Mark this comparison as active
    active_comparisons(current_comparisons)
    
    # Remove the group comparison when the remove button is clicked
    observeEvent(input[[paste0("remove_btn_", count)]], {
      removeUI(selector = paste0("#comparison_", count))
      
      # Mark this pair as inactive (set to FALSE)
      current_comparisons <- active_comparisons()
      current_comparisons[[count]] <- FALSE
      active_comparisons(current_comparisons)
    })
  })
  
  # Function to handle errors and show them in a modal
  handle_error <- function(error_message) {
    showModal(modalDialog(
      title = "Error Occurred",
      tagList(
        p("The analysis failed due to an error."),
        p("Error message:"),
        pre(error_message),
        p("Please check the error logs in the 'bulkRNAseq*.err' file for more details."),
        p("You can refresh the page and try again.")
      ),
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  }

  observeEvent(input$submit_btn, {
    shinyjs::show("spinner")  # Show spinner during processing
    disable("submit_btn")  # Disable the submit button to prevent multiple clicks
    
    count <- counter()
    comp_pair <- list()
    
    # Loop through the added pairs and extract the input values
    for (i in seq_len(count)) {
      if (active_comparisons()[[i]]) {
        test_input <- input[[paste0("group_test_", i)]]
        base_input <- input[[paste0("group_base_", i)]]
        if (!is.null(test_input) && !is.null(base_input)) {
          comp_pair <- append(comp_pair, list(c(test_input, base_input)))
        }
      }
    }
    
    shape_by_val <- if (input$shape_by == "") NULL else input$shape_by
    top_var_val <- if(input$top_var == "") NULL else as.numeric(input$top_var)
    intro_val <- gsub("\n", "<br>", input$intro)
    nf_notes_val <- gsub("\n", "<br>", input$nf_notes)
    multiqc_notes_val <- gsub("\n", "<br>", input$multiqc_notes)
    ora_msigdbr_subcategory_val <- if (input$ora_msigdbr_subcategory == "") NULL else input$ora_msigdbr_subcategory
    gsea_msigdbr_subcategory_val <- if (input$gsea_msigdbr_subcategory == "") NULL else input$gsea_msigdbr_subcategory
    
    # Collect all parameters into a list
    params <- list(
      intro = input$intro,
      nf_notes = input$nf_notes,
      multiqc_path = input$multiqc_path,
      multiqc_notes = input$multiqc_notes,
      report_out_dir = input$report_out_dir,
      nextflow_out_dir = input$nextflow_out_dir,
      gtf_file = input$gtf_file,
      protein_coding = input$protein_coding,
      metadata_file = input$metadata_file,
      sample_remove = input$sample_remove,
      color_by = input$color_by,
      shape_by = shape_by_val,
      top_var = top_var_val,
      de_method = input$de_method,
      batch_correction = input$batch_correction,
      fdr_thres = input$fdr_thres,
      fc_thres = input$fc_thres,
      comp_pair = comp_pair,
      ora_fdr_thres = input$ora_fdr_thres,
      species = input$species,
      ora_all = input$ora_all,
      ora_go = input$ora_go,
      ora_kegg = input$ora_kegg,
      ora_reactome = input$ora_reactome,
      ora_msigdb = input$ora_msigdb,
      ora_msigdbr_category = input$ora_msigdbr_category,
      ora_msigdbr_subcategory = ora_msigdbr_subcategory_val,
      gsea_msigdb = input$gsea_msigdb,
      gsea_fdr_thres = input$gsea_fdr_thres,
      gsea_msigdbr_category = input$gsea_msigdbr_category,
      gsea_msigdbr_subcategory = gsea_msigdbr_subcategory_val
    )
    
    # Wrap the rmarkdown::render process in a tryCatch block to handle errors
    tryCatch({
      output_file <- file.path(input$report_out_dir, "report.html")
      
      # Render the report, and if successful, show a success modal
      rmarkdown::render("report.Rmd", output_file = output_file, params = params, envir = new.env())
      
      shinyjs::hide("spinner")  # Hide spinner when done
      showModal(modalDialog(
        title = "Analysis Complete",
        tagList(paste("Your report has been generated at:", output_file)),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      
      enable("submit_btn")  # Re-enable the submit button after success
    }, error = function(e) {
      # Handle the error, capture the message, and show a user-friendly modal
      error_message <- conditionMessage(e)
      handle_error(error_message)
      
      shinyjs::hide("spinner")  # Hide spinner when an error occurs
      enable("submit_btn")  # Re-enable the submit button after an error
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

