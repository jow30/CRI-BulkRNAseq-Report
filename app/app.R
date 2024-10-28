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
           textInput("title", "Report Title:", placeholder = "Enter the title of the report", width = "100%"),
           textInput("author", "Report Author:", placeholder = "Enter the author of the report", width = "100%"),
           textAreaInput("intro", "Project Introduction and Experimental Design:",
                         placeholder = "Write the project introduction and the experimental design here.",
                         width = "100%", height = "200px"),
           textAreaInput("nf_notes", "A Breif Description of the Executed Pipeline:", " - We used the nf-core/rnaseq v3.16.0 pipeline for pre-processing of raw reads. \n - We used the STAR->Salmon route for read alignment and quantification. \n - We used the GRCm39 reference genome for read mapping and Gencode vM27 for gene annotation.",
                         width = "100%", height = "200px"),
           helpText("The default drop-down menu include all HTML files in your working directory. If your multiQC report file is stored somewhere else on Randi, please enter the full path to the file into the blank below."),
           selectizeInput("multiqc_path", "MultiQC report to show:", choices = list.files(path = Sys.getenv("PWD"), pattern = ".html", full.names = TRUE, recursive = FALSE), options = list(create = TRUE), width = "100%"),
           textAreaInput("multiqc_notes", "Notable Facts in the MultiQC Report:", "Key points to check in a MultiQC report:\n - Sequencing Depth: If total reads are below 25 million, inform the client as this may impact downstream analysis.\n - rRNA Contamination: Notify the client if more than 10% of reads map to rRNA, as this suggests issues with sample quality.\n - Alignment Rate: If less than 40% of reads align, flag this as a potential problem with the reference or sequencing quality.\n - Per Base Sequence Content: Significant fluctuations after the first few bases indicate potential bias, which should be communicated.\n - Per Sequence GC Content: Deviations from expected GC content patterns may suggest contamination or sequencing bias.\n - Adapter Content: High levels of adapters indicate incomplete trimming, which may require additional preprocessing.\n - Gender-Related Effects: Check the X and Y chromosome counts to ensure no unexpected imbalances due to sample mix-up or biological effects.\n - Gene Coverage Profile: A good profile should resemble a smooth “bridge.” Large fluctuations in coverage suggest technical or sample issues that may need discussion.\n - Read Distribution: For RNAseq data, a high proportion of reads mapping to CDS is expected. Unusual distributions should be reviewed.\n - Inner Distance: An inner distance less than 0 suggests that the fragment may have been sequenced twice by R1 and R2.\n - Strand-Specific Libraries: In the Infer Experiment results, confirm that the library is strand-specific (forward for sense, reverse for antisense).",
                         width = "100%", height = "200px"),
           helpText("The default drop-down menu include all existing folders in your working directory. If your would like to create a new directory to save the results, please enter the full path to the new directory into the blank below. It's recommended to create a new directory for each analysis to avoid overwriting the previous results."),
           selectizeInput("report_out_dir", "Output Directory for DE Analysis:", choices = c(Sys.getenv("PWD"), grep("\\/\\.", list.dirs(path = Sys.getenv("PWD"), full.names = TRUE, recursive = FALSE), invert = TRUE, value = TRUE)), options = list(create = TRUE), width = "100%"),
           selectInput("nextflow_out_dir", "Nextflow Output Directory:",
                       choices = grep("\\/\\.", list.dirs(path = Sys.getenv("PWD"), full.names = TRUE, recursive = FALSE), invert = TRUE, value = TRUE), width = "100%"),
           helpText("The default drop-down menu include all GTF files in /gpfs/data/referenceFiles/. If your GTF file is stored somewhere else on Randi, please enter the full path to the GTF file into the blank below."),
           selectizeInput("gtf_file", "GTF File Used in the nf-core/rnaseq Run:", choices = list.files(path = "/gpfs/data/referenceFiles", pattern = "\\.gtf$", full.names = TRUE, recursive = TRUE), options = list(create = TRUE), width = "100%"),
           helpText("The protein-coding gene filtering is based on the gene_type attribute at the 9th column of the GTF file. So, if you did not use Gencode references, the GTF format could be incompatible with this program and the protein-coding gene filtering won't work properly."),
           checkboxInput("protein_coding", "Retain Protein-Coding Genes Only", TRUE),
           helpText("The default drop-down menu include all TXT files in your working directory. If your metadata file is stored somewhere else on Randi, please enter the full path to the file into the blank below."),
           selectizeInput("metadata_file", "Metadata File:", choices = list.files(path = Sys.getenv("PWD"), pattern = ".txt", full.names = TRUE, recursive = FALSE), options = list(create = TRUE), width = "100%"),
           uiOutput("metadata_error"),
           selectizeInput("sample_remove", "Samples to Remove:", choices = NULL, multiple = TRUE, options = list(create = FALSE), width = "100%"),
           selectizeInput("color_by", "Sample Attribute for PCA Color:", choices = NULL, width = "100%"),
           selectizeInput("shape_by", "Sample Attribute for PCA Shape:", choices = NULL, selected = NULL, options = list(placeholder = "Leave empty if no secondary attribute to examine."), width = "100%"),
           textInput("top_var", "Top N Most Variable Genes for PCA:", value = 500, placeholder = "Leave empty to use all genes", width = "100%"),
           uiOutput("top_var_error"),
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
           selectInput("species", "Species:", choices = c("human", "mouse", "rat"), selected = "human", width = "100%"),
           helpText("Select the option below to merge up-/down-regulated DEGs into a single list for ORA. Otherwise, they will be analyzed separately. "),
           checkboxInput("ora_all", "Perform ORA with All DEGs", TRUE),
           checkboxInput("ora_go", "Perform ORA with GO Terms", TRUE),
           checkboxInput("ora_kegg", "Perform ORA with KEGG Pathways", TRUE),
           checkboxInput("ora_reactome", "Perform ORA with Reactome Pathways", FALSE),
           checkboxInput("ora_msigdb", "Perform ORA with MSigDB Gene Sets", FALSE),
           selectInput("ora_msigdbr_category", "MSigDB Category for ORA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_cat), selected = "H", width = "100%"),
           selectInput("ora_msigdbr_subcategory", "MSigDB Subcategory for ORA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_subcat), width = "100%"),
           checkboxInput("gsea_msigdb", "Perform GSEA with MSigDB Gene Sets:", FALSE),
           numericInput("gsea_fdr_thres", "FDR Cutoff for GSEA Results:", value = 0.05, width = "100%"),
           selectInput("gsea_msigdbr_category", "MSigDB Category for GSEA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_cat), selected = "C2", width = "100%"),
           selectInput("gsea_msigdbr_subcategory", "MSigDB Subcategory for GSEA:", 
                       choices = unique(msigdbr::msigdbr_collections()$gs_subcat), selected = "CP:KEGG", width = "100%"),
           actionButton("submit_btn", "Submit", class = "btn btn-success btn-lg")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Initialize submit button as disabled
  disable("submit_btn")
  
  # load parameters from a RData file
  observe({
    if (file.exists(file.path(Sys.getenv("PWD"), "params.RData"))) {
      params <- readRDS(file.path(Sys.getenv("PWD"), "params.RData"))
      updateTextInput(session, "title", value = params$title)
      updateTextInput(session, "author", value = params$author)
      updateTextAreaInput(session, "intro", value = params$intro)
      updateTextAreaInput(session, "nf_notes", value = params$nf_notes)
      updateSelectizeInput(session, "multiqc_path", selected = params$multiqc_path)
      updateTextAreaInput(session, "multiqc_notes", value = params$multiqc_notes)
      updateSelectizeInput(session, "report_out_dir", selected = params$report_out_dir)
      updateSelectInput(session, "nextflow_out_dir", selected = params$nextflow_out_dir)
      updateSelectizeInput(session, "gtf_file", selected = params$gtf_file)
      updateCheckboxInput(session, "protein_coding", value = params$protein_coding)
      updateSelectizeInput(session, "metadata_file", selected = params$metadata_file)
      updateSelectizeInput(session, "sample_remove", selected = params$sample_remove)
      updateSelectizeInput(session, "color_by", selected = params$color_by)
      updateSelectizeInput(session, "shape_by", selected = params$shape_by)
      updateTextInput(session, "top_var", value = params$top_var)
      updateSelectInput(session, "de_method", selected = params$de_method)
      updateCheckboxInput(session, "batch_correction", value = params$batch_correction)
      updateNumericInput(session, "fdr_thres", value = params$fdr_thres)
      updateNumericInput(session, "fc_thres", value = params$fc_thres)
      updateNumericInput(session, "ora_fdr_thres", value = params$ora_fdr_thres)
      updateSelectInput(session, "species", selected = params$species)
      updateCheckboxInput(session, "ora_all", value = params$ora_all)
      updateCheckboxInput(session, "ora_go", value = params$ora_go)
      updateCheckboxInput(session, "ora_kegg", value = params$ora_kegg)
      updateCheckboxInput(session, "ora_reactome", value = params$ora_reactome)
      updateCheckboxInput(session, "ora_msigdb", value = params$ora_msigdb)
      updateSelectInput(session, "ora_msigdbr_category", selected = params$ora_msigdbr_category)
      updateSelectInput(session, "ora_msigdbr_subcategory", selected = params$ora_msigdbr_subcategory)
      updateCheckboxInput(session, "gsea_msigdb", value = params$gsea_msigdb)
      updateNumericInput(session, "gsea_fdr_thres", value = params$gsea_fdr_thres)
      updateSelectInput(session, "gsea_msigdbr_category", selected = params$gsea_msigdbr_category)
      updateSelectInput(session, "gsea_msigdbr_subcategory", selected = params$gsea_msigdbr_subcategory)
    }
  })

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

  # Reactive to validate metadata
  validate_metadata <- reactive({
    # Ensure a file is uploaded
    if (is.null(input$metadata_file) || input$metadata_file == "") {
      output$metadata_error <- renderUI({
        div(style = "color: red;", "Please select a metadata file.")
      })
      return(FALSE)  # No file, return invalid
    }
    
    # Try reading the metadata file, handle errors
    metadata <- tryCatch({
      read.delim(input$metadata_file)
    }, error = function(e) {
      output$metadata_error <- renderUI({
        div(style = "color: red;", "Error reading the file. Please check if the file format is correct. Only TXT format is acceptable.")
      })
      return(FALSE)  # Error reading file, return invalid
    })

    # Check if 'group' column exists
    if (!"group" %in% colnames(metadata)) {
      output$metadata_error <- renderUI({
        div(style = "color: red;", "The metadata file does not contain the required 'group' column.")
      })
      return(FALSE)  # Missing group column, return invalid
    }
    
    if (!"batch" %in% colnames(metadata)) {
      updateCheckboxInput(session, "batch_correction", value = FALSE)  # Uncheck batch correction
      disable("batch_correction")  # Disable batch correction checkbox
    } else {
      enable("batch_correction")  # Enable batch correction checkbox if "batch" column exists
    }

    # If everything is valid, remove any error message
    output$metadata_error <- renderUI({ NULL })
    
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

  # Dynamically render an error message if 'top_var' is invalid
  output$top_var_error <- renderUI({
    if (!validate_top_var()) {
      div(style = "color: red;", "Please enter an integer greater than 0 or leave it empty.")
    } else {
      NULL  # No error message if input is valid
    }
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

  observeEvent(input$submit_btn, {
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
    comp_pair <- comp_pair[!sapply(comp_pair, function(x) any(x == ""))]

    shape_by_val <- if (input$shape_by == "") NULL else input$shape_by
    top_var_val <- if(input$top_var == "") NULL else as.numeric(input$top_var)
    intro_val <- gsub("\n", "<br>", input$intro)
    nf_notes_val <- gsub("\n", "<br>", input$nf_notes)
    multiqc_notes_val <- gsub("\n", "<br>", input$multiqc_notes)
    ora_msigdbr_subcategory_val <- if (input$ora_msigdbr_subcategory == "") NULL else input$ora_msigdbr_subcategory
    gsea_msigdbr_subcategory_val <- if (input$gsea_msigdbr_subcategory == "") NULL else input$gsea_msigdbr_subcategory

    if(!dir.exists(input$report_out_dir)) dir.create(input$report_out_dir, recursive = TRUE)

    # Collect all parameters into a list
    params <- list(
      title = input$title,
      author = input$author,
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
    
    # save parameters to a RData file
    saveRDS(params, file = file.path(Sys.getenv("PWD"), "params.RData"))

    showModal(modalDialog(
      title = "Generating Report",
      div(
        p("Generating report, please wait...")
      ),
      footer = NULL,  # No footer while processing
      easyClose = FALSE  # Prevent closing the modal until the process finishes
    ))
    
    # Use tryCatch to handle any potential errors
    tryCatch({
      rmarkdown::render(input = "report.Rmd", output_file = file.path(input$report_out_dir, "report.html"), params = params, envir = new.env(), quiet = FALSE)
      
      # After rendering is complete, close the modal and show success
      removeModal()
      showModal(modalDialog(
        title = "Analysis Complete",
        tagList(
          p(paste0("Please find the following results at ", input$report_out_dir, ":")), 
          p(" - report.html"),
          p(" - report.RData"),
          p(" - 2.Pre-processing_of_raw_reads/"),
          p(" - 3.Differential_expression_analysis/"),
          p(" - 4.Functional_analysis/")
        ),
        easyClose = TRUE,
        size = "l",
        footer = modalButton("OK")
      ))
      
      enable("submit_btn")
    }, error = function(e) {
      # Handle errors and display an error message in a new modal
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tagList(
          p("The analysis failed due to an error."),
          p("Error message:"),
          pre(conditionMessage(e)),
          p("Please check the error logs (in *.err file) for more details.")
        ),
        easyClose = TRUE,
        size = "l",
        footer = modalButton("OK")
      ))
      
      enable("submit_btn")
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
