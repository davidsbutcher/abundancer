
# options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)

# Functions ---------------------------------------------------------------


# Server ------------------------------------------------------------------

shinyServer(
   function(input, output, session) {

      # Initial setup -----------------------------------------------------------


      # Hide panels which are only shown conditionally



      # Disable buttons which are selectively enabled later



      # Establish params to use for shinyFiles input (local only)

      volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())

      shinyFileChoose(
         input = input,
         "fileUpload",
         session = session,
         roots = volumes
      )

      # Reactives ---------------------------------------------------------------

      # Create reactive expressions for csv path and name

      file_path <-
         reactive(
            {
               input$fileInputUpload$datapath
            }
         )

      file_name <-
         reactive(
            {
               input$fileInputUpload$name
            }
         )

      file_extension <-
         reactive(
            {
               validate(
                  need(
                     file_path(), ""
                  )
               )

               fs::path_ext(
                  file_path()
               )
            }
         )

      peaklist_reactive <-
         reactive(
            {
               validate(
                  need(file_path(), "")
               )

               if (file_extension() == "csv") {
                  readr::read_csv(
                     file_path()
                  )
               } else if (file_extension() == "xlsx") {
                  readxl::read_xlsx(
                     file_path()
                  )
               }

            }
         )

      chemical_formula <-
         reactive(
            {
               validate(
                  need(input$scoremat_sequence, "Protein sequence not provided")
               )

               # Make chemical formula ---------------------------------------------------

               chemform <-
                  OrgMassSpecR::ConvertPeptide(input$scoremat_sequence) %>%
                  unlist() %>%
                  purrr::imap(~paste0(.y, .x, collapse = '')) %>%
                  paste0(collapse = '') %>%
                  enviPat::mergeform(
                     paste0("H", input$scoremat_charge)
                  )

               # Remove C0|H0|N0|O0|P0|S0 from formula to prevent errors

               if (
                  stringr::str_detect(chemform, "C0|H0|N0|O0|P0|S0") == TRUE
               ) {
                  chemform <-
                     stringr::str_remove_all(chemform, "C0|H0|N0|O0|P0|S0")
               }

               return(chemform)
            }
         )

      is_valid <-
         reactive(
            {
               validate(
                  need(file_path(), "")
               )

               truth_vec <- vector(mode = "logical", length = 2)

               truth_vec[[1]] <-
                  {
                     if (file_extension() == "csv") {
                        TRUE
                     } else if (file_extension() == "xlsx") {
                        TRUE
                     } else {
                        FALSE
                     }
                  }

               truth_vec[[2]] <-
                  {
                     if (!is.null(peaklist_reactive())) {
                        peaklist_reactive() %>%
                           dplyr::summarise_all(class) %>%
                           tidyr::gather(variable, class) %>%
                           dplyr::pull(class) %>%
                           {. == "numeric"} %>%
                           all()
                     } else {
                        FALSE
                     }
                  }

               all(truth_vec)

            }
         )

      score_mat_reac <-
         reactive(
            {
               progressr::withProgressShiny(
                  score_matrix_out <-
                     calculate_score_matrix(
                        mz = dplyr::pull(peaklist_reactive(), 1),
                        intensity = dplyr::pull(peaklist_reactive(), 2),
                        sequence = trimws(input$scoremat_sequence),
                        charge = input$scoremat_charge,
                        start12C = input$scoremat_12C,
                        start14N = input$scoremat_14N,
                        abundStep =
                           as.double(input$scoremat_abundStep),
                        isoAbundCutoff =
                           as.double(input$scoremat_isoAbundCutoff),
                        SNR = input$peakpick_SNR,
                        compFunc = "dotproduct",
                        binSize = input$scoremat_binSize
                     )
               )

               score_matrix_out
            }
         )


      # Listeners ---------------------------------------------------------------

      listener_upload <-
         reactive(
            {
               validate(
                  need(file_path(), ""),
                  need(length(file_path()) > 0, "")
               )

               list(
                  input$fileInputUpload
               )
            }
         )


      listener_peaklist <-
         reactive(
            {
               list(
                  input$peaklist
               )
            }
         )

      # Observers ---------------------------------------------------------------

      observeEvent(
         listener_upload(),
         {
            if (is_valid() == FALSE) {
               showModal(
                  modalDialog(
                     title = "Invalid file",
                     "You must upload a CSV or XLSX with numeric columns"
                  )
               )
            }
         }
      )

      observeEvent(
         listener_upload(),
         {
            # Clear existing plot and html output on upload

            output$output_plot_MS <- NULL
            output$output_plot_scoremat <- NULL
            output$output_html_MS <- NULL

            output$isValid <-
               renderText(
                  {
                     is_valid()
                  }
               )

            output$peaklist_table <-
               renderTable(
                  {
                     validate(
                        need(
                           is_valid(), ""
                        )
                     )

                     peaklist_reactive() %>%
                        dplyr::slice_head(n = 3)
                  }
               )

            output$output_plot_MS <-
               renderPlot(
                  {
                     validate(
                        need(
                           is_valid(), ""
                        )
                     )

                     make_spectrum_peakpick(
                        mz = dplyr::pull(peaklist_reactive(), 1),
                        intensity = dplyr::pull(peaklist_reactive(), 2),
                        SNR = input$peakpick_SNR,
                        method = "MAD",
                        refineMz = "kNeighbors",
                        k = input$peakpick_k
                     )
                  }
               )

            output$output_html_MS <-
               renderUI(
                  {
                     p(
                        strong("Filename: "), file_name(),
                        br(), br(),
                        strong("Chemical formula: "), chemical_formula(),
                        br(), br(),
                        strong("Molecular mass: "),
                        OrgMassSpecR::ConvertPeptide(input$scoremat_sequence) %>%
                           OrgMassSpecR::MolecularWeight()
                     )
                  }
               )

         }
      )

      observeEvent(
         input$abundancer_start,
         {

            validate(
               need(
                  is_valid(), ""
               )
            )

            # Clear existing output on start

            output$output_plot_MS <- NULL
            output$output_plot_scoremat <- NULL
            output$output_html_MS <- NULL

            output$output_plot_scoremat <-
               renderPlot(
                  {
                     generate_score_heatmap(
                        score_matrix = score_mat_reac(),
                        fillOption = input$heatmap_fill
                     )
                  }
               )

            output$output_plot_MS <-
               renderPlot(
                  {
                     validate(
                        need(
                           is_valid(), ""
                        )
                     )

                     make_spectrum_combined(
                        mz = dplyr::pull(peaklist_reactive(), 1),
                        intensity = dplyr::pull(peaklist_reactive(), 2),
                        SNR = input$peakpick_SNR,
                        method = "MAD",
                        refineMz = "kNeighbors",
                        k = input$peakpick_k,
                        score_matrix = score_mat_reac(),
                        sequence = trimws(isolate(input$scoremat_sequence)),
                        charge = isolate(input$scoremat_charge),
                        isoAbundCutoff = 1,
                        hClust_height = 0.005
                     )

                  }
               )


            output$output_html_MS <-
               renderUI(
                  {
                     score_matrix_melt <-
                        reshape2::melt(score_mat_reac()) %>%
                        tibble::as_tibble()

                     names(score_matrix_melt) <-
                        c("14N Abundance", "12C Abundance", "Cosine\nSimilarity")

                     optimal_12c <-
                        score_matrix_melt %>%
                        dplyr::arrange(desc(`Cosine\nSimilarity`)) %>%
                        dplyr::pull(`12C Abundance`) %>%
                        dplyr::first()

                     optimal_14n <-
                        score_matrix_melt %>%
                        dplyr::arrange(desc(`Cosine\nSimilarity`)) %>%
                        dplyr::pull(`14N Abundance`) %>%
                        dplyr::first()

                     p(
                        strong("Filename: "), file_name(),
                        br(), br(),
                        strong("Chemical formula: "), chemical_formula(),
                        br(), br(),
                        strong("Molecular mass: "),
                        OrgMassSpecR::ConvertPeptide(input$scoremat_sequence) %>%
                           OrgMassSpecR::MolecularWeight(),
                        br(), br(),
                        strong("Max cosine similarity: "), round(max(score_mat_reac()), digits = 3),
                        br(), br(),
                        strong("14N abundance at max cos. sim.: "), optimal_14n,
                        br(), br(),
                        strong("12C abundance at max cos. sim.: "), optimal_12c,
                        br()
                     )
                  }
               )

         }

      )

      observeEvent(
         input$fileUpload,
         {
            output$output_test <-
               renderText(
                  if (length(file_path()) == 0) {
                     "FILE PATH HAS LENGTH ZERO"
                  } else {
                     file_path()
                  }
               )
         }
      )

      observeEvent(
         input$testButton1,
         {
            updateNumericInput(
               session = session,
               inputId = "scoremat_charge",
               value = 13
            )
         }
      )

      # Plot expressions --------------------------------------------------------

      # Plot download handlers -------------------------------------------------------

      output$downloadPDF <-
         downloadHandler(
            filename = glue::glue(
               "{format(Sys.time(), '%Y%m%d_%H%M%S')}_abundance_plot.pdf"
            ),
            content = function(file) {
               cairo_pdf(
                  filename = file,
                  width = 8,
                  height = 5,
                  bg = "transparent"
               )
               print(
                  generate_score_heatmap(
                     score_matrix = score_mat_reac(),
                     fillOption = input$heatmap_fill
                  )
               )
               dev.off()
            }
         )

      output$downloadPNG <-
         downloadHandler(
            filename = glue::glue(
               "{format(Sys.time(), '%Y%m%d_%H%M%S')}_abundance_plot.png"
            ),
            content = function(file) {
               png(
                  file = file,
                  width = 8,
                  height = 5,
                  units = "in",
                  res = 300,
                  bg = "transparent"
               )
               print(
                  generate_score_heatmap(
                     score_matrix = score_mat_reac(),
                     fillOption = input$heatmap_fill
                  )
               )
               dev.off()
            }
         )

   }
)
