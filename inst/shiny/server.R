
# options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)

# Functions ---------------------------------------------------------------


# Server ------------------------------------------------------------------

shinyServer(
   function(input, output, session) {

      # Initial setup -----------------------------------------------------------


      # Hide panels which are only shown conditionally


      # Disable buttons which are selectively enabled later


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

            progressr::withProgressShiny(
               score_matrix_out <-
                  calculate_score_matrix_dual(
                     mz = dplyr::pull(peaklist_reactive(), 1),
                     intensity = dplyr::pull(peaklist_reactive(), 2),
                     sequence = trimws(input$scoremat_sequence),
                     PTMformula = input$scoremat_PTM,
                     charge = input$scoremat_charge,
                     start12C = input$scoremat_12C,
                     start14N = input$scoremat_14N,
                     SNR = input$peakpick_SNR,
                     method = "MAD",
                     refineMz = "kNeighbors",
                     k = input$peakpick_k,
                     compFunc = input$scoremat_compfunc,
                     binSize = input$scoremat_binSize,
                     resolvingPower = input$scoremat_resolvingPower
                  )
            )

            output$output_plot_scoremat1 <-
               renderPlot(
                  {
                     generate_score_heatmap(
                        score_matrix = score_matrix_out[[1]],
                        fillOption = input$heatmap_fill,
                        scaleFillLimits =
                           c(
                              input$heatmap_scale_coarse_start,
                              input$heatmap_scale_coarse_end
                           ),
                        compFunc = isolate(input$scoremat_compfunc)
                     )
                  }
               )

            output$output_plot_scoremat2 <-
               renderPlot(
                  {
                     generate_score_heatmap(
                        score_matrix = score_matrix_out[[2]],
                        fillOption = input$heatmap_fill,
                        scaleFillLimits =
                           c(
                              input$heatmap_scale_fine_start,
                              input$heatmap_scale_fine_end
                           ),
                        compFunc = isolate(input$scoremat_compfunc)
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
                        score_matrix = score_matrix_out[[2]],
                        sequence = trimws(isolate(input$scoremat_sequence)),
                        charge = isolate(input$scoremat_charge),
                        resolvingPower = isolate(input$scoremat_resolvingPower)
                     )

                  }
               )


            output$output_html_MS <-
               renderUI(
                  {
                     optimal <-
                        get_optimal_abundances(score_matrix_out[[2]])

                     optimal_14n <-
                        optimal[[1]]

                     optimal_12c <-
                        optimal[[2]]

                     p(
                        strong("Filename: "), file_name(),
                        br(), br(),
                        strong("Chemical formula: "), chemical_formula(),
                        br(), br(),
                        strong("Molecular mass: "),
                        OrgMassSpecR::ConvertPeptide(input$scoremat_sequence) %>%
                           OrgMassSpecR::MolecularWeight(),
                        br(), br(),
                        strong("Max score: "),
                        round(max(score_matrix_out[[2]]), digits = 3),
                        br(), br(),
                        strong("14N abundance at max score: "), optimal_14n,
                        br(), br(),
                        strong("12C abundance at max score: "), optimal_12c,
                        br()
                     )
                  }
               )

            # Plot download handlers -------------------------------------------------------

            output$downloadPDFcoarse <-
               downloadHandler(
                  filename = glue::glue(
                     "{format(Sys.time(), '%Y%m%d_%H%M%S')}_abundance_plot_coarse.pdf"
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
                           score_matrix = score_matrix_out[[1]],
                           fillOption = input$heatmap_fill,
                           scaleFillLimits =
                              c(
                                 input$heatmap_scale_coarse_start,
                                 input$heatmap_scale_coarse_end
                              ),
                           compFunc = isolate(input$scoremat_compfunc)
                        )
                     )
                     dev.off()
                  }
               )

            output$downloadPDFfine <-
               downloadHandler(
                  filename = glue::glue(
                     "{format(Sys.time(), '%Y%m%d_%H%M%S')}_abundance_plot_fine.pdf"
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
                           score_matrix = score_matrix_out[[2]],
                           fillOption = input$heatmap_fill,
                           scaleFillLimits =
                              c(
                                 input$heatmap_scale_fine_start,
                                 input$heatmap_scale_fine_end
                              ),
                           compFunc = isolate(input$scoremat_compfunc)
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
                           score_matrix = score_matrix_out[[2]],
                           fillOption = input$heatmap_fill,
                           scaleFillLimits =
                              c(
                                 input$heatmap_scale_fine_start,
                                 input$heatmap_scale_fine_end
                              ),
                           compFunc = isolate(input$scoremat_compfunc)
                        )
                     )
                     dev.off()
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


   }
)
