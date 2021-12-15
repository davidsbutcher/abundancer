
library(abundancer)
library(magrittr)
library(dplyr)
library(assertthat)
library(glue)
library(shiny)
library(shinyBS)
library(tippy)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(purrr)
library(Peptides)
library(OrgMassSpecR)
library(fs)
library(enviPat)
library(reshape2)
library(progressr)

# options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)


# UI ----------------------------------------------------------------------

shinyUI(

   fixedPage(
      titlePanel("abundancer"),
      theme = "maglab_theme_old.css",
      # theme = bslib::bs_theme(version = 4),

      useShinyjs(),

      setBackgroundColor(
         color = c("#FFFFFF"),
         gradient = "linear",
         direction = "bottom"
      ),

      tags$head(
         tags$style(HTML("hr {border-top: 1px solid #A9A9A9;}")),
         tags$style(
            HTML(
               ".shiny-notification {
                        position:fixed;
                        top: calc(20%);
                        left: calc(40%);
                    }"
            )
         )
      ),

      sidebarLayout(

         # Sidebar -----------------------------------------------------------------

         sidebarPanel(
            tabsetPanel(
               id = "mainpanel",
               type = "pills",
               tabPanel(
                  "Main",
                  hr(),
                  tabsetPanel(
                     id = "guppipanel",
                     type = "tabs",
                     tabPanel(
                        "Peaklist",
                        br(),
                        fileInput(
                           "fileInputUpload",
                           tippy(
                              "Upload peaklist",
                              placement = "right",
                              arrow = TRUE,
                              allowHTML = TRUE,
                              animation = "scale",
                              duration = 250,
                              theme = "light-border",
                              tooltip =
                                 "<span style='font-size:14px;'>Upload a peaklist in CSV or XLSX format. The first column should contain m/z values and the second should contain intensities. The peaklist should contain a full isotopic distribution of a single charge state of a single protein.</span>"
                           ),
                           accept = c(".xlsx", ".csv")
                        ),
                        numericInput(
                           "scoremat_charge",
                           tippy(
                              "Charge",
                              placement = "right",
                              arrow = TRUE,
                              allowHTML = TRUE,
                              animation = "scale",
                              duration = 250,
                              theme = "light-border",
                              tooltip =
                                 "<span style='font-size:14px;'>Charge state corresponding to the isotopic distribution in the peaklist. Only positive charges are supported.</span>"
                           ),
                           value = 1,
                           min = 1,
                           max = 100,
                           step = 1
                        ),
                        textAreaInput(
                           "scoremat_sequence",
                           tippy(
                              "Sequence",
                              placement = "right",
                              arrow = TRUE,
                              allowHTML = TRUE,
                              animation = "scale",
                              duration = 250,
                              theme = "light-border",
                              tooltip =
                                 "<span style='font-size:14px;'>Sequence corresponding to the protein whose isotopic distribution is represented by the peaklist. Only single-letter codes of canonical amino acids are supported.</span>"
                           ),
                           placeholder =
                              "Single letter canonical amino acids only",
                           height = "100px"
                        ),
                        textInput(
                           "scoremat_PTM",
                           tippy(
                              "PTM formula",
                              placement = "right",
                              arrow = TRUE,
                              allowHTML = TRUE,
                              animation = "scale",
                              duration = 250,
                              theme = "light-border",
                              tooltip =
                                 "<span style='font-size:14px;'>Chemical formula (in Hill notation) of post-translational modifications to the protein sequence. The formula should be provided as a delta compared to the protein's chemical formula, i.e. for a methylation the PTM formula is CH3 and not CH4. For multiple PTMs, all formulae should be combined.</span>"
                           ),
                           placeholder = "Leave blank if unmodified"
                        ),
                        br(),
                        br(),
                        actionButton(
                           "abundancer_start",
                           "Analyze abundance"
                        ),
                        br(),
                        br(),
                        downloadButton("downloadPDFcoarse", label = "Coarse Heatmap"),
                        downloadButton("downloadPDFfine", label = "Fine Heatmap")
                     ),
                     tabPanel(
                        "Settings",
                        # h5(em("Hover parameter names for more info")),

                        bsCollapse(

                           bsCollapsePanel(

                              h5(strong("Peak picking parameters")),
                              br(),
                              splitLayout(
                                 numericInput(
                                    "peakpick_SNR",
                                    tippy(
                                       "SNR",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>From MSnbase documentation: \"A local maximum is considered a peak if its intensity is SNR times larger than the estimated noise.\" Noise is estimated using the MSnbase median absolute deviation method.</span>"
                                    ),
                                    value = 25,
                                    min = 0,
                                    max = 100,
                                    step = 1
                                 ),
                                 numericInput(
                                    "peakpick_k",
                                    tippy(
                                       "KNN",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>From MSnbase documentation: \"m/z values and intensities of the 2 * KNN closest signals to the centroid are used in the intensity weighted average calculation.\"</span>"
                                    ),
                                    value = 2,
                                    min = 1,
                                    max = 10,
                                    step = 1
                                 )
                              ),
                              br()
                           ),

                           bsCollapsePanel(

                              h5(strong("Score matrix parameters")),
                              br(),
                              splitLayout(
                                 numericInput(
                                    "scoremat_12C",
                                    tippy(
                                       "<sup>12</sup>C abundance start",
                                       placement = "top-end",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>Starting value for <sup>12</sup>C abundance to use for calculating the coarse score matrix. <sup>12</sup>C abundance will be incremented from this value to 1 by the abundance step.</span>"
                                    ),
                                    value = 0.987,
                                    min = 0.9,
                                    max = 0.999,
                                    step = 0.001
                                 ),
                                 numericInput(
                                    "scoremat_14N",
                                    tippy(
                                       "<sup>14</sup>N abundance start",
                                       placement = "top-end",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>Starting value for <sup>14</sup>N abundance to use for calculating the coarse score matrix. <sup>14</sup>N abundance will be incremented from this value to 1 by the abundance step.</span>"
                                    ),
                                    value = 0.994,
                                    min = 0.9,
                                    max = 0.999,
                                    step = 0.001
                                 )
                              ),
                              splitLayout(
                                 numericInput(
                                    "scoremat_binSize",
                                    tippy(
                                       "Bin size",
                                       placement = "top-end",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>Prior to comparing observed and theoretical spectra, spectra are binned using this bin size value. See MSnbase documentation for more information.</span>"
                                    ),
                                    value = 0.05,
                                    min = 0.01,
                                    max = 0.1,
                                    step = 0.01
                                 ),
                                 numericInput(
                                    "scoremat_resolvingPower",
                                    tippy(
                                       "Resolving power",
                                       placement = "top-end",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>Resolving power to be used in generating theoretical isotopic distributions. See enviPat::envelope documentation for more information.</span>"
                                    ),
                                    value = 300000,
                                    min = 1000,
                                    max = 1000000,
                                    step = 1000
                                 )
                              ),
                              splitLayout(
                                 numericInput(
                                    "scoremat_coarseMult",
                                    tippy(
                                       "Coarse step mult.",
                                       placement = "top-end",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>Value by which coarse step is multiplied to determine range for calculating fine matrix.</span>"
                                    ),
                                    value = 1,
                                    min = 1,
                                    max = 10,
                                    step = 1
                                 ),
                                 selectInput(
                                    "scoremat_compfunc",
                                    tippy(
                                       "Spec. comp. func.",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>Sepctrum comparison function to be used for comparing observed and theoretical isotopic envelopes.</span>"
                                    ),
                                    choices =
                                       list(
                                          "Cosine similarity" = "dotproduct",
                                          "Cor" = "cor",
                                          "ScoreMFA" = "scoremfa",
                                          "ScoreMFA (C++)" = "scoremfacpp"
                                       ),
                                    selected = "Cosine similarity"
                                 )
                              ),
                              # div(
                              #    style="display: inline-block;vertical-align:top; width: 150px;",
                              #    selectInput(
                              #       "scoremat_compfunc",
                              #       tippy(
                              #          "Spectrum comparison function",
                              #          placement = "right",
                              #          arrow = TRUE,
                              #          allowHTML = TRUE,
                              #          animation = "scale",
                              #          duration = 250,
                              #          theme = "light-border",
                              #          tooltip =
                              #             "<span style='font-size:14px;'>Function to be used for comparing observed and theoretical isotopic envelopes.</span>"
                              #       ),
                              #       choices =
                              #          list(
                              #             "Cosine similarity" = "dotproduct",
                              #             "Cor" = "cor",
                              #             "ScoreMFA" = "scoremfa",
                              #             "ScoreMFA (C++)" = "scoremfacpp"
                              #          ),
                              #       selected = "Cosine similarity"
                              #    )
                              # ),
                              br()
                           ),

                           bsCollapsePanel(

                              h5(strong("Heatmap parameters")),
                              br(),
                              div(
                                 style="display: inline-block;vertical-align:top; width: 125px;",
                                 selectInput(
                                    "heatmap_fill",
                                    tippy(
                                       "Heatmap fill",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>Color palette to use for score matrix heatmap. All palettes are from the viridis R package and are accessible to colorblind readers.</span>"
                                    ),
                                    choices =
                                       c(
                                          "magma",
                                          "inferno",
                                          "plasma",
                                          "viridis",
                                          "cividis"
                                       ),
                                    selected = "plasma"
                                 )
                              ),
                              div(
                                 style="display: inline-block;vertical-align:top; width: 125px;"
                              ),
                              div(
                                 style="display: inline-block;vertical-align:top; width: 125px;",
                                 numericInput(
                                    "heatmap_scale_coarse_start",
                                    tippy(
                                       "Coarse heatmap fill range, start",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>xxx</span>"
                                    ),
                                    value = 0,
                                    min = 0,
                                    max = 1,
                                    step = 0.01
                                 )
                              ),
                              div(
                                 style="display: inline-block;vertical-align:top; width: 125px;",
                                 numericInput(
                                    "heatmap_scale_coarse_end",
                                    tippy(
                                       "Coarse heatmap fill range, end",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>xxx</span>"
                                    ),
                                    value = 1,
                                    min = 0.01,
                                    max = 1,
                                    step = 0.01
                                 )
                              ),
                              div(
                                 style="display: inline-block;vertical-align:top; width: 125px;",
                                 numericInput(
                                    "heatmap_scale_fine_start",
                                    tippy(
                                       "Fine heatmap fill range, start",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>xxx</span>"
                                    ),
                                    value = 0,
                                    min = 0,
                                    max = 1,
                                    step = 0.001
                                 )
                              ),
                              div(
                                 style="display: inline-block;vertical-align:top; width: 125px;",
                                 numericInput(
                                    "heatmap_scale_fine_end",
                                    tippy(
                                       "Fine heatmap fill range, end",
                                       placement = "right",
                                       arrow = TRUE,
                                       allowHTML = TRUE,
                                       animation = "scale",
                                       duration = 250,
                                       theme = "light-border",
                                       tooltip =
                                          "<span style='font-size:14px;'>xxx</span>"
                                    ),
                                    value = 1,
                                    min = 0.001,
                                    max = 1,
                                    step = 0.001
                                 )
                              ),
                              br()
                           )
                        )
                     )
                  )
               ),
               tabPanel(
                  "About",
                  hr(),
                  includeMarkdown("about.md")
               )
            )
         ),


         # Main Panel --------------------------------------------------------------

         mainPanel(
            textOutput("output_test"),
            splitLayout(
               plotOutput("output_plot_MS"),
               uiOutput("output_html_MS")
            ),
            br(),
            splitLayout(
               plotOutput("output_plot_scoremat1"),
               plotOutput("output_plot_scoremat2")
            )
         ),
         fluid = FALSE
      )
   )
)

