
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
library(rhandsontable)

# options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 1000*1024^2)

# Hidden tabs -------------------------------------------------------------



# UI ----------------------------------------------------------------------

shinyUI(

   fixedPage(
      titlePanel("abundancer"),
      theme = "maglab_theme_old.css",

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
                  "Analyze",
                  hr(),
                  tabsetPanel(
                     id = "guppipanel",
                     type = "tabs",
                     tabPanel(
                        "Upload File",
                        br(),
                        fileInput(
                           "fileInputUpload",
                           "Select a CSV or XLSX file",
                           accept = c(".xlsx", ".csv")
                        ),
                        numericInput(
                           "scoremat_charge",
                           "Charge",
                           value = 1,
                           min = 1,
                           max = 100,
                           step = 1
                        ),
                        textAreaInput(
                           "scoremat_sequence",
                           "Sequence",
                           placeholder =
                              "Single letter canonical amino acids only",
                           height = "100px"
                        ),
                        textInput(
                           "scoremat_PTM",
                           "PTM chemical formula",
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
                        splitLayout(
                           downloadButton("downloadPDF", label = "Heatmap PDF"),
                           downloadButton("downloadPNG", label = "Heatmap PNG")
                        )
                     ),
                     tabPanel(
                        "Settings",
                        h5(em("Hover parameter names for more info")),
                        h5(strong("Peak picking parameters")),
                        hr(),
                        div(
                           style="display: inline-block;vertical-align:top; width: 125px;",
                           numericInput(
                              "peakpick_SNR",
                              tippy(
                                 "SNR",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "From MSnbase documentation: \"A local maximum is considered a peak if its intensity is SNR times larger than the estimated noise.\" Noise is estimated using the MSnbase median absolute deviation method."
                              ),
                              value = 10,
                              min = 0,
                              max = 100,
                              step = 1
                           )
                        ),
                        div(
                           style="display: inline-block;vertical-align:top; width: 150px;",
                           numericInput(
                              "peakpick_k",
                              tippy(
                                 "KNN",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "From MSnbase documentation: \"m/z values and intensities of the 2 * KNN closest signals to the centroid are used in the intensity weighted average calculation.\""
                              ),
                              value = 2,
                              min = 1,
                              max = 10,
                              step = 1
                           )
                        ),
                        h5(strong("Score matrix parameters")),
                        hr(),
                        div(
                           style="display: inline-block;vertical-align:top; width: 150px;",
                           numericInput(
                              "scoremat_12C",
                              tippy(
                                 "<sup>12</sup>C abundance start",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "Starting value for <sup>12</sup>C abundance to use for calculating the cosine similarity score matrix. <sup>12</sup>C abundance will be incremented from this value to 1 by the abundance step."
                              ),
                              value = 0.985,
                              min = 0.9,
                              max = 0.999,
                              step = 0.001
                           )
                        ),
                        div(
                           style="display: inline-block;vertical-align:top; width: 150px;",
                           numericInput(
                              "scoremat_14N",
                              tippy(
                                 "<sup>14</sup>N abundance start",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "Starting value for <sup>14</sup>N abundance to use for calculating the cosine similarity score matrix. <sup>14</sup>N abundance will be incremented from this value to 1 by the abundance step."
                              ),                              value = 0.985,
                              min = 0.9,
                              max = 0.999,
                              step = 0.001
                           )
                        ),
                        div(
                           style="display: inline-block;vertical-align:top; width: 150px;",
                           selectInput(
                              "scoremat_abundStep",
                              tippy(
                                 "Abundance step",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "Value by which the <sup>12</sup>C and <sup>14</sup>N abundance starts will be stepped during calculation of the cosine similarity score matrix. A smaller value will yield longer calculation times."
                              ),
                              choices =
                                 c(
                                    0.001,
                                    0.0001
                                 ),
                              selected = 0.001
                           )
                        ),
                        div(
                           style="display: inline-block;vertical-align:top; width: 150px;",
                           numericInput(
                              "scoremat_binSize",
                              tippy(
                                 "Bin size",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "Prior to comparing observed and theoretical spectra, spectra are binned using this bin size value. See MSnbase documentation for more information."
                              ),
                              value = 0.05,
                              min = 0.01,
                              max = 0.1,
                              step = 0.01
                           )
                        ),
                        div(
                           style="display: inline-block;vertical-align:top; width: 150px;",
                           numericInput(
                              "scoremat_isoAbundCutoff",
                              tippy(
                                 "Isotopic abundance cutoff",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "During generation of theoretical spectra, theoretical isotopologues whose relative abundance are lower than this value will be removed and not used for spectral comparison."
                              ),
                              value = 5,
                              min = 1,
                              max = 99,
                              step = 1
                           )
                        ),
                        h5(strong("Heatmap parameters")),
                        hr(),
                        div(
                           style="display: inline-block;vertical-align:top; width: 150px;",
                           selectInput(
                              "heatmap_fill",
                              tippy(
                                 "Heatmap fill",
                                 placement = "right",
                                 arrow = TRUE,
                                 allowHTML = TRUE,
                                 animation = "scale",
                                 duration = 250,
                                 theme = "light",
                                 tooltip =
                                    "Color palette to use for score matrix heatmap. All palettes are from the viridis R package and are accessible to colorblind readers."
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
                        )
                        # actionButton(
                        #    "testButton1",
                        #    "Test button 1"
                        # )
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
               plotOutput("output_plot_scoremat")
            )
         ),
         fluid = FALSE
      )
   )
)
