library(shiny)
library(shinythemes)
library(stats)
library(shinycssloaders)
library(shinyFeedback)

estimate_beta <- function(q_vals, sen_range = c(0.5, 100), spe_range = c(0.5, 100)) {
  best_loss <- Inf
  best_params <- c(NA, NA)
  p_vals <- c(0.025, 0.5, 0.975)
  for (a in seq(sen_range[1], sen_range[2], length.out = 1000)) {
    for (b in seq(spe_range[1], spe_range[2], length.out = 1000)) {
      q_model <- qbeta(p_vals, a, b)
      loss <- sum((q_model - q_vals)^2)
      if (loss < best_loss) {
        best_loss <- loss
        best_params <- c(a, b)
      }
    }
  }
  return(data.frame(sen = best_params[1], spe = best_params[2]))
}

correction <- function(rho_c, psi_c, phi_c, alpha, beta) {
  numerator <- rho_c * (1 - rho_c) * (psi_c - 1 + phi_c) * (alpha + beta - 1)
  denominator <- (rho_c + beta - 1) * (psi_c * rho_c * (1 - alpha) - alpha * (1 - phi_c) * (1 - rho_c))
  denominator <- ifelse(abs(denominator) < 1e-08, NA, denominator)
  
  psi_temp <- 1 - numerator / denominator
  psi <- ifelse(is.na(psi_temp) | psi_temp < 0 | psi_temp > 1, NA, psi_temp)
  
  c(psi)
}


monte_carlo <- function (n, sp, sn, n_iter,
                         rho_c, psi_c, phi_c,
                         alpha_a, alpha_b, beta_a, beta_b) {
  set.seed(123)
  n_iter <- as.numeric(n_iter)
  alpha <- rbeta(n_iter, alpha_a, alpha_b)
  beta  <- rbeta(n_iter, beta_a, beta_b)
  rho_c <- rbinom(n_iter, n, rho_c) / n
  psi_c <- rbinom(n_iter, sp, psi_c) / sp
  phi_c <- rbinom(n_iter, sn, phi_c) / sn
  
  psi_hat <- correction(rho_c, psi_c, phi_c, alpha, beta)
  
  #psi_l = quantile(psi_hat, probs = 0.25, na.rm = TRUE)
  #psi = quantile(psi_hat, probs = 0.5, na.rm = TRUE)
  #psi_u = quantile(psi_hat, probs = 0.75, na.rm = TRUE)
  
  #return(c(psi_l, psi, psi_u))
  return(psi_hat)
}

downloadButton_fixed <- function(...) {
  tag <- shiny::downloadButton(...)
  tag$attribs$download <- NULL
  tag
}

ui <- navbarPage(theme = shinytheme("superhero"),
                tags$head(
                  tags$style(HTML("
                 .shiny-input-container label {
                 font-size: 20px;
                 }
                 .justified-text {
                 text-align: justify;
                 }
                "))
                ),
                  tabPanel("Application",
                           titlePanel("ShinyPsi"),
                           h5("Adjusts asymptomatic-infection prevalence for nonspecific symptoms and test characteristics"),
                           fluidRow(
                             column(3,
                                    h5("Essential parameters"),
                                    useShinyFeedback(),
                                    numericInput("rho_c1", label = HTML(paste0("ρ", tags$sub("c"))), value = 0.028, min = 0, max = 1),
                                    useShinyFeedback(),
                                    numericInput("psi_c1", label = HTML(paste0("ψ", tags$sub("c"))), value = 0.304, min = 0, max = 1)
                                    
                             ),
                             column(3,
                                    h5("Point estimation parameters"),
                                    useShinyFeedback(),
                                    numericInput("alpha1", label = "α", value = 0.967, min = 0, max = 1),
                                    useShinyFeedback(),
                                    numericInput("beta1", label = "β", value = 0.995, min = 0, max = 1)
                             ),
                             column(3,
                                    h5("Confidence estimation set1"),
                                    numericInput("n1", label = "n", value = 13095, min = 1),
                                    selectInput("n_iter1", label = "iterations", c("select", "100000", "10000", "1000"))
                             ),
                             column(3,
                                    h5("Confidence estimation set2"),
                                    numericInput("alpha_l", label = HTML(paste0("α", " [95% CI: lower limit]")), value = 0.924, min = 0, max = 1),
                                    numericInput("alpha_u", label = HTML(paste0("α", " [95% CI: upper limit]")), value = 0.986, min = 0, max = 1)
                             )
                           ),
                           fluidRow(
                             column(3,
                                    useShinyFeedback(),
                                    numericInput("phi_c1", label = HTML(paste0("φ", tags$sub("c"))), value = 0.432, min = 0, max = 1),
                                    h4("Point Estimate"),
                                    actionButton("Calculate", label = "Calculate!", class = "btn-warning"),
                                    p(HTML(paste0("ψ", tags$sub("point")))),
                                    textOutput("point_estim")
                             ),
                             column(6,
                                    withSpinner(plotOutput("plot", width = "600px", height = "230px"), 
                                                caption = "Hold tight, press no buttons, we're doing the math now.")
                             ),
                             column(3,
                                    numericInput("beta_l", label = HTML(paste0("β", " [95% CI: lower limit]")), value = 0.987, min = 0, max = 1),
                                    numericInput("beta_u", label = HTML(paste0("β", " [95% CI: upper limit]")), value = 0.998, min = 0, max = 1),
                                    p("% of failed estimates in the virtual sample:"),
                                    textOutput("failed_estimates")
                             )
                           ),
                           fluidRow(
                             column(3,
                                    h4("Confidence Estimate"),
                                    actionButton("Simulate", label = "Simulate!", class = "btn-warning"
                                    ),
                                    p("Median (IQR)"),
                                    textOutput("conf_estim")
                             ),
                             column(6,
                                    h4("If you find ShinyPsi useful, please consider citing our manuscript:"),
                                    p("Tiwari et al. (2026). To be updated soon.")
                             ),
                             column(3,
                                    br(),
                                    actionButton("reset", label = "Reset", class = "btn-warning")
                             )
                           )
                  ),
                  tabPanel("User Manual",
                           column(4,
                                  h4("A concise description of the parameters"),
                                  h4(HTML("ρ<sub>c")),
                                  p("Unajusted prevalence of test-positive cases among sampled individuals"),
                                  p(),
                                  h4(HTML("ψ<sub>c")),
                                  p("Unadjusted prevalence of asymptomatic cases among test positive individuals"),
                                  p(),
                                  h4(HTML("φ<sub>c")),
                                  p("Unadjusted prevalence of symptomatic cases among test-negative individuals"),
                                  p(),
                                  h4(HTML("α")),
                                  p("Test sensitivity"),
                                  p(),
                                  h4(HTML("β")),
                                  p("Test specificity"),
                                  p(),
                                  h4(HTML("n")),
                                  p("Number of sampled individuals"),
                                  p(),
                                  h4(HTML("iterations")),
                                  p("Virtual sample size for generating the distribution of ψ"),
                                  p(),
                                  h4("For a more detailed discussion on the paremeters please check out our paper:"),
                                  p("Tiwari et al. (2026). To be updated soon.")
                           ),
                           column(8,
                                  h2("Usage guidelines"),
                                  p(class = "justified-text",
                                  "The following instructions are written for a user who is using the app from their laptop or desktop with any operating system (we have tested Windows, Ubuntu and MacOS) and any web browser (we have tested Chrome, Edge, Safari, Firefox). The app can be used on the phones as well, although we do not recommend it. All calculations in ShinyPsi happen locally at your browser and that way no data ever leaves your computer."),
                                  p(class = "justified-text",
                                    "The app, by default, has certain input values for each input space (measurements from a serosuvey, Espenhain et al., 2021). Thus, immediately after loading or refreshing the app page, if you press “Calculate!” and “Simulate!” – i.e., the action buttons on the left bottom of the screen, results are displayed (action buttons are orange). You have to choose the number of iterations though. These results correspond to the latter default values, which we have kept for demonstration purposes."),
                                  p(class = "justified-text",
                                    "When you begin the analysis with your own data, first press the reset button on the bottom right of the screen. This will remove all the default values from the input boxes. Note that this step (the “Reset” action) is also important if you want to download your analysis outputs at the end (see below). Reset makes sure previous input values are not used in the fresh analysis."),
                                  p(class = "justified-text",
                                    "Your first inputs should be the first column of the app, i.e., the essential parameters. We recommend that you input these fractions with as high precision as possible since these fractions are utilized further for estimating the distribution of psi."),
                                  p(class = "justified-text",
                                    "If you are attempting to get just the point estimate of psi, you need the second column of the app, i.e., the point estimation parameters. Columns 3 and 4 may remain blank. Press the “Calculate!” action button."),
                                  p(class = "justified-text",
                                    "If you are looking for both the point estimate and the distribution of ψ, you need the next two columns (see the description of parameters). After having the values of the parameters in all boxes, press “Simulate!” action button. Depending on the number of the iteration it might take few seconds to maybe minutes. All input values except n and iterations should be between 0 and 1. Any other input values or missing values will generate feedback for the user."),
                                  p(class = "justified-text",
                                    "Once the process is complete results will be displayed immediately below the action buttons. The plot that shows the distribution and the point estimate on it along with the IQRs will be generated simultaneously. Resetting will remove the plot alongside the inputs."),
                                  p(class = "justified-text",
                                    "Notice that above the reset button, the percentage of virtual sample that could not be successfully estimated is recorded. If this percentage is high (above 40%, for example) it is indicative of lower quality of the prediction of the confidence interval. The application will show outputs till this value reaches 95%. The user should be cautious in interpreting the results with high % failure rate. Beyond 95%, the app will generate an error message. Similarly, for certain combiantions of the parameter values the point estimate may not be reliable. In such cases, the app will generate a dialogue box to inform the user and will not display the unreliable point estimates."),
                                  p(class = "justified-text",
                                    "In case you want to keep a record of your analysis this is the time to back it up. Move to the Downloads tab. This tab will show the tables of point and confidence estimation. It will also show the same plot of the simulated distribution that is displayed on the application page. The tables are ready for downloading when you can see them on the top of the download buttons. Resetting removes these tables and the plot. Download the csv files as necessary. Note that if you do not keep a record of your analysis and close/reload the app, all calculations will be lost. As stated earlier, your data never leaves your computer when you use our app – so we cannot help you retrieve your past analysis. Please back up your calculations before you end your session."),
                                  p(class = "justified-text",
                                    "We welcome your feedback - please write to us at narendra@iisc.ac.in.")
                           )
                  ),
                  tabPanel("Downloads",
                  fluidRow(         
                  column(4, 
                         h2("Point estimate"),
                         tableOutput("Point_estimate_table"),
                         downloadButton_fixed("Point_estimate_download")
                  ),
                  column(4, 
                         h2("Confidence estimate"),
                         tableOutput("Confidence_estimate_table"),
                         downloadButton_fixed("Confidence_estimate_download")
                  ),
                  column(4, 
                         h2("Distribution"),
                         plotOutput("plot_copy", width = "400px", height = "230px"),
                         downloadButton_fixed("Distribution_download")
                  )
                  ),
                  fluidRow(),
                  fluidRow(
                    column(10,
                           h4("If you do not see the tables and the plot above the respective download buttons, do not press download."))
                  )
                )
                )


server <- function(input, output, session) {
  reset_state <- reactiveVal(1)
  
  observeEvent(input$reset, {
    reset_state(3)
  })
  observeEvent(input$Calculate, {
    reset_state(2)
  })
  observeEvent(input$Simulate, {
    reset_state(1)
  })
  
  # reset
  observeEvent(input$reset, {
    updateNumericInput(session, "rho_c1", value = NA)
    updateNumericInput(session, "psi_c1", value = NA)
    updateNumericInput(session, "phi_c1", value = NA)
    updateNumericInput(session, "alpha1", value = NA)
    updateNumericInput(session, "beta1", value = NA)
    updateNumericInput(session, "phi_c1", value = NA)
    updateNumericInput(session, "n1", value = NA)
    updateNumericInput(session, "alpha_l", value = NA)
    updateNumericInput(session, "alpha_u", value = NA)
    updateNumericInput(session, "beta_l", value = NA)
    updateNumericInput(session, "beta_u", value = NA)
    updateSelectInput(session, "n_iter1", choices = c("select", "100000", "10000", "1000"))
  })
  
  # input validation for point estim
  observeEvent(input$Calculate, {
    feedbackDanger("rho_c1", (input$rho_c1>1)|(input$rho_c1<0)|(is.na(input$rho_c1)), "input range 0-1")
    feedbackDanger("psi_c1", (input$psi_c1>1)|(input$psi_c1<0)|(is.na(input$psi_c1)), "input range 0-1")
    feedbackDanger("phi_c1", (input$phi_c1>1)|(input$phi_c1<0)|(is.na(input$phi_c1)), "input range 0-1")
    feedbackDanger("alpha1", (input$alpha1>1)|(input$alpha1<0)|(is.na(input$alpha1)), "input range 0-1")
    feedbackDanger("beta1", (input$beta1>1)|(input$beta1<0)|(is.na(input$beta1)), "input range 0-1")
  })
  
  observeEvent(input$Simulate, {
    feedbackDanger("rho_c1", (input$rho_c1>1)|(input$rho_c1<0)|(is.na(input$rho_c1)), "input range 0-1")
    feedbackDanger("psi_c1", (input$psi_c1>1)|(input$psi_c1<0)|(is.na(input$psi_c1)), "input range 0-1")
    feedbackDanger("phi_c1", (input$phi_c1>1)|(input$phi_c1<0)|(is.na(input$phi_c1)), "input range 0-1")
    feedbackDanger("alpha1", (input$alpha1>1)|(input$alpha1<0)|(is.na(input$alpha1)), "input range 0-1")
    feedbackDanger("beta1", (input$beta1>1)|(input$beta1<0)|(is.na(input$beta1)), "input range 0-1")
    feedbackDanger("n1", (input$n1<0)|(input$n1 %% 1 != 0)|(is.na(input$n1)), "positive integer")
    feedbackDanger("n_iter1", input$n_iter1=="select", "Please select input")
    feedbackDanger("alpha_l", (input$alpha_l>1)|(input$alpha_l<0)|(is.na(input$alpha_l)), "input range 0-1")
    feedbackDanger("alpha_u", (input$alpha_u>1)|(input$alpha_u<0)|(is.na(input$alpha_u)), "input range 0-1")
    feedbackDanger("beta_l", (input$beta_l>1)|(input$beta_l<0)|(is.na(input$beta_l)), "input range 0-1")
    feedbackDanger("beta_u", (input$beta_u>1)|(input$beta_u<0)|(is.na(input$beta_u)), "input range 0-1")
  })
  
  # point estim
  psi <- eventReactive(input$Calculate, {
    req((input$rho_c1<=1),(input$rho_c1>=0),(!is.na(input$rho_c1)))
    req((input$psi_c1<=1),(input$psi_c1>=0),(!is.na(input$psi_c1)))
    req((input$phi_c1<=1),(input$phi_c1>=0),(!is.na(input$phi_c1)))
    req((input$alpha1<=1),(input$alpha1>=0),(!is.na(input$alpha1)))
    req((input$beta1<=1),(input$beta1>=0),(!is.na(input$beta1)))
    correction(input$rho_c1, input$psi_c1, input$phi_c1, input$alpha1, input$beta1)
  })
  
  observeEvent(input$Calculate, {
    if (reset_state()%in%c(1,2)) {
      if (is.na(psi())) {
        showModal(
          modalDialog(class = "justified-text",
                      "The point estimate of ψ with the current parameter set (Q) is falling in the zones of unreliability. For more details on when and why that might happen please check figure 1 and the relevant text section of our manuscript (see at the bottom of this page for the reference). We also suggest that the user may please take a look at the simulated confidence intervals of ψ as the potential alternative.",
                      title = "Unreliable point estimate",
                      footer = tagList(
                        actionButton("cancel1", "Ok, I understand.", class = "btn-warning")
                      )
          )
        )
      }
    }
  })
  
  observeEvent(input$cancel1, {
    removeModal()
  })
  
  output$point_estim <- renderText({
    if (reset_state() == 3) {
      ""
    }
    if (reset_state()%in%c(1,2)) {
      if (is.na(psi())) {
        "unreliable point estimate"
      } else {
        as.character(round(psi(), 3))
      }
    }
  })
  
  output$Point_estimate_table <- renderTable({
    if (reset_state() == 3) {
      ""
    }
    if (reset_state()%in%c(1,2)) {
      data.frame(Params = c("rho_c", "psi_c", "phi_c", "alpha", "beta", "psi_point"), 
                 Values = c(input$rho_c1, input$psi_c1, input$phi_c1, input$alpha1, input$beta1, round(psi(), 3)))
    }
  })
  
  output$Point_estimate_download <- downloadHandler(
      filename = function() {
        paste0("point_estimate",".csv")
      },
      content = function(file) {
      write.csv(data.frame(Params = c("rho_c", "psi_c", "phi_c", "alpha", "beta", "psi_point"), 
                 Values = c(input$rho_c1, input$psi_c1, input$phi_c1, input$alpha1, input$beta1, round(psi(), 3))),
                file)
      }
  )
  
  # confidence estim
  sp <- eventReactive(input$Simulate, {round(input$rho_c1 * input$n1)})
  sn <- eventReactive(input$Simulate, {round((1 - input$rho_c1) * input$n1)})
  
  Alpha <- eventReactive(input$Simulate, {
    req((input$rho_c1<=1),(input$rho_c1>=0),(!is.na(input$rho_c1)))
    req((input$psi_c1<=1),(input$psi_c1>=0),(!is.na(input$psi_c1)))
    req((input$phi_c1<=1),(input$phi_c1>=0),(!is.na(input$phi_c1)))
    req((input$alpha1<=1),(input$alpha1>=0),(!is.na(input$alpha1)))
    req((input$beta1<=1),(input$beta1>=0),(!is.na(input$beta1)))
    req((input$alpha_l<=1),(input$alpha_l>=0),(!is.na(input$alpha_l)))
    req((input$alpha_u<=1),(input$alpha_u>=0),(!is.na(input$alpha_u)))
    req((input$beta_l<=1),(input$beta_l>=0),(!is.na(input$beta_l)))
    req((input$beta_u<=1),(input$beta_u>=0),(!is.na(input$beta_u)))
    req(input$n_iter1!="select")
    req((input$n1>0),(input$n1 %% 1 == 0),(!is.na(input$n1)))
    estimate_beta(c(input$alpha_l, input$alpha1, input$alpha_u))})
  
  Beta <- eventReactive(input$Simulate, {
    estimate_beta(c(input$beta_l, input$beta1, input$beta_u))})
  
  
  psi1 <- eventReactive(input$Simulate, {
                        monte_carlo(input$n1, sp(), sn(), input$n_iter1,
                                     input$rho_c1, input$psi_c1, input$phi_c1,
                                     Alpha()$sen, Alpha()$spe, Beta()$sen, Beta()$spe)
                        })
  
  outstring1 <- eventReactive(input$Simulate, {
                             psi_l2 = quantile(psi1(), probs = 0.25, na.rm = TRUE)
                             psi2 = quantile(psi1(), probs = 0.5, na.rm = TRUE)
                             psi_u2 = quantile(psi1(), probs = 0.75, na.rm = TRUE)
                             paste0(as.character(round(psi2,2)), ' (', as.character(round(psi_l2, 2)), ' ,', as.character(round(psi_u2, 2)), ')')
                             })
  output$conf_estim <- renderText({
    if (reset_state() == 3) {
      ""
    }
    if (reset_state() == 2) {
      "Run simulate for the distribution"
    }
    if (reset_state() == 1) {
      if ((length(which(is.na(psi1())))/as.numeric(input$n_iter1))>0.95) {
        "unreliable confidence estimate"
      } else {
        outstring1()
      }
    }
  })
  
  outstring2 <- eventReactive(input$Simulate, {
    paste0(as.character(round((length(which(is.na(psi1())))/as.numeric(input$n_iter1))*100, 2)), "%")
  })
  output$failed_estimates <- renderText({
    if (reset_state() == 3) {
      ""
    }
    if (reset_state() == 2) {
      "Simulate first"
    }
    if (reset_state() == 1) {
      outstring2()
    }
  })
  
  output$Confidence_estimate_table <- renderTable({
    if (reset_state() == 3) {
      ""
    }
    if (reset_state() == 2) {
      "Run simulate for the distribution"
    }
    if (reset_state() == 1) {
      data.frame(Params = c("rho_c", "psi_c", "phi_c", "alpha", "beta", "n", "iterations", 
                            "alpha_l", "alpha_u", "beta_l", "beta_u", "Q1", "Median", "Q3"), 
                 Values = c(input$rho_c1, input$psi_c1, input$phi_c1, input$alpha1, input$beta1, input$n1, input$n_iter1,
                            input$alpha_l, input$alpha_u, input$beta_l, input$beta_u,
                            round(quantile(psi1(), probs = 0.25, na.rm = TRUE), 2),
                            round(quantile(psi1(), probs = 0.5, na.rm = TRUE), 2),
                            round(quantile(psi1(), probs = 0.75, na.rm = TRUE), 2)
                            ))
    }
  })
  
  output$Confidence_estimate_download <- downloadHandler(
    filename = function() {
      paste0("confidence_estimate",".csv")
    },
    content = function(file) {
      write.csv(data.frame(Params = c("rho_c", "psi_c", "phi_c", "alpha", "beta", "n", "iterations", 
                                      "alpha_l", "alpha_u", "beta_l", "beta_u", "Q1", "Median", "Q3"), 
                           Values = c(input$rho_c1, input$psi_c1, input$phi_c1, input$alpha1, input$beta1, input$n1, input$n_iter1,
                                      input$alpha_l, input$alpha_u, input$beta_l, input$beta_u,
                                      round(quantile(psi1(), probs = 0.25, na.rm = TRUE), 2),
                                      round(quantile(psi1(), probs = 0.5, na.rm = TRUE), 2),
                                      round(quantile(psi1(), probs = 0.75, na.rm = TRUE), 2)
                           )),
                file)
    }
  )
  
  output$Distribution_download <- downloadHandler(
    filename = function() {
      paste0("whole_distribution",".csv")
    },
    content = function(file) {
      write.csv(data.frame(Dist = na.omit(psi1())),
                file)
    }
  )
  
  output$plot <- renderPlot({
    if (reset_state() == 3) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
           xlim = c(0, 10), ylim = c(0, 10), main = "")
      text(x = 5, y = 5, labels = "The distribution of ψ will appear here after you simulate", 
           cex = 1, col = "blue")
    }
    if (reset_state() == 2) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
           xlim = c(0, 10), ylim = c(0, 10), main = "")
      text(x = 5, y = 5, labels = "The distribution of ψ will appear here after you simulate", 
           cex = 1, col = "blue")
    }
    if (reset_state() == 1) {
      if ((length(which(is.na(psi1())))/as.numeric(input$n_iter1))>0.95) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
             xlim = c(0, 10), ylim = c(0, 10), main = "")
        text(x = 5, y = 9, labels = "Less than 5% of the virtual samples are successfully", 
             cex = 1, col = "red")
        text(x = 5, y = 6, labels = "estimated. Please revise this parameter set.", 
             cex = 1, col = "red")
        text(x = 5, y = 3, labels = "See the Figure 1 in the manuscript for details.", 
             cex = 1, col = "red")
      } else {
      par(
        mar = c(2.85, 2.7, 0.8, 1),   # bottom, left, top, right margins
        mgp = c(1.8, 0.6, 0),  # axis title, axis labels, axis line
        tcl = -0.25            # shorter tick marks
      )
      plot(
        density(psi1(), na.rm = TRUE),
        main = "",
        xlab = "ψ",
        ylab = "density",
        col = "darkgray",
        lwd = 3.5
      )
      abline(v = quantile(psi1(), probs = 0.25, na.rm = TRUE), col = "#076098", lwd = 2)
      abline(v = quantile(psi1(), probs = 0.50, na.rm = TRUE), col = "black", lwd = 2)
      abline(v = quantile(psi1(), probs = 0.75, na.rm = TRUE), col = "#076098", lwd = 2)
      box()
      legend("topleft", legend = c("Q1, Q3", "Median"), 
             col = c("#076098", "black"), pch = "—", 
             bty = "n")
      }
    }
  }, res = 96)
  
  output$plot_copy <- renderPlot({
    if (reset_state() == 3) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
           xlim = c(0, 10), ylim = c(0, 10), main = "")
      text(x = 5, y = 6.5, labels = "Simulate the distribution", 
           cex = 1, col = "blue")
      text(x = 5, y = 3.5, labels = "before downloading it.", 
           cex = 1, col = "blue")
    }
    if (reset_state() == 2) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
           xlim = c(0, 10), ylim = c(0, 10), main = "")
      text(x = 5, y = 6.5, labels = "Simulate the distribution", 
           cex = 1, col = "blue")
      text(x = 5, y = 3.5, labels = "before downloading it.", 
           cex = 1, col = "blue")
    }
    if (reset_state() == 1) {
      if ((length(which(is.na(psi1())))/as.numeric(input$n_iter1))>0.95) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
             xlim = c(0, 10), ylim = c(0, 10), main = "")
        text(x = 5, y = 6.5, labels = "Confidence estimates are unreliable.", 
             cex = 1, col = "red")
        text(x = 5, y = 3.5, labels = "Return to the application tab.", 
             cex = 1, col = "red")
      } else {
      par(
        mar = c(2.85, 2.7, 0.8, 1),   # bottom, left, top, right margins
        mgp = c(1.8, 0.6, 0),  # axis title, axis labels, axis line
        tcl = -0.25            # shorter tick marks
      )
      plot(
        density(psi1(), na.rm = TRUE),
        main = "",
        xlab = "ψ",
        ylab = "density",
        col = "darkgray",
        lwd = 3.5
      )
      abline(v = quantile(psi1(), probs = 0.25, na.rm = TRUE), col = "#076098", lwd = 2)
      abline(v = quantile(psi1(), probs = 0.50, na.rm = TRUE), col = "black", lwd = 2)
      abline(v = quantile(psi1(), probs = 0.75, na.rm = TRUE), col = "#076098", lwd = 2)
      box()
      legend("topleft", legend = c("Q1, Q3", "Median"), 
             col = c("#076098", "black"), pch = "—", 
             bty = "n")
      }
    }
  }, res = 96)
}

shinyApp(ui, server)
  