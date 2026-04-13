library(shiny)
library(shinythemes)
library(stats)

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

set.seed(123)
monte_carlo <- function (n, sp, sn, n_iter,
                         rho_c, psi_c, phi_c,
                         alpha_a, alpha_b, beta_a, beta_b) {
  
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

ui <- fluidPage(theme = shinytheme("superhero"),
                tags$head(
                  tags$style(HTML("
                  .nav-tabs {
                  display: flex;
                  justify-content: flex-end;
                  }
                  .shiny-input-container label {
                  font-size: 20px;
                  }
                  "))
                ),
                tabsetPanel(
                  tabPanel("Application",
                           titlePanel("ShinyPsi"),
                           h6("beta version"),
                           fluidRow(
                             column(3,
                                    h5("Essential parameters"),
                                    numericInput("rho_c1", label = HTML(paste0("ρ", tags$sub("c"))), value = 0.027, min = 0, max = 1),
                                    numericInput("psi_c1", label = HTML(paste0("ψ", tags$sub("c"))), value = 0.121, min = 0, max = 1)
                                    
                             ),
                             column(3,
                                    h5("point estimation parameters"),
                                    numericInput("alpha1", label = "α", value = 0.848, min = 0, max = 1),
                                    numericInput("beta1", label = "β", value = 1.0, min = 0, max = 1)
                             ),
                             column(3,
                                    h5("Confidence estimation set1"),
                                    numericInput("n1", label = "n", value = 31869, min = 1),
                                    numericInput("n_iter1", label = "iterations", value = 100000, min = 1)
                             ),
                             column(3,
                                    h5("Confidence estimation set2"),
                                    numericInput("alpha_l", label = HTML(paste0("α", " (lower bound)")), value = 0.75, min = 0, max = 1),
                                    numericInput("alpha_u", label = HTML(paste0("α", " (upper bound)")), value = 0.9, min = 0, max = 1)
                             )
                           ),
                           fluidRow(
                             column(3,
                                    numericInput("phi_c1", label = HTML(paste0("φ", tags$sub("c"))), value = 0.578, min = 0, max = 1),
                                    h4("Point Estimate"),
                                    actionButton("Calculate", label = "Calculate!", class = "btn-warning"),
                                    p("psi_point"),
                                    textOutput("point_estim")
                             ),
                             column(6,
                                    plotOutput("plot", width = "600px", height = "200px")
                             ),
                             column(3,
                                    numericInput("beta_l", label = HTML(paste0("β", " (lower bound)")), value = 0.9, min = 0, max = 1),
                                    numericInput("beta_u", label = HTML(paste0("β", " (upper bound)")), value = 1.0, min = 0, max = 1)
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
                                    actionButton("reset", label = "Reset", class = "btn-warning")
                             )
                           )
                  ),
                  tabPanel("User Manual",
                           column(4,
                                  h4("A concise description of the parameters")
                           ),
                           column(8,
                                  h2("Usage guidelines"),
                                  p("The following instructions are written for a user who is using the app from their laptop or desktop with any operating system and any web browser (well, just avoid Internet Explorer, just in case that is still installed in your computer). The app can be used on the phones as well, although we do not recommend it. However, the inputs and outputs are rendered as a long list on a phone screen to adapt to the smaller screen real estate. We recommend any ‘standard’ laptop or desktop. As we discussed in the paper, all calculations in ShinyPsi happen locally at your browser and that way no data ever leaves your computer.
The app, by default, has certain input values for each input space. Thus, immediately after loading or refreshing the app page, if you press “calculate!” and “simulate!” – i.e., the action buttons on the left bottom of the screen, results are displayed (action buttons are orange). These results correspond to the latter default values, which we have kept for demonstration purposes.
When you begin the analysis with your own data, first press the reset button on the bottom right of the screen. This will remove all the default values from the input boxes, except for the value for iteration. You can manually alter this value in specific scenarios, let’s say you want to explore the convergence of the estimates for a spectrum of number of iterations. But for most purposes the default value should be preferred. Note that this step (“Reset Blank” action) is also important if you want to download your analysis report at the end of your analysis (see below). Pressing the reset button at the very beginning of the analysis initiates the process of eventually exporting and generating a report at the end of the analysis. 
Your first inputs should be the first column of the app: the essential parameters: XX, XX and XX. We recommend that you input these fractions with as high precision as possible since these fractions are utilized further for estimating the distribution of psi.
If you are attempting to get just the point estimate of psi, i.e., psi point, you need the second column of the app: i.e., alpha and beta (see the description of the parameters). Columns 3 and 4 should remain blank. Press the “Calculate!” action button.
If you are looking for both the point and the distribution, you need the next two columns (see the description of parameters). After having the values of the parameters in all boxes, press “Simulate!” action button. Depending on the value of the iteration it might take few seconds to maybe minutes. Just hang on there.
Once the process is complete results will be displayed immediately below the action buttons. We recommend that you press the “Visualize” action button next to generate the plot that shows the distribution and the point estimate on it along with the IQRs.
In case you want to keep a record of your analysis this is the time to back it up. Press the action button “export report” and move to the tab Download Report. This tab will summarize your analysis. It will ask for the title of your report. Add a suitable title. Next download the report as txt or pdf file. Note that if you do not keep a record of your analysis and close/reload the app, everything will be lost. As stated earlier, your data never leaves your computer when you use our app – so we cannot help you retrieve your past analysis using our app – even if we intend to. Please back up your calculations before you end your session.")
                           )
                  ),
                  tabPanel("Download Report"),
                )
)

server <- function(input, output, session) {
  # reset
  observeEvent(input$reset, {
    updateNumericInput(session, "rho_c1", value = "")
    updateNumericInput(session, "psi_c1", value = "")
    updateNumericInput(session, "phi_c1", value = "")
    updateNumericInput(session, "alpha1", value = "")
    updateNumericInput(session, "beta1", value = "")
    updateNumericInput(session, "phi_c1", value = "")
    updateNumericInput(session, "n1", value = "")
    updateNumericInput(session, "alpha_l", value = "")
    updateNumericInput(session, "alpha_u", value = "")
    updateNumericInput(session, "beta_l", value = "")
    updateNumericInput(session, "beta_u", value = "")
  })
  
  # point estim
  psi <- eventReactive(input$Calculate, {
    correction(input$rho_c1, input$psi_c1, input$phi_c1, input$alpha1, input$beta1)
  })
  
  output$point_estim <- renderText({
    as.character(round(psi(), 2))
  })
  
  # confidence estim
  sp <- eventReactive(input$Simulate, {round(input$rho_c1 * input$n1)})
  sn <- eventReactive(input$Simulate, {round((1 - input$rho_c1) * input$n1)})
  
  Alpha <- eventReactive(input$Simulate, {estimate_beta(c(input$alpha_l, input$alpha1, input$alpha_u))})
  
  Beta <- eventReactive(input$Simulate, {estimate_beta(c(input$beta_l, input$beta1, input$beta_u))})
  
  
  
  psi1 <- eventReactive(input$Simulate,
                        {monte_carlo(input$n1, sp(), sn(), input$n_iter1,
                                     input$rho_c1, input$psi_c1, input$phi_c1,
                                     Alpha()$sen, Alpha()$spe, Beta()$sen, Beta()$spe)
                        })
  
  
  outstring <- eventReactive(input$Simulate,
                             {psi_l2 = quantile(psi1(), probs = 0.25, na.rm = TRUE)
                             psi2 = quantile(psi1(), probs = 0.5, na.rm = TRUE)
                             psi_u2 = quantile(psi1(), probs = 0.75, na.rm = TRUE)
                             paste0(as.character(round(psi2,2)), ' (', as.character(round(psi_l2, 2)), ' ,', as.character(round(psi_u2, 2)), ')')
                             })
  output$conf_estim <- renderText({
    outstring()
  })
  output$plot <- renderPlot({
    par(
      mar = c(2, 2, 1, 1),   # bottom, left, top, right margins
      mgp = c(1.8, 0.6, 0),  # axis title, axis labels, axis line
      tcl = -0.25            # shorter tick marks
    )
    plot(
      density(psi1(), na.rm = TRUE),
      main = "",
      xlab = "",
      ylab = "",
      lwd = 2
    )
    abline(v = quantile(psi1(), probs = 0.25, na.rm = TRUE), col = "gray")
    abline(v = quantile(psi1(), probs = 0.50, na.rm = TRUE), col = "black")
    abline(v = quantile(psi1(), probs = 0.75, na.rm = TRUE), col = "gray")
    box()
  }, res = 96)
}

shinyApp(ui, server)
  