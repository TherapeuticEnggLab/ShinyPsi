library(shiny)

ui <- fluidPage(
  titlePanel("ShinyPsi"),
    fluidRow(
      column(3,
         h5("Essential estim params"),
         numericInput("rho_c1", label = "rho_c", value = 0.027, min = 0, max = 1),
         numericInput("psi_c1", label = "psi_c", value = 0.121, min = 0, max = 1),
         numericInput("phi_c1", label = "phi_c", value = 0.578, min = 0, max = 1)
      ),
      column(3,
         h5("point estim params"),
         numericInput("alpha1", label = "alpha", value = 0.848),
         numericInput("beta1", label = "beta", value = 1.0)
      ),
      column(3,
         h5("Confidence estim params 1"),
         numericInput("n1", label = "n", value = 31869, min = 1),
         numericInput("sp1", label = "sp", value = 849, min = 1),
         numericInput("sn1", label = "sn", value = 31020, min = 1),
         numericInput("n_iter1", label = "n_iter", value = 100000, min = 1)
      ),
      column(3,
         h5("Confidence estim params 2"),
         numericInput("alpha_a1", label = "alpha_a", value = 100),
         numericInput("alpha_b1", label = "alpha_b", value = 17.3),
         numericInput("beta_a1", label = "beta_a", value = 68.3),
         numericInput("beta_b1", label = "beta_b", value = 0.2)
      )
    ),
    fluidRow(
      column(4,
         h4("Point Estimate"),
         actionButton("Calculate", label = "Calculate!"),
         p("psi_point"),
         textOutput("point_estim")
      )
    ),
  br(),
    fluidRow(
      column(4,
         h4("Confidence Estimate"),
         actionButton("Simulate", label = "Simulate!"),
         p("Median (IQR)"),
         textOutput("conf_estim")
      ),
      column(2,
      ),
      column(7,
      h4("If you find ShinyPsi useful please consider citing our manuscript:"),
      p("Tiwari, A., Chowdhury, S., James, A., Chatterjee, B., & Dixit, N. M. (2024). Adjusting for specificity of symptoms reveals higher prevalence of asymptomatic SARS-CoV-2 infections than previously estimated. medRxiv, 2024-09.")
      )
    )
)
server <- function(input, output, session) {
  # point estim
  numerator <- eventReactive(input$Calculate,
                             {input$rho_c1 * (1 - input$rho_c1) * (input$psi_c1 - 1 + input$phi_c1) * (input$alpha1 + input$beta1 - 1)
                             })
  denominator <- eventReactive(input$Calculate,
                               {(input$rho_c1 + input$beta1 - 1) * (input$psi_c1 * input$rho_c1 * (1 - input$alpha1) - input$alpha1 * (1 - input$phi_c1) * (1 - input$rho_c1))
                               })
  psi_temp <- eventReactive(input$Calculate,
                            {denominator1 <- ifelse(abs(denominator()) < 1e-08, NA, denominator())
                            1 - numerator() / denominator1
                            })
  psi <- eventReactive(input$Calculate,
                       {ifelse(is.na(psi_temp()) | psi_temp() < 0 | psi_temp() > 1, NA, psi_temp())
                       })
  output$point_estim <- renderText({
    as.character(round(psi(), 2))
  })
  
  # confidence estim
  alpha_d <- eventReactive(input$Simulate,
                           {rbeta(input$n_iter1, input$alpha_a1, input$alpha_b1)
                           })
  beta_d  <- eventReactive(input$Simulate,
                           {rbeta(input$n_iter1, input$beta_a1, input$beta_b1)
                           })
  rho_c_d <- eventReactive(input$Simulate,
                           {rbinom(input$n_iter1, input$n1, input$rho_c1) / input$n1
                           })
  psi_c_d <- eventReactive(input$Simulate,
                           {rbinom(input$n_iter1, input$sp1, input$psi_c1) / input$sp1
                           })
  phi_c_d <- eventReactive(input$Simulate,
                           {rbinom(input$n_iter1, input$sn1, input$phi_c1) / input$sn1
                           })
  
  numerator2 <- eventReactive(input$Simulate,
                              {rho_c_d() * (1 - rho_c_d()) * (psi_c_d() - 1 + phi_c_d()) * (alpha_d() + beta_d() - 1)
                              })
  
  denominator2 <- eventReactive(input$Simulate,
                                {(rho_c_d() + beta_d() - 1) * (psi_c_d() * rho_c_d() * (1 - alpha_d()) - alpha_d() * (1 - phi_c_d()) * (1 - rho_c_d()))
                                })
  
  psi_temp1 <- eventReactive(input$Simulate,
                             {denominator3 <- ifelse(abs(denominator2()) < 1e-08, NA, denominator2())
                             1 - numerator2() / denominator3
                             })
  psi1 <- eventReactive(input$Simulate,
                        {ifelse(is.na(psi_temp1()) | psi_temp1() < 0 | psi_temp1() > 1, NA, psi_temp1())
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
}

shinyApp(ui, server)