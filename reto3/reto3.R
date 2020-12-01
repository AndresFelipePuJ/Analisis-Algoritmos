
library(shiny)
library (shinydashboard) 
library(shinyjs)
library(PolynomF)
library(deSolve)
library(phaseR)
#datos 


dias<-100

p<-0.18







ui <- dashboardPage(
  dashboardHeader (), 
  dashboardSidebar (), 
  dashboardBody () 
)
header <- dashboardHeader(title = "reto3")


sidebar <- dashboardSidebar(
  
  sidebarMenu(
    
    menuItem("modeloSIRcovid19", tabName = "SIR", icon = icon("caret-right") )
  )
  
)

frow1 <- fluidRow(
  valueBoxOutput("v1")
  ,valueBoxOutput("v2")
  ,valueBoxOutput("v3")
)



TabSIR <- fluidPage(
  box(
    title = "Modelo SIR.", 
    width = 10,
    shinyjs::useShinyjs(),
    
    sliderInput("infectadosSIR", "Infectados:", 1,48000000, 48000000, step = 1000000),
    
   
    sliderInput("suceptiblesSIR", "Suceptibles:", 12,47978057 ,47978057, step = 1000000),
    
    sliderInput("diasSIR", "Dias:", 30, 100, 50, step = 1),
    
    sliderInput("aSIR", "a:", 0.1000, 0.2000, 0.1800, step = 0.0001),
    
    sliderInput("bSIR", "b:", 0.2500, 0.3500, 0.2200, step = 0.0001),
    actionButton("botonCalcularSIR", "Calcular")
  ),
  tabBox(
    title = "Metodos",
    width = 8,
    tabPanel("Adams",
             width = 400,
             plotOutput("plot2SIR", height = 380),
             plotOutput("plot2ErrorSIR", height = 400, width = 400),
             dataTableOutput("tablaErrorA2")
    ),
    tabPanel(
             "RK4",
             width = 400,
             plotOutput("plot1SIR", height = 380),
             plotOutput("plot1ErrorSIR", height = 400, width = 400),
             dataTableOutput("tablaErrorR2")
    )
  )
)


body <- dashboardBody(
  tabItems(
    tabItem("SIR",TabSIR)
  )
  
)

ui <- fluidPage(title = 'Reto 3', header, sidebar, body, theme = "bootstrap.css")

# Define server logic to show current time, update every second ----
server <- function(input, output, session) {

  
  CalcularSIR <- function(){
    
    
    output$plot1SIR <- renderPlot({
      N = 1
      init <- c(S = input$suceptiblesSIR,
                I = input$infectadosSIR,
                R = 4256)
      param <- c(a = input$aSIR,
                 b = input$bSIR)
      
      sir <- function(times, init, param) {
        with(as.list(c(init, param)), {
          dS <- (-a) * S * (I)
          dI <-  a * S * I - b * I
          dR <-  b * I
          return(list(c(dS, dI, dR)))
        })
      }
      times <- seq(0,input$diasSIR, by = 1)
      out <- ode(y = init, times = times, func = sir, parms = param, method = "rk4")
      out <- as.data.frame(out*N) 
      attach(out)
      
      plot(out$time, out$S, type="l", col="red", ylim=c(0,48000000), xlab="Dias", ylab="Infectados",main = "Runge Kutta 4")
      par(new=T)
      plot(out$time, out$I ,col="blue", type="l", xlab="Dias", ylab="Infectados", ylim = c(0,48000000))
      par(new=T)
      plot(out$time, out$R ,type = "l", col="black", xlab="Dias", ylab="Infectados", ylim = c(0, 48000000))
      legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("red", "blue", "black"), lty=rep(1, 2))
      
    })                                                                                    
    
    
    
    output$plot1ErrorSIR <- renderPlot({
      N = 1
      init <- c(S = input$suceptiblesSIR,
                I = input$infectadosSIR,
                R = 4256)
      param <- c(a = input$aSIR,
                 b = input$bSIR)
      
      sir <- function(times, init, param) {
        with(as.list(c(init, param)), {
          
          dS <- -a * S * (I) 
          dI <-  a * S * I - b * I
          dR <-                 b * I
          return(list(c(dS, dI, dR)))
        })
      }
      times <- seq(0,input$diasSIR, by = 1)
      out <- ode(y = init, times = times, func = sir, parms = param, method = "lsoda")
      out2 <- ode(y = init, times = times, func = sir, parms = param, method = "adams")
      out <- as.data.frame(out*N) 
      out2 <- as.data.frame(out2*N)
      attach(out)
      attach(out2)
      
      error <- c()
      i <- 1
      for(i in seq(1,input$diasSIR+1, by = 1)){
        error <- c(error, (abs(out$I[i] - out2$I[i]))/out$I[i])
      }
      plot(times, error, col = "purple", lwd = "2", type = "l", xlab = "Dias", ylab = "Error relativo", main = "Error", ylim = c(0, 10))
    })
    
    output$tablaErrorR2 <- renderDataTable({
      N = 1
      init <- c(S = input$suceptiblesSIR,
                I = input$infectadosSIR,
                R = 4256)
      param <- c(a = input$aSIR,
                 b = input$bSIR)
      
      sir <- function(times, init, param) {
        with(as.list(c(init, param)), {
          dS <- -a * S * (I) 
          dI <-  a * S * I - b * I
          dR <-                 b * I
          return(list(c(dS, dI, dR)))
        })
      }
      times <- seq(0,input$diasSIR, by = 1)
      out <- ode(y = init, times = times, func = sir, parms = param, method = "lsoda")
      out2 <- ode(y = init, times = times, func = sir, parms = param, method = "adams")
      out <- as.data.frame(out*N) 
      out2 <- as.data.frame(out2*N)
      attach(out)
      attach(out2)
      
      error <- c()
      i <- 1
      for(i in seq(1,input$diasSIR+1, by = 1)){
        error <- c(error, (abs(out$I[i] - out2$I[i]))/out$I[i])
      }
      
      valores <- data.frame(
        "Dia" = (1:(input$diasSIR + 1 )),
        "Real" = round(out$I, 4),
        "Aproximado" = round(out2$I, 4),
        "Error" = round(error,4)
      )
    })
    
    
    output$plot2SIR <- renderPlot({
      N = 1
      init <- c(S = input$suceptiblesSIR,
                I = input$infectadosSIR,
                R = 4256)
      param <- c(a = input$aSIR/100,
                 b = input$bSIR)
      
      sir <- function(times, init, param) {
        with(as.list(c(init, param)), {
          dS <- -a * S * (I) 
          dI <-  a * S * I - b * I
          dR <-                 b * I
          return(list(c(dS, dI, dR)))
        })
      }
      times <- seq(0,input$diasSIR, by = 1)
      out <- ode(y = init, times = times, func = sir, parms = param, method = "adams")
      out <- as.data.frame(out*N) 
      attach(out)
      
      plot(out$time, out$S, type="l", col="red", ylim=c(0,48000000), xlab="Dias", ylab="Infectados",main = "Adams")
      par(new=T)
      plot(out$time, out$I ,col="blue", type="l", xlab="Dias", ylab="Infectados", ylim = c(0,48000000))
      par(new=T)
      plot(out$time, out$R ,type = "l", col="black", xlab="Dias", ylab="Infectados", ylim = c(0, 48000000))
      legend(x = "topright", legend=c("Susceptibles", "Infectados", "Recuperados"), col=c("red", "blue", "black"), lty=rep(1, 2))
      
    })
    
    
    
    output$plot2ErrorSIR <- renderPlot({
      N = 1
      init <- c(S = input$suceptiblesSIR,
                I = input$infectadosSIR,
                R = 4256)
      param <- c(a = input$aSIR,
                 b = input$bSIR)
      
      sir <- function(times, init, param) {
        with(as.list(c(init, param)), {
          dS <- -a * S * (I)
          dI <-  a * S * I - b * I
          dR <-                 b * I
          return(list(c(dS, dI, dR)))
        })
      }
      times <- seq(0,input$diasSIR, by = 1)
      out <- ode(y = init, times = times, func = sir, parms = param, method = "adams")
      out <- as.data.frame(out*N) 
      attach(out)
      
      error <- c()
      i <- 1
      for(i in seq(1,input$diasSIR+1, by = 1)){
        error <- c(error, (abs(total[i] - out$I[i]))/total[i])
      }
      plot(times, error, col = "green", lwd = "2", type = "l", xlab = "Dias", ylab = "Error relativo", main = "Error", ylim = c(0, 10))
    })
    
    output$tablaErrorA2 <- renderDataTable({
      N = 1
      init <- c(S = input$suceptiblesSIR,
                I = input$infectadosSIR,
                R = 4256)
      param <- c(a = input$aSIR,
                 b = input$bSIR)
      
      sir <- function(times, init, param) {
        with(as.list(c(init, param)), {
          dS <- -a * S * (I)
          dI <-  a * S * I - b * I
          dR <-                 b * I
          return(list(c(dS, dI, dR)))
        })
      }
      times <- seq(0,input$diasSIR, by = 1)
      out <- ode(y = init, times = times, func = sir, parms = param, method = "adams")
      b <- lsoda
      out2 <- ode(y = init, times = times, func = sir, parms = param, method = "b")
      out <- as.data.frame(out*N) 
      out2 <- as.data.frame(out2*N)
      attach(out)
      attach(out2)
      
      error <- c()
      i <- 1
      for(i in seq(1,input$diasSIR+1, by = 1)){
        error <- c(error, (abs(out$I[i] - out2$I[i]))/out$I[i])
      }
      
      valores <- data.frame(
        "Dia" = (1:(input$diasSIR + 1 )),
        "Real" = round(out$I, 4),
        "Aproximado" = round(out2$I, 4),
        "Error" = round(error,4)
      )
    })
  }
  
  
 
  observeEvent(input$botonCalcularSIR, {
    CalcularSIR()
  })
}


# Create Shiny app ----
shinyApp(ui, server)