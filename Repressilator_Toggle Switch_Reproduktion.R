library(deSolve)

# REPRESSILATOR
# Definiere Repressilatorfkt mit Konzentration von mRNA und Protein
repressilator <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dM1 <- s * (-M1 + a/(1 + P2^m))
    dM2 <- s * (-M2 + a/(1 + P3^m))
    dM3 <- s * (-M3 + a/(1 + P1^m))
    dP1 <- s * (b * M1 - y * P1)
    dP2 <- s * (b * M2 - y * P2)
    dP3 <- s * (b * M3 - y * P3)
    return(list(c(dM1, dM2, dM3, dP1, dP2, dP3)))
  })}

# Setze Parameter in Repressilatorgleichungen ein
parameters <- c(a = 100, 
                m = 2, 
                b = 5, 
                y = 5, 
                s = 1)

# Startbedingungen
initial_state <- c(M1 = 0.5, 
                   M2 = 20.5, 
                   M3 = 20.5, 
                   P1 = 0.2, 
                   P2 = 0.45, 
                   P3 = 0.35)

# Zeitintervall
times <- seq(0, 50, by = 0.1)

################################################################################
# Löse gewöhnliche Differentialgleichung
output <- ode(y = initial_state, 
              times = times, 
              func = repressilator, 
              parms = parameters)

# Plot mRNA-Konzentration
matplot(x = output[, "time"], 
        y = output[, c("M1", "M2", "M3")], 
        type = "l", 
        lty = 1, 
        col = c("#009999", "#00FF33", "#FF66CC"), 
        xlab = "Zeit", 
        ylab = "mRNA-Konzentration", 
        main = "Repressilator")
legend("topright", 
       legend = c("M1", "M2", "M3"), 
       col = c("#009999", "#00FF33", "#FF66CC"), 
       lty = c(1, 5, 10))

# Plot Protein-Konzentration
matplot(x = output[, "time"], 
        y = output[, c("P1", "P2", "P3")], 
        type = "l", 
        lty = 1, 
        col = c("#009999", "#00FF33", "#FF66CC"),
        xlab = "Zeit", 
        ylab = "Protein-Konzentration", 
        main = "Repressilator")
legend("topright", 
       legend = c("P1", "P2", "P3"), 
       col = c("#009999", "#00FF33", "#FF66CC"), 
       lty = c(1, 5, 10))

################################################################################
 # TOGGLE SWITCH

 # Definiere Togglefkt
 toggleswitch <- function(time, state, parameters) {
   with(as.list(c(state, parameters)), {
     dX <- e/(1 + Y^n) - d * X
     dY <- e/(1 + X^n) - d * Y
     return(list(c(dX, dY)))
   })}

 # Setze Parameter in Toggle Switch-Gleichung ein
 parameters <- c(e = 2, 
                 n = 4, 
                 d = 1)

 # Startbedingungen
 initial_state_x <- c(X = 0.9, 
                      Y = 1)
 initial_state_y <- c(X = 1.1, 
                      Y = 1)

 # Zeitintervall
 times <- seq(0, 20, by = 0.1)

 # Löse gewöhnliche Differentialgleichung
 output_x <- ode(y = initial_state_x, 
                 times = times, 
                 func = toggleswitch, 
                 parms = parameters)
 output_y <- ode(y = initial_state_y, 
                 times = times, 
                 func = toggleswitch, 
                 parms = parameters)

 # Plot Konzentration von Gen X<Y
 matplot(x = output_x[, "time"],
         y = output_x[, c("X","Y")],
         type = "l",
         lty = 1,
         col = c("#009999","#00FF33"),
         xlab = "Zeit",
         ylab = "Konzentration von Gen X und Y",
         main = "Toggle Switch")
 legend("topright",
        legend = c("X", "Y"),
        col = c("#009999","#00FF33"),
        lty = c(1))

 # Plot Konzentration von Gen X>Y
 matplot(x = output_y[, "time"],
         y = output_y[, c("X","Y")],
         type = "l",
         lty = 1,
         col = c("#009999","#00FF33"),
         xlab = "Zeit",
         ylab = "Konzentration von Gen X und Y",
         main = "Toggle Switch")
 legend("topright",
        legend = c("X", "Y"),
        col = c("#009999","#00FF33"),
        lty = c(1))
 