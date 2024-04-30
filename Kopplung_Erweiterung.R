library(deSolve)

################################################################################
# KOPPLUNG 1: Toggle Switch kontrolliert Repressilator
repressilator_unter_toggle_switch <- function(time, state, parameters, signal) {
  ff1 <- signal(time)
  with(as.list(c(state, parameters)), {
    # Toggle Switch
    dX1 <- e/(1 + Y1^n) - d * X1 + ff1
    dY1 <- e/(1 + X1^n) - d * Y1 + f2
    
    # Definition von a1 als Funktion von X (Hill-Funktion)
    a1 <- A * (X1^k / (K^k + X1^k)) + B
    
    # Repressilator
    dM11 <- s * (-M11 + a1/(1 + P21^m))
    dM21 <- s * (-M21 + a/(1 + P31^m))
    dM31 <- s * (-M31 + a/(1 + P11^m))
    dP11 <- s * (b * M11 - y * P11)
    dP21 <- s * (b * M21 - y * P21)
    dP31 <- s * (b * M31 - y * P31)
    
    return(list(c(dX1, dY1, dM11, dM21, dM31, dP11, dP21, dP31), f1 = ff1))
  })
}

# Wechsel des Zustands von X als Signal 
signal <- approxfun(x = c(0, 10, 50), 
                    y = c(0, 1, 1), 
                    method = "constant", rule = 2)


# Setze Parameter für Kopplung mit a1
kopplungs_par <- c(e = 2, 
                   n = 4, 
                   d = 1, 
                   A = 100, 
                   B = 5, 
                   a = 100, 
                   m = 2, 
                   b = 5, 
                   y = 5, 
                   s = 1, 
                   K = 0.5, 
                   f1 = 0, 
                   f2 = 0, 
                   k = 3)



# Setze den Startzustand
initial_state_kopplung1 <- c(X1 = 0.1, 
                             Y1 = 2, 
                             M11 = 0.5, 
                             M21 = 20.5, 
                             M31 = 20.5, 
                             P11 = 0.2, 
                             P21 = 0.45, 
                             P31 = 0.35)

# Löse DGLS für Kopplung 1
kopplung1 <- ode(y = initial_state_kopplung1, 
                 times = seq(0, 50, by = 0.1), 
                 func = repressilator_unter_toggle_switch, 
                 parms = kopplungs_par, 
                 signal = signal)

matplot(x = kopplung1[, "time"], 
        y = kopplung1[, c("X1", "Y1")], 
        type = "l", 
        lty = 1, 
        col = c("#009999", "#00FF33", "#FF66CC"), 
        xlab = "Zeit", 
        ylab = "Konzentration von Gen X und Y", 
        main = "Erweiterung von Kopplung 1")
legend("topright", 
       legend = c("X1", "Y1"), 
       col = c("#009999", "#00FF33", "#FF66CC"), 
       lty = c(1, 5, 10))

matplot(x = kopplung1[, "time"], 
        y = kopplung1[, c("M11", "M21", "M31")], 
        type = "l", 
        lty = 1, 
        col = c("#009999", "#00FF33", "#FF66CC"), 
        xlab = "Zeit", 
        ylab = "mRNA-Konzentration", 
        main = "Erweiterung von Kopplung 1")
legend("topright", 
       legend = c("M11", "M21", "M31"), 
       col = c("#009999", "#00FF33", "#FF66CC"), 
       lty = c(1, 5, 10))
################################################################################
# KOPPLUNG 2: Repressilator kontrolliert Toggle Switch
toggle_switch_unter_repressilator <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Repressilator
    dM11 <- s * (-M11 + a/(1 + P21^m))
    dM21 <- s * (-M21 + a/(1 + P31^m))
    dM31 <- s * (-M31 + a/(1 + P10^m))
    dP10 <- s * (b * M11 - y * P10)
    dP21 <- s * (b * M21 - y * P21)
    dP31 <- s * (b * M31 - y * P31)
    
    # Definition von e1 als Funktion von P10 (Hill-Funktion)
    e1 <- A1 * (P10^l / (K^l + P10^l)) + B1
    
    # Toggle Switch
    dX1 <- e1/(1 + Y1^n) - d * X1 
    dY1 <- e/(1 + X1^n) - d * Y1
    
    return(list(c(dX1, dY1, dM11, dM21, dM31, dP10, dP21, dP31)))
  })
}

# Setze Parameter für Kopplung mit a1
kopplungs_par2 <- c(e = 2, 
                    n = 3, 
                    d = 3, 
                    A1 = 10, 
                    B1 = 0.25, 
                    a = 220, 
                    m = 2, 
                    b = 5, 
                    y = 5, 
                    s = 1, 
                    K = 60, 
                    l = 4)

# Setze den Startzustand
initial_state_kopplung2 <- c(X1 = 2, 
                             Y1 = 0.1,
                             M11 = 0.5, 
                             M21 = 40.5, 
                             M31 = 40.5, 
                             P10 = 400, 
                             P21 = 500, 
                             P31 = 400)

# Löse DGLS für Kopplung 2
kopplung2 <- ode(y = initial_state_kopplung2, 
                 times = seq(0, 50, by = 0.1), 
                 func = toggle_switch_unter_repressilator, 
                 parms = kopplungs_par2)


matplot(x = kopplung2[, "time"], 
        y = kopplung2[, c("M11", "M21", "M31")], 
        type = "l", 
        lty = 1, 
        col = c("#009999", "#00FF33", "#FF66CC"), 
        xlab = "Zeit", 
        ylab = "mRNA-Konzentration", 
        main = "Erweiterung von Kopplung 2")
legend("topright", 
       legend = c("M11", "M21", "M31"), 
       col = c("#009999", "#00FF33", "#FF66CC"), 
       lty = c(1, 5, 10))

matplot(x = kopplung2[, "time"], 
        y = kopplung2[, c("X1", "Y1")], 
        type = "l", 
        lty = 1, 
        col = c("#009999", "#00FF33", "#FF66CC"), 
        xlab = "Zeit", 
        ylab = "Konzentration von Gen X und Y", 
        main = "Erweiterung von Kopplung 2")
legend("topright", 
       legend = c("X1", "Y1"), 
       col = c("#009999", "#00FF33", "#FF66CC"), 
       lty = c(1, 5, 10))

