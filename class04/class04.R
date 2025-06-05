#' ---
#' title: "My First R Script"
#' author: "Christopher Leone (A16731724)"
#' ---

# My first R script
x <- 1:50
plot(x, sin(x))

plot(x, sin(x), typ="c", col="blue"
     , lwd="3", xlab="Silly X Axis"
     , ylab="Sensible Y Axis")
