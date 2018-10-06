# Set up the mar
par(mar = c(3, 3.1, 3, 1))

# Set the plot area
plot(0, 0, type = "n",
     xlim = c(-100, 100), ylim = c(-75, 150),
     xaxs = "i", yaxs = "i", 
     axes = F, bty = "n", ann = F
    )

# Enter in the gray rectangle around 0
rect(xleft = -100, xright = 100, ybottom = -4, ytop = 4, col = "gray")

# Enter in the forward read distribution
x_fwd = seq(-75, 25, 0.1)
y_fwd = dnorm(x_fwd, mean = -25, sd = 15)
lines(x_fwd, (y_fwd * 2000) + 4, col = "red", lwd = 2)

x_rev = seq(-25, 75, 0.1)
y_rev = dnorm(x_rev, mean = 25, sd = 15)
lines(x_rev, (y_rev * 2000) + 4, col = "darkgreen", lwd = 2)

# Make the ORC read
rect(xleft = -25, xright = 25, ybottom = 81, ytop = 89, col = "#000099")
rect(xleft = -10, xright = 10, ybottom = 75, ytop = 95, col = "darkgreen")
# plot_cdc6(x0 = 0, x_w = 15, y0 = 85, yh = 8, obj_col = "darkgreen")
text(x = 0, y = 85, labels = "ORC", col = "white", cex = 0.8)

# Insert the read sequence
segments(x0 = -25, y0 = 89, x1 = -25, y1 = 102.5, lwd = 2)
arrows(x0 = -25, y0 = 102.5, x1 = -10, y1 = 102.5, lwd = 2, angle = 35, length = 0.125)
arrows(x0 = -25, y0 = 80, x1 = -25, y1 = 60, col = "red", lwd = 2, angle = 35, length = 0.125)

# Add the vertical segments
segments(x0 = -25, y0 = 57, x1 = -25, y1 = 4, col = "red", lty = 2, lwd = 2)
segments(x0 = 25, y0 = 57, x1 = 25, y1 = 4, col = "darkgreen", lty = 2, lwd = 2)

# Make the ORC read
rect(xleft = -25, xright = 25, ybottom = -34, ytop = -26, col = "#000099")
rect(xleft = -10, xright = 10, ybottom = -40, ytop = -20, col = "darkgreen")
# plot_cdc6(x0 = 0, x_w = 15, y0 = -30, yh = 8, obj_col = "darkgreen")
text(x = 0, y = -30, labels = "ORC", col = "white", cex = 0.8)

# Insert the read sequence
segments(x0 = 25, y0 = -34, x1 = 25, y1 = -47.5, lwd = 2)
arrows(x0 = 25, y0 = -47.5, x1 = 10, y1 = -47.5, lwd = 2, angle = 35, length = 0.125)
arrows(x0 = 25, y0 = -25, x1 = 25, y1 = -5, col = "darkgreen", lwd = 2, angle = 35, length = 0.125)

# Add in the legend
legend("topleft", legend = "Fwd strand\ndensity", col = "red", lwd = 2, bty = "n")
legend("topright", legend = "Rev strand\ndensity", col = "darkgreen", lwd = 2, bty = "n")
