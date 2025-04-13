#' Enhanced quality control plots for daily precipitation

#' @description The function create a matrix plot of the enhanced quality control tests for daily precipitation
#' @param xts_obj xts of a single time series  
#' @export
#' @importFrom xts xts
#' @importFrom lattice xyplot barchart panel.barchart panel.text panel.abline
#' @importFrom grid textGrob gpar grid.draw
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#' @details
#' Six plots are built based on daily precipitation: time series, truncation, time series (threshold = 5 mm), small gaps, precision and rounding patterns and weekly cycle.
#' These provide a visual inspection (Hunziker et al., 2017) but also how the automatic application  (Huerta et al., 2020) of the enhanced quality control is applied. 
#' @references Hunziker, S., Gubler, S., Calle, J., Moreno, I., Andrade, M., Velarde, F., ... & Brönnimann, S. (2017). Identifying, attributing, and overcoming common data quality issues of manned station observations. https://doi.org/10.1002/joc.5037
#' @references Huerta, A., Serrano-Notivoli, R., & Brönnimann, S. (2024). SC-PREC4SA: A serially complete daily precipitation dataset for South America. https://doi.org/10.31223/X57D8R
#' @examples
#' \dontrun{
#' set.seed(123)
#' 
#' # creating fake daily precipitation data
#' dates_t <- seq(as.Date("1980-01-01"), as.Date("2015-12-31"), by = "day")
#' prec <- round(rnorm(length(dates_t), mean = 1.2, sd = 6), 1)
#' prec[prec<0] <- 0
#' xts_obj <- xts::xts(prec, dates_t)
#' names(xts_obj) <- "prec"
#' 
#' # enhanced qc plots
#' eqc_Plot(xts_obj)
#' 
#' # it also work if there is some empty data (but not if all is NA)
#' xts_obj["1990/2010"] <- NA
#' eqc_Plot(xts_obj)
#' }
#' 

eqc_Plot <- function(xts_obj) {
  
  variable <- NULL
  bin_test <- NULL
  dec <- NULL
  
  if (is.null(xts_obj) || all(is.na(xts_obj))) {
    stop("Error: The xts object is NULL or completely NA.")
  }
  
  ######################## theme for plots ########################
  
  theme.novpadding <-
    list(par.main.text = list(font = 1, just = "center"),
         layout.heights =
           list(main = .1,
                top.padding = 0,
                main.key.padding = 0,
                key.axis.padding = 0,
                axis.xlab.padding = -1,
                xlab.key.padding = 0,
                key.sub.padding = 0,
                bottom.padding = 3),
         layout.widths =
           list(left.padding = 0,
                key.ylab.padding = 0,
                ylab.axis.padding = 0,
                axis.key.padding = 0,
                right.padding = 1))
  
  theme.novpadding2 <-
    list(par.main.text = list(font = 1, just = "center"),
         layout.heights =
           list(main = 0.1,
                top.padding = 0,
                main.key.padding = 0,
                key.axis.padding = 0,
                axis.xlab.padding = -1,
                xlab.key.padding = 0,
                key.sub.padding = 0,
                bottom.padding = 2),
         layout.widths =
           list(left.padding = 0,
                key.ylab.padding = 0,
                ylab.axis.padding = 0,
                axis.key.padding = 0,
                right.padding = -1))
  
  
  ##################### enhanced qc plots ########################
  
  # simple time series plot
  
  xts_ts_plt <- lattice::xyplot(
    xts_obj, type = "p", cex = .1, pch = 19, col = "black",
    xlab = "", ylab = "(mm)", main = "time series",
    xlim = c(stats::time(xts_obj)[1] - 1*365,
             stats::time(xts_obj)[length(stats::time(xts_obj))] + 1*365),
    ylim = c(min(xts_obj, na.rm = TRUE), max(xts_obj, na.rm = TRUE)),
    par.settings = theme.novpadding,
    scales = list(y = list(rot = 90))
  )
  
  # truncation shape plot
  
  ## building borderline of time series
  xts_trunc <- get_ts_borderline(xts_obj = xts_obj)
  
  xts_trunc <- lattice::xyplot(
    xts_trunc, type = "l", cex = .1, pch = 19, col = "red", lwd = 3,
    xlab = "", ylab = "(mm)", main = "truncation",
    xlim = c(stats::time(xts_trunc)[1] - 1*365,
             stats::time(xts_trunc)[length(stats::time(xts_trunc))] + 1*365),
    par.settings = theme.novpadding,
    scales = list(y = list(rot = 90))
  )
  
  # simple time series plot with threshold
  
  xts_ts_thrs_plt <- lattice::xyplot(
    xts_obj, type = "p", cex = .1, pch = 19, ylim = c(0, 5), col = "black",
    xlab = "", ylab = "(mm)", main = "time series (threshold = 5 mm)",
    xlim = c(stats::time(xts_obj)[1] - 1*365,
             stats::time(xts_obj)[length(stats::time(xts_obj))] + 1*365),
    par.settings = theme.novpadding,
    scales = list(y = list(rot = 90))
  )
  
  
  # small gaps plot

  ## counting the number of values within ranges
  xts_small_gaps_plt <- get_ngaps(xts_obj = xts_obj, lmn_yday = 365 * 80 / 100)
  xts_small_gaps_plt <- xts::xts(
    xts_small_gaps_plt[, 
                   -match(c("year", "year_size"), colnames(xts_small_gaps_plt))
                   ],
    as.Date(
      paste(xts_small_gaps_plt$year, "-01-01", sep = "")
      )
  )
  
  ymax <- max(xts_small_gaps_plt) * 1.05
  ymin <- min(xts_small_gaps_plt) * 1.05
  
  if(ymax == ymin){
    
    ymax <- 1 * 1.05
    ymin <- 0 * 1.05
  }
  
  xts_small_gaps_plt <- data.frame(Date = stats::time(xts_small_gaps_plt),
                                   coredata(xts_small_gaps_plt))
  xts_small_gaps_plt <- reshape2::melt(xts_small_gaps_plt, id.vars = "Date")
  
  xts_small_gaps_plt <- lattice::xyplot(
    value ~ Date, data = xts_small_gaps_plt,
    groups = variable, ylim = c(ymin, ymax),
    xlab = "", ylab = "number of values", main = "small gaps",
    type = "l", lwd = 3,
    col = c(1, 2, 3, 4, 5),
    par.settings = c(
      theme.novpadding2,
      list(
        layout.widths = list(
          right.padding = 0,
          axis.right = 0
        )
      )
    ),
    scales = list(y = list(rot = 90, alternating = 1, tck = c(1, 0)))
  )
  
  
  # weekly cycle plots

  if(dim(xts_obj)[1] < 1) {
    
    ## handling NA time series object
    xts_wd_df <- data.frame(year = NA, dec = NA)
    
  } else {
    
    ## computing binomial test
    xts_wd_df <- get_wd_fraction(xts_obj)
    xts_wd_df$frac_wd <- xts_wd_df$frac_wd / 100
    
  }
  
  xts_wd_plt <- lattice::barchart(
    frac_wd ~ week,
    data = xts_wd_df,
    group = bin_test,
    xlab = " ",
    ylab = "wet day fraction (>= 1 mm) ",
    main = "weekly cycle",
    ylim = c(0, sum(xts_wd_df$count_wd) / sum(xts_wd_df$count) + 0.06),
    auto.key = list(space = "bottom", columns = 2,
                    padding.text = 0, between.columns = 0,
                    text.width = 4, pch = 20, between = .5, cex = .75, size = 1,
                    rectangles = TRUE, points = FALSE),
    scales = list(y = list(rot = 90, tck = c(1, 0))),
    par.settings = c(
      theme.novpadding2,
      list(superpose.polygon = list(col = c("red", "gray50")))
      ),
    panel = function(...) {
      lattice::panel.barchart(...)
      args <- list(...)
      lattice::panel.text(
        args$x, args$y + 0.02,
        xts_wd_df$count,
        pos = 4, offset = -1, alpha = .5, cex = 1
      )
      lattice::panel.abline(
        h = sum(xts_wd_df$count_wd) / sum(xts_wd_df$count),
        lty = 2, lwd = 1, col = "black"
      )
    }
  )
  
  # precision and rounding patterns plot
  
  ## counting number of decimals
  xts_ts_dec_plt <- get_ndec(xts_obj = xts_obj, lmn_yday = NA, toPlot = TRUE)
  xts_ts_dec_plt <- reshape2::melt(table(xts_ts_dec_plt))
  xts_ts_dec_plt$dec <- factor(
    xts_ts_dec_plt$dec,
    levels = c("x.0", "x.1", "x.2", "x.3", "x.4",
               "x.5", "x.6", "x.7", "x.8", "x.9")
  )
  
  ymax_n <- max(seq_along(xts_ts_dec_plt$year))
  ymin_n <- min(seq_along(xts_ts_dec_plt$year))
  ymax <- max(xts_ts_dec_plt$year)
  ymin <- min(xts_ts_dec_plt$year)
  
  xts_ts_dec_plt <- lattice::barchart(
    value ~ year, data = xts_ts_dec_plt,
    group = dec, stack = TRUE, horizontal = FALSE,
    xlab = " ", ylab = "frequency (days/year)",
    main = "precision and rounding patterns",
    lwd = .1,
    auto.key = list(space = 'bottom', columns = 10,
                    padding.text = 0, between.columns = 0, text.width = 4,
                    pch = 20, between = .5, cex = .75,
                    size = 1, lwd = .1),
    scales = list(y = list(rot = 90, tck = c(1, 0))),
    par.settings = c(
      theme.novpadding2,
      list(superpose.polygon = list(
        col = c("black", "yellow", "orange", "red", "darkslateblue",
                "darkgray", "magenta", "blue", "cyan", "darkgreen"),
        lwd = .1)
        )
    ),
    axis = function(side, ...) {
      if (side == "bottom")
        lattice::panel.axis(
          at = seq(
            ymin_n,
            ymax_n,
            by = 10
          ),
          label = seq(
            ymin,
            ymax,
            by = 10
          ),
          outside = TRUE, rot = 0, tck = 0)
      else
        lattice::axis.default(side, ...)
    }
  )
  
  
  ######################## all plots ########################
  
  res <- gridExtra::grid.arrange(
    xts_ts_plt, xts_trunc,
    xts_ts_thrs_plt, xts_small_gaps_plt,
    xts_ts_dec_plt, xts_wd_plt,
    ncol = 2,
    top = grid::textGrob(
      names(xts_obj),
      gp = grid::gpar(fontsize = 17)
    )
  )
  
  invisible(grid::grid.draw(res))
  
}