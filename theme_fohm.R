#' Fohm ggplot theme and colour palette
#'
#' A theme and colour palette (almost) adhereing to the guidelines in
#' \url{http://intranet.folkhalsomyndigheten.se/arbetssatt-och-stod/kommunicera/grafisk-manual/diagram/}
#'
#' Updated: 2020-11-03
#'
#' @examples
#' \dontrun{
#'   trig_data <- data.frame(x = rep(grid_x, 2), value = c(sin(grid_x),
#'   y_cos = cos(grid_x)), group=rep(c("sin", "cos"), each=length(grid_x)))
#'
#'   ggplot(trig_data, aes(x=x, y=value, group=group, col=group)) +
#'     geom_line(alpha = 0.75, size=1.25) +
#'     scale_colour_manual(values=fohm_colours(), labels=c("sin", "cos")) +
#'     labs(y = "", subtitle = "f(x)") +
#'     theme_fohm()
#' }
#'
#' @name theme_fohm
NULL
#> NULL

#' @param text_size Base text size, in pts.
#'
#' @export
#' @rdname theme_fohm
theme_fohm <- function(text_size=11) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  ggplot2::`%+replace%`(ggplot2::theme_bw(),
                        ggplot2::theme(
                          text = ggplot2::element_text(size = text_size),
                          panel.grid.major.x = ggplot2::element_blank(),
                          panel.grid.minor.x = ggplot2::element_blank(),
                          panel.border = ggplot2::element_blank(),
                          axis.line = ggplot2::element_line(colour = "black"),
                          axis.title.x = ggplot2::element_text(size = text_size + 1, vjust=-1),
                          axis.text.x = ggplot2::element_text(angle = 0, hjust = 1, size=text_size),
                          axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size=text_size),
                          plot.title.position = "plot",
                          plot.title = ggplot2::element_text(size = text_size + 2),
                          plot.subtitle = ggplot2::element_text(vjust=-9, size = text_size + 1),
                          legend.title = ggplot2::element_blank(),
                          legend.position="top",
                          legend.text=ggplot2::element_text(size = text_size - 1),
                          plot.background=ggplot2::element_rect(color="#BFBFBF", size=1),
                          complete = TRUE
                        ))
}


#' @export
#' @rdname theme_fohm
fohm_colours<- function(){
  c("#009F80", "#82368C", "#ED7003", "#1E3C90",
    "#E4609F", "#006D68", "#1C98D5", "#757575", "#CC1423", "#A86DAB", "#942503")
}

#' @export
#' @rdname theme_fohm
fohm_colors <- fohm_colours

