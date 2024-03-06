
# Adapted from cowplot::theme_nothing
# TODO: Add documentation
theme_blank <- function() {
    ggplot2::theme_void() %+replace%
        theme(
              # Elements in this first block aren't used directly, but are inherited
              line = ggplot2::element_blank(),
              rect = ggplot2::element_rect(),
              text = ggplot2::element_blank(),
              aspect.ratio = 1,
              axis.line =          ggplot2::element_blank(),
              axis.line.x =        NULL,
              axis.line.y =        NULL,
              axis.text =          ggplot2::element_blank(),
              axis.text.x =        NULL,
              axis.text.x.top =    NULL,
              axis.text.y =        NULL,
              axis.text.y.right =  NULL,
              axis.ticks =         ggplot2::element_blank(),
              axis.ticks.length =  unit(0, "pt"),
              axis.title =         ggplot2::element_blank(),
              axis.title.x =       NULL,
              axis.title.x.top =   NULL,
              axis.title.y =       NULL,
              axis.title.y.right = NULL,

              legend.background =  ggplot2::element_blank(),
              legend.spacing =     NULL,
              legend.spacing.x =   NULL,
              legend.spacing.y =   NULL,
              legend.margin =      margin(0, 0, 0, 0),
              legend.key =         ggplot2::element_blank(),
              legend.key.size =   NULL,
              legend.key.height =  NULL,
              legend.key.width =   NULL,
              legend.text =        ggplot2::element_blank(),
              legend.text.align =  NULL,
              legend.title =       ggplot2::element_text(hjust = 0),
              legend.title.align = NULL,
              legend.position =    "none",
              legend.direction =   NULL,
              legend.justification = "center",
              legend.box =         NULL,
              legend.box.margin =  margin(0, 0, 0, 0),
              legend.box.background = ggplot2::element_blank(),
              legend.box.spacing = unit(0, "pt"),

              panel.grid =         ggplot2::element_blank(),
              panel.grid.major =   NULL,
              panel.grid.minor =   NULL,
              panel.spacing =      unit(0, "pt"),
              panel.spacing.x =    NULL,
              panel.spacing.y =    NULL,
              panel.ontop    =     FALSE,

              strip.background =   ggplot2::element_blank(),
              strip.text =         ggplot2::element_blank(),
              strip.text.x =       NULL,
              strip.text.y =       NULL,
              strip.placement =    "inside",
              strip.placement.x =  NULL,
              strip.placement.y =  NULL,
              strip.switch.pad.grid = unit(0., "cm"),
              strip.switch.pad.wrap = unit(0., "cm"),

              plot.background =    ggplot2::element_blank(),
              plot.title =         ggplot2::element_blank(),
              plot.subtitle =      ggplot2::element_blank(),
              plot.caption =       ggplot2::element_blank(),
              plot.tag           = ggplot2::element_blank(),
              plot.margin =        margin(0, 0, 0, 0),

              panel.background = ggplot2::element_rect(fill = "#ffffffcc",
                                                       colour = "#ffffffcc"),
              panel.border = ggplot2::element_rect(colour = "black",
                                                   fill = NA),
              complete = TRUE
        )

}




# Adapted from cowplot::theme_nothing
# TODO: Add documentation
theme_blank <- function() {
    ggplot2::theme_void() %+replace%
        theme(
              # Elements in this first block aren't used directly, but are inherited
              line = ggplot2::element_blank(),
              rect = ggplot2::element_rect(),
              text = ggplot2::element_blank(),
              aspect.ratio = 1,
              axis.line =          ggplot2::element_blank(),
              axis.line.x =        NULL,
              axis.line.y =        NULL,
              axis.text =          ggplot2::element_blank(),
              axis.text.x =        NULL,
              axis.text.x.top =    NULL,
              axis.text.y =        NULL,
              axis.text.y.right =  NULL,
              axis.ticks =         ggplot2::element_blank(),
              axis.ticks.length =  unit(0, "pt"),
              axis.title =         ggplot2::element_blank(),
              axis.title.x =       NULL,
              axis.title.x.top =   NULL,
              axis.title.y =       NULL,
              axis.title.y.right = NULL,

              legend.background =  ggplot2::element_blank(),
              legend.spacing =     NULL,
              legend.spacing.x =   NULL,
              legend.spacing.y =   NULL,
              legend.margin =      margin(0, 0, 0, 0),
              legend.key =         ggplot2::element_blank(),
              legend.key.size =   NULL,
              legend.key.height =  NULL,
              legend.key.width =   NULL,
              legend.text =        ggplot2::element_blank(),
              legend.text.align =  NULL,
              legend.title =       ggplot2::element_text(hjust = 0),
              legend.title.align = NULL,
              legend.position =    "none",
              legend.direction =   NULL,
              legend.justification = "center",
              legend.box =         NULL,
              legend.box.margin =  margin(0, 0, 0, 0),
              legend.box.background = ggplot2::element_blank(),
              legend.box.spacing = unit(0, "pt"),

              panel.grid =         ggplot2::element_blank(),
              panel.grid.major =   NULL,
              panel.grid.minor =   NULL,
              panel.spacing =      unit(0, "pt"),
              panel.spacing.x =    NULL,
              panel.spacing.y =    NULL,
              panel.ontop    =     FALSE,

              strip.background =   ggplot2::element_blank(),
              strip.text =         ggplot2::element_blank(),
              strip.text.x =       NULL,
              strip.text.y =       NULL,
              strip.placement =    "inside",
              strip.placement.x =  NULL,
              strip.placement.y =  NULL,
              strip.switch.pad.grid = unit(0., "cm"),
              strip.switch.pad.wrap = unit(0., "cm"),

              plot.background =    ggplot2::element_blank(),
              plot.title =         ggplot2::element_blank(),
              plot.subtitle =      ggplot2::element_blank(),
              plot.caption =       ggplot2::element_blank(),
              plot.tag           = ggplot2::element_blank(),
              plot.margin =        margin(0, 0, 0, 0),

              panel.background = ggplot2::element_rect(fill = "#ffffffcc",
                                                       colour = "#ffffffcc"),
              panel.border = ggplot2::element_rect(colour = "black",
                                                   fill = NA),
              complete = TRUE
        )

}


# FIXME: Improve color palette
# TODO: Add documentation
scale_color_mpimg <- function(name = "mpimg", ...) {

  mpi_colors <- c(
    "#006c66", # MPG-CD-Grün
    "#777777", # MPG-Dunkelgrau
    "#a7a7a8", # MPG-Grau
    "#c6d325", # MPG Hellgrün
    "#29485d", # MPG Dunkelblau
    "#00b1ea", # MPG Hellblau
    "#ef7c00" # MPG Orange
  )
  # mpi colors extended.
  mpimg_colors <- c(
    "#006C66",
    "#29485D",
    "#009ACD",
    "#777777",
    "#C6D325",
    "#21AE2D",
    "#EF7C00",
    "#F5B742",
    "#57219D",
    "#950980"
  )

  if (name == "mpimg") {

    ggplot2::discrete_scale(
      scale_name = "mpimg",
      aesthetics = "color",
      palette = scales::manual_pal(values = mpimg_colors),
      ...
    )

  } else if (name == "mpi") {

    ggplot2::discrete_scale(
      scale_name = "mpi",
      aesthetics = "color",
      palette = scales::manual_pal(values = mpi_colors),
      ...
    )

  }
}

