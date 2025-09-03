
# Adapted from cowplot::theme_nothing
#' @importFrom ggplot2 '%+replace%'
# ggplot2 theme that strips all elements from a plot.
theme_blank <- function(title = ggplot2::element_blank(),
                        text = ggplot2::element_blank()) {
    ggplot2::theme_void() %+replace%
        ggplot2::theme(
            # Elements in this first block aren't used directly, but are inherited
            line = ggplot2::element_blank(),
            rect = ggplot2::element_rect(),
            text = text,
            aspect.ratio = 1,
            axis.line = ggplot2::element_blank(),
            axis.line.x = NULL,
            axis.line.y = NULL,
            axis.text = ggplot2::element_blank(),
            axis.text.x = NULL,
            axis.text.x.top = NULL,
            axis.text.y = NULL,
            axis.text.y.right = NULL,
            axis.ticks = ggplot2::element_blank(),
            axis.ticks.length = ggplot2::unit(0, "pt"),
            axis.title = ggplot2::element_blank(),
            axis.title.x = NULL,
            axis.title.x.top = NULL,
            axis.title.y = NULL,
            axis.title.y.right = NULL,
            legend.background = ggplot2::element_blank(),
            legend.spacing = NULL,
            legend.spacing.x = NULL,
            legend.spacing.y = NULL,
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.key = ggplot2::element_blank(),
            legend.key.size = NULL,
            legend.key.height = NULL,
            legend.key.width = NULL,
            legend.text = ggplot2::element_blank(),
            legend.text.align = NULL,
            legend.title = ggplot2::element_text(hjust = 0),
            legend.title.align = NULL,
            legend.position = "none",
            legend.direction = NULL,
            legend.justification = "center",
            legend.box = NULL,
            legend.box.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.background = ggplot2::element_blank(),
            legend.box.spacing = ggplot2::unit(0, "pt"),
            panel.grid = ggplot2::element_blank(),
            panel.grid.major = NULL,
            panel.grid.minor = NULL,
            panel.spacing = ggplot2::unit(0, "pt"),
            panel.spacing.x = NULL,
            panel.spacing.y = NULL,
            panel.ontop = FALSE,
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_blank(),
            strip.text.x = NULL,
            strip.text.y = NULL,
            strip.placement = "inside",
            strip.placement.x = NULL,
            strip.placement.y = NULL,
            strip.switch.pad.grid = ggplot2::unit(0., "cm"),
            strip.switch.pad.wrap = ggplot2::unit(0., "cm"),
            plot.background = ggplot2::element_blank(),
            plot.title = title,
            plot.subtitle = ggplot2::element_blank(),
            plot.caption = ggplot2::element_blank(),
            plot.tag = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0),
            panel.background = ggplot2::element_rect(
                fill = "#ffffffcc",
                colour = "#ffffffcc"
            ),
            panel.border = ggplot2::element_rect(
                colour = "black",
                fill = NA
            ),
            complete = TRUE
        )
}

#' @importFrom ggplot2 '%+replace%'
# ggplot2 theme that strips all elements from a plot, but leaves axis ticks and elements.
theme_axis_only <- function(title = ggplot2::element_blank(),
                        text = ggplot2::element_blank()) {
    ggplot2::theme_classic() %+replace%
        ggplot2::theme(
            # Elements in this first block aren't used directly, but are inherited
            # line = ggplot2::element_blank(),
            rect = ggplot2::element_rect(),
            text = text,
            aspect.ratio = 1,
            axis.line = ggplot2::element_blank(),
            axis.line.x = NULL,
            axis.line.y = NULL,
            axis.title = ggplot2::element_blank(),
            axis.title.x = NULL,
            axis.title.x.top = NULL,
            axis.title.y = NULL,
            axis.title.y.right = NULL,
            legend.background = ggplot2::element_blank(),
            legend.spacing = NULL,
            legend.spacing.x = NULL,
            legend.spacing.y = NULL,
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.key = ggplot2::element_blank(),
            legend.key.size = NULL,
            legend.key.height = NULL,
            legend.key.width = NULL,
            legend.text = ggplot2::element_blank(),
            legend.text.align = NULL,
            legend.title = ggplot2::element_text(hjust = 0),
            legend.title.align = NULL,
            legend.position = "none",
            legend.direction = NULL,
            legend.justification = "center",
            legend.box = NULL,
            legend.box.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.background = ggplot2::element_blank(),
            legend.box.spacing = ggplot2::unit(0, "pt"),
            panel.grid = ggplot2::element_blank(),
            panel.grid.major = NULL,
            panel.grid.minor = NULL,
            panel.spacing = ggplot2::unit(0, "pt"),
            panel.spacing.x = NULL,
            panel.spacing.y = NULL,
            panel.ontop = FALSE,
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_blank(),
            strip.text.x = NULL,
            strip.text.y = NULL,
            strip.placement = "inside",
            strip.placement.x = NULL,
            strip.placement.y = NULL,
            strip.switch.pad.grid = ggplot2::unit(0., "cm"),
            strip.switch.pad.wrap = ggplot2::unit(0., "cm"),
            plot.background = ggplot2::element_blank(),
            plot.title = title,
            plot.subtitle = ggplot2::element_blank(),
            plot.caption = ggplot2::element_blank(),
            plot.tag = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0),
            panel.background = ggplot2::element_rect(
                fill = "#ffffffcc",
                colour = "#ffffffcc"
            ),
            panel.border = ggplot2::element_rect(
                colour = "black",
                fill = NA
            ),
            complete = TRUE
        )
}


#' ggplot2 scale for the MPI colors and extended swatches.
#' @param name The name of the scale. Either "mpi" or "mpimg".
#' @param ... Further arguments to ggplot2::discrete_scale.
#' @export
scale_color_mpimg <- function(name = "mpimg", ...) {
    if (name == "mpimg") {
        ggplot2::discrete_scale(
            scale_name = "mpimg",
            aesthetics = "color",
            palette = mpimg_pal(),
            ...
        )
    } else if (name == "mpi_extend") {
        ggplot2::discrete_scale(
            scale_name = "mpi_extend",
            aesthetics = "color",
            palette = mpi_extend_pal(),
            ...
        )
    }
}

#' ggplot2 scale for the MPI colors and extended swatches.
#' @param name The name of the scale. Either "mpi" or "mpimg".
#' @param ... Further arguments to ggplot2::discrete_scale.
#' @export
scale_fill_mpimg <- function(name = "mpimg", ...) {
    if (name == "mpimg") {
        ggplot2::discrete_scale(
            scale_name = "mpimg",
            aesthetics = "fill",
            palette = mpimg_pal(),
            ...
        )
    } else if (name == "mpi_extend") {
        ggplot2::discrete_scale(
            scale_name = "mpi_extend",
            aesthetics = "fill",
            palette = mpi_extend_pal(),
            ...
        )
    }
}


scale_fill_gradient_mpimg <- function(name = "orange", ...) {
    if (name == "orange") {
        ggplot2::scale_fill_gradientn(
            colours = c("#29485d", "#006c66", "#c6d325", "#ef7c00"),
            ...
        )
    } else if (name == "green") {
        ggplot2::scale_fill_gradientn(
            colours = c("#29485d", "#006c66", "#c6d325"),
            ...
        )
    } else {
        rlang::abort("Pick a correct name!")
    }
}

#' MPIMG color palette.
#' @returns A function that can be used to generate colors.
#' @export
mpimg_pal <- function() {
    mpi_colors <- c(
        "#006c66", # MPG-CD-Grün
        "#c6d325", # MPG Hellgrün
        "#ef7c00", # MPG Orange
        "#29485d", # MPG Dunkelblau
        "#00b1ea", # MPG Hellblau
        "#777777", # MPG-Dunkelgrau
        "#a7a7a8" # MPG-Grau
    )

    scales::manual_pal(values = mpi_colors)
}

#' MPI color palette
#' @returns A function that can be used to generate colors.
#' @export
mpi_extend_pal <- function() {
    # mpi colors extended.
    mpi_extend_colors <- c(
        "#006c66", # MPG-CD-Grün
        "#c6d325", # MPG Hellgrün
        "#ef7c00", # MPG Orange
        "#29485d", # MPG Dunkelblau
        "#00b1ea", # MPG Hellblau
        "#777777", # MPG-Dunkelgrau
        "#d44a3d", # Warm Red
        "#8c5fa8", # Soft Purple
        "#a7a7a8", # MPG-Grau
        "#f4c542", # Golden Yellow
        "#d95276", # Rich Pink
        "#a25b43", # Warm Brown
        "#00CED1" # Soft Cyan
    )

    scales::manual_pal(values = mpi_extend_colors)
}

#' @export
scale_colour_rnd <- function(...) {
    require(randomcoloR)

    ggplot2::discrete_scale(
        scale_name = "mpi_extend",
        aesthetics = "color",
        palette = randomcoloR::distinctColorPalette,
        ...
    )
}
