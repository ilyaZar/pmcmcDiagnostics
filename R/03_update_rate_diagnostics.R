#' Generates update rate plots
#'
#' If latent states are provided only!
#'
#' @param states latent state values
#' @param settings_urs settings used for viewing/saving update rate plots
#'
#' @return side effect function; plotting and/or saving
#' @export
get_update_rates <- function(states, settings_urs, settings_plots) {
  urs <- compute_states_ur(trajectories = states, states_in_cols = TRUE)
  if (settings_urs$ur_view) {
    graphics::par(mfrow = c(1, 1))
    graphics::matplot(urs , type = "l")
  }
  if (settings_urs$ur_save) {
    browser()
    current_plot_name <- file.path(settings_urs$ur_path,
                                   paste0("00_", settings_urs$ur_name, ".pdf"))
    grDevices::setEPS()
    grDevices::postscript(current_plot_name, width = 9, height = 5.25)
    graphics::matplot(urs , type = "l")
    grDevices::dev.off()
    print(paste("Saved update rate plot in: ", current_plot_name))
  }
  return(urs)
}
