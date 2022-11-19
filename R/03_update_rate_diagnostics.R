#' Generates update rate plots
#'
#' If latent states are provided only!
#'
#' @param states latent state values
#' @param settings_urs settings used for viewing/saving update rate plots
#' @param settings_plots general plot settings
#'
#' @return side effect function; plotting and/or saving
#' @export
get_update_rates <- function(states, settings_urs, settings_plots) {
  if (settings_urs$ur_view) {
    graphics::par(mfrow = c(1, 1))
    analyse_states_ur(trajectories = states)
  }
  if (settings_urs$ur_save) {
    current_plot_name <- file.path(settings_plots$plot_path,
                                   paste0("00_", settings_urs$ur_name, ".pdf"))
    grDevices::setEPS()
    grDevices::postscript(current_plot_name, width = 9, height = 5.25)
    analyse_states_ur(trajectories = states)
    grDevices::dev.off()
    print(paste("Saved update rate plot in: ", current_plot_name))
  }
}
