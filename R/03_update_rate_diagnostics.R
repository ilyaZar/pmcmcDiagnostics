#' Generates update rate plots
#'
#' If latent states are provided only!
#'
#' @param states latent state values
#' @inheritParams analyse_states_convergence
#'
#' @return side effect function; plotting and/or saving
#' @export
get_update_rates <- function(states, settings_urs, dim_list, WITH_CHECKS) {
  urs <- compute_states_urs(trajectories = states, dim_list, WITH_CHECKS)
  if (settings_urs$ur_view) {
    graphics::par(mfrow = c(1, 1))
    graphics::matplot(urs, type = "l")
  }
  if (settings_urs$ur_save) {
    current_plot_name <- file.path(settings_urs$ur_path,
                                   paste0("00_", settings_urs$ur_name, ".pdf"))
    grDevices::setEPS()
    grDevices::postscript(current_plot_name, width = 9, height = 5.25)
    graphics::matplot(urs, type = "l")
    grDevices::dev.off()
    print(paste("Saved update rate plot in: ", current_plot_name))
  }
  return(urs)
}
