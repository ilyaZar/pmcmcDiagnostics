#' Generates update rate plots
#'
#' If latent states are provided only!
#'
#' @param states latent state values
#' @param type_output_device output device type ("pdf" [default] or "eps")
#' @inheritParams analyse_states_convergence
#'
#' @return side effect function; plotting and/or saving
#' @export
get_update_rates <- function(states, settings_urs, dim_list,
                             WITH_CHECKS, type_output_device = "pdf") {
  urs <- compute_states_urs(trajectories = states, dim_list, WITH_CHECKS)

  if (settings_urs$ur_view) {
    graphics::par(mfrow = c(1, 1))
    graphics::matplot(urs, type = "l", ylim = c(0, 1))
  }

  if (settings_urs$ur_save) {
    # build filename based on device
    ext <- ifelse(type_output_device == "eps", "eps", "pdf")
    current_plot_name <- file.path(
      settings_urs$ur_path,
      paste0("00_", settings_urs$ur_name, ".", ext)
    )

    # open appropriate graphics device
    if (type_output_device == "eps") {
      grDevices::setEPS()
      grDevices::postscript(current_plot_name, width = 9, height = 5.25)
    } else {
      grDevices::pdf(current_plot_name, width = 9, height = 5.25)
    }

    graphics::matplot(urs, type = "l", ylim = c(0, 1))
    grDevices::dev.off()

    message("Saved update rate plot in: ", current_plot_name)
  }

  return(urs)
}
