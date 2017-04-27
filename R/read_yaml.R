#' Read YAML file
#'
#' @param yaml_file YAML file
#'
#' @return List of parsed YAML data if file found, otherwise \code{NULL}
#' @export
read_yaml <- function(yaml_file) {
    if (file.exists(yaml_file)) {
        yaml <- yaml.load_file(yaml_file)
    } else {
        yaml <- NULL
    }
    return(yaml)
}
