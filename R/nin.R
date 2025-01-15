#' @export `%nin%`
`%nin%` = Negate(`%in%`) # a negation of %in% operator, needed for samples removal to work (apparently lessR .()'s 'rows' argument can't begin with a '!')
