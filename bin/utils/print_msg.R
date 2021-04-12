print_msg <- function(...){

    vector = c(...)
    msg = paste(vector, collapse="")
    msg = paste(Sys.time(), msg, "\n", collapse="\t")
    cat(msg)
}
