# obtain exec file from package path
obtain_exec = function(caller = c("cnv", "mutation", "indel"), out = getwd()) {
    caller = match.arg(caller)
    if (caller == "cnv") {
        caller_path = system.file("exec/cnv_pipe.R", package = "VSHunter")
        file.copy(from = caller_path, to = out)
    } else {
        stop ("Not support for now.")
    }
}
