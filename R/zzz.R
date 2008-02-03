
.First.lib <- function(lib, pkg)
{
  library.dynam("gRc", package = pkg, lib.loc = lib)  
  return(invisible(0))
}
