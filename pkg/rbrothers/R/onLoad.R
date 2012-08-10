.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  Sys.setenv(CLASSPATH=paste(system.file("java/DualBrothers.jar",package="dualbrothers"),":",system.file("java/colt.jar",package="dualbrothers"),sep=""))
}