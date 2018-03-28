runTarMet <- function(){
  appdir <- system.file('app', package = 'TarMet')
  runApp(appdir, display.mode = 'normal')
}