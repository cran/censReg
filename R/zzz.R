.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'censReg' package as:\n",
         "Henningsen, Arne (2017). ",
         "censReg: Censored Regression (Tobit) Models. ",
         "R package version 0.5. ",
         "http://CRAN.R-Project.org/package=censReg.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'censReg' package, ",
         "please use a forum or 'tracker' at the R-Forge site ",
         "of the 'sampleSelection' project:\n",
         "https://r-forge.r-project.org/projects/sampleselection/"),
      domain = NULL,  appendLF = TRUE )
}
