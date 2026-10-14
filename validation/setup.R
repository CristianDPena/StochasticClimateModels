# Run from repository root. No changes to the user's global R library.
dir.create('validation/.library',recursive=TRUE,showWarnings=FALSE)
install.packages('evir',repos='https://cloud.r-project.org',lib='validation/.library',type='source')
.libPaths(c(normalizePath('validation/.library'),.libPaths()))
stopifnot(as.character(packageVersion('evir'))=='1.7.4')
