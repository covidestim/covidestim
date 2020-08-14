dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
if (!file.exists(M)) file.create(M)
cat(if( grepl("^darwin", R.version$os)) "\nCXX14FLAGS=-O3 -march=native -mtune=native -arch x86_64 -ftemplate-depth-256" else 
      if (.Platform$OS.type == "windows") "\nCXX14FLAGS=-O3 -mtune=native -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2" else
            "CXX14FLAGS = -fPIC -O3",
          file = M, sep = "\n", append = TRUE)
