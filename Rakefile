desc "Build files for packaging"
task :default do
    sh 'Rscript -e "devtools::document()"'
    sh 'Rscript -e "knitr::knit(\"README.Rmd\")"'
    sh 'pandoc README.md -o README.html'
    sh 'Rscript -e "pkgdown::build_site()"'
end
