repo = "https://cloud.r-project.org"
install.packages("devtools", repos=repo)
install.packages("FME", repos=repo)
install.packages("bigleaf", repos=repo)
devtools::install_github("git@github.com:oscarperezpriego/ETpartitioning.git") 
