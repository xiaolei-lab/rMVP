devtools::document()
file.remove(file.path("man", dir(path = "man", pattern = "^FarmCPU")))
file.remove(file.path("man", dir(path = "man", pattern = "_")))

devtools::build()

