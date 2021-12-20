using Pkg
路径 = pwd()
println(路径)
Pkg.activate(路径)
Pkg.instantiate()