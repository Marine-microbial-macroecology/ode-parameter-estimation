# julia -t auto make.jl 

using Weave

cd("/Users/airwin/Dropbox/Julia/ode-parameter-estimation/")
files = readdir()
# [ i.match for i in match.(r".jmd$", files) ]

z = filter( x -> !isnothing(x), match.(r".*.jmd$", files))
jmd_files = [ i.match for i in z ]
# jmd_stat = stat.(jmd_files)
# html_files = replace.(jmd_files, "jmd" => "html")
# html_mtime = stat.(html_files)

for f in jmd_files
    if stat(f).mtime > stat(replace("build" * f, "jmd" => "html")).mtime
      print("Starting " * f * "\n")
      weave(f; doctype = "md2html", out_path = "build", cache = :off) # cache options :all, :off, :refresh
    end
end

# I don't cache using weave since I only recompute if the jmd is newer than the html.
