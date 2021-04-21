# julia -t auto make.jl 

using Weave, Dates

cd("/Users/airwin/Dropbox/Julia/ode-parameter-estimation/")
files = readdir()
# [ i.match for i in match.(r".jmd$", files) ]

z = filter( x -> !isnothing(x), match.(r".*.jmd$", files))
jmd_files = [ i.match for i in z ]
# jmd_stat = stat.(jmd_files)
# html_files = replace.(jmd_files, "jmd" => "html")
# html_mtime = stat.(html_files)

for f in jmd_files
    jmd_info = stat(f)
    html_info = stat(replace("docs/" * f, "jmd" => "html"))

    if jmd_info.mtime > html_info.mtime
      print("Starting " * f * "\n")
      weave(f; doctype = "md2html", out_path = "docs", cache = :refresh) # cache options :all, :off, :refresh
      weave(f; doctype = "github", out_path = "docs", cache = :all) 
    else
      print("Skipping " * f * ": last modified " * string(unix2datetime(jmd_info.mtime)) * 
             " older than html modified at " * string(unix2datetime(html_info.mtime)) * "\n")
    end
end

# I cache using weave to enable output in two formats; only recompute if the jmd is newer than the html.
