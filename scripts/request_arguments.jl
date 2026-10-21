# Shared by presentation and homology experiments. Parsing and validation finish
# before a runner creates output or starts an arithmetic computation.
function request_arguments(args;default_discriminants,default_output)
    root=default_output;resume=false;requests=Dict{String,Any}[]
    for arg in args
        if startswith(arg,"--output=")
            root=abspath(split(arg,'=';limit=2)[2])
        elseif arg=="--resume"
            resume=true
        elseif startswith(arg,"--requests=")
            path=split(arg,'=';limit=2)[2]
            data=JSON3.read(read(path,String))
            data isa AbstractVector || throw(ArgumentError("request file must contain a JSON array of records"))
            append!(requests,normalize_voronoi_request.(data))
        elseif startswith(arg,"--")
            throw(ArgumentError("unknown option $arg"))
        else
            push!(requests,voronoi_request(parse(Int,arg)))
        end
    end
    if isempty(requests)
        any(arg->startswith(arg,"--requests="),args) && throw(ArgumentError("request files contain no requests"))
        append!(requests,voronoi_request.(default_discriminants))
    end
    unique!(requests)
    (;root,resume,requests)
end

function matching_request(entry,request)
    haskey(entry,"request") && entry["request"]==request
end
