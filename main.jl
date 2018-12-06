using DelimitedFiles

mutable struct CCMatrix
    p::Array{UInt128}
    i::Array{UInt128}
    x::Array{Int8}
    nz::UInt128
end

function makeCCSparseFromFile(x::String)

    data = readdlm("email-Enron.txt",Int128, header=false, skipblanks=true, comments=true, comment_char='#')##the filename will be changed to x afterwards
    data = reshape(data',1,length(data)) ##all data is on a line 
    p = Array{Int128}(undef,length(data))
    t = Array{Int8}(undef,length(data))
    for i::UInt128 = 1:length(data)
        p[i] = (i*2-1)
        if i%2==0
            t[i]=1
        else
            t[i]=-1
        end
    end
    sparseMatrix = CCMatrix(p,data,t,length(data))
    return sparseMatrix
end

function ReadArray(path)
    using DelimitedFiles
    Array = readdlm(path, comments=true, comment_char='#')
    Array_i = Array[:, 1]
    Array_v = Array[:, 2]
end