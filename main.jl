using DelimitedFiles

mutable struct CCMatrix
    p::Array{UInt128}
    i::Array{UInt128}
    x::Array{Int8}
    nz::UInt128
end

mutable struct CCArray
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

function makeCCArray(path)
    data = readdlm(path, comments=true, comment_char='#')
    i= Array[:, 1]
    x = Array[:, 2]
    sparseArray = CCArray(i, x, length(data))
    return sparseArray
end

function permColumn!(sparseMatrix, sparseArray, column1, column2)
    #remove the value from struct and put them into temporal array
    j::UInt8 = sparseMatrix.p[column1]
    k::UInt8 = sparseMatrix.p[column1 + 1]
    iTmp1::Array{UInt128} = sparseMatrix.i[j:k-1]
    xTmp1::Array{Int8} = sparseMatrix.x[j:k-1]
    while j < k
        deleteat!(sparseMatrix.i, j)
        deleteat!(sparseMatrix.x, j)
        j += 1
    end
    j = sparseMatrix.p[column2]
    k = sparseMatrix.p[column2 + 1]
    iTmp2::Array{UInt128} = sparseMatrix.i[j:k-1]
    xTmp2::Array{Int8} = sparseMatrix.x[j:k-1]
    while j < k
        deleteat!(sparseMatrix.i, j)
        deleteat!(sparseMatrix.x, j)
        j += 1
    end
    #insert value at their new place
    j = sparseMatrix.p[column2]
    while j < sparseMatrix.p[column1]
        insert!(sparseMatrix.i, iTmp1[j])
        insert!(sparseMatrix.x, xTmp1[j])
        j += 1
    end
    j = sparseMatrix.p[column1]
    while j < sparseMatrix.p[column2]
        insert!(sparseMatrix.i, iTmp2[j])
        insert!(sparseMatrix.x, xTmp2[j])
        j += 1
    end
    #update p
    j = sparseMatrix.p[column1 + 1] - sparseMatrix.p[column1]
    k = sparseMatrix.p[column1 + 1]
    while k < sparseMatrix.p[column2]
        sparseMatrix.p[k] +=j
        k += 1
    end
    j = sparseMatrix.p[column2 + 1] - sparseMatrix.p[column2]
    k = sparseMatrix.p[column2 + 1]
    while k < sparseMatrix.p[column2]
        sparseMatrix.p[k] +=j
        k += 1
    end
    j = sparseX[column1]
    sparseX[column1] = sparseX[column2]
    sparseX[column2] = j
end