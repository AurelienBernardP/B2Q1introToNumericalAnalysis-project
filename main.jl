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

function permColumn!(A, X, column1, column2)
    #remove the value from struct and put them into temporal array
    j::UInt = A.p[column1]
    k::UInt = A.p[column1 + 1]
    if column1 == lastindex(A.p)
        k += 1
    end
    iTmp1::Array{UInt128} = A.i[j:k-1]
    xTmp1::Array{Int8} = A.x[j:k-1]
    while j < k
        deleteat!(A.i, j)
        deleteat!(A.x, j)
        j += 1
    end
    j = A.p[column2]
    k = A.p[column2 + 1]
    if column2 == lastindex(A.p)
        k += 1
    end
    iTmp2::Array{UInt128} = A.i[j:k-1]
    xTmp2::Array{Int8} = A.x[j:k-1]
    while j < k
        deleteat!(A.i, j)
        deleteat!(A.x, j)
        j += 1
    end
    #insert value at their new place
    j = A.p[column2]
    i = column2
    k = 1
    while i < column2 + 1
        insert!(A.i, j, iTmp1[k])
        insert!(A.x, j, xTmp1[K])
        j += 1
        k += 1
        i += 1
    end
    j = A.p[column1]
    i = column1
    k = 1
    while i < column2 + 1
        insert!(A.i, j, iTmp2[k])
        insert!(A.x, j, xTmp2[k])
        j += 1
        k += 1
        i += 1
    end
    #update p
    j = A.p[column1 + 1] - A.p[column1]
    k = column1 + 1
    while k < column2
        A.p[k] += j
        k += 1
    end
    j = A.p[column2 + 1] - A.p[column2]
    k = column2 + 1
    while k < lastindex(A.p)
        A.p[k] += j
        k += 1
    end
    j = X[column1]
    X[column1] = X[column2]
    X[column2] = j
end

function Substitution(A, B, X)
    if A.i[1] == 1
        if B.i[1] == 1
            X[1] = B.x[1] / A.x[1]
        end
    end
    else
        X[1] = 0
    end
    line::UInt = 2
    over::UInt = B.i[length(B.i)]
    a::Int = 0
    b::Int = 0
    while line <= over
        if A.x[line] == 0
            X[line] = 0
            line += 1
        end
        else
            for k::UInt=1:line-1
                if A.i[k] == line
                    a += A.x[k] + X[k]
                end
            end
            if B.i[line] == line
                b = B.x[line]
            end
            b -= a
            X[line] = b/A.x[line]
            line ++
        end
    end
end

function lineSubstraction!(A, pivot, column, line)
    a::UInt = A.p[column] #start at the pivot's column
    mul::Int = A.x[A.p[column]]/pivot
    b::Uint
    nbrLines = maximum(A.i)
    if column == lastindex(A.p)
        b = A.nz + 1
    else
        b = A.p[column+1]
    end
    gap::Uint = b-a
    if (a+1 < b) && (A.i[a+1] != lineS)
        #the first non zero element of the line should be under the pivot
        #if not, we do nothing 
        return
    end
    i::Uint = column
    currentColumn::UInt = column
    j::Uint
    while i <= A.p[lastindex(A_p)]
        if (A.i[i] == lineS) && (A.[i - lineS + 1] == lineS) #both line are not empty
            A.x[i] += mul * A.[i - lineS + 1] # value in the pivot's line in the current column
            if A.x[i] == 0
                deleteat!(A.i, i)
                deleteat!(A.x, i)
                #update p and nz
                for k::Uint=currentColumn+1:lastindex(A_p)
                    A.p[k] -= 1
                end
                A.nz -= 1
            end

        elseif (A.i[i] != lineS) && (A.[i - lineS + 1] == lineS) #modified line is empty
            newVal = mul * A.[i - lineS + 1]
            insert!(A.x, i, newVal)
            insert!(A.i, i, lineS);
            #update p and nz
            for k::Uint=currentColumn+1:lastindex(A_p)
                A.p[k] += 1
            end
                A.nz += 1

        #if pivot's line is empty, we do nothing
        end
        i += 1
        if (A.p[i] == currentColumn + 1) && (currentColumn < lastindex(A.p))
            currentColumn += 1
        end
    end
    #update B
    i = 0
    j = 0
    i = findfirst(isequal(lineS), B.i)
    j = findfirst(isequal(line), B.i)
    if i && j
        B.x[i] += mul * B.x[j]
    elseif !i && j
        insert!(B.i, i)
        insert!(B.x,  mul * B.x[j])
        B.nz += 1
    end
end
