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

mutable struct PivotTupple
    l::UInt128##line of the pivot
    c::UInt128##column of the pivot
    val::Int128##value of pivot

end

function makeCCSparseFromFile(filename::String)

    data = readdlm("email-Enron.txt",Int128, header=false, skipblanks=true, comments=true, comment_char='#')##the filename will be d to "filename" afterwards
    data = reshape(data',1,length(data)) ##all data is on a line 
    p = Array{Int128}(undef, convert(UInt128,length(data)/2+1))## lenght of the number of columns in the matrix plus one, the number of non 0 elements in the matrix
    t = Array{Int8}(undef,length(data))##lenght of all non 0 elements in matrix

    for i in 1:length(p)
        p[i] = (i*2-1)
    end
    p[length(p)] = p[length(p)] - 1 ## last element of p is the number of non 0 elements in the sparse array
    
    for i in 1:length(t)
        if i%2==0
            t[i]=1
        else
            t[i]=-1
        end
    end

    for i in 1:length(data)
        data[i] = data[i] + 1 ##reformat data to fit julia arrays(data started at index 0 julia arrays start at 1)
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
end

function swapLines!(matrix, firstLine, secondLine)
    # scan the sub-arrays defined by A.p[k ... k+1]
    for k in 1:lastindex(matrix.p)-1
        # idx of the currect sub-array
        beginning::UInt = matrix.p[k]
        ending::UInt = matrix.p[k + 1] - 1

        if (k == lastindex(p) - 1)
            ending += 1
        end

        firstIdx::Int = -1
        secondIdx::Int = -1
        counter::UInt = 0

        j::UInt = beginning
        while j <= ending
            # Switch the "number" of the line(s)
            if (matrix.i[j] == firstLine)
                matrix.i[j] = secondLine
                firstIdx = j
                counter += 1
            elseif (matrix.i[j] == secondLine)
                matrix.i[j] = firstLine
                secondIdx = j
                counter += 1
            end
            j += 1
        end

        # reorganize if necessary
        if (counter > 0 && beginning < ending)
            sortByColumns(matrix, beginning, ending, firstIdx, secondIdx, counter)
        end
    end
end

function sortByColumns(matrix, beginning, ending, firstIdx, secondIdx, counter)
    if (counter == 2)
        # Just need to swap in matrix.i and matrix.x
        matrix.i[firstIdx], matrix.i[secondIdx] = matrix.i[secondIdx], matrix.i[firstIdx]
        matrix.x[firstIdx], matrix.x[secondIdx] = matrix.x[secondIdx], matrix.x[firstIdx]
        return
    end

    # if counter == 1
    # Need to shift elements to reorganize

    mainIdx::UInt = max(firstIdx, secondIdx)
    if (mainIdx == beginning)
        shiftLeft(matrix, ending, mainIdx)
        return
    end

    if (mainIdx == ending)
        shiftRight(matrix, beginning, mainIdx)
        return
    end

    if (matrix.i[mainIdx] > matrix.i[mainIdx + 1])
        shiftLeft(matrix, ending, mainIdx)
        return
    end

    if (matrix.i[mainIdx] < matrix.i[mainIdx - 1])
        shiftRight(matrix, beginning, mainIdx)
        return
    end

end

function shiftRight(matrix, beginning, mainIdx)
    k::UInt = mainIdx - 1
    while k >= beginning
        if (matrix.i[k] > matrix.i[k + 1])
            matrix.i[k], matrix.i[k + 1] = matrix.i[k + 1], matrix.i[k]
        else
            break
        end
        k -= 1
    end

    j::UInt = mainIdx - 1
    while j > k
        matrix.x[j], matrix.x[j + 1] = matrix.x[j + 1], matrix.x[j]
        j -= 1
    end
end

function shiftLeft(matrix, ending, mainIdx)
    k::UInt = mainIdx + 1
    while k <= ending
        if (matrix.i[k] < matrix.i[k - 1])
            matrix.i[k], matrix.i[k - 1] = matrix.i[k - 1], matrix.i[k]
        else
            break
        end
        k += 1
    end

    j::UInt = mainIdx + 1
    while j < k
        matrix.x[j], matrix.x[j - 1] = matrix.x[j - 1], matrix.x[j]
        j += 1
    end
end

function lineSubstraction!(A, pivot, column, line, lineS)
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



function findNon0Pivot(A:: CCMatrix, pivot:: PivotTupple)
    if (A.p[pivot.c+1] - A.p[pivot.c]) <= 0
        ## column is empty
        return -1
    else
        ##column not empty
        
        for k:A.p[pivot.c] to A.p[pivot.c+1]
            if A.x[k] < pivot.l
                #there is a non 0 pivot available
                return k
            else
                #all the non zero elements of the column are in already triangulized nbEmptyColumns 
                #nothing to be done; continue to next column with Gauss GaussElimination
                return 0
            end
        end

    end
end

function changePivot!(A::CCMatrix,X::CCArray, B::CCArray, pivot::PivotTupple, nbEmptyColumns)
    nonZeroRow = findNon0Pivot
    if findNon0Pivot == -1
        ##column is empty, swap to end of matrix 
        permColumn!(A, X, pivot.c, (length(A.p)-1)-nbEmptyColumns)
        nbEmptyColumns += 1
        changePivot!(A, X, B, pivot, nbEmptyColumns)
        return
    elseif findNon0Pivot == 0
        ## there are no non 0 values in the column above the pivot 
        ## elimination of this column is complete
        return -1
    else
        # partialPivot #swap lines
        # a non zero pivot is available 
        swapLines!(A, nonZeroRow, pivot.l)
        swapLines!(B, nonZeroRow, pivot.l)
        pivot.val = pivotValue(A,pivot)
        return 0
    end

end

function pivotValue(A,pivot)
    for n in A.p[pivot.c]: A.p[pivot.c+1]
        if A.i[n] == pivot.l
            return A.x[n]
        end
    end
    return 0

end

####################################################
# makeZeroColumn!
# Arguments:
#       A: Matrix coded in CCsparse form || coefficients of equations
#       B: Column Matrix coded in a modified CCsparse || results of equations
#       X: Column Matrix coded in a modified CCsparse || unknowns of equations
#       Pivot: Pivot structure represeententing the pivot in a Gauss eliminations
# Description:
#       Eliminates all the non 0 elements above the pivot for a gauss elimination in an effort to make a Lower triangle matrix
# Output:
#       All the elements above the pivot are equal to 0 in the matrix A, also the Matrix B has beeen modified acordingly
#
####################################################
function makeZeroColumn!(A:: CCMatrix,X::CCArray,B::CCArray, pivot::PivotTupple)
    n = A.p[pivot.c -1 ]
    while n < A.p[pivot.c]
        if(A.i[n] < pivot.l )
            lineSubstraction!(A,B, pivot.val, pivot.c, pivot.l, A.i[n])     
        n+=1
    end   
end

####################################################
# GaussElimination!
# Arguments:
#       A: Matrix coded in CCsparse form || coefficients of equations
#       B: Column Matrix coded in a modified CCsparse || results of equations
#       X: Column Matrix coded in a modified CCsparse || unknowns of equations
# 
# Description:
#       Gauss elimination in an effort to make a Lower triangle matrix to solve A·X = B
# Output:
#       Gauss elimination has been completed and the matrix A X and B have been modified acordingly
#
####################################################
function GaussElimination!(A:: CCMatrix, X::CCArray, B::CCArray)
    if length(A.p <= 1)
        println(The input matrix is too small or indequately initialized)
        return
    end
    nbEmptyColumns = 0
    i = findmax(A.i)

    while i >= 2 && (nbEmptyColumns - length(A.p)) <= i
        pivot = PivotTupple(i, i, 0)
        pivot.val = pivotValue(A,pivot)
        if pivot.val == 0
           case = changePivot!(A, pivot)
           if case = -1##column is already eliminated
                i -= 1
                continue
           end
        end
        makeZeroColumn!(A,X,B,pivot)

        i -= 1
    end
    return
end

####################################################
# solveMatrixEquation!
# Arguments:
#       A: Matrix coded in CCsparse form || coefficients of equations
#       B: Column Matrix coded in a modified CCsparse || results of equations
# 
# Description:
#       Uses Gauss elimination solve A·X = B where A and B are sparse matricies and X is the matrix of unknowns to be found
# Output:
#       Prints a txt file with the unknowns X of the equations
#
####################################################
function solveMatrixEquation(A:: CCMatrix, B::CCArray)
    X = Array{Int128}(undef, convert(UInt128,length(A.p)-1)### to be changed !!!!!!!!!!!!!!!!!!! this is suposed to be a sparse array
    GaussElimination!(A:: CCMatrix, X::CCArray, B::CCArray)
    Substitution(A, B, X)


end