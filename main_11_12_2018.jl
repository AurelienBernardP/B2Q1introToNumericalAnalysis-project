module MatrixAG
    using DelimitedFiles

    mutable struct CCMatrix
        p::Array{Int128}
        i::Array{Int128}
        x::Array{Int64}
        nz::UInt128
    end


    mutable struct PivotTupple
        l::Int128 ##line of the pivot
        c::Int128 ##column of the pivot
        val::Int128 ##value of pivot

    end

    function makeCCSparseFromFile(path::String)

        data = readdlm(path,Int128, header=false, skipblanks=true, comments=true, comment_char='#')##the filename will be d to "filename" afterwards
        data = reshape(data',1,length(data)) ##all data is on a line 
        p = Array{Int128}(undef, convert(UInt128,length(data)/2+1))## length of the number of columns in the matrix plus one, the number of non 0 elements in the matrix
        t = Array{Int64}(undef,length(data))##length of all non 0 elements in matrix

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
        sparseMatrix.i = vec(sparseMatrix.i)
        sparseMatrix.x = vec(sparseMatrix.x)
        return sparseMatrix
    end

    function makeCCMatrix(path)
        data = readdlm(path, comments=true, comment_char='#')
        i = data[:, 1]
        for n in 1: length(i)
            i[n] += 1
        end
        x = data[:, 2]
        p = [1, length(x)]
        sparseArray = CCMatrix(p, i, x, length(data))
        return sparseArray
    end

    #########################################################
    # permColumn!
    # Arguments:
    #       sparseMatrix: Matrix coded in CCsparse form || coefficients of equations
    #       sparseArray: Array coded in CCsparse form || independent members of equations
    #       column1 and column2: the columns of the matrix to swap
    # Description:
    #       Swap two column of a matrix
    # Output:
    #       The modified Matrix
    #
    #########################################################
    function permColumn!(sparseMatrix::CCMatrix, sparseArray::CCMatrix, column1, column2)
        #remove the value from struct and put them into temporal array
        j::UInt8 = sparseMatrix.p[column1]
        k::UInt8 = sparseMatrix.p[column1 + 1]
        iTmp1::Array{UInt128} = sparseMatrix.i[j:k-1]
        xTmp1::Array{Int64} = sparseMatrix.x[j:k-1]
        while j < k
            deleteat!(sparseMatrix.i, j)
            deleteat!(sparseMatrix.x, j)
            j += 1
        end
        j = sparseMatrix.p[column2]
        k = sparseMatrix.p[column2 + 1]
        iTmp2::Array{UInt128} = sparseMatrix.i[j:k-1]
        xTmp2::Array{Int64} = sparseMatrix.x[j:k-1]
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
            sparseMatrix.p[k] += j
            k += 1
        end
        j = sparseMatrix.p[column2 + 1] - sparseMatrix.p[column2]
        k = sparseMatrix.p[column2 + 1]
        while k < sparseMatrix.p[column2]
            sparseMatrix.p[k] += j
            k += 1
        end
    end

    function swapLines!(matrix, firstLine, secondLine)
        # scan the sub-arrays defined by A.p[k ... k+1]
        for k in 1:lastindex(matrix.p)-1
            # idx of the currect sub-array
            beginning::UInt = matrix.p[k]
            ending::UInt = matrix.p[k + 1] - 1

            if (k == lastindex(matrix.p) - 1)
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
                sortByColumns!(matrix, beginning, ending, firstIdx, secondIdx, counter)
            end
        end
    end

    function sortByColumns!(matrix, beginning, ending, firstIdx, secondIdx, counter)
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
            shiftLeft!(matrix, ending, mainIdx)
            return
        end

        if (mainIdx == ending)
            shiftRight!(matrix, beginning, mainIdx)
            return
        end

        if (matrix.i[mainIdx] > matrix.i[mainIdx + 1])
            shiftLeft!(matrix, ending, mainIdx)
            return
        end

        if (matrix.i[mainIdx] < matrix.i[mainIdx - 1])
            shiftRight!(matrix, beginning, mainIdx)
            return
        end

    end

    function shiftRight!(matrix, beginning, mainIdx)
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

    function shiftLeft!(matrix, ending, mainIdx)
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

    #########################################################
    # lineSubstraction!
    # Arguments:
    #       A: Matrix coded in CCsparse form || coefficients of equations
    #       pivot: the element of the matrix used to do the substraction
    #       column: the column where the pivot is
    #       line: the line where the pivot is
    #       lineS: the line of the matrix to be modified
    # 
    # Description:
    #       Modify a line of a CCsparse coded Matrix with a linear combination of another line
    # Output:
    #       The modified Matrix
    #
    #########################################################
    function lineSubstraction!(A::CCMatrix, B::CCMatrix, pivot::Int128, column, line, linesub, nblines::UInt128)
        if pivot == 0
            return -1
        end
        mul::Int128 = convert(Int128, A.x[A.p[column]]/pivot)
        i::Int128 = A.p[column] 
        LastColumn::Int128 = 0
        j::Int128 = 1
        k::Int128 = 1
        while i < lastindex(A.p)
                k = A.i[A.p[i]]
                if i + 1 == lastindex(A.p)
                    LastColumn = 1
                end
                while (A.i[k] != linesub) && (k - LastColumn < A.p[i+1])
                    k += 1 # find if column i has a non zero element in linesub 
                end

                if k - LastColumn == A.p[i+1] # 0 in linesub
                    j = A.i[A.p[i]]
                else
                    j = A.i[k] # not 0 in line sub, element in line come after in the sparse structure                                                            
                end

                while (A.i[j] != line) && (j - LastColumn < A.p[i+1])
                    j += 1 # find if column i has a non zero element in line
                end

                if (A.i[j] == line) && (A.i[k] == linesub) #both line in column i are not empty
                    A.x[k] += mul * A.i[j]
                    if A.x[k] == 0
                        deleteat!(A.i, k)
                        deleteat!(A.x, k)
                        #update p and nz
                        for k::UInt=currentColumn+1:lastindex(A_p)
                            A.p[k] -= 1
                        end
                        A.nz -= 1
                    end
                end
                if (A.i[j] == line) && (k - LastColumn == A.p[i+1]) #column i has a 0 element in linesub
                    newVal = mul * A.i[j]
                    insert!(A.x, k, newVal)
                    insert!(A.i, k, linesub)
                    #update p and nz
                    for k::Unt=currentColumn+1:lastindex(A.p)
                        A.p[k] += 1
                    end
                    A.nz += 1
                # if column i has a 0 element in line, we do nothing
                end
            i += 1
        end
        #update B
        
        d = findfirst(isequal(linesub), B.i)
        f = findfirst(isequal(line), B.i)

        if (d != nothing) && (f != nothing)
            B.x[d] += mul * B.x[f]
        elseif (d == nothing) && (f != nothing)
            insert!(B.i, d)
            insert!(B.x,  mul * B.x[f])
            B.nz += 1
        end
    end

    function Substitution!(A::CCMatrix, B::CCMatrix, X::CCMatrix, nblines)
        if A.i[1] == 1
            if B.i[1] == 1
                X.x[1] = B.x[1] / A.x[1]
            end
        else
            X.x[1] = 0
        end
        i::UInt = lastindex(A.p)
        while (i > nblines) && (i > 1)
            if A.p[i-1] != A.p[i]
                splice!(A.i, A.p[i-1]:A.p[i]-1)
                splice!(A.x, A.p[i-1]:A.p[i]-1)
            end
            A.p[i] = A.nz
            i -= 1
        end
        line::UInt = 2
        p::Int = 2
        a::Int = 0
        b::Int = 0
        lineInB = 1
        endColumn::Int = 1
        while p <= nblines
            #see if diagonal element is null
            if p + 1 == lastindex(A.p)
                #if we're at the last column
                endColumn = 0
            end
            if A.i[A.p[p+1]-endColumn] == line
                i = 1
                #going column by column
                for j::Int=1:A.p[p]-1 #case last column?
                    if A.i[j] == line
                        a += A.x[j] * X.x[i]
                        # j see each non null element of matrix, i is the current column
                        i += 1
                    elseif j == A.p[i+1]
                        i += 1
                        #element in column i is null at line line. j is at the first non null element in the next column
                    end
                end
                b = 0
                #does B have a non zero element at line ?
                if B.i[lineInB] == line
                    b = B.x[lineInB]
                    lineInB += 1
                end
                if A.x[A.p[p+1]-endColumn] == 0
                    #should never happen
                    println("a = 0")
                else                    
                    X.x[line] = (b - a) / (A.x[A.p[p+1]-endColumn]) 
                end
            else
                X.x[line] = 0
            end
        line += 1
        p += 1
        end
    end

    function findNon0Pivot(A:: CCMatrix, pivot:: PivotTupple)
        if (A.p[pivot.c+1] - A.p[pivot.c]) <= 0
            ## column is empty
            return -1
        else
            ##column not empty
            for k in A.p[pivot.c] : A.p[pivot.c+1]
                if A.i[k] < pivot.l
                    #there is a non 0 pivot available
                    return A.i[k]
                else
                    #all the non zero elements of the column are in already triangulized nbEmptyColumns 
                    #nothing to be done; continue to next column with Gauss GaussElimination
                    return 0
                end
            end
        end
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
    function changePivot!(A::CCMatrix, X::CCMatrix, B::CCMatrix, pivot::PivotTupple, nbEmptyColumns)
        nonZeroRow = findNon0Pivot(A, pivot)
        #println(nonZeroRow)
        if nonZeroRow == -1
            ##column is empty, swap to end of matrix 
            permColumn!(A, X, pivot.c, (length(A.p)-1)-nbEmptyColumns)
            nbEmptyColumns += 1
            pivot.val = pivotValue(A, pivot)
            if pivot.val == 0 
                changePivot!(A, X, B, pivot, nbEmptyColumns)
            end
            return
        elseif nonZeroRow == 0
            ## there are no non 0 values in the column above the pivot 
            ## elimination of this column is complete
            return -1
        else
            # partialPivot #swap lines
            # a non zero pivot is available
            #println(A.i) 
            swapLines!(A, nonZeroRow, pivot.l)
            swapLines!(B, nonZeroRow, pivot.l)
            pivot.val = pivotValue(A, pivot)##update pivot value
            #println("lines were swapped")
            #println(A.i)
            #println(pivot.val)
            return 0
        end
    end

    ####################################################
    # pivotValue
    # Arguments:
    #       A: Matrix coded in CCsparse form || coefficients of equations
    #       Pivot: Pivot structure represeententing the pivot in a Gauss eliminations
    # Description:
    #       Determins the value of the pivot from its position in matrix A
    # Output:
    #       The value at the pivot position in the Sparse matrix A
    #
    ####################################################
    function pivotValue(A,pivot)
        for n in A.p[pivot.c]: A.p[pivot.c+1]
            if A.i[n] == pivot.l
                #println("piv non 0")
                return A.x[n]
            end
        end
        #println("piv 0")
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
    function makeZeroColumn!(A:: CCMatrix,X::CCMatrix,B::CCMatrix, pivot::PivotTupple, nblines::UInt128)
        n = A.p[pivot.c]
        linesToEliminate = Int128[]
        while n < A.p[pivot.c + 1]
            if (A.i[n] < pivot.l)
                push!(linesToEliminate, A.i[n])
            end
            n += 1
        end
      #  if length(linesToEliminate) != 0 && length(linesToEliminate) != 1
      #      println(length(linesToEliminate))
      #  end
        for i in 1:length(linesToEliminate)
            lineSubstraction!(A, B, pivot.val, pivot.c, pivot.l, linesToEliminate[i], nblines)
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
    function GaussElimination!(A::CCMatrix, X::CCMatrix, B::CCMatrix, nblines::UInt128)
        if length(A.p) <= 1
            #println("The input matrix is too small or indequately initialized")
            return
        end
        nbEmptyColumns = 0
        i = nblines
        if ((length(A.p)-1) - nbEmptyColumns) < nblines
           println("too many empty columns")
        end
        while i >= 2 && (((length(A.p)-1) - nbEmptyColumns) > nblines)
            pivot = PivotTupple(i, i, 0)
            pivot.val = pivotValue(A,pivot)

            case = 0
            if pivot.val == 0
                case = changePivot!(A, X, B, pivot, nbEmptyColumns)
                pivot.val = pivotValue(A,pivot)##### pivot value is not necesary we already update it in changePivot!()
            end
            if case == -1 ##column is already eliminated
                i -= 1
                continue
            end
            makeZeroColumn!(A,X,B,pivot, nblines)
            i -= 1
            print
        end
        return
    end

    #########################################################
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
    #########################################################
    function solveMatrixEquation(A::CCMatrix, B::CCMatrix)
        
        X = CCMatrix([1 length(A.p)-1], Array{Int128}(undef,length(A.p)-1), Array{Int64}(undef,length(A.p)-1), length(A.p)-1)
        for n in 1:length(X.i)
            X.i[n] = n
            X.x[n] = 0
        end
        nblines::UInt128 = maximum(A.i)
    #    println("*****before elimination*****")
    #    println(A)
    #    println(B)
    #    println("****before elimination*****") 
        isTriag(A, 0)
        GaussElimination!(A, X, B, nblines)
        isTriag(A, 0)
        for i in 1:length(A.x)
            if A.x[i] == 0
                println("hi! I am a zero and I shouldn't be here")
            end
        end
        Substitution!(A, B, X, nblines)
        
        return X
    end

    function isTriag(A::CCMatrix, print)
        nblines::Int128 = maximum(A.i)
        counter = 0
        for i in 1:nblines
            for j in A.p[i]:A.p[i+1]
                if A.i[j] < i
                    counter +=1
                end
                if A.x[j] == 0
                    println("hi! I am a zero and I shouldn't be here")
                end
                if print == 1
                    print("  line =", A.i[j])
                    print("  Column =", i)
                    print("  Value =", A.x[j])
                    print("\n")
                end
            end
        end
        println("#########################")
        println(counter)
        println("#########################")
    end

end

