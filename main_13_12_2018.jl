module MatrixAG
    using DelimitedFiles



    #########################################################
    # CCMatrix Structure
    # Description:
    #       Structure for a Column compressed sparse matrix representation
    # Fields:
    #       p: Index in fields i and x of where each column starts
    #       i: Line index of the non zero elements in the matrix
    #       x: Value of the non zero elements in the matrix
    #       nz: Number of non zero elements in the matrix
    #
    #########################################################
    mutable struct CCMatrix
        p::Array{Int128}
        i::Array{Int128}
        x::Array{Int64}
        nz::UInt128
    end

    #########################################################
    # PivotTupple Structure
    # Description:
    #       Structure to store a pivot and its information
    # Fields:
    #       l: index of the line where the pivot is found
    #       c: index of the column where the pivot is found
    #       val: Value of the pivot
    #
    #########################################################
    mutable struct PivotTupple
        l::Int128 ##line of the pivot
        c::Int128 ##column of the pivot
        val::Int128 ##value of pivot

    end

    #########################################################
    # makeCCSparseFromFile
    # Arguments:
    #       path: a string containing the path to a delimited txt file with information about the arcs of a graph
    # Description:
    #       loads the contaigence matrix of the graph described by "path" in its Column compressed format
    # Output:
    #       A column compressed Sparse Matrix decribing the graph of path
    #
    #########################################################
    function makeCCSparseFromFile(path::String)

        data = readdlm(path,Int128, header=false, skipblanks=true, comments=true, comment_char='#')
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

        i::Int = 1 ## reorganizing the vectors
        while i <= length(sparseMatrix.i)
            if sparseMatrix.i[i] > sparseMatrix.i[i+1]
                sparseMatrix.i[i], sparseMatrix.i[i+1]  = sparseMatrix.i[i+1], sparseMatrix.i[i]
                sparseMatrix.x[i], sparseMatrix.x[i+1]  = sparseMatrix.x[i+1], sparseMatrix.x[i]
            end
            i += 2
        end
        return sparseMatrix
    end

    #########################################################
    # makeCCMatrix
    # Arguments:
    #       path: a string containing the path to a delimited txt file with information about the demands on the nodes of a graph
    # Description:
    #       loads the contaigence matrix of the demands described by "path" in its Column compressed format
    # Output:
    #       A column compressed Sparse Matrix decribing the demands of path
    #
    #########################################################
    function makeCCMatrix(path)
        data = readdlm(path, comments=true, comment_char='#')
        i = data[:, 1]
        for n in 1: length(i)
            i[n] += 1
        end
        x = data[:, 2]
        p = [1, length(x)]
        sparseArray = CCMatrix(p, i, x, length(data))
        sparseArray.i = vec(sparseArray.i)
        sparseArray.x = vec(sparseArray.x)

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
        gap1::UInt = sparseMatrix.p[column1 + 1] - sparseMatrix.p[column1]
        gap2::UInt = sparseMatrix.p[column2 + 1] - sparseMatrix.p[column2]
        iTmp1::Array{UInt128} = sparseMatrix.i[j:k-1]
        xTmp1::Array{Int64} = sparseMatrix.x[j:k-1]
        while j < k
            deleteat!(sparseMatrix.i, j)
            deleteat!(sparseMatrix.x, j)
            j += 1
        end
        j = sparseMatrix.p[column2]
        k = sparseMatrix.p[column2 + 1]
        if column2 == lastindex(A.p)
            k += 1
        end
        iTmp2::Array{UInt128} = sparseMatrix.i[j:k-1]
        xTmp2::Array{Int64} = sparseMatrix.x[j:k-1]
        while j < k
            deleteat!(sparseMatrix.i, j)
            deleteat!(sparseMatrix.x, j)
            j += 1
        end
        #update p
        k = lastindex(sparseMatrix.p)
        while k > sparseMatrix.p[column2]
            sparseMatrix.p[k] -= gap2
            k -= 1
        end
        k = lastindex(sparseMatrix.p)
        while k >  sparseMatrix.p[column1]
            sparseMatrix.p[k] -= gap1
            k -= 1
        end
        #insert value at their new place
        j = sparseMatrix.p[column2]

        while j < sparseMatrix.p[column1]
            insert!(sparseMatrix.i, iTmp1[j])
            insert!(sparseMatrix.x, xTmp1[j])
            j += 1
        end
        j = sparseMatrix.p[column1 + 1]
        while j < sparseMatrix.p[column2]
            insert!(sparseMatrix.i, iTmp2[j])
            insert!(sparseMatrix.x, xTmp2[j])
            j += 1
        end
        #update p
        k = sparseMatrix.p[column1 +1 ]
        while k < lastindex(sparseMatrix.p)
            sparseMatrix.p[k] += gap1
            k += 1
        end
        k = sparseMatrix.p[column2 + 1]
        while k < lastindex(sparseMatrix.p)
            sparseMatrix.p[k] += gap2
            k += 1
        end
    end

    #########################################################
    # swapLines!
    # Arguments:
    #       matrix: Matrix coded in CCsparse form || coefficients of equations
    #       firstLine: first line to swap
    #       secondLine: second line to swap
    # Description:
    #       Swap two lines of a matrix
    # Output:
    #       The modified Matrix
    #
    #########################################################
    function swapLines!(matrix, firstLine, secondLine)
        # scan the sub-arrays defined by A.p[k ... k+1]
        for k in 1:lastindex(matrix.p)-1
            # idx of the currect sub-array
            beginning::UInt = matrix.p[k]
            ending::UInt = matrix.p[k + 1] - 1

            # Last element in p equals nz 
            # Then we have to create ending = nz                                                                        
            if (k == lastindex(matrix.p) - 1)
                ending += 1
            end

            # -1 is a symbolic value. Idx is valid if != -1                        
            firstIdx::Int = -1
            secondIdx::Int = -1

            # Counts the number of element != 0 to permute in a column            
            counter::UInt = 0

            # Aurelien's (horrible) trick            
            j::UInt = beginning
            if ending > length(matrix.i)
                ending = length(matrix.i)
            end
            
            # Scans the whole current sub-array                                    
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

    #########################################################
    # sortByColumns!
    # Arguments:
    #       matrix: Matrix coded in CCsparse form || coefficients of equations  
    #       beginning: beginning's index of the sub-array                        
    #       ending: ending's index of the sub-array
    #       firstLine: first line to swap
    #       secondLine: second line to swap
    #       counter: number of non-zero element to swap in a column        
    # Description:
    #       Reorganize the vector i and x
    # Output:
    #       The modified Matrix
    #
    #########################################################
    function sortByColumns!(matrix, beginning, ending, firstIdx, secondIdx, counter)
        if (counter == 2)  && (firstIdx >= 1) && (secondIdx >= 1)
            # Just need to swap in matrix.i and matrix.x
            matrix.i[firstIdx], matrix.i[secondIdx] = matrix.i[secondIdx], matrix.i[firstIdx]
            matrix.x[firstIdx], matrix.x[secondIdx] = matrix.x[secondIdx], matrix.x[firstIdx]
            return
        end

        # if counter == 1
        # Need to shift elements to reorganize

        mainIdx::UInt = max(firstIdx, secondIdx)

        # Special case to avoid Segmentation Fault: index is the first value in the sub-array
        if (mainIdx == beginning)
            shiftLeft!(matrix, ending, mainIdx)
            return
        end
        
        # Special case to avoid Segmentation Fault: index is the last value in the sub-array                
        if (mainIdx == ending)
            shiftRight!(matrix, beginning, mainIdx)
            return
        end

        # No Segmentation Fault here because of the two previous conditions
        if (matrix.i[mainIdx] > matrix.i[mainIdx + 1])
            shiftLeft!(matrix, ending, mainIdx)
            return
        end

        if (matrix.i[mainIdx] < matrix.i[mainIdx - 1])
            shiftRight!(matrix, beginning, mainIdx)
            return
        end

    end

    #########################################################
    # shiftRight!
    # Arguments:
    #       matrix: Matrix coded in CCsparse form || coefficients of equations  
    #       beginning: beginning's index of the sub-array
    #       mainIdx: index of the value to move
    # Description:
    #       Shifts the elements in i and x
    # Output:
    #       The modified Matrix
    #
    #########################################################
    function shiftRight!(matrix, beginning, mainIdx)
        # Only one element is misplaced. The sub-array must be sorted
        # Then we use the property of a sorted array                                    
        k::UInt = mainIdx - 1
        while k >= beginning
            if (matrix.i[k] > matrix.i[k + 1])
                matrix.i[k], matrix.i[k + 1] = matrix.i[k + 1], matrix.i[k]
            else
                break
            end
            k -= 1
        end

        # k is now the index of the destination
        # i and x must correspond
        j::UInt = mainIdx - 1
        while j > k
            matrix.x[j], matrix.x[j + 1] = matrix.x[j + 1], matrix.x[j]
            j -= 1
        end
    end

    #########################################################
    # shiftLeft!
    # Arguments:
    #       matrix: Matrix coded in CCsparse form || coefficients of equations  
    #       ending: ending's index of the sub-array
    #       mainIdx: index of the value to move
    # Description:
    #       Shifts the elements in i and x
    # Output:
    #       The modified Matrix
    #
    #########################################################
    function shiftLeft!(matrix, ending, mainIdx)
        # Only one element is misplaced. The sub-array must be sorted
        # Then we use the property of a sorted array        
        k::UInt = mainIdx + 1
        while k <= ending
            if (matrix.i[k] < matrix.i[k - 1])
                matrix.i[k], matrix.i[k - 1] = matrix.i[k - 1], matrix.i[k]
            else
                break
            end
            k += 1
        end

        # k is now the idx of the destination
        # i and x must correspond
        j::UInt = mainIdx + 1
        while j < k
            matrix.x[j], matrix.x[j - 1] = matrix.x[j - 1], matrix.x[j]
            j += 1
        end
    end

    #########################################################
    # Substitution!
    # Arguments:
    #       A: Matrix coded in CCsparse form || coefficients of equations
    #       B: Column Matrix coded in a modified CCsparse || results of equations
    #       X: Column Matrix coded in a modified CCsparse || unknowns of equations
    #       nblines: number of lines in matrix A
    # Description:
    #       Performs top to bottom substitution on a lower left triangularized CCsparse Matrix
    #       to find the values in the vector X
    # Output:
    #       The values of the unknons in Sparse vector X are found and placed in it
    #
    #########################################################
    function Substitution!(A::CCMatrix, B::CCMatrix, X::CCMatrix, nblines)
    if A.i[1] == 1
        X.x[1] = B.x[1] / A.x[1]
    else
        X.x[1] = 0
    end
    i::UInt = lastindex(A.p)

    line::UInt = 2
    p::Int = 2
    a::Int = 0
    b::Int = 0
    lineInB = 2
        while p <= nblines 
            #see if diagonal element is null
        if A.i[A.p[p]] == line
            i = 1
            #going column by column
            for j::Int=1:A.p[p]-1 #case last column only happend if input matrix is square
                if A.i[j] == line
                    a += A.x[j] * X.x[i]
                    # j see each non null element of matrix, i is the current column
                    i += 1
                elseif j >= A.p[i+1]
                    i += 1
                    #element in column i is null at line line. j is at the first non null element in the next column
                end
            end
            b = 0
            #does B have a non zero element at line ?
            while B.i[lineInB] < line
                lineInB +=1
            end
            if B.i[lineInB] == line
                b = B.x[lineInB]
                lineInB += 1
            end

            if A.x[A.p[p]] == 0
                #should never happen
                println("a = 0")
            else
                X.x[p] = (b - a) / (A.x[A.p[p]])
            end
            else
                X.x[p] = 0
            end
            line += 1
            p += 1
        end
    end

    #########################################################
    # findNon0Pivot
    # Arguments:
    #       A: Matrix coded in CCsparse form || coefficients of equations
    #       pivot: Pivot structure represeententing the pivot in a Gauss eliminations
    # Description:
    #       finds a non 0 pivot abouve the pivot if there are any in the sparse matrix
    # Output:
    #       -1 : if the column in which we are looking for the pivot is empty
    #       0 : if the column in which weare looking for the pivot is already correctly triangulized
    #       strictly positive integer: the line in which a suitable non zero pivot is 
    #########################################################
    function findNon0Pivot(A:: CCMatrix, pivot:: PivotTupple)
        if (A.p[pivot.c+1] - A.p[pivot.c]) <= 0
            ## column is empty
            return -1
        else
            ##column not empty
            for k in A.p[pivot.c] : A.p[pivot.c+1]-1
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
    # changePivot!!
    # Arguments:
    #       A: Matrix coded in CCsparse form || coefficients of equations
    #       B: Column Matrix coded in a modified CCsparse || results of equations
    #       X: Column Matrix coded in a modified CCsparse || unknowns of equations
    #       Pivot: Pivot structure represeententing the pivot in a Gauss eliminations
    #       nbEmptyColumns: number of empty columns in the matrix
    # Description:
    #       Swaps columns and lines of the matrix to have a non 0 pivot for a gauss elimination 
    # Output:
    #       -1: if there are no non zero elements above the pivot and therefore the elimination of that line is completed
    #       0: if a partial pivoting took place and now the pivot is non 0
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

            #for i = 1:length(B.i)
            #    if B.i[i] == pivot.l || B.i[i] == nonZeroRow
                    swapLines!(B, nonZeroRow, pivot.l)
            #    end
           # end

            pivot.val = pivotValue(A, pivot)##update pivot value
 
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
                return A.x[n]
            end
        end
        ##if no value with the column and line index of the pivot is found, return 0
        return 0
    end

    ####################################################
    # makeZeroColumn!
    # Arguments:
    #       A: Matrix coded in CCsparse form || coefficients of equations
    #       B: Column Matrix coded in a modified CCsparse || results of equations
    #       X: Column Matrix coded in a modified CCsparse || unknowns of equations
    #       Pivot: Pivot structure represeententing the pivot in a Gauss eliminations
    #       nblines: the number of lines in the matrix
    # Description:
    #       Eliminates all the non 0 elements above the pivot for a gauss elimination in an effort to make a Lower triangle matrix
    # Output:
    #       All the elements above the pivot are equal to 0 in the matrix A, also the Matrix B has beeen modified acordingly
    #
    ####################################################
    function makeZeroColumn!(A:: CCMatrix,X::CCMatrix,B::CCMatrix, pivot::PivotTupple, nblines::UInt128)
        n = A.p[pivot.c]
        linesToEliminate = Int128[]
        while n < A.p[pivot.c + 1] - 1
            if (A.i[n] < pivot.l)
                push!(linesToEliminate, A.i[n])
            end
            n += 1
        end
      #  if length(linesToEliminate) != 0 && length(linesToEliminate) != 1
      #      println(length(linesToEliminate))
      #  end
        for i in 1:length(linesToEliminate)
            #print(pivot.l, linesToEliminate[i])
            SubstractLine!(A, B, pivot, linesToEliminate[i])
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
                println(i)
                continue
            end
            makeZeroColumn!(A,X,B,pivot, nblines)
            i -= 1
            println(i)
        end
        #println(nbEmptyColumns)
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
    function isTriag(A::CCMatrix, print)
            nblines::Int128 = maximum(A.i)
            counter = 0
            for i in 1:nblines
                for j in A.p[i]:A.p[i+1]-1
                    if A.i[j] < i
                        counter +=1
                        if print == 1
                            print("  line =", A.i[j])
                            print("  Column =", i)
                            print("  Value =", A.x[j])
                            print("\n")
                        end
                    end
                    if A.x[j] == 0
                        println("hi! I am a zero and I shouldn't be here")
                    end
                end
            end
            println("#######nb of not triag elements ########")
            println(counter)
            println("size of a.i", length(A.i))
            println("#########################")
            counter = 0
            for i in 1:nblines
                for j in A.p[i]:A.p[i+1]
                    if A.i[j] >= i
                        counter +=1
                    end
                    if A.x[j] == 0
                        println("hi! I am a zero and I shouldn't be here")
                    end
                end
            end
            println("####### nb of well triag elements ########")
            println(counter)
            println("#########################")

    end

    function SubstractLine!(A::CCMatrix, B::CCMatrix, pivot::PivotTupple, linesub)
        lineVal = 0
        for i in A.p[pivot.c]:A.p[pivot.c+1]-1 ### find the non zero element to eliminate
            if A.i[i] == linesub
                lineVal = A.x[i]
            end
        end

        if lineVal == 0
            println("error")
        end
        multiplicator = lineVal / pivot.val
        SubstractionInB(B, pivot, linesub, multiplicator)

        elemOfLineIsNonZero = false
        elemOfPivotLineIsNonZero = false

        for i = 1:(length(A.p)-1) ### iterate over all non zero elements
            elemOfLineIsNonZero = false
            elemOfPivotLineIsNonZero = false

            valueOfLineElement::Int = 0
            valueOfPivotLineElement::Int = 0

            indexOfLineElement::Int = 0

            j = A.p[i]
            endSubarray = A.p[i+1]
            if endSubarray > length(A.i)
               ## println("################ Corrected #####################")
                endSubarray = length(A.i)
            end

            while j < endSubarray

                if A.i[j] == pivot.l
                     elemOfPivotLineIsNonZero = true
                     valueOfPivotLineElement = A.x[j]

                end
                if A.i[j] == linesub
                    elemOfLineIsNonZero = true
                    valueOfLineElement = A.x[j]
                    indexOfLineElement = j
                end
                if A.i[j] > pivot.l## as the columns are always sorted it can break the loop as it will not find the values
                   break
                end
                j += 1
            end

            if elemOfPivotLineIsNonZero == true || elemOfLineIsNonZero == true
                ##one of the elements in non zero
                if elemOfPivotLineIsNonZero == true && elemOfLineIsNonZero == false
                    ## the pivot line element is non zero and the element of the line is zero
                    for n in A.p[i]:A.p[i+1]-1
                        if A.i[n] >= linesub
                            insert!(A.x, n, -(valueOfPivotLineElement*multiplicator))
                            insert!(A.i, n, linesub)
                            A.nz += 1
                            for k in i+1:length(A.p)
                                A.p[k]+=1
                            end
                            break### it has found the place of the place
                        end
                    end
                end
                if elemOfPivotLineIsNonZero == false && elemOfLineIsNonZero == true
                    ## the pivot element is zero but he element of the line is not zero
                    ##A.x[indexOfLineElement] = -A.x[indexOfLineElement]

                end
                if elemOfPivotLineIsNonZero == true && elemOfLineIsNonZero == true
                ## both elements are non zero
                    newValue::Int = valueOfLineElement - valueOfPivotLineElement*multiplicator
                    if newValue == 0

                        deleteat!(A.x, indexOfLineElement)
                        deleteat!(A.i, indexOfLineElement)
                        A.nz -= 1
                        for k in i+1:length(A.p)
                            A.p[k] -= 1
                        end
                    else
                        A.x[indexOfLineElement] = newValue
                    end
                end
            end
        end
    end


    function SubstractionInB(B, pivot, linesub, multiplicator)
    
        valueInPivotLine::Int = 0
        valueInLinesub::Int = 0
        indexInBIofLinesub::Int = 0
        for i in 1:length(B.i)
            if B.i[i] == pivot.l
                valueInPivotLine = B.x[i]
            end
            if B.i[i] == linesub
                valueInLinesub = B.x[i]
                indexInBIofLinesub = i
            end
        end

        newValue::Int = valueInPivotLine * multiplicator - valueInLinesub
        if indexInBIofLinesub == 0
            for i in 1:length(B.i)-1
                if B.i[i] >= linesub
                    insert!(B.i, i, linesub)
                    insert!(B.x, i, newValue)
                    break
                end
            end
            B.p[2] += 1
            B.nz += 1
        else
            B.x[indexInBIofLinesub] = newValue
        end
    end### end substractionInB

end### end of module
