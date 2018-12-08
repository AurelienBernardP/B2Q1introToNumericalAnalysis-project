using DelimitedFiles

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
