
function ReadArray(path)
    using DelimitedFiles
    Array = readdlm(path, comments=true, comment_char='#')
    Array_i = Array[:, 1]
    Array_v = Array[:, 2]
end