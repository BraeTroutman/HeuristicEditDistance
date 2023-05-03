module Edist
import FastaIO
export Full, Bounded, Hirschberg, get_fasta, align, score

include("Full.jl")
include("Bounded.jl")
include("Hirschberg.jl")

"""
    align(moduleName, sequence, query)
align the two sequences using the edit distance algorithm implemented in the specified module--
Full, Bounded, or Hirschberg

hides the underlying implementation parameters to make calculations conform to one unified interface
"""
function align(moduleName, sequence, query, prnt::Bool = false)
    if moduleName == Hirschberg
        (sal, qal) = Hirschberg.alignment(sequence, query)

        prnt && print("$(sal)\n$(qal)\n")

        score = 0
        for (a, b) in zip(sal, qal)
            if a == '-' || b == '-'
                score -= 2
            elseif a == b
                score += 1
            else
                score -= 1
            end
        end

        return (score = score, seq_alignment = sal, query_alignment = qal, memory_used = 0)
    end
    res = moduleName.construct(sequence, query)
    (sal, qal, score) = (
        (moduleName == Bounded) ? moduleName.alignment(res..., sequence, query, prnt) :
        moduleName.alignment(res, sequence, query, prnt)
    )
    return (
        score = score,
        seq_alignment = sal,
        query_alignment = qal,
        memory_used = Base.summarysize(res),
    )
end


"""
    score(moduleName, sequence, query)

backend-agnostic function to calculate score for alignment
"""
function score(moduleName, sequence, query)
    return moduleName.score(sequence, query)
end


""" get_fasta(filename)

return the sequences encoded in the file `filename`
"""
function get_fasta(filename)
    return getindex.(FastaIO.readfasta(filename), 2)
end

end # module Edist
