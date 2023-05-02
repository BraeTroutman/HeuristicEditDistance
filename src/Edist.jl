module Edist
import FastaIO
export Full, Bounded, Hirschberg, get_fasta, align, score

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

        return (score = score, seq_alignment = sal, query_alignment = qal)
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

module Full
import ..Edist, Plots

"""
    construct(sequence, query)

build the dynamic programming matrix for alignment of two strings `sequence` and `query`

scores with a gap penalty of `-2` and a match/mismatch of `+-1`
"""
function construct(sequence, query)
    M = length(sequence)
    N = length(query)

    align_mtx = zeros(Int, M + 1, N + 1)
    align_mtx[1, :] .= (0:N) * -2
    align_mtx[:, 1] .= (0:M) * -2

    # iterate over each column
    for j = 2:(N+1)
        # iterate over the elements of the current column
        for i = 2:(M+1)
            match = sequence[i-1] == query[j-1] ? 1 : -1
            align_mtx[i, j] = max(
                align_mtx[i-1, j-1] + match,
                align_mtx[i-1, j] - 2,
                align_mtx[i, j-1] - 2,
            )
        end
    end
    return align_mtx
end

"""
    score(sequence, query)

calculate the optimal score for the alignment of the two given sequences
"""
function score(sequence, query)
    res = construct(sequence, query)
    return (score = res[end, end], memory_used = Base.summarysize(res))
end

"""
    alignment(score_mtx, s, q; prnt)

return the score and alignment strings of the two sequences scored with the given score matrix

if prnt is passed as true the alignment will also be printed to the console
"""
function alignment(score_mtx, s, q, prnt::Bool = false)
    sequence = ""
    query = ""

    M, N = size(score_mtx)
    i = M
    j = N

    while i != 1 || j != 1
        match =
            (i > 1 && j > 1) &&
            (score_mtx[i-1, j-1] + 1) == score_mtx[i, j] &&
            s[i-1] == q[j-1]
        mismatch =
            (i > 1 && j > 1) &&
            (score_mtx[i-1, j-1] - 1) == score_mtx[i, j] &&
            s[i-1] != q[j-1]
        query_gap = (i > 1) && (score_mtx[i-1, j] - 2) == score_mtx[i, j]
        seq_gap = (j > 1) && (score_mtx[i, j-1] - 2) == score_mtx[i, j]

        if match || mismatch
            sequence = s[i-1] * sequence
            query = q[j-1] * query
            i -= 1
            j -= 1
        elseif query_gap
            sequence = s[i-1] * sequence
            query = '-' * query
            i -= 1
        else
            sequence = '-' * sequence
            query = q[j-1] * query
            j -= 1
        end
    end

    prnt && println(sequence)
    prnt && println(query)

    return (sequence, query, score_mtx[end, end])
end

"""
    traceback(score_mtx)

generate a matrix the same shape of the scoring matrix passed in, with matrix
entries along the traceback path set to 1 and all others set to 0
"""
function traceback(score_mtx, s, q)
    al_mtx = fill(Int16(0), size(score_mtx))
    traceback(score_mtx, al_mtx, s, q)
    return al_mtx
end

function traceback(score_mtx, al_mtx, s, q)
    M, N = size(score_mtx)

    al_mtx[M, N] = 1

    match = false
    if M > 1 && N > 1
        match = (s[M-1] == q[N-1]) ? 1 : -1
        match = (M > 1 && N > 1) && (score_mtx[M-1, N-1] + match) == score_mtx[M, N]
    end
    query_gap = (M > 1) && (score_mtx[M-1, N] - 2) == score_mtx[M, N]
    seq_gap = (N > 1) && (score_mtx[M, N-1] - 2) == score_mtx[M, N]

    if match
        traceback(view(score_mtx, 1:(M-1), 1:(N-1)), view(al_mtx, 1:(M-1), 1:(N-1)), s, q)
    elseif query_gap
        traceback(view(score_mtx, 1:(M-1), :), view(al_mtx, 1:(M-1), :), s, q)
    elseif seq_gap
        traceback(view(score_mtx, :, 1:(N-1)), view(al_mtx, :, 1:(N-1)), s, q)
    end

    return
end

"""
    traceback_all(score_mtx)

traceback every possible path for the given scoring matrix--

note that this can be potentially very computationally intensive for large matrices, as the
branching factor for each entry is 3
"""
function traceback_all(score_mtx, s, q)
    al_mtx = fill(Int16(0), size(score_mtx))
    traceback_all(score_mtx, al_mtx, s, q)
    return al_mtx
end

function traceback_all(score_mtx, al_mtx, s, q)
    M, N = size(score_mtx)

    al_mtx[M, N] = 1

    match = false
    if M > 1 && N > 1
        match = (s[M-1] == q[N-1]) ? 1 : -1
        match = (M > 1 && N > 1) && (score_mtx[M-1, N-1] + match) == score_mtx[M, N]
    end

    query_gap = (M > 1) && (score_mtx[M-1, N] - 2) == score_mtx[M, N]
    seq_gap = (N > 1) && (score_mtx[M, N-1] - 2) == score_mtx[M, N]

    if match
        traceback(view(score_mtx, 1:(M-1), 1:(N-1)), view(al_mtx, 1:(M-1), 1:(N-1)), s, q)
    end
    if query_gap
        traceback(view(score_mtx, 1:(M-1), :), view(al_mtx, 1:(M-1), :), s, q)
    end
    if seq_gap
        traceback(view(score_mtx, :, 1:(N-1)), view(al_mtx, :, 1:(N-1)), s, q)
    end

    return
end

"""
    visualize(fasta_filename)

take a file in FASTA format and align the first sequence within with every other--
then produce a heatmap demonstrating where the optimal alignments fall in the dynamic
programming matrix
"""
function visualize(fasta_filename)
    total = fill(0, 100, 100)

    sequences = Edist.get_fasta(fasta_filename)

    base = sequences[1]
    for query in sequences[2:end]
        score = construct(base[1:99], query[1:99])
        total += traceback_all(score, base, query)
    end

    Plots.heatmap(reverse(total, dims = 1))
end

end # module Full

module Hirschberg
import ..Edist

function score(s, q)
    sal, qal = alignment(s, q)
    score = 0
    for (a, b) in zip(sal, qal)
        if a == b
            score += 1
        elseif a == '-' || b == '-'
            score -= 2
        else
            score -= 1
        end
    end

    return (score = score)
end

function needleman_wunsch_score(s, q)
    M = length(s)
    N = length(q)

    scores = zeros(Int, 2, N + 1)
    scores[1, :] .= (0:N) * -2

    for i = 2:M+1
        for j = 1:N+1
            if j == 1
                scores[2, j] = scores[1, j] - 2
            else
                match = s[i-1] == q[j-1] ? 1 : -1
                scores[2, j] =
                    max(scores[1, j-1] + match, scores[1, j] - 2, scores[2, j-1] - 2)
            end
        end
        scores[1, :] .= @view scores[2, :]
    end

    return scores[2, :]
end

function alignment(sequence, query)
    M = length(sequence)
    N = length(query)

    sal = ""
    qal = ""

    if M == 0
        for i = 1:N
            sal *= '-'
            qal *= query[i]
        end
    elseif N == 0
        for i = 1:M
            sal *= sequence[i]
            qal *= '-'
        end
    elseif M == 1 || N == 1
        score, sal, qal = Edist.align(Edist.Full, sequence, query)
    else
        leftSequence = sequence[1:div(end, 2)]
        rightSequence = sequence[div(end, 2)+1:end]

        leftScore = needleman_wunsch_score(leftSequence, query)
        rightScore = needleman_wunsch_score(reverse(rightSequence), reverse(query))

        mid = argmax(leftScore + reverse(rightScore)) - 1

        leftSal, leftQal = alignment(leftSequence, query[1:mid])
        rightSal, rightQal = alignment(rightSequence, query[mid+1:end])
        sal = leftSal * rightSal
        qal = leftQal * rightQal
    end
    return sal, qal
end

end # module Hirschberg

module Bounded
import ..Full

@enum Front Bottom Right

"""
    score(sequence, query; k)

return the score of the given sequence and query calculated using a spatially bounded dp matrix--
that is, limiting the search space to a smaller subset of entries around the diagonal of the matrix

`k` specifies the number of elements to the left of the diagonal to consider
"""
function score(sequence::String, query::String)
    res = construct(sequence, query)
    return (score = res[2][end, end], memory_used = Base.summarysize(res))
end

"""
    construct(sequence, query)

calculate the scoring matrices for the given sequence and query. 

returns the top left portion of the full scoring matrix, plus the full set of frontiers from across iterations.
"""
function construct(sequence::String, query::String)
    M = length(sequence)
    N = length(query)

    d = abs(M - N)
    k = d > 4 ? div(d, 2) : 3

    if M > N
        return construct(query, sequence, k, d)
    else
        return construct(sequence, query, k, d)
    end
end

function construct(sequence::String, query::String, k::Int)
    M = length(sequence)
    N = length(query)

    d = abs(M - N)

    if M > N
        return construct(query, sequence, k, d)
    else
        return construct(sequence, query, k, d)
    end
end

function construct(sequence::String, query::String, k::Int, d::Int)
    # assume that the query is longer than the sequence
    M = length(sequence) + 1
    N = length(query) + 1

    top_left = Full.construct(sequence[1:k-1], query[1:k+d-1])
    right_frontier = zeros(Int64, M - k, k - 1)
    bottm_frontier = zeros(Int64, M - k, d + k)

    # calculate the first right frontier
    match = (query[k+d] == sequence[1]) ? 1 : -1
    right_frontier[1, 1] = max(top_left[1, k+d] + match, top_left[2, k+d] - 2)
    for i = 2:k-1
        match = (query[k+d] == sequence[i]) ? 1 : -1
        right_frontier[1, i] = max(
            top_left[i+1, k+d] - 2,
            right_frontier[1, i-1] - 2,
            top_left[i, k+d] + match,
        )
    end

    # calculate the first bottom frontier
    match = (query[1] == sequence[k]) ? 1 : -1
    bottm_frontier[1, 1] = max(top_left[k, 1] + match, top_left[k, 2] - 2)
    for i = 2:d+k-1
        match = (query[i] == sequence[k]) ? 1 : -1
        bottm_frontier[1, i] =
            max(top_left[k, i] + match, top_left[k, i+1] - 2, bottm_frontier[1, i-1] - 2)
    end
    match = (query[d+k] == sequence[k]) ? 1 : -1
    bottm_frontier[1, end] = max(
        top_left[end, end] + match,
        bottm_frontier[1, end-1] - 2,
        right_frontier[1, end] - 2,
    )

    for f = 2:M-k
        # right frontier 
        # special case-- first element in frontier has no element directly above
        match = (query[d+k+f-1] == sequence[f]) ? 1 : -1
        right_frontier[f, 1] =
            max(right_frontier[f-1, 1] + match, right_frontier[f-1, 2] - 2)
        # generic case
        for i = 2:k-2
            match = (query[d+k+f-1] == sequence[f+i-1]) ? 1 : -1
            right_frontier[f, i] = max(
                right_frontier[f-1, i] + match,
                right_frontier[f-1, i+1] - 2,
                right_frontier[f, i-1] - 2,
            )
        end
        # special case-- last element in frontier must access last element in bottom frontier
        match = (query[d+k+f-1] == sequence[f+k-2]) ? 1 : -1
        right_frontier[f, end] = max(
            right_frontier[f-1, end] + match,
            bottm_frontier[f-1, end] - 2,
            right_frontier[f, end-1] - 2,
        )
        # bottom frontier
        # special case-- first element in frontier has no element directly to left
        match = (query[f] == sequence[k+f-1]) ? 1 : -1
        bottm_frontier[f, 1] =
            max(bottm_frontier[f-1, 1] + match, bottm_frontier[f-1, 2] - 2)

        # generic case
        for i = 2:d+k-1
            match = (query[f+i-1] == sequence[k+f-1]) ? 1 : -1
            bottm_frontier[f, i] = max(
                bottm_frontier[f-1, i] + match,
                bottm_frontier[f-1, i+1] - 2,
                bottm_frontier[f, i-1] - 2,
            )
        end

        # special case-- last element in frontier must access last element of right_frontier
        match = (query[f+d+k-1] == sequence[k+f-1]) ? 1 : -1
        bottm_frontier[f, end] = max(
            bottm_frontier[f-1, end] + match,
            right_frontier[f, end] - 2,
            bottm_frontier[f, end-1] - 2,
        )
    end
    return top_left, bottm_frontier, right_frontier
end

"""
    alignment(top_left, frontier_b, frontier_r, s, q; prnt)

return the alignment strings for the two sequences `s` and `q`, calculated from the arguments
`top_left`, `frontier_b`, `frontier_r` which are the results of the `Bounded.construct` function

# Example
```julia-repl
julia> score = Bounded.construct("hello", "hello world");
julia> Bounded.alignment(score..., "hello", "hello world")
("hello", "hello world", -7)
```
"""
function alignment(top_left, frontier_b, frontier_r, s, q, prnt::Bool = false)
    m = length(s)
    n = length(q)
    sequence = m < n ? s : q
    query = m < n ? q : s

    M = length(sequence)
    N = length(query)

    (k, d) = size(top_left)
    d -= k

    score_mtx = fill(typemin(Int), M + 1, N + 1)
    score_mtx[1:k, 1:k+d] .= top_left

    for f = 1:M-k+1
        score_mtx[k+f, f+1:f+d+k] = @view frontier_b[f, :]
        score_mtx[f+1:f+k-1, d+k+f] = @view frontier_r[f, :]
    end

    sal, qal, score = Full.alignment(score_mtx, sequence, query, prnt)

    return m < n ? (sal, qal, score) : (qal, sal, score)
end

end # module Bounded

end # module Edist
