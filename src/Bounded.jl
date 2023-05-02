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

    if M > N
        k = d < div(N, 2) && d != 0 ? d : 2
        k = d == 1 ? 2 : k
        return construct(query, sequence, k, d)
    else
        k = d < div(M, 2) && d != 0 ? d : 2 
        k = d == 1 ? 2 : k
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
            max(right_frontier[f-1, 1] + match, k > 2 ? right_frontier[f-1, 2] - 2 : typemin(Int))
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
            k > 2 ? right_frontier[f, end-1] - 2 : typemin(Int)
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
