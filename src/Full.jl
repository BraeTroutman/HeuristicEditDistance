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

function traceback_total(fasta_filename)
    total = fill(0, 100, 100)
    sequences = Edist.get_fasta(fasta_filename)

    base = sequences[1]
    for query in sequences[2:end]
        score = construct(base[1:99], query[1:99])
        total += traceback_all(score, base, query)
    end

    return total
end

function traceback_total(fasta_filename, m, n)
    total = fill(0, m, n)
    sequences = Edist.get_fasta(fasta_filename)

    base = sequences[1]
    for query in sequences[2:end]
        score = construct(base[1:m-1], query[1:n-1])
        total += traceback_all(score, base, query)
    end

    return total
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

function visualize(fasta_filename, m, n)
    total = fill(0, m, n)
    sequences = Edist.get_fasta(fasta_filename)

    base = sequences[1]
    for query in sequences[2:end]
        score = construct(base[1:m-1], query[1:n-1])
        total += traceback_all(score, base, query)
    end

    Plots.heatmap(reverse(total[10:end,:], dims = 1))
end

function visualize(total::Matrix{Int})
    Plots.heatmap(reverse(total, dims = 1))
end

end # module Full

