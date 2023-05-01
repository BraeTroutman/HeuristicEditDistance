module Edist
import FastaIO, Plots
export Full, get_fasta

""" get_fasta(filename)

return the sequences encoded in the file `filename`
"""
function get_fasta(filename)
    return getindex.(FastaIO.readfasta(filename), 2)
end

module Full
	"""
	    construct(sequence, query)
	
	build the dynamic programming matrix for alignment of two strings `sequence` and `query`
	
	scores with a gap penalty of `-2` and a match/mismatch of `+-1`
	"""
	function construct(sequence, query)
	    M = length(sequence)
	    N = length(query)
	
	    align_mtx = zeros(Int32, M+1, N+1)
	    align_mtx[1,:] .= (0:N) * -2
	    align_mtx[:,1] .= (0:M) * -2
	
	    # iterate over each column
	    for j in 2:(N+1)
	        # iterate over the elements of the current column
	        for i in 2:(M+1)
	            match = sequence[i-1] == query[j-1] ? 1 : -1
	            align_mtx[i,j] = max(
	                                 align_mtx[i-1,j-1] + match, 
	                                 align_mtx[i-1,j]   - 2,
	                                 align_mtx[i,j-1]   - 2
	                                )
	        end
	    end
	    return align_mtx
	end
    
    function score(sequence, query)
        return construct(sequence, query)[end,end]
    end

    function alignment(score_mtx, s, q, prnt::Bool=false)
        sequence = ""
        query    = ""

        M, N = size(score_mtx)
        i = M
        j = N

        while i != 1 || j != 1
    	    match       = (i > 1 && j > 1) && (score_mtx[i-1, j-1] + 1) == score_mtx[i, j]
	        mismatch    = (i > 1 && j > 1) && (score_mtx[i-1, j-1] - 1) == score_mtx[i, j]
	        query_gap   = (i > 1)          && (score_mtx[i-1, j]   - 2) == score_mtx[i, j]
	        seq_gap     = (j > 1)          && (score_mtx[i, j-1]   - 2) == score_mtx[i, j]
            
            if match || mismatch
                sequence = s[i-1] * sequence 
                query    = q[j-1] * query
                i -= 1
                j -= 1
            elseif query_gap
                sequence = s[i-1] * sequence
                query    = '-'    * query
                i -= 1
            else
                sequence = '-'    * sequence
                query    = q[j-1] * query
                j -= 1
            end
        end
        
        prnt && println(sequence)
        prnt && println(query)

        return (sequence, query) 
    end

	"""
	    traceback(score_mtx)
	
	generate a matrix the same shape of the scoring matrix passed in, with matrix
	entries along the traceback path set to 1 and all others set to 0
	"""
	function traceback(score_mtx)
	    al_mtx = fill(Int16(0), size(score_mtx))
	    traceback(score_mtx, al_mtx)
	    return al_mtx
	end
	
	function traceback(score_mtx, al_mtx)
	    M, N = size(score_mtx)
	    
	    al_mtx[M,N] = 1
	
	    match       = (M > 1 && N > 1) && (score_mtx[M-1, N-1] + 1) == score_mtx[M, N]
	    mismatch    = (M > 1 && N > 1) && (score_mtx[M-1, N-1] - 1) == score_mtx[M, N]
	    query_gap   = (M > 1)          && (score_mtx[M-1, N]   - 2) == score_mtx[M, N]
	    seq_gap     = (N > 1)          && (score_mtx[M, N-1]   - 2) == score_mtx[M, N]
	
	    if match
	        traceback(view(score_mtx, 1:(M-1),1:(N-1)), view(al_mtx, 1:(M-1),1:(N-1)))
	    elseif mismatch
	        traceback(view(score_mtx, 1:(M-1), 1:(N-1)), view(al_mtx, 1:(M-1),1:(N-1)))
	    elseif query_gap
	        traceback(view(score_mtx, 1:(M-1), :), view(al_mtx, 1:(M-1), :))
	    elseif seq_gap
	        traceback(view(score_mtx, :, 1:(N-1)), view(al_mtx, :, 1:(N-1)))
	    end
	
	    return
	end
	
	"""
	    traceback_all(score_mtx)
	
	traceback every possible path for the given scoring matrix--
	
	note that this can be potentially very computationally intensive for large matrices, as the
	branching factor for each entry is 3
	"""
	function traceback_all(score_mtx)
	    al_mtx = fill(Int16(0), size(score_mtx))
	    traceback_all(score_mtx, al_mtx)
	    return al_mtx
	end
	
	function traceback_all(score_mtx, al_mtx)
	    M, N = size(score_mtx)
	    
	    al_mtx[M,N] = 1
	
	    match       = (M > 1 && N > 1) && (score_mtx[M-1, N-1] + 1) == score_mtx[M, N]
	    mismatch    = (M > 1 && N > 1) && (score_mtx[M-1, N-1] - 1) == score_mtx[M, N]
	    query_gap   = (M > 1)          && (score_mtx[M-1, N]   - 2) == score_mtx[M, N]
	    seq_gap     = (N > 1)          && (score_mtx[M, N-1]   - 2) == score_mtx[M, N]
	
	    if match
	        traceback(view(score_mtx, 1:(M-1),1:(N-1)), view(al_mtx, 1:(M-1),1:(N-1)))
	    end
	    if mismatch
	        traceback(view(score_mtx, 1:(M-1), 1:(N-1)), view(al_mtx, 1:(M-1),1:(N-1)))
	    end
	    if query_gap
	        traceback(view(score_mtx, 1:(M-1), :), view(al_mtx, 1:(M-1), :))
	    end
	    if seq_gap
	        traceback(view(score_mtx, :, 1:(N-1)), view(al_mtx, :, 1:(N-1)))
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
	
	    sequences = get_fasta(fasta_filename)
	
	    base = sequences[1]
	    for query in sequences[2:end]
	        score = construct(base[1:99], query[1:99])
	        total += traceback_all(score)
	    end
	
	    Plots.heatmap(reverse(total, dims=1))
	end

end # module Full

module Hirschberg

end # module Hirschberg

module Bounded
    import ..Full
    
    function score(sequence::String, query::String)
        return construct(sequence, query)[3][end,end]
    end

    function construct(sequence::String, query::String)
        M = length(sequence)
        N = length(query)
        
        k = 3
        d = abs(M - N)

        if M > N
            return construct(query, sequence, k, d)  
        else
            return construct(sequence, query, k, d)
        end
    end

    function construct(sequence::String, query::String, k::Int, d::Int)
        # assume that the query is longer than the sequence
        M = length(sequence)+1
        N = length(query)+1
        
        top_left = Full.construct(sequence[1:k-1], query[1:k+d-1])
        right_frontier = zeros(Int64, M-k, k-1)
        bottm_frontier = zeros(Int64, M-k, d+k)

        # calculate the first right frontier
        match = (query[k+d] == sequence[1]) ? 1 : -1
        right_frontier[1,1] = max(
                                   top_left[1, k+d] + match,
                                   top_left[2, k+d] - 2
                                  )
        for i in 2:k-1
            match = (query[k+d] == sequence[i]) ? 1 : -1
            right_frontier[1,i] = max(
                                      top_left[i+1, k+d] - 2,
                                      right_frontier[1,i-1] -2,
                                      top_left[i, k+d] + match
                                     )
        end

        # calculate the first bottom frontier
        match = (query[1] == sequence[k]) ? 1 : -1
        bottm_frontier[1,1] = max(
                                  top_left[k,1] + match,
                                  top_left[k,2] - 2
                                 )
        for i in 2:d+k-1
            match = (query[i] == sequence[k]) ? 1 : -1
            bottm_frontier[1,i] = max(
                                      top_left[k, i] + match,
                                      top_left[k, i+1] - 2,
                                      bottm_frontier[1, i-1] - 2
                                     )
        end
        match = (query[d+k] == sequence[k]) ? 1 : -1
        bottm_frontier[1,end] = max(
                                    top_left[end,end] + match,
                                    bottm_frontier[1,end-1] - 2,
                                    right_frontier[1,end] - 2
                                   )

        for f in 2:M-k
            # right frontier 
            # special case-- first element in frontier has no element directly above
            match = (query[d+k+f-1] == sequence[f]) ? 1 : -1
            right_frontier[f, 1] = max(
                                              right_frontier[f-1, 1] + match,
                                              right_frontier[f-1, 2] - 2
                                             )
            # generic case
            for i in 2:k-2
                match = (query[d+k+f-1] == sequence[f+i-1]) ? 1 : -1
                right_frontier[f, i] = max(
                                                  right_frontier[f-1, i] + match,
                                                  right_frontier[f-1, i+1] - 2,
                                                  right_frontier[f, i-1] - 2
                                                 )
            end
            # special case-- last element in frontier must access last element in bottom frontier
            match = (query[d+k+f-1] == sequence[f+k-2]) ? 1 : -1
            right_frontier[f, end] = max(
                                                right_frontier[f-1, end] + match,
                                                bottm_frontier[f-1, end] - 2,
                                                right_frontier[f, end-1] - 2
                                               )
            # bottom frontier
            # special case-- first element in frontier has no element directly to left
            match = (query[f] == sequence[k+f-1]) ? 1 : -1
            bottm_frontier[f, 1] = max(
                                       bottm_frontier[f-1, 1] + match,
                                       bottm_frontier[f-1, 2] - 2
                                      )

            # generic case
            for i in 2:d+k-1
                match = (query[f+i-1] == sequence[k+f-1]) ? 1 : -1
                bottm_frontier[f, i] = max(
                                           bottm_frontier[f-1, i] + match,
                                           bottm_frontier[f-1, i+1] - 2,
                                           bottm_frontier[f, i-1] - 2
                                          )
            end

            # special case-- last element in frontier must access last element of right_frontier
            match = (query[f+d+k-1] == sequence[k+f-1]) ? 1 : -1
            bottm_frontier[f,end] = max(
                                        bottm_frontier[f-1, end] + match,
                                        right_frontier[f, end] - 2,
                                        bottm_frontier[f, end-1] - 2
                                       )
        end
        return top_left, right_frontier, bottm_frontier
    end

end # module Bounded

end # module Edist

