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

