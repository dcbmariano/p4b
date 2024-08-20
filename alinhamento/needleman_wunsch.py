def needleman_wunsch(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-1):
    # Inicialização da matriz de pontuações
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0] * (m+1) for _ in range(n+1)]
    
    # Inicialização da matriz de direções
    traceback_matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    
    # Preenchimento da primeira coluna com penalidades de lacuna
    for i in range(1, len(seq1) + 1):
        score_matrix[i][0] = i * gap_penalty
        traceback_matrix[i][0] = "up"
    
    # Preenchimento da primeira linha com penalidades de lacuna
    for j in range(1, len(seq2) + 1):
        score_matrix[0][j] = j * gap_penalty
        traceback_matrix[0][j] = "left"
    
    # Preenchimento das outras células
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            if score_matrix[i][j] == match:
                traceback_matrix[i][j] = "diagonal"
            elif score_matrix[i][j] == delete:
                traceback_matrix[i][j] = "up"
            else:
                traceback_matrix[i][j] = "left"
    
    # Reconstrução do alinhamento
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if traceback_matrix[i][j] == "diagonal":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == "up":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    return score_matrix[len(seq1)][len(seq2)], aligned_seq1, aligned_seq2

# Exemplo de uso
seq1 = "AAGCTAA"
seq2 = "ACGTGCG"
score, aligned_seq1, aligned_seq2 = needleman_wunsch(seq1, seq2)

print("Seq1:", aligned_seq1)
print("Seq2:", aligned_seq2)
print("Score:", score)
