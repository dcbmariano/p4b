def smith_waterman(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-1):
    # Inicialização da matriz de pontuações
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0] * (m+1) for _ in range(n+1)]
    
    # Preenchimento da matriz de pontuações
    max_score = 0
    max_pos = None
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(0, match, delete, insert)
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)
    
    # Reconstrução do alinhamento
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_pos
    while score_matrix[i][j] != 0:
        if score_matrix[i][j] == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif score_matrix[i][j] == score_matrix[i-1][j] + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        elif score_matrix[i][j] == score_matrix[i][j-1] + gap_penalty:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    return max_score, aligned_seq1, aligned_seq2

# Exemplo de uso
seq1 = "ACAGTGT"
seq2 = "ACGCCGT"
score, aligned_seq1, aligned_seq2 = smith_waterman(seq1, seq2)


print("Seq1:", aligned_seq1)
print("Seq2:", aligned_seq2)
print("Score:", score)

