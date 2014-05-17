function score = calculate_greedy(M)

% calculates the greedy mapping score
score = sum(max(M))/sum(sum(M));
