function randomized=shuffincol(matrix)

[m, n] = size(matrix);
randomized = zeros(m, n);  % 결과 행렬 초기화

max_tries = 100;
tries = 0;

while tries < max_tries
    for col = 1:n
        vals = matrix(:, col);           % 해당 열의 전체 값
        shuffled = vals(randperm(m));      % 행 안에서 무작위 섞기
        randomized(:, col) = shuffled; % 새로운 열에 넣기
    end
    
    % 모든 행에 하나 이상의 0이 아닌 값이 있는지 확인
    row_nonzero_counts = sum(randomized ~= 0, 2);
    
    if all(row_nonzero_counts > 0)
        break;  % 조건 만족 시 종료
    end
    
    tries = tries + 1;
end
