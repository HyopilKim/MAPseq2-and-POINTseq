function [corrMatrix, pvalMatrix] = nonzero_pairwise_corr(X,type)
    % X: (n_cells x n_genes) matrix
    % Output:
    %   corrMatrix: (n_genes x n_genes) pairwise Pearson correlation matrix
    %   pvalMatrix: (n_genes x n_genes) corresponding p-values

    [n_rows, n_genes] = size(X);
    corrMatrix = nan(n_genes);  % 상관계수 행렬 초기화
    pvalMatrix = nan(n_genes);  % p-value 행렬 초기화

    for i = 1:n_genes
        xi = X(:, i);
        for j = i:n_genes  % j = i부터 시작 (symmetric)
            xj = X(:, j);
            
            % 둘 다 0인 셀 제외
            valid_idx = ~(xi == 0 & xj == 0);
            xi_valid = xi(valid_idx);
            xj_valid = xj(valid_idx);

            % 유효한 샘플 수가 너무 적으면 계산 생략
            if sum(valid_idx) >= 3
                [r, p] = corr(xi_valid, xj_valid, 'Type', type);
                corrMatrix(i, j) = r;
                pvalMatrix(i, j) = p;
                
                % symmetric matrix니까 대칭 위치도 채워줌
                if i ~= j
                    corrMatrix(j, i) = r;
                    pvalMatrix(j, i) = p;
                end
            end
        end
    end
end