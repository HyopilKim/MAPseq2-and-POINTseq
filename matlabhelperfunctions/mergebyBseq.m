function merged_mat = mergebyBseq(mat1, seq1, mat2, seq2)

    % 1) Create the union of all sequence identifiers
    all_seq = unique([seq1; seq2], 'stable','rows');

    % 2) Initialize result matrices filled with 0
    % Rows correspond to all_seq
    res1 = zeros(length(all_seq), size(mat1, 2));
    res2 = zeros(length(all_seq), size(mat2, 2));

    % 3) Find row indices mapping original seq to all_seq
    [~, idx1] = ismember(seq1, all_seq,'rows');
    [~, idx2] = ismember(seq2, all_seq,'rows');

    % 4) Fill in values for mat1
    res1(idx1, :) = mat1;

    % 5) Fill in values for mat2
    res2(idx2, :) = mat2;

    % 6) Concatenate matrices column-wise
    merged_mat = [res1, res2];

end