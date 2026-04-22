function seqs = multiplexedlibrary(n, l, m)
% make_dna_sequences  n개의 DNA 시퀀스를 만든다.
%   - 각 시퀀스 길이는 l
%   - 앞 4개 위치에 ATCG가 한 번씩 등장 (순열)
%   - 모든 쌍의 Hamming distance >= m
%
% 사용 예:
%   seqs = make_dna_sequences(10, 20, 5);

    if l < 4
        error('시퀀스 길이 l은 최소 4 이상이어야 합니다.');
    end

    alphabet = 'ATCG';
    numAlphabet = numel(alphabet);

    % 결과 저장용 (char 배열: n x l)
    seqs = repmat('N', n, l);  % 초기값 'N'(dummy)

    maxTrials = 1e5;  % 무한 루프 방지용

    for i = 1:n
        trial = 0;
        while true
            trial = trial + 1;
            if trial > maxTrials
                error('조건을 만족하는 시퀀스를 만들 수 없습니다. n, l, m을 다시 확인하세요.');
            end

            % 1. 앞 4개는 ATCG 순열
            prefix = alphabet(randperm(4));

            % 2. 나머지는 랜덤으로 ATCG 중 하나
            if l > 4
                rest = alphabet(randi(numAlphabet, 1, l-4));
                cand = [prefix rest];
            else
                cand = prefix;
            end

            % 3. Hamming distance 조건 검사
            if i == 1
                % 첫 번째 시퀀스는 무조건 채택
                seqs(i, :) = cand;
                break;
            else
                % 이전까지 만들어진 시퀀스들과의 Hamming distance 계산
                prevSeqs = seqs(1:i-1, :);

                % 각 행과 cand의 Hamming distance: 열 방향 비교 후 합산
                dists = sum(prevSeqs ~= cand, 2);  % (i-1)x1 벡터

                if all(dists >= m)
                    seqs(i, :) = cand;
                    break;
                end
            end

        end % while
    end % for
end