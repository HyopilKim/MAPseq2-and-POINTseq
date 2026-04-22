function nearestNeighborData = findnearestDNA(dnaSequences)
%FINDNEARESTDNA 각 DNA 시퀀스에 대해 가장 가까운 다른 시퀀스의 인덱스와 해밍 거리를 계산합니다.
%
% 입력:
%   dnaSequences - DNA 시퀀스를 포함하는 문자형 배열 (char matrix) 또는 셀 배열 (cell array of char vectors).
%                  문자형 배열인 경우 각 행은 하나의 시퀀스여야 하며 길이는 모두 같아야 합니다.
%
% 출력:
%   nearestNeighborData - 두 컬럼으로 구성된 행렬.
%                         첫 번째 컬럼: 가장 가까운 다른 시퀀스의 행 인덱스.
%                         두 번째 컬럼: 해당 시퀀스와의 해밍 거리.

    % 입력 시퀀스의 형식 확인 및 변환
    if ischar(dnaSequences)
        % 문자형 배열인 경우
        X = dnaSequences;
        numSequences = size(X, 1);
        sequenceLength = size(X, 2);
    elseif iscell(dnaSequences)
        % 셀 배열인 경우 (길이가 같다고 가정)
        X = char(dnaSequences); % knnsearch를 위해 char matrix로 변환
        numSequences = size(X, 1);
        sequenceLength = size(X, 2);
    else
        error('입력은 문자형 배열 또는 문자형 벡터의 셀 배열이어야 합니다.');
    end
    
    % KNNSearch를 위한 입력 데이터는 행이 관측치(시퀀스)이고 열이 변수(위치)여야 합니다.
    % 이미 char matrix X가 이 형태를 따르고 있습니다.
    
    % knnsearch를 사용하여 가장 가까운 이웃을 찾습니다.
    % 'K', 2: 자기 자신(K=1)을 제외한 가장 가까운 이웃 하나(K=2)를 찾습니다.
    % 'Distance', 'hamming': 해밍 거리를 사용합니다. 해밍 거리는 일치하지 않는 좌표의 비율로 반환됩니다.
    % 'Exclude', 'self': 쿼리 데이터와 모델 데이터가 같으므로 자기 자신을 제외합니다.
    %                   이 옵션은 k=1로 설정하고 사용해야 하지만, 여기서는 k=2로 설정하고 
    %                   가장 가까운 2개 중 자기 자신을 제외한 1개를 선택하는 방식으로 처리합니다.
    %                   knnsearch는 기본적으로 자기 자신을 가장 가까운 이웃으로 간주합니다 (거리가 0이므로).
    
    % 자기 자신을 제외한 가장 가까운 이웃 하나(K=1)를 찾기 위해 K=2를 사용합니다.
    K = 2;
    [Idx, D] = knnsearch(X, X, 'K', K, 'Distance', 'hamming');
    
    % 결과 처리:
    % Idx의 첫 번째 컬럼(Idx(:,1))은 항상 자기 자신(인덱스)이고 거리는 0입니다 (D(:,1)).
    % 따라서 가장 가까운 *다른* 시퀀스는 두 번째 컬럼(Idx(:,2))에 있으며, 거리는 D(:,2)입니다.
    
    % 각 행에 대해 가장 가까운 이웃의 인덱스와 거리를 추출합니다.
    % 1. 인덱스 (가장 가까운 다른 시퀀스의 행 번호)
    nearestIndex = Idx(:, 2);
    
    % 2. 해밍 거리 (불일치하는 요소의 비율)
    % 해밍 거리는 불일치하는 문자의 개수로 표시되기를 원한다면, 비율에 시퀀스 길이를 곱합니다.
    hammingDistance = D(:, 2) * sequenceLength;
    
    % 결과를 두 컬럼 행렬로 반환합니다.
    nearestNeighborData = [nearestIndex, hammingDistance];
end