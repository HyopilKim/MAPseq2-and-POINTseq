
% SHC 기반 계층적 클러스터링 + leafOrder 정렬 버전

%% 1. 데이터 준비 및 linkage 생성
Z = linkage(data, 'ward', 'euclidean');
leafOrder = optimalleaforder(Z, pdist(data));
sortedData = data(leafOrder, :);

%% 2. SHC 기반 클러스터링 수행 (사용자 정의 함수)
% 입력: Z, shc_pval, alpha, k, data
% 출력: T (cluster labels), cluster_roots (internal node IDs)
[T, cluster_roots] = hierarchical_withSHC(Z, shc_pval, alpha, k, data);

%% 3. leafOrder 기반으로 클러스터 결과 정렬
T_sorted = T(leafOrder);  % leafOrder 순서대로 클러스터 번호 정렬

%% 4. leafOrder에서 등장 순서대로 클러스터 번호 재정의
numClusters = max(T_sorted);
firstAppearance = zeros(numClusters, 1);

for i = 1:numClusters
    firstAppearance(i) = find(T_sorted == i, 1, 'first');
end

[~, newOrder] = sort(firstAppearance);  % leafOrder 기준 정렬된 클러스터 순서

T_reordered = zeros(size(T_sorted));
for i = 1:numClusters
    T_reordered(T_sorted == newOrder(i)) = i;
end

%% 5. cluster_roots도 동일한 순서로 정렬
cluster_roots_sorted = cluster_roots(newOrder);

%% 결과 변수
% sortedData         : leafOrder 기준으로 정렬된 데이터
% T_reordered        : sortedData 순서에 맞춰진 클러스터 번호 (1 ~ numClusters)
% cluster_roots_sorted : 각 클러스터를 대표하는 노드 번호 (정렬된 순서대로)
