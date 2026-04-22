name="DAT-5";
processed=DAT5;
regions_raw=[targets "VTA-Ip" "SNc-Ip"];

%% processed
% 1. 숫자 행렬을 테이블로 변환하면서 열 이름 지정
BCtable = array2table(processed, 'VariableNames', targets);
Barcode_sequence=char(Bseq);
% 3. seq를 테이블로 만들어서 앞에 붙이기
T = [table(Barcode_sequence), BCtable];
% 4. CSV로 저장
writetable(T, name+'_processedmatrix.csv');

%%raw
% 1. 숫자 행렬을 테이블로 변환하면서 열 이름 지정
BCtable = array2table(barcodematrix, 'VariableNames', regions_raw);
Barcode_sequence=char(refbarcodes);
% 3. seq를 테이블로 만들어서 앞에 붙이기
T = [table(Barcode_sequence), BCtable];
% 4. CSV로 저장
writetable(T, name+'_rawmatrix.csv');