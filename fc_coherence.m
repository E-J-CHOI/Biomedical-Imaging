%% Load Images
fmri_nii = nifti('filtered_func_data_tempfilt_m24_wmcsf_drift.nii'); %Preprocessed fMRI 데이터를 로드한다
fmri = double(fmri_nii.dat);    % 데이터를 fmri에 실수형태로 변환한다
atlas_nii = nifti('aal2example_func.nii');       %   Volumetric region-of-interest (ROI) 데이터를 로드한다
atlas = double(atlas_nii.dat);       % 데이터를 atlas에 실수형태로 변환한다
roi = load('aal_values.txt');              %Index value of region-of-interest (ROI) volume 데이터를 로드한다

%% Compute average BOLD signals
n_time = size(fmri,4);   % fmri의 네번째 dimension의 길이를 n_time으로 지정한다
n_region = length(roi);        % 가장 긴 차원의 roi의 개수를 n_region으로 지정한다.
ts = zeros(n_time,n_region);    % n_time X n_region의 모두 0으로 이루어진 행렬을 ts로 지정한다.
for i = 1: 90                      %i를 1부터 90까지 반복하여 실행한다.
    idx = bsxfun(@plus, 0:numel(atlas):numel(fmri)-numel(atlas), find(atlas==roi(i)));  % fmri개수를 atlas개수 단위로 증분시킨 값과 atlas에서 roi(i)값에 해당하는 인덱스를 더해준다. 
    ts(:,i) = mean(fmri(idx),1)';  % idx에 해당하는 fmri 데이터의 각 열에 있는 평균값 행을 전치하여 ts의 i열에 입력한다.
end

%% Functional connectivity based on correlation
fc_correlation = corr(ts);   % ts에서의 roi의 correlation 구하기
 
%% Functional connectivity based on coherence
fc_coherence = zeros(n_region, n_region);    % fc_coherence로 n_region X n_region 제로 행렬을 생성한다  
for i = 1:90               %i가 1부터 90까지 반복한다
    for j = 1:90                   %j에 1부터 90까지 반복한다.
        Cxy = mscohere(ts(:,i), ts(:,j), [], [], 128, 1/3);    %mscohere 함수를 이용하여 비교해준다 % 1/3 은 샘플 레이트로 단위 시간당 샘플 개수를 뜻한다. repetition time(TR)이 3000ms 이므로 1/3Hz이다.    
        fc_coherence(i,j) = mean(Cxy);        % fc_coherence(i,j) 자리에 Cxy의 평균값을 넣어준다.
    end
end
fc_coherence = fc_coherence + fc_coherence' + eye(n_region);    % fc_coherence, fc_coherence의 전치,단위행렬을 더해준다.


%% Plot
close all;                    % 모두 닫기(종료)
figure('Color', 'w');    %배경색 : white
subplot(1,2,1);            %  1열 2행의 그래프 중 첫번째 그래프 plot
imagesc(fc_correlation);        % fc_correlation 그래프를 이미지로 불러오기
axis equal off;                   % axis line 보이지 않도록 설정
title('FC based on correlation');           %그래프 제목 :'FC based on correlation'
subplot(1,2,2);           % 1열 2행의 그래프 중 두번째 그래프 plot
imagesc(fc_coherence);        % fc_coherence 그래프를 이미지로 불러오기  
axis equal off;                % axis line 보이지 않도록 설정
title('FC based on coherence');       %그래프 제목 :'FC based on coherence'


