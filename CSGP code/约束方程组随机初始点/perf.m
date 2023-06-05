 function perf(T,logplot)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
%%% Failures on a given problem are represented by a NaN.   ���ܹؼ���
% The optional argument logplot is used to produce a 
% log (base 2) performance plot.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004
% clear all
% clc
% T=rand(50,4)
% nargin=2
% logplot=1;
if (nargin < 2) logplot = 0; end

%plot LineSpec����
%colors  = ['m' 'b' 'r' 'g' 'k' 'c' 'y'];% 'm'Ʒ��ɫ 'w'Ϊ��ɫû����
colors  = ['m' 'b' 'r' 'g' 'k' 'c' 'y' 'm' 'b' 'r' 'g' 'k' 'c' 'y' 'b' 'r' 'g' 'k'];% �޸ĵģ�����һ����ɫ��7*2=14
%lines   = [':' '-' {'-.'} {'--'}];%�����������������������޸ģ�����4����������
%lines   = [':' '-' {'-.'} {'--'} ':'];     ����ԭ���ģ�ֻ��5���ߣ�4+1=5��
lines   = [':' '-' {'-.'} {'--'} ':' '-' {'-.'} {'--'} ':' '-' {'-.'} {'--'} ':' '-' {'-.'} {'--'} ':'];      %����޸ĵģ���12���ߣ�4*3=12��
%markers = ['x' '*' 's' 'd' 'v' '^' 'o'];% �ο�http://wiki.ilovematlab.cn/doc-view-54.html
markers = ['x' '*' 's' 'd' 'v' '^' 'o' '+' 'h' 'p' '>' '<' '.' 'x' '*' 's' 'd' 'v'];        %�޸ĵģ�����6��marker����13��

[np,ns] = size(T); 
T(find(T==0))=10^(-17);   %��T�����е�0�滻Ϊ10^(-17)

% Minimal performance per solver

minperf = min(T,[],2); %ÿ���е���С�����һ��

% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p);
end

if (logplot) r = log2(r); end

max_ratio = max(max(r));%ȡ��r��������Ԫ���е������ֵ

% Replace all NaN's with twice the max_ratio and sort.

r(find(isnan(r))) = 2*max_ratio;%������r�е�NaN�ֵΪ���r�������ֵ��2��
r = sort(r);

% Plot stair graphs with markers.

%clf;
for s = 1: ns
 [xs,ys] = stairs(r(:,s),[1:np]/np);
 option = [lines(s) colors(s) markers(s)];
 plot(xs,ys,cell2mat(option),'MarkerSize',3);
 hold on;
end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.
%set(gca,'xlim',[0,3]) ���� x�� �ķ�Χ

axis([ 0 1.1*max_ratio 0 1 ]);

% Legends and title should be added.



