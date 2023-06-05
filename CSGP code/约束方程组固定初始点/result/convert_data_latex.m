clear all
clc;
format long
format compact
%%
data = textread('data.txt');
[row, col] = size(data);
prob=unique(data(:,1));
for j = 1:length(prob)
    fid = fopen(['data_latex',num2str(j),'.txt'],'w');
    for i = 1: row
        tem1 = data(i,:);
        if tem1(1) == prob(j)
            inti_n = ['x',num2str(tem1(3)),'(',num2str(tem1(2)/1000),')'];
            pgp = tem1(4:7);
            zyl = tem1(8:11);
            mscg = tem1(12:col);
        fprintf(fid,'%s & %.2e/%d/%d/%.2e & %.2e/%d/%d/%.2e & %.2e/%d/%d/%.2e %s\n',inti_n,pgp,zyl,mscg, '\\');
        end
    end
    fclose(fid);
end