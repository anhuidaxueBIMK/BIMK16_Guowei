function write_out_run_more(output,i,j,k,istance_Number,pt,name,dataset,trainfilepath)
%   output = output + 1e-10*randn(length(output),1);  % Break ties at random
fname1 = ['Result/' dataset '/outputWX/' name '.fold' num2str(i) '_' num2str(j) '_' num2str(k) '_' num2str(istance_Number) '_' num2str(pt)];
floder = ['fold' num2str(i)];
if ~exist(['Result/' dataset '/MulitiRANK/' floder],'file')
    mkdir(['Result/' dataset '/MulitiRANK' ],[floder]);
end
  fname2 = ['Result/' dataset '/MulitiRANK/fold' num2str(i) '/' name '.fold' num2str(i) '_' num2str(j) '_' num2str(k)  '_' num2str(istance_Number) '_' num2str(pt)];
  save(fname1,'output','-ascii');
  % Either copy the evaluation script in the current directory or
  % change the line below with the correct path 
%   %跑3.0数据集
  system(['perl Eval-Score-3.0.pl ' trainfilepath  name ...
           '.txt ' fname1 ' ' fname2 '.metric 0']);
%跑4.0数据集
%  system(['perl Eval-Score-4.0.pl ' trainfilepath  name ...
%            '.txt ' fname1 ' ' fname2 '.metric 0']);