files=dir('P*_GND.mat');
for i = 1:length(files)
    n = files(i).name;
    load(n);
    dlmwrite([n(1:8),'_gt.txt'],line_gnd,' ');
end
