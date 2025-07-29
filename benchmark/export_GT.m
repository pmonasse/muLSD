folders = dir('P*');
for i = 1:length(folders)
  n = folders(i).name;
  load([n,'/',n,'LinesAndVP.mat']);
  r = lines';
  r = reshape(r(:),4,[])';
  dlmwrite([n,'_gt.txt'],r,' ');
end
