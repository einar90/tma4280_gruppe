clear all; close all; clc;
%%%%%%%%%%%%%%%
% With OpenMP %
%%%%%%%%%%%%%%%
filename = 'output';
fid = fopen(filename, 'r');
ln_nr = 0;
data_nr = 0;
dataset = zeros(3*12, 4);
ompcpu = 0;
mpicpu = 0;

while true
  ln = fgetl(fid);
  ln_nr = ln_nr+1;
  if ln == -1
    break;
  end
  elems = strsplit(ln, '\t');
  nums = str2double(elems);

  if mod(ln_nr,13) == 1
    mpicpu = nums(1);
  else
    n = nums(1);
    t = nums(2);
    e = nums(4);
    data_nr = data_nr + 1;
    dataset(data_nr,:) = [mpicpu n t e];
  end
end

p2 = dataset(1:12,:);
p4 = dataset(13:24,:);
p8 = dataset(25:36,:)

fclose(fid);


%%%%%%%%%%%%%
% No OpenMP %
%%%%%%%%%%%%%
filename = 'output_no_omp';
fid = fopen(filename, 'r');
ln_nr = 0;
data_nr = 0;
dataset = zeros(3*12, 4);
ompcpu = 0;
mpicpu = 0;

while true
  ln = fgetl(fid);
  ln_nr = ln_nr+1;
  if ln == -1
    break;
  end
  elems = strsplit(ln, '\t');
  nums = str2double(elems);

  if mod(ln_nr,13) == 1
    mpicpu = nums(1);
  else
    n = nums(1);
    t = nums(2);
    e = nums(4);
    data_nr = data_nr + 1;
    dataset(data_nr,:) = [mpicpu n t e];
  end
end

np2 = dataset(1:12,:);
np4 = dataset(13:24,:);
np8 = dataset(25:36,:)


%%%%%%%%%%%%%%%%%%%%
% Plotting results %
%%%%%%%%%%%%%%%%%%%%
fig_err = figure('name', 'error vs n');
loglog(np8(:,2),np8(:,4),'b');
hold on;
loglog(np2(:,2),np2(:,4),'r');
legend(['P=8'; 'P=2']);
xlabel('n'); ylabel('error');
title('Error vs n, no omp');

fig_time = figure('name', 'time vs n');
hold on;
plot(p2(:,2),p2(:,3),'r');
plot(p4(:,2),p4(:,3),'g');
plot(p8(:,2),p8(:,3),'b');

plot(np2(:,2),np2(:,3),'r--');
plot(np4(:,2),np4(:,3),'g--');
plot(np8(:,2),np8(:,3),'b--');

xlabel('n'); ylabel('t');
title('Time vs. n');
legend(['def omp, np=2';
        'def omp, np=4';
        'def omp, np=8';
        'omp off, np=2';
        'omp off, np=4';
        'omp off, np=8'])
legend('Location', 'EastOutside');


fclose(fid);

