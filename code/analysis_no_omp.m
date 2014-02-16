clear all; clc; close all;
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

p2 = dataset(1:12,:);
p4 = dataset(13:24,:);
p8 = dataset(25:36,:)

fig_time = figure('name', 'time vs np, no omp');
hold on;
plot(p2(:,2),p2(:,3),'r');
plot(p4(:,2),p4(:,3),'g');
plot(p8(:,2),p8(:,3),'b');


fclose(fid);
