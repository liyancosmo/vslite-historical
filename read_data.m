function [ years,RW,T,P,D ] = read_data( filename )
[data,txt] = xlsread(filename);
nrows = size(txt,1);
ncols = size(data,2);
years = NaN(1,ncols);
RW = NaN(1,ncols);
T = NaN(12,ncols);
P = NaN(12,ncols);
D = NaN(1,ncols);
for i = 1:nrows
    rowhead = txt{i};
    if strcmpi(rowhead, 'YEAR')
        years(1,:) = data(i,:);
    elseif strcmpi(rowhead, 'TRW')
        RW(1,:) = data(i,:);
    elseif strcmpi(rowhead, 'D')
        D(1,:) = data(i,:);
    elseif strcmpi(rowhead(1), 'T')
        imonth = str2double(rowhead(2:end));
        T(imonth,:) = data(i,:);
    elseif strcmpi(rowhead(1), 'P')
        imonth = str2double(rowhead(2:end));
        P(imonth,:) = data(i,:);
    end 
end
end

