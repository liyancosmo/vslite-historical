function [ years,RW,phi,T,P,D ] = read_data( filename )
[data,txt] = xlsread(filename);
nrows = size(txt,1);
ncols = size(data,2);
years = NaN(1,ncols);
RW = NaN(1,ncols);
T = NaN(12,ncols);
P = NaN(12,ncols);
D = NaN(1,ncols);
phi = NaN;
for i = 1:nrows
    rowhead = txt{i};
    if strcmpi(rowhead, 'Year')
        years(1,:) = data(i,:);
    elseif strcmpi(rowhead, 'TRW')
        RW(1,:) = data(i,:);
    elseif strcmpi(rowhead, 'Phi')
        phi = data(i,:);
        phi = phi(~isnan(phi));
        phi = phi(1);
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

