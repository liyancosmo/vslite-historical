function [] = write_data( filename,years,RW,RWhat,phi,T,P,D )
nrows = 29;
ncols = length(years);
outdata = cell(nrows,ncols+1);
outarray = [years;RW;RWhat;T;P;D;[phi,NaN(1,ncols-1)]];
outheads = cell(nrows,1);
irow = 0;
irow=irow+1; outheads{irow} = 'Year';
irow=irow+1; outheads{irow} = 'TRW';
irow=irow+1; outheads{irow} = 'TRW(predicted)';
for i=1:12; irow=irow+1; outheads{irow} = sprintf('T%02d',i); end
for i=1:12; irow=irow+1; outheads{irow} = sprintf('P%02d',i); end
irow=irow+1; outheads{irow} = 'D';
irow=irow+1; outheads{irow} = 'phi';
for i = 1:nrows
    outdata{i,1} = outheads{i};
    for j = 1:ncols
        val = outarray(i,j);
        if isnan(val); val = ''; end
        outdata{i,j+1} = val;
    end
end
xlswrite(filename,outdata);
end


