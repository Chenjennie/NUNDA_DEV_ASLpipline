function [GMdata WMdata CSFdata]=pv_Asllani(data,Pgm,Pwm,Pcsf,mask,ksize,isM0)

pve = cat(4,Pgm,Pwm);
GMdata = zeros(size(data));
WMdata = GMdata;

if isM0==1
    pve=cat(4,pve,Pcsf);
    CSFdata=GMdata;
end

if length(ksize)==1
    nsel=ksize;nzsel=ksize;
elseif length(ksize)==2
    nsel=ksize(1);nzsel=ksize(2);
else
    error('Size must be scalar or have two entries (xy dimension and z dimension)');
end

xsize = size(mask,1);ysize=size(mask,2);zsize=size(mask,3);
count=1;
for i=1:size(mask,1)
    for j=1:size(mask,2)
        for k=1:size(mask,3)
            if mask(i,j,k)>0
                submask=mask(max(i-nsel,1):min(i+nsel,xsize),max(j-nsel,1):min(j+nsel,ysize),max(k-nzsel,1):min(k+nzsel,zsize));
                % calculate the sum of all elements in submask
                if sum(sum(sum(submask)))
                    %if sum3(submask)>5
                    subdata = vols2matrix(data(max(i-nsel,1):min(i+nsel,xsize),max(j-nsel,1):min(j+nsel,ysize),max(k-nzsel,1):min(k+nzsel,zsize),:),submask);
                    subpve = vols2matrix(pve(max(i-nsel,1):min(i+nsel,xsize),max(j-nsel,1):min(j+nsel,ysize),max(k-nzsel,1):min(k+nzsel,zsize),:),submask);
                    pveinv = pinv(subpve);
                    pveprop = sum(subpve)/size(subpve,1);
                    for ti=1:size(data,4)
                        subd = pveinv*subdata(:,ti);
                        GMdata(i,j,k,ti) = subd(1);
                        WMdata(i,j,k,ti) = subd(2);
                        if isM0==1
                            CSFdata(i,j,k,ti)=subd(3);
                        end
                    end
                    %deal with cases where there is very little GM or WM
                    %within the sample volume
                    if (pveprop(1)<0.01), GMdata(i,j,k,:) = 0;end
                    if (pveprop(2)<0.01), WMdata(i,j,k,:) = 0;end
                    if isM0==1
                        if (pveprop(3)<0.01), CSFdata(i,j,k,:) = 0;end
                    end
                    
                end
            end
        end
        count=count+1;
%         if rem(count,100)==0
%             disp(count);
%         end
    end
end
end

