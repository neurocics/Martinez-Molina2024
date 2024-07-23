function MATr = reduce_cluster(MATo, nele, conn)

if nargin<3
    conn= 26;
end

MATr = MATo;
MAT = ~isnan(MATo);
clusig = bwlabeln(MAT, conn);
n=0;
if max(clusig(:))>0
    for nc = 1:max(clusig(:))
        if sum(clusig(:)==nc)<nele
            MATr(clusig==nc) = NaN;
            n = n+1;
        end

    end
end
n
end