
clearvars -except matSnaps;
clc;

nSta    =   1;              % start snap number
nEnd    =   384;            % end snap number
nnode   =   12181;
nvar    =   1;
ndim    =   nnode*nvar;
nSnap   =   nEnd - nSta + 1;

snaps(1:ndim,1:nSnap)           =   0.0;

for isnap                       =   nSta : nEnd
    
    disp ( ['ReadFlowSnap    '  ,   int2str(isnap)] )
    fname                       =   GetFName ( isnap, 1           );
    snap                        =   ReadSnap ( fname, nnode, nvar );
        k                       =   0;
    for ivar                    =   1:nvar
    for inode                   =   1:nnode

        k                       =   k + 1;
        snaps( k, isnap-nSta+1 )=   snap     ( inode, ivar        );

    end
    end

end

%user specified
%load matSnaps.mat      % MATLAB restart
matSnaps(:,nSta:nEnd)  =   snaps(:,1:nSnap);

save ('matSnaps.mat', 'matSnaps', '-v7.3');

disp(['maximum value = ', num2str(max(max(matSnaps)))])
disp(['mimimum value = ', num2str(min(min(matSnaps)))])

