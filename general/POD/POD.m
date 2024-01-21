function [ eigs, rt, rt_cumu, coef, probesRed ]=   POD ( nvar,    nnode,   nmode,               ...
                                                         nRedOrd, lAddAve, nRedField, RedField, ... 
                                                         lProbe,  nProbe,  probeInd                 )
                                          
    global snaps;
    ndim                            =   size( snaps, 1   );
    nsnap                           =   size( snaps, 2   );
    snaps_ave                       =   mean( snaps, 2   );
    
    snaps_fluc(1:ndim,1:nsnap)      =   0.0;
    for isnap                       =   1 : nsnap
        snaps_fluc( :, isnap )      =   snaps(:,isnap) - snaps_ave(:);
    end
    clearvars snaps;

    disp('svd...');
    [ u, s, v ]                     =   svd    ( snaps_fluc, 'econ'  );
    modes                           =   u      ( :         , 1:nsnap );
    eigs                            =   diag   ( s                   );
    eigs                            =   reshape( eigs, [ nsnap, 1 ]  );
    clearvars u s v;
    
    disp('coef...');
    coef(1:nsnap,1:nRedOrd)         =   0;
    for isnap                       =   1 : nsnap
    for imode                       =   1 : nRedOrd
       
        coef(isnap,imode)           =   dot( snaps_fluc(:,isnap), modes(:,imode) );
    
    end
    end
    
    clearvars snaps_fluc;
    disp('snaps_red...');
        snaps_red(1:ndim,1:nsnap)   =   0.0;
    for isnap                       =   1 : nsnap
        for imode                   =   1 : nRedOrd
            snaps_red(:,isnap)      =   snaps_red(:,isnap) + ...
                                        coef(isnap,imode) * modes(:,imode);
        end
        if( lAddAve == 1 )
            snaps_red(:,isnap)      =   snaps_red(:,isnap) +  snaps_ave(:);
        end
    end
    
    
    rt                              =   eigs/sum(eigs);                    % rt  -> ratio
    rt_cumu( 1:nsnap )              =   0.0;
    for ieig    =   1 : nsnap
        rt_cumu(ieig)               =   sum( eigs(1:ieig) ) / sum( eigs ); % cum -> cumulate
    end
    rt_cumu                         =   reshape( rt_cumu, [ nsnap, 1 ]  );
    
    probesRed( 1:nsnap, 1:nProbe*nvar )   =   0.0;
    if( lProbe == 1 )
        for ivar                    =   1 : nvar 
        for iProbe                  =   1 : nProbe
            probesRed( 1:nsnap, iProbe + (ivar-1)*nProbe ) ...
        =   snaps_red ( probeInd(iProbe)+(ivar-1)*nnode, 1:nsnap );
        end
        end
    end
    
    disp('WritePODModes...');
    for imode                       =   0 : nmode
        
        disp ( [ 'PODMode      ', int2str(imode) ] )
        fname                       =   GetFName( imode, 2 );
        
        if( imode == 0 ) 
            WriteField ( fname, nvar, nnode, snaps_ave      );
        else
            WriteField ( fname, nvar, nnode, modes(:,imode) );
        end
        
    end
    clearvars modes;
    
    disp('WriteFieldRed...');
    for iRedField                   =   1 : nRedField
        isnap                       =   RedField(iRedField);     
        disp ( [ 'PODFieldRed  ', int2str( isnap ) ] )
        fname                       =   GetFName  ( isnap, 3   );
        WriteField ( fname, nvar, nnode, snaps_red(:,isnap)    );
    end
    clearvars snaps_red;
end


