
function [ eigs, amps, grs, freqs]  =   DMD ( nvar, nnode, nmode, dt )
    
    global  snaps;
    nsnap                           =   size( snaps    , 2         );
    snaps_m                         =   snaps( :       , 1:nsnap-1 );       % _m -> minus
    snaps_p                         =   snaps( :       , 2:nsnap   );       % _p -> plus                                    
    snaps_ave                       =   mean ( snaps   , 2         );
    
    [ u, s, v ]                     =   svd  ( snaps_m , 'econ'    );
    u                               =   u    ( :       , 1:nmode   );       % unit
    s                               =   s    ( 1:nmode , 1:nmode   );       % fix
    v                               =   v    ( :       , 1:nmode   );       % unit
    Atilde                          =   u'*snaps_p*v*inv(s);
    
    [ vecs, eigs ]                  =   eig  ( Atilde              );
    eigs                            =   diag ( eigs                );
    modes                           =   u*vecs;
    
    for imode                   =   0 : nmode
        
        disp ( [ 'DMDMode    ', int2str(imode) ] )
        fname1                  =   GetFName( imode, 4 );
        fname2                  =   GetFName( imode, 5 );
        if( imode == 0 ) 
            WriteField ( fname1, nvar, nnode, snaps_ave );               % average snap
        else
            WriteField ( fname1, nvar, nnode, real(modes(:,imode)) );    % real part
            WriteField ( fname2, nvar, nnode, imag(modes(:,imode)) );    % imag part
        end
        
    end
    
    % mode amplitudes (amps), growth rates (grs), frequencies (freqs)
    amps                        =   inv(vecs) * u' * snaps_m(:,1);
    grs                         =   real   ( log10(eigs) )   / dt;
    freqs                       =   imag   ( log10(eigs) )   / dt;
    eigs                        =   reshape( eigs , [ nmode, 1 ] );
    amps                        =   reshape( amps , [ nmode, 1 ] );    
    grs                         =   reshape( grs  , [ nmode, 1 ] );
    freqs                       =   reshape( freqs, [ nmode, 1 ] );
    
end

