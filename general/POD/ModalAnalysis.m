clear all;
clc;

ipod        =   1;              % ipod = 1 -> use POD tool
idmd        =   0;              % idmd = 1 -> use DMD tool
nSta        =   1;
nEnd        =   384;
nnode       =   12181;
nvar        =   1;
nmode       =   10;              % number of POD modes ( control output)/ DMD Modes
dt          =   0.1;             % required, DMD
ndim        =   nnode*nvar;
nsnap       =   nEnd - nSta + 1;
nRedOrd     =   10;                 
lAddAve     =   0;
RedField    =   [1:1];
nRedField   =   length( RedField );
lRedProbe   =   0;
nProbe      =   0;
probeInd    =   [];


if( lRedProbe == 1 ) 
    fname               =   '../ProbeLocs.dat';
    funit               =   fopen ( fname               );
    line                =   fgets ( funit               );
    nProbe              =   sscanf( line, '%d'          );
    probeInd(1:nProbe)  =   0;   
    for i               =   1 : nProbe
        line            =   fgets ( funit               );
        [ var,count ]   =   sscanf( line, '%d %d %d %f' );
        probeInd(i)     =   var   ( 3                   );
    end
end


%user specified
disp('load...')
load matSnaps.mat
global snaps;
snaps(1:ndim,1:nsnap)       =   matSnaps(1:ndim,nSta:nEnd);
clearvars matSnaps;


%addition operations on variables
snaps = ( snaps - 1.0/1.4 ) / ( 0.5*0.17^2 );


% eigs( 1: nmode ), coef( 1:nmode, 1:nsnap )
% amps( 1: nmode ), grs ( 1:nmode          ), freqs( 1:nmode )
if( ipod == 1 )     
    [ eigs_pod, rt   ,rt_cumu ,coef, probesRed  ]  =   POD( nvar,      nnode,   nmode,              ...
                                                            nRedOrd,   lAddAve, nRedField, RedField,...
                                                            lRedProbe, nProbe,  probeInd                );
end
if( idmd == 1 )     
    [ eigs_dmd, amps ,grs     ,freqs            ]  =   DMD( nvar,      nnode,   nmode,     dt           );
end

if( ipod == 1 ) 
    WriteCoef( '../PODCoef/TimeCoef.dat'     , nsnap, nRedOrd+1, [ [1:nsnap]'    , coef(:,1:nRedOrd)  ] );
    WriteCoef( '../PODCoef/Eigs.dat'         , nsnap, 2        , [ [1:nsnap]'    , eigs_pod           ] );
    WriteCoef( '../PODCoef/CumEigs.dat'      , nsnap, 2        , [ [1:nsnap]'    , rt_cumu            ] );
    
    if( lRedProbe == 1 )
        postFix(1:6)    =   'V1V2V3';
        for ivar        =   1 : nvar
            fname       =   [ '../PODRedProbe/PODRedProbe.', postFix(2*ivar-1:2*ivar) ];
            
            WriteCoef( fname, nsnap, nProbe+1, ...
                        [ [1:nsnap]', probesRed( :, [1:nProbe]+(ivar-1)*nProbe ) ] );
        end
    end 
end

if( idmd == 1 ) 
    WriteCoef( '../DMDCoef/Eigs.dat'         , nmode, 2      , [ real(eigs_dmd), imag(eigs_dmd)] );
    WriteCoef( '../DMDCoef/GrowthAndFreq.dat', nmode, 2      , [ grs           , freqs         ] );
    WriteCoef( '../DMDCoef/Amps.dat'         , nmode, 2      , [ [1:nmode]'    , amps          ] );
end


