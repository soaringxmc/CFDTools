
function fname  =   GetFName( isnap, opt )
    
    ctmp        =   int2str ( isnap );
    if    ( opt == 1      ) 
        prefix  =   '../FlowSnap/FlowSnap_';
        %prefix  =   '../../rime_20deg/03_WR/FlowSnap_2560000-2619900/FlowSnap_';
    elseif( opt == 2      )
        prefix  =   '../PODMode/PODMode_';
    elseif( opt == 3      )
        prefix  =   '../PODFieldRed/PODFieldRed_';
    elseif( opt == 4      )
        prefix  =   '../DMDMode/ModeImag/ModeReal_';
    elseif( opt == 5      )
        prefix  =   '../DMDMode/ModeReal/ModeImag_';
    end
    
    if    ( isnap < 10    ) 
         fname  =   [ prefix, '000', ctmp, '.dat' ];
    elseif( isnap < 100   ) 
         fname  =   [ prefix, '00' , ctmp, '.dat' ];
    elseif( isnap < 1000  )  
         fname  =   [ prefix, '0'  , ctmp, '.dat' ];
    elseif( isnap < 10000 )  
         fname  =   [ prefix, ''   , ctmp, '.dat' ];
    end
        
end