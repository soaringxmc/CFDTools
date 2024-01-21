function WriteField ( fname, nvar, nnode, mode )
    
    funit           =   fopen( fname, 'wb'      );
    ofs             =   [];
    for ivar        =   1 : nvar
        ofs         =   [ ofs, nnode*(ivar-1)   ];  % ofs -> offset
    end
    
    modeOut         =   reshape( mode, [nnode,nvar] );
    fwrite ( funit, modeOut(:,:), 'double' );
    
    fclose( funit );
    
end


% function WriteField ( fname, nvar, nnode, mode )
%     
%     funit           =   fopen( fname, 'wt'      );
%     fmt             =   ' ';
%     ofs             =   []; 
%     for ivar        =   1  : nvar
%         fmt         =   strcat( fmt , ' %14f'   );
%         ofs         =   [ ofs, nnode*(ivar-1)   ];    % ofs -> offset
%     end
%     fmt             =   [ fmt        , '\n'     ];
%     
%     for inode       =   1 : nnode   
%         fprintf( funit, fmt, mode( inode + ofs ) );
%     end
%     fclose( funit );
%     
% end