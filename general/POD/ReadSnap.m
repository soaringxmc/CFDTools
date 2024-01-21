function snap               =   ReadSnap ( fname, nnode, nvar )

    funit                   =   fopen( fname, 'rb' );    
    snap( 1:nnode, 1:nvar ) =   0.0;
    snap                    =   fread(funit,[nnode,nvar],'double');
    fclose(funit);
    
end



% 
% function snap   = ReadSnap ( fname, nnode, nvar )
%     
%     funit       =   fopen( fname );
%     for i       =   1 : 4
%         line    =   fgets( funit );
%     end
% 
%     fmt         =   ' ';
%     for ivar    =   1  : nvar
%         fmt     =   strcat( fmt, ' %f');
%     end
%     
%     snap( 1:nnode, 1:nvar )         =   0.0;  
%     for inode                       =   1  : nnode
%         
%         if ( mod(inode,10000) == 0 ) 
%             disp(['inode = ', int2str(inode)]);
%         end
%         
%         line                        =   fgets( funit );
%         [ var,count ]               =   sscanf( line, fmt );
%         if (count == nvar)
%             snap( inode, 1:nvar )   =   var( 1:nvar ); 
%         else
%             disp( 'count ~= nvar' );
%         end
%     end
%     
% end