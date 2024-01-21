

function    WriteCoef   ( fname, nrow, ncol, coef )

    funit           =   fopen ( fname, 'wt' );

    fmt             =   [];
    for icol        =   1 : ncol                % col -> column
        fmt         =   [ fmt, ' %16f' ];
    end
    fmt             =   [ fmt, '\n'    ];  
    
%         fprintf ( funit, '#variables = V1 V2 \n');
    for irow        =   1 : nrow                % row -> row
        fprintf ( funit, fmt, coef( irow,: ) );
    end
    
    fclose( funit );
    
end