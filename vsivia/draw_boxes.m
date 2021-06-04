function draw_boxes(S, E, N, varargin)

    switch numel(varargin)
        
        case 2
            
            x = varargin{1} ;
            y = varargin{2} ;
            
            for i = 1:size(x,2)
                
                figure ;
                
                hold on ;
   
                if ~isempty(N)
                draw(N(:,x(1,i)), N(:,y(1,i)), 'b') ;
                else 'No Non-solution Boxes Found'
                end
                if ~isempty(E)
                draw(E(:,x(1,i)), E(:,y(1,i)), 'w') ;
                else 'No Undefined Boxes Found'
                end                
                if ~isempty(S)
                draw(S(:,x(1,i)), S(:,y(1,i)), 'r') ;
                else 'No Solution Boxes Found'
                end

            end
            
        case 3
            
            x = varargin{1} ;
            y = varargin{2} ;
            z = varargin{3} ;
            
            for i = 1:size(x,2)
                
                figure ;
                
                hold on ;
                    
                if ~isempty(S)
                S2 = mid(S) ;
                plot3(S2(:,x(1,i)),S2(:,y(1,i)),S2(:,z(1,i)),'o','MarkerSize',4,'MarkerEdgeColor',[.25 0 0],'MarkerFaceColor',[1 0 0]) ;
                else 'No Solution Boxes Found'
                end
                if ~isempty(E)
                E2 = mid(E) ;
                plot3(E2(:,x(1,i)),E2(:,y(1,i)),E2(:,z(1,i)),'o','MarkerSize',2,'MarkerEdgeColor',[.25 .25 0],'MarkerFaceColor',[1 1 0]) ;
                else 'No Undefined Boxes Found'
                end
                if ~isempty(N)
                N2 = mid(N) ;
                plot3(N2(:,x(1,i)),N2(:,y(1,i)),N2(:,z(1,i)),'o','MarkerSize',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]) ;
                else 'No Non-solution Boxes Found'
                end
            end
            
    end

end
