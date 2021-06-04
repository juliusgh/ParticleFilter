%%  classe contact_roue_sol
% Domaine défini par l'intersection 
%                     de la boule creuse 
%                                        (Rayon_min)²<||(M-C)||²<Rayon²
%                     du cone defini par 
%                                    (M-C),F(:,i))<0 ; 0<i<Nf+1
%                     de la boite initiale U0 definie dans la classe  
%                                    (M-C),F(:,i))<0 ; 0<i<Nf+1
%
% We set the initial box to [-10 10]x[-10 10]x[-10 10].
%
% Cette classe est une instance de vsivia_parameters dont la notice est 
%                dans <../../html/vsivia_parameters.html vsivia_parameters>

    classdef contact_roue_sol_a < vsivia_parameters

    properties 
       
        algorithm = 'inversion' ;
        % pavé initial de la recherche
        U0 = [-10 10 ; -10 10 ; -10 10] ;
        % Domaine admissible dans l'espace image
        %Rayon_min=7.8; Rayon=8;
        Y0 = interval([7.8^2 8^2 ; -inf 0 ; -inf 0   ; -inf 0    ; -inf 0  ])' ;
        % epaisseur de la frontière
        epsilon = (8-7.8)/2;
        
        C;
        
        F;
        

    end % properties
        
    methods 
        
        
        function obj = contact_roue_sol_a
        
            obj.C = [ 0; 0; 10];

            C1 = [ -1; +1; +0.5];                            
            C2 = [ -1; -1; +1.0];                                 
            C3 = [ +1; -1; -0.5];                                 
            C4 = [ +1; +1; +1.0]; 

            % normales exterieures aux faces
            n1 = C1-obj.C; n1 = n1/norm(n1);         
            n2 = C2-obj.C; n2 = n2/norm(n2); 
            n3 = C3-obj.C; n3 = n3/norm(n3); 
            n4 = C4-obj.C; n4 = n4/norm(n4); 
            
            F1 = cross(n1,(n2-n1));      
            F2 = cross(n2,(n3-n2));      
            F3 = cross(n3,(n4-n3));      
            F4 = cross(n4,(n1-n4)); 
            obj.F = [F1, F2, F3, F4];
        
        end      
        
        
        
        function w = compute(this,x,y,z)
            
            %global C F
            %n = F(:,1);
            
            C_ = this.C;
            F_ = this.F;
            
            w = [(x-C_(1))^2+(y-C_(2))^2+(z-C_(3))^2, ...
                (x-C_(1))*F_(1,1)+(y-C_(2))*F_(2,1)+(z-C_(3))*F_(3,1), ...
                (x-C_(1))*F_(1,2)+(y-C_(2))*F_(2,2)+(z-C_(3))*F_(3,2), ...
                (x-C_(1))*F_(1,3)+(y-C_(2))*F_(2,3)+(z-C_(3))*F_(3,3), ...
                (x-C_(1))*F_(1,4)+(y-C_(2))*F_(2,4)+(z-C_(3))*F_(3,4)] ;         
        end
        
    end % methods
    
    end


