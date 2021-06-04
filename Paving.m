classdef Paving
    %PAVING Paving of many HyperRects
    %   Not in use
    
    properties
        Kin % contained in the set
        Ki % intersected with the set
    end
    
    methods
        function obj = Paving(Kin,Ki)
            %PAVING Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0
                obj.Kin = Kin;
                obj.Ki = Ki;
            else
                obj.Kin = HyperRect.empty();
                obj.Ki = HyperRect.empty();
            end
        end
        
        function plot(obj, color)
            if nargin < 2
                color = 'g';
            end
            arrayfun(@(x) x.plot(color), obj.Kin);
            arrayfun(@(x) x.plot(color), obj.Ki);
        end
        
        function p = polygon(obj)
            p = polyshape();
            for i = 1:numel(obj.Ki)
                pi = obj.Ki(i).polygon();
                p = union(p,pi);
            end
            for i = 1:numel(obj.Kin)
                pin = obj.Kin(i).polygon();
                p = union(p,pin);
            end
        end
        
        function p = polygonIn(obj)
            p = polyshape();
            for i = 1:numel(obj.Kin)
                pin = obj.Kin(i).polygon();
                p = union(p,pin);
            end
        end
        
        function c = convhull(obj)
            c = convhull(obj.polygon());
        end
        
        function obj = apply(obj,f)
            obj.Kin = arrayfun(f, obj.Kin);
            obj.Ki = arrayfun(f, obj.Ki);
        end
        
        function p = intersect(obj,p2)
            % very inefficient, better: intersection of polyshapes
            p = Paving();
            for i = 1:numel(obj.Ki)
               for j = 1:numel(p2.Ki)
                   if obj.Ki(i).intersects(p2.Ki(j))
                       if obj.Ki(i).width > p2.Ki(j).width
                            p.Ki(end+1) = p2.Ki(j);
                       else
                            p.Ki(end+1) = obj.Ki(i);
                       end
                   end
               end
               for j = 1:numel(p2.Kin)
                   if obj.Ki(i).intersects(p2.Kin(j))
                       if obj.Ki(i).width > p2.Kin(j).width
                            p.Ki(end+1) = p2.Kin(j);
                       else
                            p.Ki(end+1) = obj.Ki(i);
                       end
                   end
               end
            end
            for i = 1:numel(obj.Kin)
               for j = 1:numel(p2.Ki)
                   if obj.Kin(i).intersects(p2.Ki(j))
                       if obj.Kin(i).width > p2.Ki(j).width
                            p.Ki(end+1) = p2.Ki(j);
                       else
                            p.Ki(end+1) = obj.Kin(i);
                       end
                   end
               end
               for j = 1:numel(p2.Kin)
                   if obj.Kin(i).intersects(p2.Kin(j))
                       if obj.Kin(i).width > p2.Kin(j).width
                            p.Kin(end+1) = p2.Kin(j);
                       else
                            p.Kin(end+1) = obj.Kin(i);
                       end
                   end
               end
            end
        end
    end
end

