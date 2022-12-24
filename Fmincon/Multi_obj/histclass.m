classdef histclass < handle
    properties 
        x= [];
        fval= [];
    end
    methods 
        function add(obj,newx,newfval)
            obj.x= [obj.x; newx];
            obj.fval= [obj.fval  newfval];
        end
    end
end
