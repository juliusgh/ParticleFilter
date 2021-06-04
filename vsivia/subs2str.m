function s = subs2str(subs)

 s = cellfun(@process, subs, 'UniformOutput', false) ;
 
 s = [ '{' cell2mat(s) '}' ] ;

end


function s = process(s)
    if isa(s, 'numeric')
        s = [ ' [' num2str(s) '] ' ] ;
    else
        s = [ ' ''' s ''' '] ;
    end
end