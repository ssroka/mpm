% remove the paths I added around my code

maxit = 100;
Iter = 1;
gp = path;
start_char = 1;
str_ptn = '/home/ssroka/toEngaging/';
stop_ptn = '/cm/shared';

str_ptn_length = length(str_ptn);
stop_ptn_length = length(stop_ptn);
while ~strcmp(gp(start_char:start_char+stop_ptn_length-1),stop_ptn) &&  ~any(strfind(gp(start_char:start_char+stop_ptn_length-1),stop_ptn(2:end))) && Iter < maxit 
    
    % fprintf('%s\n',gp(start_char:start_char+str_ptn_length))
    
    if strcmp(gp(start_char:start_char+str_ptn_length-1),str_ptn)
        rm_str = sprintf('%s',gp(start_char:start_char+find(gp(start_char:end)==':',1)-2));
        fprintf('Removing: %s\n',rm_str)
        rmpath(rm_str)
        gp = path;
    else
        fprintf('Keeping: %s\n',gp(start_char:start_char+find(gp(start_char:end)==':',1)-1))
        start_char = start_char + find(gp(start_char:end)==':',1) ;
    end
    
    
    Iter = Iter +1;
end





