



vars = who;

for ii = 1:length(vars)
    eval(sprintf('ic.%s = %s;',vars{ii},vars{ii}))
end


