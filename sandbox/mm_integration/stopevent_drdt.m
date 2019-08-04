function [value,isterminal,direction] = stopevent_drdt(t,r,m_s)


value = r(1)-(3*m_s/(4*pi*2165))^(1/3); % is equal to zero when r becomes equal to just the salt
isterminal = 1;
direction  = -1;


end