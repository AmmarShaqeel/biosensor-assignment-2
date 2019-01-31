function [t, p] = model(p0, t_seg, I)

    C = 1 ;                 %% unit: uF/cm^2

    v0 = -60.045;           %% unit: mV
    vNa = 55.17;            %% unit: mV
    vK = -72.14 ;           %% unit: mV
    vL = -49.42 ;           %% unit: mV


    gNa = 120 ;             %% unit: mS/cm^2
    gK = 36 ;               %% unit: mS/cm^2
    gL = 0.3 ;              %% unit: mS/cm^2

    idx_seg = 0;
    func = @(t,p) [ (1/C)*(I(idx_seg)-(gK*p(2).^4.*(p(1) - vK)) - (gNa*p(3).^3.*p(4).*(p(1)-vNa))-(gL.*(p(1)-vL)));
        ( (an(p(1)).*(1-p(2))) - (bn(p(1)).*p(2)) );
        ( (am(p(1)).*(1-p(3))) - (bm(p(1)).*p(3)) ); 
        ( (ah(p(1)).*(1-p(4))) - (bh(p(1)).*p(4)) ) ] ;
    
    t_start = 0;
    t = t_start;
    p = p0;
    while idx_seg < length(t_seg)
        idx_seg = idx_seg + 1;
        t_end = t_seg(idx_seg);
        [t_sol, p_sol] = ode15s(func, [t_start, t_end], p(end, :));
        t = [t; t_sol(2 : end)];
        p = [p; p_sol(2 : end, :)];
        t_start = t_end;
    end
end