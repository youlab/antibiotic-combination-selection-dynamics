% `dimODEs` implements the system's ODEs
% works with partial private inhibition mediated by c

function dydt = nondimODEsvariable(t, y, p)
    nr = y(1);
    ns = y(2);        
    b = y(3);
    a = y(4);
    s = y(5); 
 
    xi=p(1); % nutrient release
    kappab= p(2); % extracellular bla degradation of antibiotic
    phimax = p(3); % cell degradation of antibiotic
    da = p(4); % natural degradation of antibiotic
    ha = p(5); % hill coefficient of antibiotic
    gamma = p(6); % lysis coefficient
    db = p(7); % natural degradation of extracellular bla
    hi = p(8); % hill coefficient of inhibitor
    alpha = p(9);
    betamin = p(10); % private good
    c = p(11);
    inh = p(12); % inhibitor concentration

    beta = betamin + c*((1 - betamin)*inh^hi)/(1 + inh^hi);
    phi = phimax*(1 - c*(inh^hi/(1 + inh^hi)));
    
    g = s/(1+s);
    l = gamma * a^ha/(1+a^ha) * g;

    dnrdt = (alpha*g-beta*l)*nr;
    dnsdt = (g-l)*ns;
    dbdt = beta*l*nr - ((db * inh^hi) / (1 + inh^hi))*b;
    dadt =-kappab*b*a-phi*nr*a - da*a;
    dsdt = (xi*l - g) * ns + (xi*beta*l - alpha*g) * nr;        

    dydt = [dnrdt; dnsdt; dbdt; dadt; dsdt];
end
