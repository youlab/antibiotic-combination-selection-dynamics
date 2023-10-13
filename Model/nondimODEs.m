% `dimODEs` implements the system's ODEs
% works with partial private inhibition mediated by c

function dydt = nondimODEs(t, y, p)
    nr = y(1); % resistant cells
    ns = y(2); % sensitive cells
    b = y(3);  % extracellular Bla
    a = y(4);  % antibiotic 
    s = y(5);  % nutrient
 
    xi=p(1); % nutrient release
    kappab= p(2); % extracellular bla degradation of antibiotic
    phimax = p(3); % resistant cell degradation of antibiotic
    da = p(4); % natural degradation of antibiotic
    ha = p(5); % hill coefficient of antibiotic
    gamma = p(6); % lysis coefficient
    betamin = p(7); % private benefit
    alpha = p(8); % burden
    db = p(9); % natural degradation of extracellular bla
    hi = p(10); % hill coefficient of inhibitor
    inh = p(11); % inhibitor concentration
    c = p(12); % intracellular inhibition

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
