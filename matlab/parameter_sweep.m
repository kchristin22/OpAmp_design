% last two digits of the AEM
x = 94;

% prerequisites
L = 1e-6;
CL = (2 + 0.01 * x) * 1e-12;
SR = (18 + 0.01 * x) / 1e-6 + 1e-6;
Vdd = (1.8 + 0.003 * x);
Vss = -Vdd;
GB = (7 + 0.01 * x) * 1e+6;

ln = 0.04;
lp = 0.05;
l2 = ln;
l3 = lp;
l6 = lp;
l7 = ln;

Vinmax = 100e-3;
Vinmin = -100e-3;
Kp = 2.9352e-5;
Kn = 9.6379e-5;
VT03max = -0.9056;
VT01min = 0.7860;
VT01max = VT01min;
Cox = 2.47e-3;
Voutmax = Vdd*0.9; % my choice, almost rail-to-rail

K3 = Kp;
K4 = Kp;
K6 = Kp;
K1 = Kn;
K2 = Kn;
K5 = Kn;

% equations that get altered depending the conditions
Cc = 0.22 * CL + 0.1e-12;

I5 = SR * Cc;

S3 = I5 / (K3 * (Vdd - Vinmax - abs(VT03max) + VT01min)^2);

% conditions initialization
cond1 = Cc > 0.22*CL;
cond2 = S3 >= 1;
cond3 = 1;
cond4 = 1;
cond5 = 1;
cond6 = 1;

while((cond1 && cond2 && cond3 && cond4 && cond5 && cond6)==0)
    if(cond2 == 0)
        S3 = 1;
        I5 = S3 * (K3 * (Vdd - Vinmax - abs(VT03max) + VT01min)^2);
        Cc = I5 / SR;
        cond1 = Cc > 0.22*CL;
    end
    if(cond1 == 0)
        Cc = Cc + 0.1e-12;
        I5 = SR * Cc;
        S3 = I5 / (K3 * (Vdd - Vinmax - abs(VT03max) + VT01min)^2);
    end
    S4 = S3;

    I3 = I5 / 2;
    W3 = S3 * L;
    p3 = sqrt(2*K3*S3*I3) / (2*0.667*W3*L*Cox);
    cond3 = p3 > 10*GB;
    if (cond3 == 0)
        p3 = 10 * GB + 1;
        S3 = 2*K3*I3 / (p3 * (2*0.667*L^2*Cox))^2;
        cond2 = S3 >= 1;
        if (cond2 == 0)
            continue;
        end
    end
    
    gm1 = GB * Cc * 2 * pi;
    S2 = gm1^2 / (K2*I5);
    S1 = S2;
    
    b1 = K1 * S1;
    VDS5sat = Vinmin - Vss - sqrt(I5/b1) - VT01max;
    S5 = 2 * I5 / (K5 * VDS5sat^2);
    cond4 = VDS5sat >= 100e-3;
    if (cond4 == 0)
        VDS5sat = 100e-3;
        S5 = 2 * I5 / (K5 * VDS5sat^2);
    end
    
    gm6 = 2.2 * gm1 * CL / Cc;
    I4 = I3;
    gm4 = sqrt(K4*S4*I4);
    S6phase = S4 * gm6 / gm4;
    
    VDS6sat = Vdd - Voutmax;
    S6sat = gm6 / (K6*VDS6sat);
    
%     S6 = max(S6phase,S6sat);
    S6 = S6phase;

    I6 = (gm6^2) / (2*K6*S6);
    
    S7 = I6*S5/I5;
    
    A1 = - 2*gm1 / (I5 * (l2 + l3));
    A2 = - gm6 / (I6 * (l6 + l7));
    A = A1*A2;
    Pdiss = (I5+I6) * (Vdd + abs(Vss));
    cond5 = 20*log10(A) > (20 + 0.01 * x);
    cond6 = Pdiss < (50 + 0.01 * x) * 1e-3;
    if(cond5 == 0 || cond6 == 0)
        Cc = Cc - 0.1e-12;
        I5 = SR * Cc;
        S3 = I5 / (K3 * (Vdd - Vinmax - abs(VT03max) + VT01min)^2);
    end

    cond1 = Cc > 0.22*CL;
    cond2 = S3 >= 1;
    cond3 = p3 > 10*GB;
    cond4 = VDS5sat >= 100e-3;
end

W1 = S1 * L
W2 = S2 * L
W3 = S3 * L
W4 = S4 * L
W5 = S5 * L
W6 = S6 * L
W7 = S7 * L



