%-------Paramater for Prototype Lowpass Fillter-------%
N = 2;
nfz = 0;
RL = 22;
syms omega s f;

TZs = [inf inf];
%TZs = [1 -1];
%TZs = [1.3217 1.8082 inf inf];
%TZs = [-3.7431 -1.8051 1.5699 6.1910];
%TZs = [0.9218i-0.1546 -0.9218i-0.1546 1.2576 inf inf inf inf];

loadFlag = false;
filename = 'test.mat';

%-------Paramater for Frequancey Mapping-------%
f0 = 2.45*10^9;
omega0 = 2*pi*f0;
BW = 100*10^6;

%-------Design Prototype Lowpass Fillter-------%
%calcurate Fs and Vs
if loadFlag == true
    load(filename)
else
    
    Upre = omega - 1/TZs(1);
    Vpre = sqrt(omega.^2-1)*sqrt(1-1/TZs(1).^2);
    
    for i = 2:N
        Unext = omega*Upre-Upre/TZs(i)+sqrt((omega.^2)-1)*sqrt(1-1/(TZs(i).^2))*Vpre;
        Vnext = omega*Vpre-Vpre/TZs(i)+sqrt((omega.^2)-1)*sqrt(1-1/(TZs(i).^2))*Upre;
        Upre = Unext;
        Vpre = Vnext;
    end
    
    Us = subs(Unext,omega,s/1i);
    Vs = subs(Vnext,omega,s/1i);
    Uroots = AllRoots(Us);
    Fs = sym(1);
    
    for v = Uroots
        Fs = Fs*(s - v);
        %Fs = Fs*(s - imag(v)*1i);
    end
    
    Fs = expand(Fs);
    Fomega = subs(Fs,s,1i*omega);
    %vpa(Fs)
    
    %calcurate Ps
    Ps = sym(1);
    
    for i = 1:nfz
        Ps = Ps*(s - TZs(i)*1i);
    end
    
    Ps = expand(Ps);
    Pomega = subs(Ps,s,1i*omega);
    %vpa(Ps);
    
    if rem(N - nfz,2) == 0
        Ps = Ps*1i;
    end
    
    %calcurate Es
    e = (1/sqrt(10.^(RL/10)-1))*abs(subs(Ps,s,1i)/subs(Fs,s,1i));
    
    if N == nfz
        er = e/sqrt(e.^2-1);
    else
        er = 1;
    end
    
    Ecal = Pomega/e - 1i*Fomega/er;
    Eroots = AllRoots(Ecal);
    Eomega = sym(1);
    
    for v = Eroots
        if vpa(imag(v)) < 0
            v = conj(v);
        end
        Eomega = Eomega*(omega-v);
    end
    
    Es = subs(Eomega,omega,s/1i);
    EsC = coeffs(Es,'All');
    Es = Es/EsC(1);
    Es = expand(Es);
    %-----------------------------------------------------------%
    %S11 = vpa(20*log10(abs(Fomega / (Eomega*er))));
    %S21 = vpa(20*log10(abs(Pomega / (Eomega*e))));
    %fplot([S11 S21])
    
    
    %-------Frequancey Mapping-------%
    %omegaMap = 1/((f0/BW)*(omega/omega0-omega0/omega));
    omegaMap = (f0/BW)*(omega/omega0-omega0/omega);
    %omegaMap = omega/omega0;
    EomegaMap = subs(Eomega,omega,omegaMap);
    FomegaMap = subs(Fomega,omega,omegaMap);
    PomegaMap = subs(Pomega,omega,omegaMap);
    %S11 = vpa(20*log10(abs(FomegaMap / (EomegaMap*er))))
    %S21 = vpa(20*log10(abs(PomegaMap / (EomegaMap*e))))
    %fplot([S11 S21],[2*pi*2*10^9 2*pi*3*10^9])
    
    
    %-------------F Paramater-------------%
    Fcoeffs = coeffs(Fs,'All');
    Ecoeffs = coeffs(Es,'All');
    
    ns = sym(0);
    ms = sym(0);
    
    for i = 0:N
        if rem(i,2) == 0
            ms = ms + real(Ecoeffs(i+1)+Fcoeffs(i+1))*s^(N-i);
            ns = ns + 1i*imag(Ecoeffs(i+1)+Fcoeffs(i+1))*s^(N-i);
        else
            ms = ms + 1i*imag(Ecoeffs(i+1)+Fcoeffs(i+1))*s^(N-i);
            ns = ns + real(Ecoeffs(i+1)+Fcoeffs(i+1))*s^(N-i);
        end
    end
    
    if rem(N,2) == 0 %even
        yd = ms;
        y21n = Ps/e;
        y22n = ns;
    else %odd
        yd = ns;
        y21n = Ps/e;
        y22n = ms;
    end
    
    yd = expand(yd);
    ydC = coeffs(yd,'All');
    yd = yd/ydC(1);
    y21n = y21n/ydC(1);
    y22n = y22n/ydC(1);
    y21 = y21n/yd;
    y22 = y22n/yd;
    
    Msl = limit(y21,s,sym(inf));
    y21n_exM = y21n - Msl*yd;
    
    lambda_k = AllRoots(yd)/1i;
    
    yd_diff = diff(yd);
    
    r21k = zeros(1,N);
    r22k = zeros(1,N);
    T_NK = zeros(1,N);
    T_1K = zeros(1,N);
    
    for i = 1:N
        r21k(1,i) = subs(y21n/yd_diff,s,1i*lambda_k(i));
        r22k(1,i) = subs(y22n/yd_diff,s,1i*lambda_k(i));
        T_NK(1,i) = sqrt(r22k(1,i));
        T_1K(1,i) = r21k(1,i)/sqrt(r22k(1,i));
    end
    
    M = zeros(N+2);
    
    for i = 1:N+2
        if 1<i&&i<N+2
            M(i,i) = -1*lambda_k(i-1);
            M(i,1) = T_1K(1,i-1);
            M(1,i) = T_1K(1,i-1);
            M(N+2,i) = T_NK(1,i-1);
            M(i,N+2) = T_NK(1,i-1);
        end
    end
    
    M(1,N+2) = Msl/1i;
    M(N+2,1) = Msl/1i;
    %save(filename,'M')
    %M;
    
    %-------Frequancey Mapping-------%
    f0i = zeros(1,N);
    Qi = zeros(1,N);
    for i = 1:N+2
        if 1<i&&i<N+2
            M(i,1) = sqrt((BW/f0)*M(i,1)^2);
            M(1,i) = sqrt((BW/f0)*M(1,i)^2);
            M(N+2,i) = sqrt((BW/f0)*M(N+2,i)^2);
            M(i,N+2) = sqrt((BW/f0)*M(i,N+2)^2);
        end
    end
    M
    
    for i=2:N+1
        f0i(1,i-1) = (-BW*M(i,i)+sqrt(BW^2*M(i,i)^2+4*f0^2))/2;
        Qi(1,i-1) = 1/M(i,1)^2;
    end
    f0i
    Qi
end

%-------Function-------%
function fRoots = AllRoots(f)
f = expand(f);
fC = coeffs(f,'All');
%fRoots = roots(fC).';
fRoots = vpa(roots(fC)).';
end





