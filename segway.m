% Luka MIlovanov RA105/2022
% Filip Čonić RA126/2022
% Aleksa Balšić RA111/2022
% Doroća Apro RA118/2022

m = 85; 
M = 3.5; 
l = 1.7; 
g = 9.81;
R = 0.05; 
Iw = 0.07; 
Ir = 68.98; 
K = 1;
k1 = (M + m) * R + Iw/R;
k2 = m*l*R;
k3 = m*l;
k4 = Ir + m*l*l;
k5 = m*l*g;
k6 = m*l*R;


Im = Iw + (M+m)*R*R;
I = Ir + m*l*l;

A = [0, 1, 0, 0;
     k5*k1/(k4*k1 - k2*k3), 0, 0, 0;
     0, 0, 0, 1;
    - k2*k5/(k4*k1 - k2*k3), 0, 0, 0];
B = [0;
     (-K*k1 - k3)/(k4*k1 - k2*k3);
     0;
     1/k1 + ((K*k1+k3)/(k4*k1-k2*k3))*(k2/k1)];
C = [1, 0, 0, 0;
     0, 1, 0, 0;
     0, 0, 1, 0;
     0, 0, 0, 1];

D = 0;

mmups = ss(A, B, C, D);
funkcija_prenosa = tf(mmups);

% ***** "ručno" računanje konstanti *****

syms pk1 pk2 pk3 pk4 s

K = [pk1, pk2, pk3, pk4];
A1 = A - B*K;
S = [s, 0, 0, 0;
     0, s, 0, 0;
     0, 0, s, 0;
     0, 0, 0, s];
zeljeno = det(S - A1);
koeficijenti = coeffs(zeljeno, s);
koeficijenti_zeljeni = coeffs((s+1)^4);
pk3 = solve(koeficijenti(1) == koeficijenti_zeljeni(1), pk3);
pk4 = solve(koeficijenti(2) == koeficijenti_zeljeni(2), pk4);
pk1 = solve(koeficijenti(3) == koeficijenti_zeljeni(3), pk1);
pk2 = solve(koeficijenti(4) == koeficijenti_zeljeni(4), pk2);

% ***** računanje konstanti preko acker-a, svi polovi su u p = -1 *****
P = [-1; -1; -1; -1];
preko_funkcije_acker = acker(A, B, P);
probni1 = preko_funkcije_acker(1);
probni2 = preko_funkcije_acker(2);
probni3 = preko_funkcije_acker(3);
probni4 = preko_funkcije_acker(4);
