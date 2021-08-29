a = 10;
b = 14.5;
c = 50;
d = 60;
e = 100;
f = 620;
g = 710;
l = 720;

Rpol = 26.145;
Ft = 295.48;
Sut = 470; %MPa
Rx1 = 299.68;
Rx2 = 223.74;
Ry1 = 73.19;
Ry2 =4.43;
Fsx = 590.96;
Fsy = 77.62;
Fat = 515.025;
dab = 15;   %12,7
dbd = 17;   %15
dde = 23;
def = 30; %15
dfl = 15;
z = (0:0.1:l);  % Para gráfico
pos = (0:0.1:720);
E = 205000;
I1 = (pi*dab^4)/64;
I2 = (pi*dbd^4)/64;
I3 = (pi*dde^4)/64;
I4 = (pi*def^4)/64;
I5 = (pi*dfl^4)/64;
%Constantes para eixo x
B1 = (Rx1/2)*(b-a)^2;
B2 = ((Rx1/2)*(d-a)^2 - (Fsx/2)*(d-c)^2);
B3 = ((Rx1/2)*(e-a)^2 - (Fsx/2)*(e-c)^2);
B4 = (((Rx1/2)*(f-a)^2) - ((Fsx/2)*(f-c)^2) + (Fat/(520*6))*(f-e)^3);

A2 = ((Rx1/(6))*(b-a)^3);
A3 = (((Rx1/(6))*(d-a)^3) - ((Fsx/6)*(d-c)^3));
A4 = (((Rx1/(6))*(e-a)^3) - ((Fsx/6)*(e-c)^3));
A5 = (((Rx1/(6))*(f-a)^3) - ((Fsx/6)*(f-c)^3) + ((Fat/(520*24))*(f-e)^4));
A6 = (((Rx1/(6))*(g-a)^3) - ((Fsx/6)*(g-c)^3) + ((Fat/(520*24))*(g-e)^4) - ((Fat/(520*24))*(g-f)^4));

MAAx = [I2 -I1 0 0 0 0 0 0;I2*(b-a) -I1*b 0 0 0 -I1 0 0;0 I3 -I2 0 0 0 0 0;0 I3*d -I2*d 0 0 I3 -I2 0;0 0 I4 -I3 0 0 0 0;0 0 I4*e -I3*e 0 0 I4 -I3;0 0 0 I5 -I4 0 0 0;0 0 0 I5*f -I4*(f-g) 0 0 I5];
MRRx = [B1*(I1-I2);A2*(I1-I2);B2*(I2-I3);A3*(I2-I3);B3*(I3-I4);A4*(I3-I4);B4*(I4-I5);(A5*(I4-I5))-(A6*I4)];
MX = MAAx\MRRx;
c1=MX(1,1);
c2=MX(2,1);
c3=MX(3,1);
c4=MX(4,1);
c5=MX(5,1);
k2=MX(6,1);
k3=MX(7,1);
k4=MX(8,1);
k1=-c1*a;
k5=(-c5*g)-A6;

%Constantes para eixo y
E1=(-Ry1/2)*(b-a)^2;
E2=(-Ry1/2)*(d-a)^2 + (Fsy/2)*(d-c)^2;
E3=(-Ry1/2)*(e-a)^2 + (Fsy/2)*(e-c)^2;
E4=(-Ry1/2)*(f-a)^2 + (Fsy/2)*(f-c)^2;

D2=((-Ry1/6)*(b-a)^3);
D3=((-Ry1/6)*(d-a)^3) + (Fsy/6)*(d-c)^3;
D4=((-Ry1/6)*(e-a)^3) + (Fsy/6)*(e-c)^3;
D5=((-Ry1/6)*(f-a)^3) + (Fsy/6)*(f-c)^3;
D6=((-Ry1/6)*(g-a)^3) + (Fsy/6)*(g-c)^3;

MAAy = [I2 -I1 0 0 0 0 0 0;I2*(b-a) -I1*b 0 0 0 -I1 0 0;0 I3 -I2 0 0 0 0 0;0 I3*d -I2*d 0 0 I3 -I2 0;0 0 I4 -I3 0 0 0 0;0 0 I4*e -I3*e 0 0 I4 -I3;0 0 0 I5 -I4 0 0 0;0 0 0 I5*f -I4*(f-g) 0 0 I5];
MRRy = [E1*(I1-I2);D2*(I1-I2);E2*(I2-I3);D3*(I2-I3);E3*(I3-I4);D4*(I3-I4);E4*(I4-I5);(D5*(I4-I5))-(D6*I4)];
MY = MAAx\MRRy;

f1=MY(1,1);
f2=MY(2,1);
f3=MY(3,1);
f4=MY(4,1);
f5=MY(5,1);
u2=MY(6,1);
u3=MY(7,1);
u4=MY(8,1);
u1=-f1*a;
u5=(-f5*g)-D6;
n = 7201;
tetax = zeros (1,n);
tetay = zeros (1,n);
yx = zeros (1, n);
yy = zeros (1, n);
Mx = zeros (1, n);
My = zeros (1, n);
Vx = zeros (1, n);
Vy = zeros (1, n);
T = zeros (1, n); 
for k=1:1:7201
    if (z(k) < b)
        tetax(k) = (1/(E*I1))*((Rx1/2)*(z(k)>=a).*(z(k)-a).^2 - (Fsx/2)*(z(k)>=c).*(z(k)-c).^2 + (Fat/(520*2*3))*(z(k)>=e).*(z(k)-e).^3 - (Fat/(520*2*3))*(z(k)>=f).*(z(k)-f).^3 - (Rx2/2)*(z(k)>=g).*(z(k)-g).^2 + c1);
        tetay(k) = (1/(E*I1))*((-Ry1/2)*(z(k)>=a).*(z(k)-a).^2 + (Fsy/2)*(z(k)>=c).*(z(k)-c).^2 - (Ry2/2)*(z(k)>=g).*(z(k)-g).^2 + f1 );
        yx(k) = (1/(E*I1))*(((Rx1/(6))*(z(k)>=a).*(z(k)-a).^3) - ((Fsx/6)*(z(k)>=c).*(z(k)-c).^3) + ((Fat/(520*24))*(z(k)>=e).*(z(k)-e).^4) - ((Fat/(520*24))*(z(k)>=f).*(z(k)-f).^4)-((Rx2/6)*(z(k)>=g).*(z(k)-g).^2) + c1.*z(k) + k1);
        yy(k) = (1/(E*I1))*(((-Ry1/6)*(z(k)>=a).*(z(k)-a).^3) + ((Fsy/6)*(z(k)>=c).*(z(k)-c).^3) - ((Ry2/6)*(z(k)>=g).*(z(k)-g).^3) + (f1.*z(k)) + u1);
 
    elseif (z(k) < d)
        tetax(k) = (1/(E*I2))*((Rx1/2)*(z(k)>=a).*(z(k)-a).^2 - (Fsx/2)*(z(k)>=c).*(z(k)-c).^2 + (Fat/(520*2*3))*(z(k)>=e).*(z(k)-e).^3 - (Fat/(520*2*3))*(z(k)>=f).*(z(k)-f).^3 - (Rx2/2)*(z(k)>=g).*(z(k)-g).^2 + c2);
        tetay(k) = (1/(E*I2))*((-Ry1/2)*(z(k)>=a).*(z(k)-a).^2 + (Fsy/2)*(z(k)>=c).*(z(k)-c).^2 - (Ry2/2)*(z(k)>=g).*(z(k)-g).^2 + f2);
        yx(k) = (1/(E*I2))*(((Rx1/(6))*(z(k)>=a).*(z(k)-a).^3) - ((Fsx/6)*(z(k)>=c).*(z(k)-c).^3) + ((Fat/(520*24))*(z(k)>=e).*(z(k)-e).^4) - ((Fat/(520*24))*(z(k)>=f).*(z(k)-f).^4)-((Rx2/6)*(z(k)>=g).*(z(k)-g).^2) + c2.*z(k) + k2);
        yy(k) = (1/(E*I2))*(((-Ry1/6)*(z(k)>=a).*(z(k)-a).^3) + ((Fsy/6)*(z(k)>=c).*(z(k)-c).^3) - ((Ry2/6)*(z(k)>=g).*(z(k)-g).^3) + (f2.*z(k)) + u2);
 
    elseif (z(k) < e)
        tetax(k) = (1/(E*I3))*((Rx1/2)*(z(k)>=a).*(z(k)-a).^2 - (Fsx/2)*(z(k)>=c).*(z(k)-c).^2 + (Fat/(520*2*3))*(z(k)>=e).*(z(k)-e).^3 - (Fat/(520*2*3))*(z(k)>=f).*(z(k)-f).^3 - (Rx2/2)*(z(k)>=g).*(z(k)-g).^2 + c3);
        tetay(k) = (1/(E*I3))*((-Ry1/2)*(z(k)>=a).*(z(k)-a).^2 + (Fsy/2)*(z(k)>=c).*(z(k)-c).^2 - (Ry2/2)*(z(k)>=g).*(z(k)-g).^2 + f3);
        yx(k) = (1/(E*I3))*(((Rx1/(6))*(z(k)>=a).*(z(k)-a).^3) - ((Fsx/6)*(z(k)>=c).*(z(k)-c).^3) + ((Fat/(520*24))*(z(k)>=e).*(z(k)-e).^4) - ((Fat/(520*24))*(z(k)>=f).*(z(k)-f).^4)-((Rx2/6)*(z(k)>=g).*(z(k)-g).^2) + c3.*z(k) + k3);
        yy(k) = (1/(E*I3))*(((-Ry1/6)*(z(k)>=a).*(z(k)-a).^3) + ((Fsy/6)*(z(k)>=c).*(z(k)-c).^3) - ((Ry2/6)*(z(k)>=g).*(z(k)-g).^3) + (f3.*z(k)) + u3);
 
    elseif (z(k) < f)
        tetax(k) = (1/(E*I4))*((Rx1/2)*(z(k)>=a).*(z(k)-a).^2 - (Fsx/2)*(z(k)>=c).*(z(k)-c).^2 + (Fat/(520*2*3))*(z(k)>=e).*(z(k)-e).^3 - (Fat/(520*2*3))*(z(k)>=f).*(z(k)-f).^3 - (Rx2/2)*(z(k)>=g).*(z(k)-g).^2 + c4);
        tetay(k) = (1/(E*I4))*((-Ry1/2)*(z(k)>=a).*(z(k)-a).^2 + (Fsy/2)*(z(k)>=c).*(z(k)-c).^2 - (Ry2/2)*(z(k)>=g).*(z(k)-g).^2 + f4);
        yx(k) = (1/(E*I4))*(((Rx1/(6))*(z(k)>=a).*(z(k)-a).^3) - ((Fsx/6)*(z(k)>=c).*(z(k)-c).^3) + ((Fat/(520*24))*(z(k)>=e).*(z(k)-e).^4) - ((Fat/(520*24))*(z(k)>=f).*(z(k)-f).^4)-((Rx2/6)*(z(k)>=g).*(z(k)-g).^2) + c4.*z(k) + k4);
        yy(k) = (1/(E*I4))*(((-Ry1/6)*(z(k)>=a).*(z(k)-a).^3) + ((Fsy/6)*(z(k)>=c).*(z(k)-c).^3) - ((Ry2/6)*(z(k)>=g).*(z(k)-g).^3) + (f4.*z(k)) + u4);

    
    elseif (z(k) <= l)
        tetax(k) = (1/(E*I5))*((Rx1/2)*(z(k)>=a).*(z(k)-a).^2 - (Fsx/2)*(z(k)>=c).*(z(k)-c).^2 + (Fat/(520*2*3))*(z(k)>=e).*(z(k)-e).^3 - (Fat/(520*2*3))*(z(k)>=f).*(z(k)-f).^3 - ((Rx2/2)*(z(k)>=g).*(z(k)-g).^2) + c5);
        tetay(k) = (1/(E*I5))*((-Ry1/2)*(z(k)>=a).*(z(k)-a).^2 + (Fsy/2)*(z(k)>=c).*(z(k)-c).^2 - (Ry2/2)*(z(k)>=g).*(z(k)-g).^2 + f5);
        yx(k) = (1/(E*I5))*(((Rx1/(6))*(z(k)>=a).*(z(k)-a).^3) - ((Fsx/6)*(z(k)>=c).*(z(k)-c).^3) + ((Fat/(520*24))*(z(k)>=e).*(z(k)-e).^4) - ((Fat/(520*24))*(z(k)>=f).*(z(k)-f).^4)-((Rx2/6)*(z(k)>=g).*(z(k)-g).^2) + c5.*z(k) + k5);
        yy(k) = (1/(E*I5))*(((-Ry1/6)*(z(k)>=a).*(z(k)-a).^3) + ((Fsy/6)*(z(k)>=c).*(z(k)-c).^3) - ((Ry2/6)*(z(k)>=g).*(z(k)-g).^3) + (f5.*z(k)) + u5);

     else
        disp("ALGO DE ERRADO NÃO ESTA CERTO");
    end
    
    Mx(k) = (Rx1)*(z(k)>=a).*(z(k)-a).^1 - (Fsx)*(z(k)>=c).*(z(k)-c).^1 + (Fat/(520*2))*(z(k)>=e).*(z(k)-e).^2 - (Fat/(520*2))*(z(k)>=f).*(z(k)-f).^2 - ((Rx2)*(z(k)>=g).*(z(k)-g).^1);
    My(k) = (-Ry1)*(z(k)>=a).*(z(k)-a).^1 + (Fsy)*(z(k)>=c).*(z(k)-c).^1 - (Ry2)*(z(k)>=g).*(z(k)-g).^1;
    Vx(k) = (Rx1)*(z(k)>=a).*(z(k)-a).^0 - (Fsx)*(z(k)>=c).*(z(k)-c).^0 + (Fat/(520))*(z(k)>=e).*(z(k)-e).^1 - (Fat/(520))*(z(k)>=f).*(z(k)-f).^1 - ((Rx2)*(z(k)>=g).*(z(k)-g).^0);
    Vy(k) = (-Ry1)*(z(k)>=a).*(z(k)-a).^0 + (Fsy)*(z(k)>=c).*(z(k)-c).^0 - (Ry2)*(z(k)>=g).*(z(k)-g).^0;
    T (k) = Ft*Rpol*(z(k) >= c).*(z(k)-c).^0 - (Fat/520)*(def/2)*(z(k)>=e).*(z(k)-e).^1 + (Fat/520)*(def/2)*(z(k)>=f).*(z(k)-f).^1;
end

% tetax = tetax1 + tetax2 + tetax3 + tetax4+ tetax5;
% tetay = tetay1 + tetay2 + tetay3 + tetay4+ tetay5;
% yy = yy_1 + yy_2 + yy_3 + yy_4+ yy_5;
% yx = yx_1 + yx_2 + yx_3 + yx_4+ yx_5;
tetar = (tetax.^2 + tetay.^2).^0.5;
Mr = (Mx.^2 + My.^2).^0.5;
Vr = (Vx.^2 + Vy.^2).^0.5;
yr = (yx.^2 + yy.^2).^0.5;
% figure ("Name", "Inclinação")
% subplot (2,2,1);
% xlabel('mm');
% ylabel('rad');
% plot(z,tetax);
% title('Inclinação em x');
% 
% subplot (2,2,2);
% plot(z,tetay);
% xlabel('mm');
% ylabel('rad');
% title('Inclinação em Y');
% 
% subplot (2,2,3);
% plot(z,tetar);
% xlabel('mm');
% ylabel('rad');
% title('Inclinação Resultante');
% 
% 
% figure ("Name", "Deflexão")
% subplot (2,2,1);
% plot(z,yx);
% xlabel('mm');
% ylabel('mm');
% title('Deflexão em X');
% 
% subplot (2,2,2);
% plot(z,yy);
% xlabel('mm');
% ylabel('mm');
% title('Deflexão em Y');
% 
% subplot (2,2,3);
% plot(z,yr);
% xlabel('mm');
% ylabel('mm');
% title('Deflexão Resultante');
% 
% 
% figure ("Name", "Momento")
% subplot (2,2,1);
% xlabel('mm');
% ylabel('N.mm');
% plot(z,Mx);
% title('Momento em x');
% 
% subplot (2,2,2);
% plot(z,My);
% xlabel('mm');
% ylabel('N.mm');
% title('Momento em Y');
% 
% subplot (2,2,3);
% plot(z,Mr);
% xlabel('mm');
% ylabel('N.mm');
% title('Momento Resultante');
% 
% figure ("Name", "Esforço")
% subplot (2,2,1);
% plot(z,Vx);
% xlabel('mm');
% ylabel('N');
% title('Esforço em X');
% 
% subplot (2,2,2);
% plot(z,Vy);
% xlabel('mm');
% ylabel('N');
% title('Esforço em Y');
% 
% subplot (2,2,3);
% plot(z,Vr);
% xlabel('mm');
% ylabel('N');
% title('Esforço Resultante');

%_________________________________________________
%Concentradores de Tensão
r_rola = 0.3; %raio da ranhura
r_esc = 0.3;  %raio da diferença de diametros

% dab = 15;
% dbc = 17;   
% dcd = 23;
% ded = 30; 
% ddl =15;


%concentrador da chaveta Escalonado A
Aa_a = 0.97098;
bb_a = (-1)*0.21796;
kt_a = Aa_a*(r_rola/dab)^bb_a;
Aas_a = 0.83425;
bbs_a = (-1)*0.21649;
kts_a = Aas_a*(r_rola/dab)^bbs_a;

%concentrador da chaveta
kt_cha = 2.8;
kts_cha = 3.6;

%concentrador escalonado B
Aa_b = 0.962523;
bb_b = (-1)*0.22823;
kt_b = Aa_b*(r_esc/dbd)^bb_b;
Aas_b = 0.83425;
bbs_b = (-1)*0.21649;
kts_b = Aas_b*(r_esc/dbd)^bbs_b;

%concentrador escalonado C
Aa_c = 0.95963;
bb_c = (-1)*0.23174;
kt_c = Aa_c*(r_esc/dde)^bb_c;
Aas_c = 0.84607;
bbs_c = (-1)*0.22863;
kts_c = Aas_c*(r_esc/dde)^bbs_c;

%concentrador escalonado D
Aa_d = 0.90879;
bb_d = (-1)*0.28598;
kt_d = Aa_d*(r_esc/dde)^bb_d;
Aas_d = 0.86331;
bbs_d = (-1)*0.23865;
kts_d = Aas_d*(r_esc/dde)^bbs_d;

%concentrador Anel Elastico
Aa_an = 0.98770;
bb_an = (-1)*0.23997;
kt_an = Aa_an*(r_esc/dde)^bb_an;
Aas_an = 0.93959;
bbs_an = (-1)*0.16789;
kts_an = Aas_an*(r_esc/dde)^bbs_an;

%%_____________________________________________________
ra = 0.093; %Raiz de a  torção
q = 1/(1+ra/((r_esc/25.4)^0.5));

ras = 0.070; %Raiz de a flexão
qs = 1/(1+ras/((r_esc/25.4)^0.5));
sig = zeros (1, n);
tal = zeros (1, n);
diametro = zeros (1, n);
Cs = zeros(1,n);
for k = 1:1:n
    if (pos(k) <= b)
        diametro(k) = dab;
        if(pos(k) == b)
           kt = kt_a;
           kts = kts_a;
        else
            kt = 1;
            kts = 1;
        end
    elseif(pos(k) <= d)
        diametro(k) = dbd;
        if(pos(k) == d)
               kt = kt_b;
               kts = kts_b;  
        elseif(pos(k) > (c-10)&& pos(k)< (c+9) ) 
               kt = kt_cha;
               kts = kts_cha;
        else
            kt = 1;
            kts = 1;
        end
    elseif (pos(k) <= e)
        diametro(k) = dde;
        if(pos(k) == e)
           kt = kt_c;
           kts = kts_c;
        else
            kt = 1;
            kts = 1;
        end
    elseif (pos(k) <= f)
        diametro(k) = def;
        if(pos(k) == f)
           kt = kt_d;
           kts = kts_d;
        else
            kt = 1;
            kts = 1;
        end
    elseif (pos(k) <= l)
        diametro(k) = dfl;
        if(pos(k) > (l-15.6)&& pos(k)< (l-14.5))
           kt = kt_an;
           kts = kts_an;
        else
            kt = 1;
            kts = 1;
        end
    end
    kf = 1+q*(kt-1);
    kfs = 1+qs*(kts-1);
    sig(k) = kf*(32*Mr(k))/(pi*diametro(k)^3);
    tal(k) = kfs*(16*T(k))/(pi*diametro(k)^3);
    Cs (k) = 1.189*diametro(k)^(-0.097);
end

sig1 = sig;
sig2 = -sig;
tal1 = tal;
tal2 = tal;

siga = abs(sig1-sig2)/2; %Tensão alternada
sigm = (sig1+sig2)/2;    %Tensão média

tala = abs(tal1-tal2)/2; %cisalhamento alternado
talm = (tal1+tal2)/2;    %cisalhamento médio

sigmv = (sigm.^2 + 3*talm.^2).^(1/2); %Tensão média equivalente de von mises
sigav = (siga.^2 + 3*tala.^2).^(1/2); %Tensão alternada equivalente de von mises 

  % Condições para Cs
  % diametros entre 8 e 250, calculados anteriomente na variavel Cs
        

%----------------------------------------------------    
    % Calculo Cf para Mpa
    % usinado
    Ai = 4.51;
    bi = -0.265;
    Cf= Ai*(Sut)^bi;
    if Cf > 1
        Cf = 1;
    end
%----------------------------------------------------    
    %Calculo Cr para confiabilidades
    Cr = 0.814; %confiabilidade = 99%

    %Calculo Cl para carga
    %não há carga axial
    Cl = 1;
    
    %Calculo Ce para ambiente
    %ambiente não salubre, corrosivo
    Ce = 1;
    
    %Calculo Ce para temperatura
    %temperatura ambiente
    Ct = 1;
%----------------------------------------------------
% Calculo tensão corrigida fadiga 
Se = 237*Cs*Cf*Cr*Cl*Ce*Ct;

Nf = (Se/sigav)*(1-(sigmv/Sut));
plot(z,Nf);
xlabel('mm');
ylabel('CS');
title('Coeficiente de Segurança');

subplot (2,2,3);
plot(z,sig);
xlabel('mm');
ylabel('MPa');
title('Tensão Resultante');

subplot (1,1,1);
plot(z,sig);
xlabel('mm');
ylabel('MPa');
title('Tensão Resultante');

subplot (1,2,2);
plot(z,T);
xlabel('mm');
ylabel('N.mm');
title('Torque Resultante');
