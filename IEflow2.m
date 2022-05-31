p0 = 100000;
pinf = 100000;
uinf = 20;
vinf = 0;
beta = 20;
CFL = 0.5;
c = 1;
rho = 1.225;
R0 = 0;
R = 0;
iter = 1;
fileID = fopen('bumpgrid.dat', 'r');
formatSpec = '%f';
sizeA = [2 Inf];
A = fscanf(fileID, formatSpec, sizeA);
[rows,columns] = size(A);
%yyx = zeros(rows, columns-2);

%for i = 1:rows
 %   for j = 2:(columns-1)
  %      if i == 1
   %         yyx(i,j-1) = -A(i,j+1) + A(i,j);
    %    else
     %       yyx(i,j-1) = A(i,j+1) - A(i,j);
      %  end
   % end
%end

for k = 1:columns
    if A(1,k) == -A(1,2)
        break
    end
end
matrix = zeros(k-1,(columns-1)/(k-1) + 1);
xnx = zeros(k-1,(columns-1)/(k-1) - 1);
xny = zeros(k-1,(columns-1)/(k-1) - 1);
ynx = zeros(k-1-1, (columns-1)/(k-1));
yny = zeros(k-1-1, (columns-1)/(k-1));
Volume = zeros(k-1-1, (columns-1)/(k-1)-1);
AreaX = zeros(k-1,(columns-1)/(k-1) - 1);
AreaY = zeros(k-1-1, (columns-1)/(k-1));
u = zeros(k, (columns-1)/(k-1)+1);
v = zeros(k, (columns-1)/(k-1)+1);
p = zeros(k, (columns-1)/(k-1)+1);
UAx1 = zeros(k-1, (columns-1)/(k-1)-1);
UAy1 = zeros(k-1-1, (columns-1)/(k-1));
UAx2 = zeros(k-1, (columns-1)/(k-1)-1);
UAy2 = zeros(k-1-1, (columns-1)/(k-1));
l1max = zeros(k-1, (columns-1)/(k-1)-1);
l2max = zeros(k-1, (columns-1)/(k-1)-1);
l3max = zeros(k-1-1, (columns-1)/(k-1));
l4max = zeros(k-1-1, (columns-1)/(k-1));
delT = zeros(k-1-1, (columns-1)/(k-1)-1);
Uiplushplus = zeros(k-1-1, (columns-1)/(k-1)-1);
Uiplushminus = zeros(k-1-1, (columns-1)/(k-1)-1);
Ujplushplus = zeros(k-1-1, (columns-1)/(k-1)-1);
Ujplushminus = zeros(k-1-1, (columns-1)/(k-1)-1);
Uiminushplus = zeros(k-1-1, (columns-1)/(k-1)-1);
Uiminushminus = zeros(k-1-1, (columns-1)/(k-1)-1);
Ujminushplus = zeros(k-1-1, (columns-1)/(k-1)-1);
Ujminushminus = zeros(k-1-1, (columns-1)/(k-1)-1);
Fi1 = zeros(k, (columns-1)/(k-1)+1);
Fj1 = zeros(k, (columns-1)/(k-1)+1);
Rij1 = zeros(k, (columns-1)/(k-1)+1);
Fi2 = zeros(k, (columns-1)/(k-1)+1);
Fj2 = zeros(k, (columns-1)/(k-1)+1);
Rij2 = zeros(k, (columns-1)/(k-1)+1);
Fi3 = zeros(k, (columns-1)/(k-1)+1);
Fj3 = zeros(k, (columns-1)/(k-1)+1);
Rij3 = zeros(k, (columns-1)/(k-1)+1);

for  i = 2:k
    for m = 1:(columns-1)/(k-1)
        matrix(i-1,m+1) = A(2, i + (m-1) * (k-1));
    end
    matrix(i-1,1) = A(1,i);
end

for i = 1 : k-1
    for j = 1 : (columns-1)/(k-1) - 1
        xnx(i,j) = matrix(i,j+1) - matrix(i,j+2);
        xny(i,j) = -matrix(i,1) + matrix(i,1);
    end
end

for j = 2 : (columns-1)/(k-1) + 1
        for i = 1 : k-1-1
            yny(i,j-1) = -matrix(i+1,1) + matrix(i,1);
            ynx(i,j-1) = matrix(i+1,j) - matrix(i,j);
        end
end


for  i = 1:k-1-1
    for j = 1:(columns-1)/(k-1)-1
        Volume(i,j) = -( 0.5 * (-yny(i,j) * xnx(i,j) - (-xny(i,j)) * ynx(i,j)) + 0.5 * (-yny(i,j+1) * xnx(i+1,j) - (-xny(i+1,j)) * ynx(i,j+1)));
    end
end

for i = 1 : k
    for j = 1:(columns-1)/(k-1)+1
        u(i,j) = uinf;
        v(i,j) = vinf;
        p(i,j) = pinf-p0;
    end 
end

for i = 1: k-1
    for j = 1:(columns-1)/(k-1) - 1
        AreaX(i,j) = sqrt(xnx(i,j)^2 + xny(i,j)^2);
    end
end

for i = 1: k-1-1
    for j = 1:(columns-1)/(k-1)
        AreaY(i,j) = sqrt(ynx(i,j)^2 + yny(i,j)^2);
    end
end

while iter >= 1
%BC
for i = 1
    for j = 2 : (columns-1)/(k-1)
        p(i,j) = p(i+1,j);
        u(i,j) = uinf;
        v(i,j) = vinf; 
    end
end
for i = k
    for j = 2 : (columns-1)/(k-1)
        p(i,j) = pinf-p0;
        u(i,j) = u(i-1,j);
        v(i,j) = v(i-1,j);
    end
end
for i = 2 : k-1
    for j = 1
        p(i,j) = p(i,j+1);
        u(i,j) = u(i,j+1) - 2 * u(i,j+1) * ynx(i-1,j)^2/(AreaY(i-1,j)^2) - 2 * v(i,j+1) * ynx(i-1,j) * yny(i-1,j)/(AreaY(i-1,j)^2);
        v(i,j) = v(i,j+1) - 2 * v(i,j+1) * yny(i-1,j)^2/(AreaY(i-1,j)^2) - 2 * u(i,j+1) * ynx(i-1,j) * yny(i-1,j)/(AreaY(i-1,j)^2);
    end
end
for i = 2 : k-1
    for j = (columns-1)/(k-1)+1
        p(i,j) = p(i,j-1);
        u(i,j) = u(i,j-1) - 2 * u(i,j-1) * ynx(i-1,j-1)^2/(AreaY(i-1,j-1)^2) - 2 * v(i,j-1) * ynx(i-1,j-1) * yny(i-1,j-1)/(AreaY(i-1,j-1)^2);
        v(i,j) = v(i,j-1) - 2 * v(i,j-1) * yny(i-1,j-1)^2/(AreaY(i-1,j-1)^2) - 2 * u(i,j-1) * ynx(i-1,j-1) * yny(i-1,j-1)/(AreaY(i-1,j-1)^2);
    end
end
%wave speed
for i = 2: k -1 
    for j = 1:(columns-1)/(k-1) - 1
             UAx2(1,j) = (u(1,j) * (-xnx(1,j)) + v(1,j) * (-xny(1,j)))/sqrt(xnx(1,j)^2 + xny(1,j)^2);
             UAx1(k-1,j) = (u(k,j) * (-xnx(k-1,j)) + v(k,j) * (-xny(k-1,j)))/sqrt(xnx(k-1,j)^2 + xny(k-1,j)^2);
             UAx1(i-1,j) = (u(i,j+1) * (-xnx(i-1,j)) + v(i,j+1) * (-xny(i-1,j)))/sqrt(xnx(i-1,j)^2 + xny(i-1,j)^2);
             UAx2(i,j) = (u(i,j+1) * (-xnx(i,j)) + v(i,j+1) * (-xny(i,j)))/sqrt(xnx(i,j)^2 + xny(i,j)^2);
             l1max(i-1,j) = abs(UAx1(i-1,j))/2 + sqrt(UAx1(i-1,j)^2 + 4 * beta^2)/2;
             l1max(k-1,j) = abs(UAx1(k-1,j))/2 + sqrt(UAx1(k-1,j)^2 + 4 * beta^2)/2;
             l2max(i,j) = abs(UAx2(i,j))/2 + sqrt(UAx2(i,j)^2 + 4 * beta^2)/2;
             l2max(1,j) = abs(UAx2(1,j))/2 + sqrt(UAx2(1,j)^2 + 4 * beta^2)/2;
    end
end
for i = 1: k -1 -1
    for j = 2:(columns-1)/(k-1) 
             UAy2(i,1) = (u(i,1) * (-ynx(i,1)) + v(i,1) * (-yny(i,1)))/sqrt(ynx(i,1)^2 + yny(i,1)^2);
             UAy1(i,(columns-1)/(k-1)) = (u(i,(columns-1)/(k-1)+1) * (-ynx(i,(columns-1)/(k-1))) + v(i,(columns-1)/(k-1)+1) * (-yny(i,(columns-1)/(k-1))))/sqrt(ynx(i,(columns-1)/(k-1))^2 + yny(i,(columns-1)/(k-1))^2);
             UAy1(i,j-1) = (u(i+1,j) * (-ynx(i,j-1)) + v(i+1,j) * (-yny(i,j-1)))/sqrt(ynx(i,j-1)^2 + yny(i,j-1)^2);
             UAy2(i,j) = (u(i+1,j) * (-ynx(i,j)) + v(i+1,j) * (-yny(i,j)))/sqrt(ynx(i,j)^2 + yny(i,j)^2);
             l3max(i,j-1) = abs(UAy1(i,j-1))/2 + sqrt(UAy1(i,j-1)^2 + 4 * beta^2)/2;
             l3max(i,(columns-1)/(k-1)) = abs(UAy1(i,(columns-1)/(k-1)))/2 + sqrt(UAy1(i,(columns-1)/(k-1))^2 + 4 * beta^2)/2;
             l4max(i,j) = abs(UAy2(i,j))/2 + sqrt(UAy2(i,j)^2 + 4 * beta^2)/2;
             l4max(i,1) = abs(UAy2(i,1))/2 + sqrt(UAy2(i,1)^2 + 4 * beta^2)/2;
    end
end

for i = 1: k-1-1
    for j = 1:(columns-1)/(k-1) - 1
        delT(i,j) = Volume(i,j)/((1/CFL) * (l1max(i,j) * AreaX(i,j) + l3max(i,j) * AreaY(i,j) + l2max(i,j) * AreaX(i+1,j) + l4max(i,j) * AreaY(i,j+1))); 
    end
end

for i = 1: k -1 -1
    for j = 1:(columns-1)/(k-1) - 1 
        Uiplushplus(i,j) = max([0 UAx2(i+1,j)]);
        Uiminushplus(i,j) = max([0 UAx2(i,j)]);

        Ujplushplus(i,j) = max([0 UAy2(i,j+1)]);
        Ujminushplus(i,j) = max([0 UAy2(i,j)]);

        Uiplushminus(i,j) = min([0 UAx1(i+1,j)]);
        Uiminushminus(i,j) = min([0 UAx1(i,j)]);
        
        Ujplushminus(i,j) = min([0 UAy1(i,j+1)]);
        Ujminushminus(i,j) = min([0 UAy1(i,j)]);
        
        %{
        if Ujminushplus(i,j) > 0 && Ujminushminus(i,j) > 0
            Ujminushminus(i,j) = 0;
        end
        
        if Uiminushplus(i,j) > 0 && Uiminushminus(i,j) > 0
            Uiminushminus(i,j) = 0;
        end
        
        if Ujplushplus(i,j) > 0 && Ujplushminus(i,j) > 0
            Ujplushminus(i,j) = 0;
        end
        
        if Uiplushplus(i,j) > 0 && Uiplushminus(i,j) > 0
            Uiplushminus(i,j) = 0;
        end
        %}
    end
end

%mass flux
for i = 1: k-1-1
    for j = 1:(columns-1)/(k-1) -1 
      Fi1(i+1,j+1) =  rho * (AreaX(i+1,j) * Uiplushplus(i,j) + AreaX(i+1,j) * Uiplushminus(i,j) - AreaX(i,j) *  Uiminushplus(i,j) - AreaX(i,j) * Uiminushminus(i,j)) - (c/(l2max(i+1,j)) * AreaX(i+1,j) * p(i+2,j+1) - c/(l2max(i+1,j)) * AreaX(i+1,j) * p(i+1,j+1) + c/(l1max(i,j)) * AreaX(i,j) * p(i,j+1) - c/(l1max(i,j)) * AreaX(i,j) * p(i+1,j+1));
      Fi1(1,j+1) = rho * AreaX(1,j) * (Uiminushplus(1,j) + Uiminushminus(1,j));
      Fi1(k,j+1) = rho * AreaX(k-1,j) * (Uiplushplus(k-1-1,j) + Uiplushminus(k-1-1,j));
      Fj1(i+1,j+1) =  rho * (AreaY(i,j+1) * Ujplushplus(i,j) + AreaY(i,j+1) * Ujplushminus(i,j) - AreaY(i,j) *  Ujminushplus(i,j) - AreaY(i,j) * Ujminushminus(i,j)) - (c/(l4max(i,j+1)) * AreaY(i,j+1) * p(i+1,j+2) - c/(l4max(i,j+1)) * AreaY(i,j+1) * p(i+1,j+1) + c/(l3max(i,j)) * AreaY(i,j) * p(i+1,j) - c/(l3max(i,j)) * AreaY(i,j) * p(i+1,j+1));

    end
end
      Rij1 = Fi1 + Fj1;
      
%momentum fluxes
for i = 1: k-1-1
    for j = 1:(columns-1)/(k-1)-1
        Fi2(i+1,j+1) = AreaX(i+1,j) * (rho * (Uiplushplus(i,j) * u(i+1,j+1) + Uiplushminus(i,j) *  u(i+2,j+1)) + 1/(2*l2max(i+1,j)) * (p(i+1,j+1) + p(i+2,j+1)) * (-xnx(i+1,j)/AreaX(i+1,j))) - AreaX(i,j) * rho * Uiminushplus(i,j) * u(i,j+1) - AreaX(i,j) * rho * Uiminushminus(i,j) *  u(i+1,j+1) - AreaX(i,j) * 1/(2*l1max(i,j))  * (p(i,j+1) + p(i+1,j+1)) * (-xnx(i,j)/AreaX(i,j));
        Fi2(1,j+1) = rho * AreaX(1,j) * (Uiminushplus(1,j) * u(1,j+1) + Uiminushminus(1,j) * u(2,j+1)) + 0.5/l2max(1,j) * (p(1,j+1) + p(2,j+1)) * -xnx(1,j);
        Fi2(k,j+1) = rho * AreaX(k-1,j) * (Uiplushplus(k-1-1,j) * u(k-1,j+1) + Uiplushminus(k-1-1,j) * u(k,j+1)) + 0.5/l1max(k-1,j) * (p(k-1,j+1) + p(k,j+1)) * -xnx(k-1,j);
        Fj2(i+1,j+1) = AreaY(i,j+1) * (rho * (Ujplushplus(i,j) * u(i+1,j+1) + Ujplushminus(i,j) *  u(i+1,j+2)) + 1/(2*l4max(i,j+1)) * (p(i+1,j+1) + p(i+1,j+2)) * (-ynx(i,j+1)/AreaY(i,j+1))) - AreaY(i,j) * rho * Ujminushplus(i,j) * u(i+1,j) - AreaY(i,j) * rho * Ujminushminus(i,j) *  u(i+1,j+1) - AreaY(i,j) * 1/(2*l3max(i,j)) * (p(i+1,j) + p(i+1,j+1)) * (-ynx(i,j)/AreaY(i,j));
    end
    Fj2(i+1,1) = 0.5/l4max(i,1) * (p(i+1,1) + p(i+1,2)) * -ynx(i,1);
    Fj2(i+1,(columns-1)/(k-1)+1) = 0.5/l3max(i,(columns-1)/(k-1)) * (p(i+1,(columns-1)/(k-1)) + p(i+1,(columns-1)/(k-1)+1)) * -ynx(i,(columns-1)/(k-1));
end
        Rij2 = Fi2 + Fj2;

for i = 1: k-1-1
    for j = 1:(columns-1)/(k-1)-1
        Fi3(i+1,j+1) = AreaX(i+1,j) * (rho * (Uiplushplus(i,j) * v(i+1,j+1) + Uiplushminus(i,j) *  v(i+2,j+1)) + 1/(2*l2max(i+1,j)) * (p(i+1,j+1) + p(i+2,j+1)) * (-xny(i+1,j)/AreaX(i+1,j))) - AreaX(i,j) * rho * Uiminushplus(i,j) * v(i,j+1) - AreaX(i,j) * rho * Uiminushminus(i,j) *  v(i+1,j+1) - AreaX(i,j) * 1/(2*l1max(i,j))  * (p(i,j+1) + p(i+1,j+1)) * (-xny(i,j)/AreaX(i,j));
        Fi3(1,j+1) = rho * AreaX(1,j) * (Uiminushplus(1,j) * v(1,j+1) + Uiminushminus(1,j) * v(2,j+1)) + 0.5/l2max(1,j) * (p(1,j+1) + p(2,j+1)) * -xny(1,j);
        Fi3(k,j+1) = rho * AreaX(k-1,j) * (Uiplushplus(k-1-1,j) * v(k-1,j+1) + Uiplushminus(k-1-1,j) * v(k,j+1)) + 0.5/l1max(k-1,j) * (p(k-1,j+1) + p(k,j+1)) * -xny(k-1,j);
        Fj3(i+1,j+1) = AreaY(i,j+1) * (rho * (Ujplushplus(i,j) * v(i+1,j+1) + Ujplushminus(i,j) *  v(i+1,j+2)) + 1/(2*l4max(i,j+1)) * (p(i+1,j+1) + p(i+1,j+2)) * (-yny(i,j+1)/AreaY(i,j+1))) - AreaY(i,j) * rho * Ujminushplus(i,j) * v(i+1,j) - AreaY(i,j) * rho * Ujminushminus(i,j) *  v(i+1,j+1) - AreaY(i,j) * 1/(2*l3max(i,j)) * (p(i+1,j) + p(i+1,j+1)) * (-yny(i,j)/AreaY(i,j));
    end
    Fj3(i+1,1) = 0.5/l4max(i,1) * (p(i+1,1) + p(i+1,2)) * -yny(i,1);
    Fj3(i+1,(columns-1)/(k-1)+1) = 0.5/l3max(i,(columns-1)/(k-1)) * (p(i+1,(columns-1)/(k-1)) + p(i+1,(columns-1)/(k-1)+1)) * -yny(i,(columns-1)/(k-1));
end
        Rij3 = Fi3 + Fj3;

%convergence
for i = 1: k-1-1
    for j = 1:(columns-1)/(k-1)-1
        p(i+1,j+1) = p(i+1,j+1) - beta^2 * delT(i,j)/Volume(i,j) * Rij1(i+1,j+1);
        u(i+1,j+1) = u(i+1,j+1) - delT(i,j)/Volume(i,j) * (Rij2(i+1,j+1)/rho - Rij1(i+1,j+1) * u(i+1,j+1)/rho);
        v(i+1,j+1) = v(i+1,j+1) - delT(i,j)/Volume(i,j) * (Rij3(i+1,j+1)/rho - Rij1(i+1,j+1) * v(i+1,j+1)/rho);
    end
end
for i = 1: k-1-1
    for j = 1:(columns-1)/(k-1) - 1
        R = R + sqrt((Rij1(i+1,j+1)/(rho * uinf))^2 + (Rij2(i+1,j+1)/(rho * uinf * uinf))^2 + (Rij3(i+1,j+1)/(rho * uinf * uinf))^2);
    end
end
M = R;
if iter == 1
    R0 = R;
end

%figure(1);
%plot(iter,R,'o');
%hold on;
if R/R0<0.001
    break
else
    iter = iter + 1;
    R=0;
end
end



U=zeros(97,49)
for i=1:97
    for j=1:49
        U(i,j)=0.25*(u(i,j)+u(i,j+1)+u(i+1,j)+u(i+1,j+1));
    end
end

V=zeros(97,49)
for i=1:97
    for j=1:49
        V(i,j)=0.25*(v(i,j)+v(i,j+1)+v(i+1,j)+v(i+1,j+1));
    end
end

P=zeros(97,49)
for i=1:97
    for j=1:49
        P(i,j)=0.25*(p(i,j)+p(i,j+1)+p(i+1,j)+p(i+1,j+1))+p0;
    end
end

X=matrix(:,1)
Y=matrix(:,2:end)
Z=0

fid = fopen( 'output.csv', 'w' );
fprintf(fid,'%6s %6s %6s %6s %6s %6s\r\n','x','y','z','u','v','p');
for j = 1:49
    for i=1:97
    fprintf( fid,'%6.2f %12.8f %6.2f %12.8f %12.8f %12.8f\n',X(i),Y(i,j),Z,U(i,j),V(i,j),P(i,j));
    end
end
fclose( fid );