function [c emin emax pc pemin pemax] = compute_cn(L,pfun)
%clf
N = size(L, 1);
Vmax = randn(size(L,1),1);
Vmin = randn(size(L,1),1);
VVmin = randn(size(L,1),1);
pVmax = randn(size(L,1),1);
pVmin = randn(size(L,1),1);
pVVmin = randn(size(L,1),1);

%gen_ev_plot(L, pfun);

for k=1:300
    Vmax = L * Vmax ;
    Vmax = Vmax - mean(Vmax);
    Vmax = Vmax / norm(Vmax) ;
    
    pVmax = pfun(L * pVmax) ;
    pVmax = pVmax - mean(pVmax);
    pVmax = pVmax / norm(pVmax) ;
    
    emax(k) = Vmax'*L*Vmax ;
    pemax(k) = pVmax'*pfun(L*pVmax) ;
end

%figure(202); clf
figure; clf;
subplot(2,2,1), plot(emax)
subplot(2,2,2), plot(pemax)

emax = Vmax'*L*Vmax ;
pemax = pVmax'*pfun(L*pVmax) ;

opemax = 0.9 / pemax;
oemax = 0.9 / emax;

tx = rand(size(L,1),1) ;
tx = tx - mean(tx) ;
b = L *tx ;
b = b - mean(b) ;
x = zeros(size(L,1),1) ;

for k=1:300
    %pVmin = pVmin + opemax * (pVmin - pfun(L*pVmin)) ;
    %if(k < 50)
    pVmin = pVmin - opemax * pfun(L*pVmin) ;
    pVmin = pVmin - mean(pVmin);
    pVmin = pVmin / norm(pVmin) ;
    pemin(k) = pVmin'*pfun(L*pVmin) ;
    %end
    
    %Vmin = Vmin + oemax * (Vmin - L*Vmin) ;
    Vmin = Vmin - oemax * L*Vmin ;
    Vmin = Vmin - mean(Vmin);
    Vmin = Vmin / norm(Vmin) ; 
    emin(k) = Vmin'*L*Vmin ;
end


subplot(2,2,3), plot(emin), title('org')
hold on

if (0)
  clear emin;
  for k=1:30
    VVmin = solve(VVmin, L, pfun) ;
    VVmin = VVmin - mean(VVmin);
    VVmin = VVmin / norm(VVmin) ;
    
    emin(k) = VVmin'*L*VVmin ;
  end
  subplot(2,2,3), plot(emin,'r'), title('org')
end;

subplot(2,2,4), plot(pemin), title('pre')

n = sqrt(N);
%figure; imagesc(reshape(abs(pVmin), n, n)); title('Min eigenvector (Ours)');

emin = Vmin'*L*Vmin ;
emin = emin(end);
pemin = pVmin'*pfun(L*pVmin) ;

emax = eigs(L, 1, 'LM');
emin = eigs(L, 2, 'SM');
emin = sort(emin, 'ascend');

if (min(abs(emin)) < 1e-5)
	emin = max(emin);
else
	emin = min(emin);
end;

c = emax / emin;
pc = pemax / pemin;

function x = solve(b,L,pfun)
x = zeros(size(b)) ;
for k=1:100
    x = x + 0.5 * (pfun(b-L*x)) ;
end

function gen_ev_plot(L, pfun)
  d = eig(full(L));
  figure; plot(sort(d)); title('Original EV');
  N = size(L, 1);
  B = zeros(N, N);
  for i=1:N
    ee = zeros(N, 1);
    ee(i) = 1;
    B(:, i) = pfun(ee);
  end;
  d = eig(L*B);
  figure; plot(sort(abs(d))); title('Precond. EV');
  
  
