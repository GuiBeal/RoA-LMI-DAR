%#ok<*NBRAK,*SAGROW>

yalmip('clear');
sdpOpts = sdpsettings;
sdpOpts.solver = 'lmilab';
sdpOpts.verbose = false;

%% Dyamics
f = @(t,X) [-X(2) ; X(1) - X(2).*(1 - X(1).^2 + 0.1.*X(1).^4) ] %#ok<NOPTS>

%% Dimensions
n = 2;
nz = 3;

%% Iterate over increasing polytope
a0 = 0.1;
af = 2;
as = 0.1;

fig = figure;
grid on;
hold all;
cm = jet(ceil((af - a0)/as));
c = 1;

for a = a0:as:af
  a %#ok<NOPTS>

  % Polytope
  H = 1/a*[ 1 ,  1 , -1 , -1 ; ...
            1 , -1 , -1 ,  1 ];

  V = [ a ,  0 , -a ,  0 ; ...
        0 , -a ,  0 ,  a ];

  ne = size(H,2);

  % LMIs
  constraints = [];

  P = sdpvar(n,n,'symmetric');
  constraints = [constraints] + [P >= 0];

  Sigma = [ P           , zeros(n,n)  , zeros(n,nz)  ; ...
            zeros(n,n)  , zeros(n,n)  , zeros(n,nz)  ; ...
            zeros(nz,n) , zeros(nz,n) , zeros(nz,nz) ];

  Lambda = [ zeros(n,n)  , P           , zeros(n,nz)  ; ...
             P           , zeros(n,n)  , zeros(n,nz)  ; ...
             zeros(nz,n) , zeros(nz,n) , zeros(nz,nz) ];

  nu = sdpvar(1,1);
  constraints = [constraints] + [nu >= 0];

  L = sdpvar(2*n+nz,n+nz,'full');
  W = sdpvar(2*n+nz,n+nz,'full');
  for k = 1:ne
    R{k} = sdpvar(2*n+nz+1,2*n+nz,'full');

    h{k} = H(:,k);
    omega{k} = [ h{k}        ; ...
                 zeros(n,1)  ; ...
                 zeros(nz,1) ];

    Pi{k} = [ Sigma         , -nu*omega{k} ; ...
              -nu*omega{k}' , 2*nu - 1     ];
  end

  for j = 1:size(V,2)
    x = V(:,j);

    A1 = [ 0 , -1 ; ...
           1 , -1 ];
    A2 = [  0    , 0 , 0        ; ...
           -x(2) , 0 , 0.1*x(2) ];
    Omega1 = [ x(1) , 0 ; ...
               0    , 0 ; ...
               0    , 0 ];
    Omega2 = [ -1   ,  0   ,  0 ; ...
               x(1) , -1   ,  0 ; ...
                0   , x(1) , -1 ];

    M = [ A1     , -eye(n)     , A2     ; ...
          Omega1 , zeros(nz,n) , Omega2 ];

    N = [ M                                    , zeros(n+nz,1) ; ...
        [ -eye(n) , zeros(n,n) , zeros(n,nz) ] , x             ];

    constraints = [constraints] + [Sigma + L*M + M.'*L.' >= 0];
    constraints = [constraints] + [-Lambda + W*M + M.'*W.' >= 0];

    for k = 1:ne
      constraints = [constraints] + [Pi{k} + R{k}*N + N.'*R{k}.' >= 0];
    end
  end

  objective = trace(P);
  mSol = optimize(constraints, objective, sdpOpts);

  if mSol.problem
    break;
  end

  P = double(P) %#ok<NOPTS>

  % Plot
  step = a/20;
  x1 = -a:step:a;
  x2 = -a:step:a;
  [x1,x2] = meshgrid(x1,x2);

  r = P(1,1).*x1.^2 + P(1,2).*x1.*x2 + P(2,1).*x2.*x1 + P(2,2).*x2.^2 < 1;

  plot([ V(1,:) , V(1,1) ], [ V(2,:) , V(2,1) ], 'LineWidth', 2, 'Color', cm(c,:));
  for i = 1:size(x1,1)
    for j = 1:size(x1,2)
      if r(i,j)
        scatter(x1(i,j), x2(i,j), 'MarkerEdgeColor', cm(c,:), 'MarkerFaceColor', cm(c,:));
      end
    end
  end
  c = c+1;
end

%% Phase Portrait
step = a/20;
x1 = -a:step:a;
x2 = -a:step:a;
[x1,x2] = meshgrid(x1,x2);

x1d = -x2;
x2d = x1 - x2.*(1 - x1.^2 + 0.1.*x1.^4);

quiver(x1, x2, x1d, x2d, 'Color', 'black');

xlim(a*[-1 1]);
ylim(a*[-1 1]);
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');

saveas(fig, 'fig/exemplo2_iter.eps', 'epsc');

%% Repeat for last a
a = a -as;

%%
yalmip('clear');
sdpOpts = sdpsettings;
sdpOpts.solver = 'lmilab';
sdpOpts.verbose = false;

fig = figure;
grid on;
hold all;

% Polytope
H = 1/a*[ 1 ,  1 , -1 , -1 ; ...
          1 , -1 , -1 ,  1 ];

V = [ a ,  0 , -a ,  0 ; ...
      0 , -a ,  0 ,  a ];

ne = size(H,2);

% LMIs
constraints = [];

P = sdpvar(n,n,'symmetric');
constraints = [constraints] + [P >= 0];

Sigma = [ P           , zeros(n,n)  , zeros(n,nz)  ; ...
          zeros(n,n)  , zeros(n,n)  , zeros(n,nz)  ; ...
          zeros(nz,n) , zeros(nz,n) , zeros(nz,nz) ];

Lambda = [ zeros(n,n)  , P           , zeros(n,nz)  ; ...
           P           , zeros(n,n)  , zeros(n,nz)  ; ...
           zeros(nz,n) , zeros(nz,n) , zeros(nz,nz) ];

nu = sdpvar(1,1);
constraints = [constraints] + [nu >= 0];

L = sdpvar(2*n+nz,n+nz,'full');
W = sdpvar(2*n+nz,n+nz,'full');
for k = 1:ne
  R{k} = sdpvar(2*n+nz+1,2*n+nz,'full');

  h{k} = H(:,k);
  omega{k} = [ h{k}        ; ...
               zeros(n,1)  ; ...
               zeros(nz,1) ];

  Pi{k} = [ Sigma         , -nu*omega{k} ; ...
            -nu*omega{k}' , 2*nu - 1     ];
end

for j = 1:size(V,2)
  x = V(:,j);

  A1 = [ 0 , -1 ; ...
         1 , -1 ];
  A2 = [  0    , 0 , 0        ; ...
         -x(2) , 0 , 0.1*x(2) ];
  Omega1 = [ x(1) , 0 ; ...
             0    , 0 ; ...
             0    , 0 ];
  Omega2 = [ -1   ,  0   ,  0 ; ...
             x(1) , -1   ,  0 ; ...
              0   , x(1) , -1 ];

  M = [ A1     , -eye(n)     , A2     ; ...
        Omega1 , zeros(nz,n) , Omega2 ];

  N = [ M                                    , zeros(n+nz,1) ; ...
      [ -eye(n) , zeros(n,n) , zeros(n,nz) ] , x             ];

  constraints = [constraints] + [Sigma + L*M + M.'*L.' >= 0];
  constraints = [constraints] + [-Lambda + W*M + M.'*W.' >= 0];

  for k = 1:ne
    constraints = [constraints] + [Pi{k} + R{k}*N + N.'*R{k}.' >= 0];
  end
end

objective = trace(P);
mSol = optimize(constraints, objective, sdpOpts);

if mSol.problem
  break;
end

P = double(P) %#ok<NOPTS>

% Plot
step = a/100;
x1 = -a:step:a;
x2 = -a:step:a;
[x1,x2] = meshgrid(x1,x2);

r = P(1,1).*x1.^2 + P(1,2).*x1.*x2 + P(2,1).*x2.*x1 + P(2,2).*x2.^2 < 1;

plot([ V(1,:) , V(1,1) ], [ V(2,:) , V(2,1) ], 'LineWidth', 2, 'Color', 'black');
for i = 1:size(x1,1)
  for j = 1:size(x1,2)
    if r(i,j)
      scatter(x1(i,j), x2(i,j), 'MarkerEdgeColor', 0.75*[1 1 1], 'MarkerFaceColor', 0.75*[1 1 1]);
    end
  end
end

% Phase Portrait
d = 3;
step = d/15;
x1 = -d:step:d;
x2 = -d:step:d;
[x1,x2] = meshgrid(x1,x2);

x1d = -x2;
x2d = x1 - x2.*(1 - x1.^2 + 0.1.*x1.^4);

quiver(x1, x2, x1d, x2d, 'Color', 'black');

xlim(d*[-1 1]);
ylim(d*[-1 1]);
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex')

saveas(fig, 'fig/exemplo2.eps', 'epsc');