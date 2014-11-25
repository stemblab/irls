function [un,iter,epsilon] = irls(Phi,b,varargin)

% $\epsilon$-regularized IRLS solution to $\Phi u =b$ 
% See http://stemblab.github.io/irls
 
% $\LaTeX$ definitions
% brackets: $$\def\lr#1{\left(#1\right)}$$

% algorithm defaults

p=0; % $||\cdot||^p$
thresh=0.1; % set $u_i$ to zero if less than thresh
max_eps = 1e-8; % accuracy stopping condition
max_iter = 200; % iteration stopping condition
eps_scale = 10; % factor to scale epsilon if eps_thresh condition met
eps_thresh = 100;

% function call defaults

if nargin > 2, p = varargin{1};end
if nargin > 3, thresh = varargin{2};end
if nargin > 4, max_eps = varargin{3};end
if nargin > 5, max_iter = varargin{4};end
if nargin > 6, eps_scale = varargin{5};end
if nargin > 7, eps_thresh = varargin{6};end

% initialization

epsilon=1; % $\epsilon$ 
ul=Phi\b; % "last $u$", $u^{(n-1)}$ = 2-norm solution of $\Phi u=b$
m=size(ul,1);
Qn=zeros(m,m); % $Q_n$

% iteration

for iter=1:max_iter

    % stop if error is small
    
    if epsilon<max_eps; break; end
    
    % main algorithm

    % $w_i=\lr{\lr{u_i^{(n-1)}}^2+\epsilon}^{p/2-1}$ 
    
    
    % $Q_i=1/w_i$
    
    Qn(logical(eye(m)))= 1./((abs(ul).^2+epsilon).^(p/2-1));
    
    % $u^{(n)}=Q_n\Phi^T \lr{\Phi Q_n \Phi^T}^{-1}b$ 
    
    
    un=Qn*Phi'*inv(Phi*Qn*Phi')*b; % $u^{(n)}$
    
    % set small $u^{(n)}$ to zero
    
    un(find(abs(un)<thresh))=0;

    % scale $\epsilon$

    if norm(un-ul)/norm(un)<sqrt(epsilon)/eps_thresh
        epsilon=epsilon/eps_scale;
    end

    ul=un; % last u = this u: $u^{(n-1)}=u^{(n)}$

end

%!end (3)

