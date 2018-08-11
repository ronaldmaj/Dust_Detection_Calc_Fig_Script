function w = faddeeva1(z,N)
% FADDEEVA   Faddeeva function
%   W = FADDEEVA(Z) is the Faddeeva function, aka the plasma dispersion
%   function, for each element of Z. The Faddeeva function is defined as:
%
%     w(z) = exp(-z^2) * erfc(-j*z)
%
%   where erfc(x) is the complex complementary error function.
%
%   W = FADDEEVA(Z,N) can be used to explicitly specify the number of terms
%   to truncate the expansion (see (13) in [1]). N = 16 is used as default.
%
%   Example:
%       x = linspace(-10,10,1001); [X,Y] = meshgrid(x,x); 
%       W = faddeeva(complex(X,Y)); 
%       figure; 
%       subplot(121); imagesc(x,x,real(W)); axis xy square; caxis([-1 1]); 
%       title('re(faddeeva(z))'); xlabel('re(z)'); ylabel('im(z)'); 
%       subplot(122); imagesc(x,x,imag(W)); axis xy square; caxis([-1 1]);
%       title('im(faddeeva(z))'); xlabel('re(z)'); ylabel('im(z)'); 
%
%   Reference:
%   [1] J.A.C. Weideman, "Computation of the Complex Error Function," SIAM
%       J. Numerical Analysis, pp. 1497-1518, No. 5, Vol. 31, Oct., 1994 
%       Available Online: http://www.jstor.org/stable/2158232

if nargin<2, N = []; end
if isempty(N), N = 16; end

w = zeros(size(z)); % initialize output

M = 2*N;
M2 = 2*M;
k = (-M+1:1:M-1)'; % M2 = no. of sampling points.
L = sqrt(N/sqrt(2)); % Optimal choice of L.

theta = k*pi/M;
t = L*tan(theta/2); % Variables theta and t.
f = exp(-t.^2).*(L^2+t.^2);
f = [0; f]; % Function to be transformed.
a = real(fft(fftshift(f)))/M2; % Coefficients of transform.
a = flipud(a(2:N+1)); % Reorder coefficients.

Z = (L+1i*z)./(L-1i*z);
p = polyval(a,Z); % Polynomial evaluation.
w = 2*p./(L-1i*z).^2 + (1/sqrt(pi))./(L-1i*z); % Evaluate w(z).
