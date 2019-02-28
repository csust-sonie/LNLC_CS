% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%


%function reconstructed_image = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols,block_size)
function reconstructed_image = DCT_Decoder(y, A,At, num_rows, num_cols,block_size)

%[M N] = size(Phi);

%M=round(32*rate);
%A = @(z) A_DCTB(z,Phi);
%At = @(z) At_DCTB(z,Phi);

[NN MM]=size(y);
Psi = DCT2D_Matrix(block_size);

lambda = 6;
TOL = 0.0001;
D_prev = 0;

num_factor = 0;
max_iterations = 200;

%x = Phi' * y;
x=At(y);
%for i =1:MM
%  x(:,i)=At(y(:,i));
%end;

for i = 1:max_iterations
 % [x D] = SPLIteration(y, x, Phi, Psi, ...
 %     block_size, num_rows, num_cols, lambda);
  [x D] = SPLIteration(y, x, A,At, Psi, ...
      block_size, num_rows, num_cols, lambda);
  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
    if num_factor == 4;
      break;
    end
    lambda = lambda * 0.6;
    num_factor = num_factor + 1;
  end
  D_prev = D;
end

reconstructed_image = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distict');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%function [x D] = SPLIteration(y, x, Phi, Psi, ...
%    block_size, num_rows, num_cols, lambda)
function [x D] = SPLIteration(y, x, A, At, Psi,...
    block_size, num_rows, num_cols, lambda)


 
x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct'); 
%imshow(x,[0 255]);
[NN MM]=size(y);
x_hat = wiener2(x, [3, 3]);
x_hat = im2col(x_hat, [block_size block_size], 'distinct');

%x_hat = x_hat + Phi' * (y - Phi * x_hat);
x_hat = x_hat + At(y - A(x_hat));

x1 = col2im(x_hat, [block_size block_size], ...
    [num_rows num_cols], 'distinct');

x_check = Psi' * x_hat;
threshold = lambda * sqrt(2 * log(num_rows * num_cols)) * ...
    (median(abs(x_check(:))) / 0.6745);
x_check(abs(x_check) < threshold) = 0;

x_bar = Psi * x_check;
%x = x_bar + Phi' * (y - Phi * x_bar);
x=x_bar+At(y-A(x_bar));
%x = x_bar + At (y - A(x_bar));
x2 = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct');

D = RMS(x1, x2);