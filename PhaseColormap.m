function cmap = PhaseColormap(varargin)
% Calculates and applies a colormap for displaying phase images
% blue for negative values, red for positive values
% odd length due to symmetrical distribution [-pi ... pi] of phase values
% Default length: 257
% Optional input arguments: N, type
% Type 1: following Kim's DHM review : fully saturated colours in [-pi; pi],
% colour brightness increasing when the absolute phase value is going towards 0 (going towards white)
% Type 2: same as type 1, but colour brightness decreasing when the absolute
% phase value is going towards 0 (going towards black)
% Type 3: same as type 2, but half saturation in [-pi; pi]
% Type 4: continuous phase jump (white), 0 phase black, blue or red colour
% enforced by gamma adjustement of RGB curves

    if nargin>2
        error('myApp:argChk', 'Too many input arguments')
    end
    if nargin>=1
        N = varargin{1};
    else
        N = 256;
    end
    if nargin==2
        type = varargin{2};
    else
        type = 1;
    end
    if mod(N,2)
        N = (N + 1)/2;
    else
        N = N/2 + 1;
    end
    switch type
        case 1
            c = [linspace(0,1,N);linspace(0,1,N);ones(1,N)].';
        case 2
            c = [zeros(1,N);zeros(1,N);linspace(1,0,N)].';
        case 3
            c = [linspace(0.5,0,N);linspace(0.5,0,N);linspace(1,0,N)].';
        case 4
            gamma = 2;
            c = [linspace(1,0,N).^gamma;linspace(1,0,N).^gamma;linspace(1,0,N)].';            
    end
    cmap = [c;rot90(c(1:end-1,:),2)];
end

