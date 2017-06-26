function rotR=rotz(phi)

% this routine is to generate rotation matrix
% rotating phi around z axis (counter-clockwise active, which is different from Matlab's definition)
% phi : unit in degree 
%
% Yiming Hu, 2011

phi = phi*pi/180;%convert into rad 
rotR=[	cos(phi) 	sin(phi)	0; 
	-sin(phi) 	cos(phi)	0
	0		0		1];
return
