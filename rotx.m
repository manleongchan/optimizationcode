function rotR=rotx(phi)

% this routine is to generate rotation matrix
% rotating phi around x axis (counter-clockwise active, which is different from Matlab's definition)
% phi : unit in degree 
%
% Yiming Hu, 2011

phi = phi*pi/180;%convert into rad 
rotR=[	1 		0 		0;
	0 		cos(phi) 	sin(phi); 
	0 		-sin(phi) 	cos(phi)];
return
