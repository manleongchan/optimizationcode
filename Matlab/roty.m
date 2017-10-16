function rotR=roty(phi)

% this routine is to generate rotation matrix
% rotating phi around y axis (counter-clockwise active, which is different from Matlab's definition)
% phi : unit in degrees 
%
% Yiming Hu, 2011

phi = phi*pi/180;%convert into rad 
rotR=[	cos(phi)	0 		-sin(phi);
	0 	 	1		0; 
	sin(phi)	0 		cos(phi)];
return
