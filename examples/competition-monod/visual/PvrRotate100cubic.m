function [ POV_RAY ] = PvrRotate100cubic(image, rotate)

%==================================



% WRITE import.pov FILE FOR POVRAY
fid=fopen('import.pov','w');

fprintf(fid,'#include "shapes.inc"\n');

fprintf(fid,'// Right-handed coordinate system in which the z-axis points upwards\n');
fprintf(fid,'camera {\n');
fprintf(fid,['location <0, 60e-2, 10e-2>\n']);
fprintf(fid,'sky z\n');
fprintf(fid,'right 0.24*x*image_width/image_height\n');
fprintf(fid,'up 0.24*y\n');
fprintf(fid,'look_at <0, 0, 0>\n');
fprintf(fid,'look_at 0\n')
fprintf(fid,'rotate <0, 0, 20>\n', rotate);
fprintf(fid,'}\n');
fprintf(fid,'// Create simualtion domain\n');
fprintf(fid,'#declare b_x=0.05;\n');
fprintf(fid,'#declare b_y=0.02;\n');
fprintf(fid,'#declare b_z=0.05;\n');
fprintf(fid,'object{ \n');
fprintf(fid,['Wire_Box(<b_x, b_y, b_z>, <-b_x, -b_y, -b_z>, 0.0002, 0) \n']);
fprintf(fid,[' texture{ pigment{ color rgb<1,1,1>}}\n']);
fprintf(fid,'}\n');
fprintf(fid,'// White background\n');
fprintf(fid,'background{rgb 1}\n');
fprintf(fid,'// Two lights with slightly different colors\n');
fprintf(fid,'light_source{<-10e-2,80e-2,0e-2> color rgb <0.77,0.75,0.75>}\n');
fprintf(fid,'light_source{<10e-2,80e-2,0e-2> color rgb <0.57,0.55,0.55>}\n');
fprintf(fid,'// Radius of the Voronoi cell network\n');
fprintf(fid,'#declare r=0.05;\n');
fprintf(fid,'// Radius of the particles\n');
fprintf(fid,'#declare s=0.6;\n');
fprintf(fid,'// Particles\n');
fprintf(fid,'union{\n');
fprintf(fid,'#include "data_p.pov"\n');
fprintf(fid,'}\n');



fclose(fid);

POV_RAY=1


end





